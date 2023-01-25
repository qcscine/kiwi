/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/InitialGuess/SAD.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <Utils/DataStructures/IntegralSpecifier.h>
#include <Utils/DataStructures/Shell.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Scf/ConvergenceAccelerators/FockDiis.h>
#include <Eigen/Eigenvalues>
#include <utility>

namespace Scine {
namespace Kiwi {

SADGuess::SADGuess(Utils::AtomCollection geometry, bool isRestricted, const Utils::Integrals::BasisSet& basis)
  : geometry_(std::move(geometry)), isRestricted_(isRestricted), basis_(basis), dim_(basis_.nbf()) {
  D_.getRestrictedMatrix().resize(dim_, dim_);
  D_.getRestrictedMatrix().setZero();
  C_.restrictedMatrix().resize(dim_, dim_);
  C_.restrictedMatrix().setZero();
  energies_.resize(dim_);

  computeOccupationNumberVector();
}

auto SADGuess::getDensity() const -> Utils::DensityMatrix {
  return D_;
}

auto SADGuess::getMOs() const -> Utils::MolecularOrbitals {
  return C_;
}

auto SADGuess::getEnergies() const -> Eigen::VectorXd {
  return energies_;
}

auto SADGuess::setVerbose(bool verbose) -> void {
  verbose_ = verbose;
}

auto SADGuess::computeInitialDensityMatrices() -> void {
  for (const auto& elem : occNumMap_) {
    computeInitialDensityMatrixAt(elem.first);
  }
}

auto SADGuess::runScf() -> void {
  for (const auto& elem : atomIndexMap_) {
    auto type = elem.first;
    if (type == Utils::ElementType::H || type == Utils::ElementType::D || type == Utils::ElementType::T ||
        type == Utils::ElementType::H1) {
      evaluateProton(type);
    }
    else {
      atomicScfAt(type);
    }
  }

  finish();
}

auto SADGuess::finish() -> void {
  if (isRestricted_) {
    C_ = Utils::MolecularOrbitals::createFromRestrictedCoefficients(C_.restrictedMatrix());
  }
  else {
    C_ = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients(C_.restrictedMatrix(), C_.restrictedMatrix());
    D_.getAlphaMatrix() = 0.5 * D_.restrictedMatrix();
    D_.getBetaMatrix() = 0.5 * D_.restrictedMatrix();
  }
}

auto SADGuess::numElectronsToOccupation(Utils::ElementType type) -> void {
  auto numElectrons = Utils::ElementInfo::Z(type);

  // Aufbau principle
  // 1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p
  std::vector<std::size_t> sizes{1, 1, 3, 1, 3, 1, 5, 3, 1, 5, 3, 1, 7, 5, 3};

  for (auto size : sizes) {
    // Every orbital can be doubly occupied
    unsigned numElectronsToFillShell = 2 * size;
    if (numElectronsToFillShell <= numElectrons) {
      // If there are enough electrons to fill the shell, do that:
      for (auto j = 0UL; j < size; ++j) {
        occNumMap_.at(type).push_back(2);
      }
      numElectrons -= numElectronsToFillShell;
    }
    else {
      // If not: partially fill the shells to achieve spherical symmetry
      // Suppose we have Boron, the result will be
      // 1s^2 2s^2 2p_x^1/2 2p_y^1/2 2p_z^1/3
      for (auto j = 0UL; j < size; ++j) {
        occNumMap_.at(type).push_back(double(numElectrons) / double(size));
      }
      numElectrons = 0;
    }

    if (numElectrons == 0) {
      break;
    }
  }
}

auto SADGuess::computeOccupationNumberVector() -> void {
  const auto& atoms = geometry_;

  shellIndices_ = basis_.atomToShell(atoms);

  size_ = atoms.size();
  offest_.resize(size_);
  std::vector<int> orbitalsPerAtomVec(size_);

  // Loop over atoms:
  for (int i = 0; i < atoms.size(); ++i) {
    auto type = atomIndexToParticleType(i);

    if (atomIndexMap_.find(type) == atomIndexMap_.end()) {
      atomIndexMap_[type] = {};
    }
    atomIndexMap_[type].push_back(i);

    orbitalsPerAtomMap_[type] = 0;
    orbitalsPerAtomVec[i] = 0;

    //
    // Loop over shells that belong to the atom.
    //
    for (const auto& index : shellIndices_.at(i)) {
      auto size = basis_.at(index).size();
      orbitalsPerAtomMap_[type] += size;
      orbitalsPerAtomVec[i] += size;
    }

    // Compute the offset
    offest_[i] = 0;
    for (int j = 0; j < i; ++j) {
      offest_[i] += orbitalsPerAtomVec[j];
    }

    // Now make occupation number matrix:
    if (occNumMat_.find(type) == occNumMat_.end()) {
      basisMap_[type] = Utils::Integrals::BasisSet();

      for (auto const& index : shellIndices_.at(i)) {
        basisMap_[type].push_back(basis_.at(index));
      }

      occNumMap_[type] = std::vector<double>();
      occNumMat_[type] = Eigen::MatrixXd();
      occNumMat_[type].resize(orbitalsPerAtomMap_[type], orbitalsPerAtomMap_[type]);
      occNumMat_[type].setZero();

      numElectronsToOccupation(type);

      for (auto j = 0UL; j < occNumMap_[type].size(); ++j) {
        occNumMat_[type](j, j) = occNumMap_[type][j];
      }
    }
  }
}

auto SADGuess::computeInitialDensityMatrixAt(Utils::ElementType type) -> void {
  Utils::AtomCollection atoms;
  Utils::Atom atom(type, basisMap_.at(type).at(0).getShift());
  atoms.push_back(atom);

  Integrals::LibintIntegrals::generateShellPairs(basisMap_[type]);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Kinetic;
  auto resultMapKin = Integrals::LibintIntegrals::evaluate(specifier, basisMap_[type], basisMap_[type]);
  const auto& Kin = resultMapKin[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  specifier.op = Utils::Integrals::Operator::PointCharges;
  specifier.atoms = atoms;
  auto resultMapPC = Integrals::LibintIntegrals::evaluate(specifier, basisMap_[type], basisMap_[type]);
  const auto& pointChrgs = resultMapPC[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  specifier.op = Utils::Integrals::Operator::Overlap;
  auto resultMapS = Integrals::LibintIntegrals::evaluate(specifier, basisMap_[type], basisMap_[type]);

  overlapMap_[type] = resultMapS[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  coreMatrixMap_[type] = Kin + pointChrgs;

  Eigen::VectorXd redundant;

  updateDensity(coreMatrixMap_[type], overlapMap_[type], occNumMat_[type],
                densityMatrixMap_[type].getRestrictedMatrix(), mosMap_[type].restrictedMatrix(), redundant);
}

auto SADGuess::evaluateProton(Utils::ElementType type) -> void {
  auto dim = orbitalsPerAtomMap_[type];

  auto& D = densityMatrixMap_[type].getRestrictedMatrix();
  auto& C = mosMap_[type].restrictedMatrix();
  const auto& S = overlapMap_[type];
  const auto& H = coreMatrixMap_[type];
  const auto& occupation = occNumMat_[type];
  Eigen::VectorXd orbitalEnergies(dim);
  updateDensity(H, S, occupation, D, C, orbitalEnergies);

  for (const auto& atomIndex : atomIndexMap_[type]) {
    D_.getRestrictedMatrix().block(offest_[atomIndex], offest_[atomIndex], orbitalsPerAtomMap_[type],
                                   orbitalsPerAtomMap_[type]) = D;
    C_.restrictedMatrix().block(offest_[atomIndex], offest_[atomIndex], orbitalsPerAtomMap_[type],
                                orbitalsPerAtomMap_[type]) = C;
    energies_.middleRows(offest_[atomIndex], orbitalsPerAtomMap_[type]) = orbitalEnergies;
  }
}

auto SADGuess::atomicScfAt(Utils::ElementType type) -> void {
  auto dim = orbitalsPerAtomMap_[type];

  auto coulomb = makeCoulomb(type);
  auto exchange = makeExchange(coulomb);
  auto& D = densityMatrixMap_[type].getRestrictedMatrix();
  auto& C = mosMap_[type].restrictedMatrix();
  const auto& S = overlapMap_[type];
  const auto& H = coreMatrixMap_[type];
  const auto& occupation = occNumMat_[type];
  Eigen::MatrixXd J(dim, dim);
  Eigen::MatrixXd K(dim, dim);
  Eigen::VectorXd orbitalEnergies(dim);
  Utils::SpinAdaptedMatrix Fock;
  Fock.resize(dim);
  auto& F = Fock.restrictedMatrix();

  double currentEnergy = 1e10;
  // initial fock matrix building
  buildFockMatrix(F, J, K, coulomb, exchange, H, D);
  double lastEnergy = computeEnergy(H, F, D);

  int maxIter = 100;
  double thresh = 1e-7;
  bool converged = false;

  Utils::FockDiis diis;
  diis.setNAOs(dim);
  diis.setOrthogonal(false);
  diis.setOverlapMatrix(S);
  diis.setSubspaceSize(5);
  diis.setUnrestricted(false);

  int it = 0;
  for (; it < maxIter; ++it) {
    updateDensity(F, S, occupation, D, C, orbitalEnergies);
    buildFockMatrix(F, J, K, coulomb, exchange, H, D);
    currentEnergy = computeEnergy(H, F, D);

    if (std::abs(lastEnergy - currentEnergy) < thresh) {
      converged = true;
      break;
    }

    lastEnergy = currentEnergy;

    diis.addMatrices(Fock, densityMatrixMap_[type]);
    Fock = diis.getMixedFockMatrix();
  }

  if (!converged) {
    std::cout << "Atomic scf not converged!" << std::endl;
  }

  for (const auto& atomIndex : atomIndexMap_[type]) {
    D_.getRestrictedMatrix().block(offest_[atomIndex], offest_[atomIndex], orbitalsPerAtomMap_[type],
                                   orbitalsPerAtomMap_[type]) = D;
    C_.restrictedMatrix().block(offest_[atomIndex], offest_[atomIndex], orbitalsPerAtomMap_[type],
                                orbitalsPerAtomMap_[type]) = C;
    energies_.middleRows(offest_[atomIndex], orbitalsPerAtomMap_[type]) = orbitalEnergies;
  }
}

auto SADGuess::makeCoulomb(Utils::ElementType type) -> Eigen::MatrixXd {
  Utils::Integrals::IntegralSpecifier specifier;

  specifier.op = Utils::Integrals::Operator::Coulomb;
  auto resultMap = Integrals::LibintIntegrals::evaluate(specifier, basisMap_[type], basisMap_[type]);

  return resultMap[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
}

auto SADGuess::makeExchange(Eigen::MatrixXd& coulomb) -> Eigen::MatrixXd {
  Eigen::MatrixXd exchange = Matrix::Zero(coulomb.rows(), coulomb.cols());

  auto dim = static_cast<int>(std::sqrt(coulomb.rows()));

  for (auto p = 0; p < dim; ++p) {
    for (auto q = 0; q < dim; ++q) {
      exchange.block(p * dim, q * dim, dim, dim) = coulomb.block(p * dim, q * dim, dim, dim).transpose();
    }
  }

  return exchange;
}

auto SADGuess::buildFockMatrix(Eigen::MatrixXd& F, Eigen::MatrixXd& J, Eigen::MatrixXd& K, const Eigen::MatrixXd& coulomb,
                               const Eigen::MatrixXd& exchange, const Eigen::MatrixXd& H, Eigen::MatrixXd& D) -> void {
  Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(D.data(), D.rows() * D.cols()));
  {
    Eigen::VectorXd tmpVec = coulomb * coefficientVector;
    J = Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), J.rows(), J.cols());
  }
  {
    Eigen::VectorXd tmpVec = exchange * coefficientVector;
    K = Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), K.rows(), K.cols());
  }

  F = H + J - 0.5 * K;
}

auto SADGuess::computeEnergy(const Eigen::MatrixXd& H, const Eigen::MatrixXd& F, const Eigen::MatrixXd& D) -> double {
  return 0.5 * ((D * (H + F)).trace());
}

auto SADGuess::updateDensity(const Eigen::MatrixXd& F, const Eigen::MatrixXd& S, const Eigen::MatrixXd& occupation,
                             Eigen::MatrixXd& D, Eigen::MatrixXd& C, Eigen::VectorXd& energies) -> void {
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(F, S);
  C = ges.eigenvectors();
  energies = ges.eigenvalues();
  D = C * occupation * C.transpose();
}

} // namespace Kiwi
} // namespace Scine
