/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/HartreeFock/InitialGuess/SND.h>
#include <Kiwi/HartreeFock/Projector/Projector.h>
#include <Kiwi/HartreeFock/SecondOrder.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace Kiwi {

SNDGuess::SNDGuess(Utils::AtomCollection positions, const Utils::Integrals::BasisSet& nuclearBasis,
                   const Utils::Integrals::BasisSet& electronicBasis, int na, int nb)
  : positions_(std::move(positions)), size_(positions_.size()), type_(positions_[0].getElementType()), n_alpha(na), n_beta(nb) {
  makeBasis(nuclearBasis, electronicBasis);

  offest_.resize(size_);
  for (int i = 0; i < int(size_); ++i) {
    offest_[i] = i * dim_;
  }

  D_.getAlphaMatrix().resize(size_ * dim_, size_ * dim_);
  D_.getAlphaMatrix().setZero();
  C_.alphaMatrix().resize(size_ * dim_, size_ * dim_);
  C_.alphaMatrix().setZero();
  D_.getBetaMatrix().resize(size_ * dim_, size_ * dim_);
  D_.getBetaMatrix().setZero();
  C_.betaMatrix().resize(size_ * dim_, size_ * dim_);
  C_.betaMatrix().setZero();
}

auto SNDGuess::makeBasis(const Utils::Integrals::BasisSet& nuclearBasis, const Utils::Integrals::BasisSet& electronicBasis)
    -> void {
  atom_ = Utils::AtomCollection();
  atom_.push_back(Utils::Atom(positions_[0].getElementType(), {0, 0, 0}));

  auto nuclearShellIndices = nuclearBasis.atomToShell(positions_);
  auto electronicShellIndices = electronicBasis.atomToShell(positions_);

  nuclearBasis_ = Utils::Integrals::BasisSet(atom_);
  electronicBasis_ = Utils::Integrals::BasisSet(atom_);

  // All shells must be the same
  for (const auto& index : nuclearShellIndices.at(0)) {
    nuclearBasis_.push_back(nuclearBasis.at(index));
  }
  for (const auto& index : electronicShellIndices.at(0)) {
    electronicBasis_.push_back(electronicBasis.at(index));
  }

  if (electronicBasis_.size() == 1 || nuclearBasis_.size() == 1) {
    throw std::runtime_error("SND guess not possible in minimal basis.");
  }

  nuclearBasis_.move({0, 0, 0});
  nuclearBasis_.setPureSpherical(nuclearBasis.isPureSpherical());
  Integrals::LibintIntegrals::generateShellPairs(nuclearBasis_);

  electronicBasis_.move({0, 0, 0});
  electronicBasis_.setPureSpherical(electronicBasis.isPureSpherical());
  Integrals::LibintIntegrals::generateShellPairs(electronicBasis_);

  dim_ = nuclearBasis_.nbf();
}

auto SNDGuess::getDensity() const -> Utils::DensityMatrix {
  return D_;
}

auto SNDGuess::getMOs() const -> Utils::MolecularOrbitals {
  return C_;
}

auto SNDGuess::setVerbose(bool verbose) -> void {
  verbose_ = verbose;
}

auto SNDGuess::runScf() -> void {
  std::unordered_map<Utils::ElementType, Utils::Integrals::BasisSet> basisSetMap;
  basisSetMap[Utils::ElementType::E] = electronicBasis_;
  basisSetMap[type_] = nuclearBasis_;

  std::pair<Utils::AtomCollection, std::vector<bool>> structureBoolVectorPair;
  structureBoolVectorPair.first = atom_;
  structureBoolVectorPair.second = {true};
  MoleculeSettings molSettings;
  molSettings.charge = 0;
  molSettings.multiplicity = Utils::ElementInfo::Z(type_) % 2 + 1;
  molSettings.useHighSpinApproximation = true;
  molSettings.isRestricted = false;

  auto mol = std::make_shared<Kiwi::Molecule>(basisSetMap, structureBoolVectorPair, molSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);
  data->oneBodyIntegrals(false);
  data->makeLoewdinOrtho(10e-10, false);
  data->twoBodyIntegrals(false);
  data->makeExchange(false);

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::Core;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::Core;
  settings.neScfType = Kiwi::NuclearElectronicSCF::TRAH;

  HartreeFockMain::makeGuess(mol, data, settings, false);

  Kiwi::SecondOrder<Kiwi::SymmetryType::None> scf(mol, data, settings, ProjectionParameters(), false);
  scf.run();

  for (auto i = 0; i < int(size_); ++i) {
    if (n_beta == 0) {
      D_.getAlphaMatrix().block(offest_[i], offest_[i], dim_, dim_) = data->D.at(type_).alphaMatrix();
      C_.alphaMatrix().block(offest_[i], offest_[i], dim_, dim_) = data->C.at(type_).alphaMatrix();
    }
    else {
      if (i % 2 == 0) {
        D_.getAlphaMatrix().block(offest_[i], offest_[i], dim_, dim_) = data->D.at(type_).alphaMatrix();
        C_.alphaMatrix().block(offest_[i], offest_[i], dim_, dim_) = data->C.at(type_).alphaMatrix();
      }
      else {
        D_.getBetaMatrix().block(offest_[i], offest_[i], dim_, dim_) = data->D.at(type_).alphaMatrix();
        C_.betaMatrix().block(offest_[i], offest_[i], dim_, dim_) = data->C.at(type_).alphaMatrix();
      }
    }
  }

  finish();
}

auto SNDGuess::finish() -> void {
  D_.getRestrictedMatrix() = D_.getAlphaMatrix() + D_.getBetaMatrix();
  C_ = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients(C_.alphaMatrix(), C_.betaMatrix());
}

} // namespace Kiwi
} // namespace Scine
