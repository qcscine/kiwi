/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_COULOMBEXCHANGEMATRICES_H
#define KIWI_COULOMBEXCHANGEMATRICES_H

#include <Kiwi/KiwiUtils/Data.h>
#include <LibintIntegrals/LibintIntegrals.h>

namespace Scine {
namespace Kiwi {
namespace FockMatrix {

//! In-memory evaluation of the J and K matrices.
inline auto buildJKmatrices(const Utils::DensityMatrix& D, const std::shared_ptr<Molecule>& molecule,
                            const std::shared_ptr<Data>& data, Utils::ElementType type)
    -> std::array<Utils::SpinAdaptedMatrix, 2> {
  Utils::SpinAdaptedMatrix J;
  Utils::SpinAdaptedMatrix K;

  if (molecule->at(type).isRestricted) {
    Eigen::MatrixXd Dt = D.restrictedMatrix().transpose();

    Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dt.data(), Dt.rows() * Dt.cols()));
    {
      Eigen::VectorXd tmpVec = data->Coulomb[Data::getUnique(type, type)] * coefficientVector;
      J.setRestrictedMatrix(Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule->at(type).LAO, molecule->at(type).LAO));
    }
    {
      Eigen::VectorXd tmpVec = data->Exchange[type] * coefficientVector;
      K.setRestrictedMatrix(Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule->at(type).LAO, molecule->at(type).LAO));
    }
  }

  else {
    if (molecule->at(type).N < 2) {
      auto dim = molecule->at(type).LAO;
      J.setAlphaMatrix(Eigen::MatrixXd::Zero(dim, dim));
      J.setBetaMatrix(Eigen::MatrixXd::Zero(dim, dim));
      K.setAlphaMatrix(Eigen::MatrixXd::Zero(dim, dim));
      K.setBetaMatrix(Eigen::MatrixXd::Zero(dim, dim));
    }

    Eigen::MatrixXd Dalpha = D.alphaMatrix().transpose();
    Eigen::MatrixXd Dbeta = D.betaMatrix().transpose();

    // Coulomb alpha
    {
      Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dalpha.data(), Dalpha.rows() * Dalpha.cols()));
      Eigen::VectorXd tmpVec = data->Coulomb[Data::getUnique(type, type)] * coefficientVector;
      J.setAlphaMatrix(Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule->at(type).LAO, molecule->at(type).LAO));
    }
    // Coulomb beta
    if (molecule->at(type).msVector[1] > 0) {
      Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dbeta.data(), Dbeta.rows() * Dbeta.cols()));
      Eigen::VectorXd tmpVec = data->Coulomb[Data::getUnique(type, type)] * coefficientVector;
      J.setBetaMatrix(Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule->at(type).LAO, molecule->at(type).LAO));
    }
    // Alpha exchange
    if (molecule->at(type).msVector[0] > 0) {
      Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dalpha.data(), Dalpha.rows() * Dalpha.cols()));
      Eigen::VectorXd tmpVec = data->Exchange[type] * coefficientVector;
      K.setAlphaMatrix(Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule->at(type).LAO, molecule->at(type).LAO));
    }
    // Beta exchange
    if (molecule->at(type).msVector[1] > 0) {
      Eigen::VectorXd coefficientVector(Eigen::Map<Eigen::VectorXd>(Dbeta.data(), Dbeta.rows() * Dbeta.cols()));
      Eigen::VectorXd tmpVec = data->Exchange[type] * coefficientVector;
      K.setBetaMatrix(Eigen::Map<Eigen::MatrixXd>(tmpVec.data(), molecule->at(type).LAO, molecule->at(type).LAO));
    }
  }

  return {J, K};
}

//! Integral-direct evaluation of the J and K matrices.
inline auto buildJKmatricesDirect(const Utils::DensityMatrix& D, const std::shared_ptr<Molecule>& molecule,
                                  Utils::ElementType type, bool formIncremental = false, std::size_t numRebuilds = 0)
    -> std::array<Utils::SpinAdaptedMatrix, 2> {
  Utils::SpinAdaptedMatrix J;
  Utils::SpinAdaptedMatrix K;

  const auto& densityMatrix = D;

  if (D.restricted() && (D.restrictedMatrix().cwiseAbs().maxCoeff() < 1e-16)) {
    J.restrictedMatrix().resizeLike(D.restrictedMatrix());
    J.restrictedMatrix().setZero();
    K.restrictedMatrix().resizeLike(D.restrictedMatrix());
    K.restrictedMatrix().setZero();
    return {J, K};
  }

  if (!D.restricted() && (D.restrictedMatrix().cwiseAbs().maxCoeff() < 1e-16)) {
    J.alphaMatrix().resizeLike(D.alphaMatrix());
    J.alphaMatrix().setZero();
    K.alphaMatrix().resizeLike(D.alphaMatrix());
    K.alphaMatrix().setZero();
    J.betaMatrix().resizeLike(D.betaMatrix());
    J.betaMatrix().setZero();
    K.betaMatrix().resizeLike(D.betaMatrix());
    K.betaMatrix().setZero();
    return {J, K};
  }

  auto& basis = molecule->at(type).basisSet;

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Coulomb;

  double prescreeningThreshold = 1e-8 / (3.0 * double(basis.nbf()));

  if (formIncremental) {
    for (auto i = 0UL; i < numRebuilds; ++i) {
      prescreeningThreshold /= 100;
    }
  }

  auto result =
      Integrals::LibintIntegrals::evaluateTwoBodyDirectBo(specifier, basis, basis, densityMatrix, prescreeningThreshold);

  return {result.first, result.second};
}

} // namespace FockMatrix
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_COULOMBEXCHANGEMATRICES_H
