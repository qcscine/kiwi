/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/SimultaneousDiis/SimultaneousDiis.h>
#include <Eigen/QR>
#include <algorithm>
#include <iostream>

namespace Scine {
namespace Kiwi {

SimultaneousDiis::SimultaneousDiis() {
  setSubspaceSize(5);
}

void SimultaneousDiis::setSubspaceSize(int n) {
  bool resizeNeeded = n != subspaceSize_;

  subspaceSize_ = n;

  if (resizeNeeded) {
    resizeMembers();
  }
}

void SimultaneousDiis::resizeMembers() {
  fockMatrices.resize(subspaceSize_);
  errorMatrices.resize(subspaceSize_);

  B = Eigen::MatrixXd::Ones(subspaceSize_ + 1, subspaceSize_ + 1) * (-1);
  B(0, 0) = 0;

  rhs = Eigen::VectorXd::Zero(subspaceSize_ + 1);
  rhs(0) = -1;

  restart();
}

void SimultaneousDiis::restart() {
  C = Eigen::VectorXd::Zero(subspaceSize_ + 1);
  iterationNo_ = 0;
  index_ = 0;
}

void SimultaneousDiis::addMatrices() {
  iterationNo_++;
  lastAdded_ = index_;

  if (useOrthonormalBasis) {
    fockMatrices[index_] = data_->hartreeFockData->F_OAO;
  }
  else {
    fockMatrices[index_] = data_->hartreeFockData->F;
  }

  evalauteError(index_);

  updateBMatrix();

  index_ = (index_ + 1) % subspaceSize_;
}

void SimultaneousDiis::updateBMatrix() {
  int activeSize = iterationNo_ > subspaceSize_ ? subspaceSize_ : iterationNo_;

  // Bii element
  B(lastAdded_ + 1, lastAdded_ + 1) = evalauteBMatrixElement(lastAdded_, lastAdded_);

  // Bij elements
  for (int i = 1; i < activeSize + 1; i++) {
    if (i == lastAdded_ + 1) {
      continue;
    }
    B(lastAdded_ + 1, i) = evalauteBMatrixElement(lastAdded_, i - 1);
    B(i, lastAdded_ + 1) = B(lastAdded_ + 1, i);
  }
}

auto SimultaneousDiis::getMixedFockMatrix() -> std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> {
  if (iterationNo_ > subspaceSize_) {
    iterationNo_ = subspaceSize_;
  }

  // If we have only one Simultaneous matrix
  if (iterationNo_ < 2) {
    return fockMatrices[0];
  }

  C.head(iterationNo_ + 1) =
      B.block(0, 0, iterationNo_ + 1, iterationNo_ + 1).colPivHouseholderQr().solve(rhs.head(iterationNo_ + 1));

  return calculateLinearCombination();
}

Eigen::MatrixXd SimultaneousDiis::calculateErrorMatrixOrthonormal(const Eigen::MatrixXd& D, const Eigen::MatrixXd& F) {
  Eigen::MatrixXd DF = D * F;
  return DF - DF.transpose();
}

Eigen::MatrixXd SimultaneousDiis::calculateErrorMatrix(const Eigen::MatrixXd& D, const Eigen::MatrixXd& F,
                                                       const Eigen::MatrixXd& S) {
  Eigen::MatrixXd SDF = S * D * F;
  return SDF - SDF.transpose();
}

auto SimultaneousDiis::evalauteError(int idx) -> void {
  for (const auto& elem : *molecule_) {
    const auto type = elem.first;
    if (elem.second.isRestricted) {
      if (useOrthonormalBasis) {
        errorMatrices[idx][type].restrictedMatrix() = calculateErrorMatrixOrthonormal(
            data_->D_OAO[type].restrictedMatrix(), data_->hartreeFockData->F_OAO[type].restrictedMatrix());
      }
      else {
        errorMatrices[idx][type].restrictedMatrix() = calculateErrorMatrix(
            data_->D_OAO[type].restrictedMatrix(), data_->hartreeFockData->F_OAO[type].restrictedMatrix(), data_->S[type]);
      }
    }
    else {
      if (useOrthonormalBasis) {
        errorMatrices[idx][type].alphaMatrix() = calculateErrorMatrixOrthonormal(
            data_->D_OAO[type].alphaMatrix(), data_->hartreeFockData->F_OAO[type].alphaMatrix());
        if (elem.second.msVector[1] > 0) {
          errorMatrices[idx][type].betaMatrix() = calculateErrorMatrixOrthonormal(
              data_->D_OAO[type].betaMatrix(), data_->hartreeFockData->F_OAO[type].betaMatrix());
        }
      }
      else {
        errorMatrices[idx][type].alphaMatrix() = calculateErrorMatrix(
            data_->D_OAO[type].alphaMatrix(), data_->hartreeFockData->F_OAO[type].alphaMatrix(), data_->S[type]);
        if (elem.second.msVector[1] > 0) {
          errorMatrices[idx][type].betaMatrix() = calculateErrorMatrix(
              data_->D_OAO[type].betaMatrix(), data_->hartreeFockData->F_OAO[type].betaMatrix(), data_->S[type]);
        }
      }
    }
  }
}

auto SimultaneousDiis::evalauteBMatrixElement(int i, int j) -> double {
  double ret = 0;

  for (const auto& elem : *molecule_) {
    const auto type = elem.first;
    if (elem.second.isRestricted) {
      ret += errorMatrices[i][type].restrictedMatrix().cwiseProduct(errorMatrices[j][type].restrictedMatrix()).sum();
    }
    else {
      ret += errorMatrices[i][type].alphaMatrix().cwiseProduct(errorMatrices[j][type].alphaMatrix()).sum();
      if (elem.second.msVector[1] > 0) {
        ret += errorMatrices[i][type].betaMatrix().cwiseProduct(errorMatrices[j][type].betaMatrix()).sum();
      }
    }
  }

  return ret;
}

Utils::SpinAdaptedMatrix SimultaneousDiis::calculateLinearCombination(const Utils::ElementType type) {
  if (!molecule_->at(type).isRestricted) {
    Eigen::MatrixXd F_alpha;
    F_alpha.resizeLike(fockMatrices[0][type].alphaMatrix());
    F_alpha.setZero();
    Eigen::MatrixXd F_beta;
    F_beta.resizeLike(fockMatrices[0][type].betaMatrix());
    F_beta.setZero();
    for (int i = 0; i < iterationNo_; i++) {
      F_alpha += C(i + 1) * fockMatrices[i][type].alphaMatrix();
    }
    if (molecule_->at(type).msVector[1] > 0) {
      for (int i = 0; i < iterationNo_; i++) {
        F_beta += C(i + 1) * fockMatrices[i][type].betaMatrix();
      }
    }
    return Utils::SpinAdaptedMatrix::createUnrestricted(std::move(F_alpha), std::move(F_beta));
  }

  Eigen::MatrixXd F_restricted;
  F_restricted.resizeLike(fockMatrices[0][type].restrictedMatrix());
  F_restricted.setZero();
  for (int i = 0; i < iterationNo_; i++) {
    F_restricted += C(i + 1) * fockMatrices[i][type].restrictedMatrix();
  }
  return Utils::SpinAdaptedMatrix::createRestricted(std::move(F_restricted));
}

auto SimultaneousDiis::calculateLinearCombination() -> std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> {
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> ret;

  for (const auto& elem : *molecule_) {
    ret[elem.first] = calculateLinearCombination(elem.first);
  }

  return ret;
}

auto SimultaneousDiis::setMixedFockMatrix() -> void {
  if (useOrthonormalBasis) {
    data_->hartreeFockData->F_OAO = getMixedFockMatrix();
  }
  else {
    data_->hartreeFockData->F = getMixedFockMatrix();
  }
}

} // namespace Kiwi
} // namespace Scine
