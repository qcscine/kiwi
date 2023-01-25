/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_SIGMAVECTOREVALUATOR_H
#define KIWI_SIGMAVECTOREVALUATOR_H

#include <Kiwi/HartreeFock/SecondOrder/ArhInterface.h>
#include <Kiwi/KiwiOpt/Davidson.h>
#include <iomanip>
#include <iostream>
#include <utility>

namespace Scine {
namespace Kiwi {
namespace Arh {

class SigmaVectorEvaluator {
 private:
  std::shared_ptr<ArhInterface> ptrArhWrapper;
  const std::vector<int> dim;
  // Number of different types and spins
  const int numElements;

  std::vector<int> offset;

  Eigen::MatrixXd AV;

  Eigen::VectorXd gradient;

  double gradientNorm;

 public:
  double getGradientNorm() const {
    return gradientNorm;
  }

 public:
  int maxRedSpaceDim;
  int redSpaceDim;
  int trialVecDim;
  // Dimension of all bases combined
  int fullDim;
  // Maximum fullDimension. Required for reserving enough space.
  // Trial vectors b_i
  // Attention: the first trial vector (1, 0, 0, ... ) is implicit.
  Eigen::MatrixXd trialVectors;
  // Sigma vectors H b_i
  Eigen::MatrixXd sigmaVectors;
  // TRAH in reduces space
  // A(0,0) = 0
  // A(0,1) = ||grad||
  // A(i,j) = b_i sigma_j
  Eigen::MatrixXd redSpaceTRAH;

  // Preconditioner
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> precond;
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> eye;

  const Eigen::VectorXd& getGradient() const {
    return gradient;
  }

  /*
   * Recommended settings:
   *  -- In general: 3 trial vectors without noise
   *  -- For difficult unrestricted/multireference cases try
   *        1st 2 trial vectors, second with noise
   *        2nd 3 trial vectors, third with noise
   *        Homo-lumo (or more) mixing with combinations thereof
   */

  SigmaVectorEvaluator(std::shared_ptr<ArhInterface> arhWrapper, int maxDim, bool secondGuessVector = true,
                       bool secondGuessVectorNoise = false, bool thirdGuessVector = false, bool thirdGuessVectorNoise = false)
    : ptrArhWrapper(std::move(arhWrapper)),
      dim(ptrArhWrapper->getDim()),
      numElements(dim.size()),
      offset(ptrArhWrapper->getOffset()),
      maxRedSpaceDim(maxDim),
      redSpaceDim(2),
      trialVecDim(1) {
    fullDim = std::accumulate(dim.begin(), dim.end(), 0);

    gradient.resize(fullDim);

    populateGradient();
    gradientNorm = gradient.lpNorm<2>();
    // 1st reserve enough space
    trialVectors.resize(fullDim, maxRedSpaceDim - 1);
    sigmaVectors.resize(fullDim, maxRedSpaceDim - 1);
    redSpaceTRAH.resize(maxRedSpaceDim, maxRedSpaceDim);
    // Now that we have reserved enough space, we reduce the matrices to their starting sizes.
    resize();

    addFirstTrialVector();

    if (secondGuessVector) {
      addSecondTrialVector(secondGuessVectorNoise);
    }

    if (thirdGuessVector) {
      addThirdTrialVector(thirdGuessVectorNoise);
    }
  }

  auto addTrialVector(const Eigen::VectorXd& newDirection) -> void {
    redSpaceDim += 1;
    trialVecDim += 1;

    int indexTrialBasis = redSpaceDim - 2;
    int indexRedBasis = redSpaceDim - 1;

    if (redSpaceDim > maxRedSpaceDim) {
      throw std::runtime_error("Attention: discrepancy between maximum subspace fullDimensions in TRAH procedure.");
    }

    resize();

    trialVectors.col(indexTrialBasis) = newDirection;

    for (int i = 0; i < indexTrialBasis; ++i) {
      trialVectors.col(indexTrialBasis) -=
          trialVectors.col(i) * (trialVectors.col(i).adjoint() * trialVectors.col(indexTrialBasis));
    }
    trialVectors.col(indexTrialBasis).normalize();

    construct(indexTrialBasis, indexRedBasis);
  }

  auto addRandomTrialVector() -> void {
    Eigen::VectorXd randomDirection = Eigen::VectorXd::Random(fullDim);
    addTrialVector(randomDirection);
  }

  auto addTrialVectorAndRestart(const Eigen::VectorXd& newDirection) -> void {
    redSpaceDim = 4;
    trialVecDim = 3;

    int indexTrialBasis = redSpaceDim - 2;
    int indexRedBasis = redSpaceDim - 1;

    resize();

    trialVectors.col(indexTrialBasis) = newDirection;

    trialVectors = Optimization::Davidson::modifiedGramSchmidt(trialVectors);

    construct(indexTrialBasis, indexRedBasis);
  }

  auto setResult(const Eigen::VectorXd& fullX) -> void {
    ptrArhWrapper->updateX(fullX);
  }

  auto setLevelShift(double levelShift) -> void {
    ptrArhWrapper->setLevelShift(levelShift);
  }

  auto removeLevelShift() -> void {
    ptrArhWrapper->removeLevelShift();
  }

  auto evaluate(const Eigen::VectorXd& X) -> Eigen::VectorXd {
    Eigen::VectorXd HX(fullDim);

    ptrArhWrapper->updateX(X);
    ptrArhWrapper->makeContributionsDependentOnX();

    for (int i = 0; i < numElements; ++i) {
      HX.segment(offset[i], dim[i]) = ptrArhWrapper->evaluate(i);
    }

    return HX;
  }

  auto evaluate(const Eigen::MatrixXd& X) -> Eigen::MatrixXd {
    Eigen::MatrixXd HX(fullDim, X.cols());

    for (int col = 0; col < X.cols(); ++col) {
      ptrArhWrapper->updateX(X.col(col));
      ptrArhWrapper->makeContributionsDependentOnX();

      for (int i = 0; i < numElements; ++i) {
        HX.col(col).segment(offset[i], dim[i]) = ptrArhWrapper->evaluate(i);
      }
    }
    return HX;
  }

  auto evaluatePreconditioner() -> void {
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> ret(fullDim);
    precond.resize(fullDim);
    eye.resize(fullDim);
    eye.setIdentity();

    for (int i = 0; i < numElements; ++i) {
      // F_ii - F_aa
      precond.diagonal().segment(offset[i], dim[i]) = ptrArhWrapper->getRhDiagonal(i).diagonal();
    }
  }

  auto getRawPreconditioner() -> Eigen::DiagonalMatrix<double, Eigen::Dynamic> {
    return precond;
  }

  auto getPreconditioner(bool useLvelshift = false, double mu = 0.) -> Eigen::DiagonalMatrix<double, Eigen::Dynamic> {
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> ret(fullDim);

    ret = precond;

    if (useLvelshift) {
      ret.diagonal() -= mu * eye.diagonal();
    }

    ret = ret.inverse();

    return ret;
  }

 private:
  auto populateGradient() -> void {
    for (int i = 0; i < numElements; ++i) {
      gradient.segment(offset[i], dim[i]) = ptrArhWrapper->getGradient(i);
    }
  }

  auto populateSigmaVector(int at) -> void {
    ptrArhWrapper->updateX(trialVectors.col(at));
    ptrArhWrapper->makeContributionsDependentOnX();

    for (int i = 0; i < numElements; ++i) {
      sigmaVectors.col(at).segment(offset[i], dim[i]) = ptrArhWrapper->evaluate(i);
    }
  }

  auto resize() -> void {
    trialVectors.conservativeResize(fullDim, trialVecDim);
    sigmaVectors.conservativeResize(fullDim, trialVecDim);
    redSpaceTRAH.conservativeResize(redSpaceDim, redSpaceDim);
  }

  /**
   * This method generates b_1, which is B_2 in the full fullDimension.
   */
  auto addFirstTrialVector() -> void {
    trialVectors.col(0) = gradient; // / gradientNorm;

    trialVectors.col(0).normalize();
    // Here a copy is made. This has to be eliminated in the future.
    populateSigmaVector(0);

    redSpaceTRAH(0, 0) = 0;
    redSpaceTRAH(1, 0) = gradientNorm;
    redSpaceTRAH(0, 1) = gradientNorm;
    redSpaceTRAH(1, 1) = trialVectors.col(0).dot(sigmaVectors.col(0));
  }

  /**
   * This method generates b_1, which is B_2 in the full fullDimension.
   */
  auto addSecondTrialVector(bool addNoise = false) -> void {
    redSpaceDim += 1;
    trialVecDim += 1;

    int indexTrialBasis = redSpaceDim - 2;
    int indexRedBasis = redSpaceDim - 1;

    if (redSpaceDim > maxRedSpaceDim) {
      throw std::runtime_error("Attention: discrepancy between maximum subspace fullDimensions in TRAH procedure.");
    }

    resize();

    trialVectors.col(indexTrialBasis) = gradient;

    ptrArhWrapper->updateX(gradient);
    ptrArhWrapper->makeContributionsDependentOnX();

    for (int i = 0; i < numElements; ++i) {
      trialVectors.col(indexTrialBasis).segment(offset[i], dim[i]) += ptrArhWrapper->evaluate(i);
    }

    if (addNoise) {
      Eigen::VectorXd noise;
      noise.resizeLike(trialVectors.col(indexTrialBasis));
      noise.setRandom();
      noise *= 0.01;
      trialVectors.col(indexTrialBasis) += noise;
    }

    // Orthonormalize with Gram-Schmidt:
    for (int i = 0; i < indexTrialBasis; ++i) {
      trialVectors.col(indexTrialBasis) -=
          trialVectors.col(i) * (trialVectors.col(i).adjoint() * trialVectors.col(indexTrialBasis));
    }
    trialVectors.col(indexTrialBasis).normalize();

    construct(indexTrialBasis, indexRedBasis);
  }

  /**
   * This method generates b_3, which is B_4 in the full fullDimension.
   */
  auto addThirdTrialVector(bool addNoise = false) -> void {
    redSpaceDim += 1;
    trialVecDim += 1;

    int indexTrialBasis = redSpaceDim - 2;
    int indexRedBasis = redSpaceDim - 1;

    if (redSpaceDim > maxRedSpaceDim) {
      throw std::runtime_error("Attention: discrepancy between maximum subspace fullDimensions in TRAH procedure.");
    }

    resize();

    for (int i = 0; i < numElements; ++i) {
      trialVectors.col(indexTrialBasis).segment(offset[i], dim[i]) = ptrArhWrapper->getThirdTrialVector(i);
    }

    if (addNoise) {
      Eigen::VectorXd noise;
      noise.resizeLike(trialVectors.col(indexTrialBasis));
      noise.setRandom();
      noise *= 0.01;
      trialVectors.col(indexTrialBasis) += noise;
    }

    // Orthonormalize with Gram-Schmidt:
    for (int i = 0; i < indexTrialBasis; ++i) {
      trialVectors.col(indexTrialBasis) -=
          trialVectors.col(i) * (trialVectors.col(i).adjoint() * trialVectors.col(indexTrialBasis));
    }
    trialVectors.col(indexTrialBasis).normalize();

    construct(indexTrialBasis, indexRedBasis);
  }

  auto construct(int indexTrialBasis, int indexRedBasis) -> void {
    // Evaluate Sigma vector
    populateSigmaVector(indexTrialBasis);

    redSpaceTRAH(0, indexRedBasis) = 0;
    redSpaceTRAH(indexRedBasis, 0) = 0;

    for (int i = 1; i <= indexRedBasis; ++i) {
      redSpaceTRAH(i, indexRedBasis) = trialVectors.col(i - 1).dot(sigmaVectors.col(indexTrialBasis));
      if (i != indexRedBasis) {
        redSpaceTRAH(indexRedBasis, i) = trialVectors.col(indexTrialBasis).dot(sigmaVectors.col(i - 1));
      }
    }
  }
};

} // namespace Arh
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_SIGMAVECTOREVALUATOR_H
