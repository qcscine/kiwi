/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_DAVIDSON_H
#define KIWI_DAVIDSON_H

#include <Eigen/Eigenvalues>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {
namespace Optimization {

/**
 * @class Davidson
 * @brief Implementation of the Davidson diagonalization algorithm after the following referece:
 * https://epubs.siam.org/doi/abs/10.1137/0915004
 * SIAM J. Sci. Comput., 15(1), 62–76. (15 pages)
 */
class Davidson {
  using Matrix = Eigen::MatrixXd;
  using Vector = Eigen::VectorXd;

 private:
  int maxIterations = 100;
  int iterations;
  double thresh = 1e-10;
  double error;
  bool successful = false;
  bool useDiagonalPreconditioner = false;

  Matrix projector;
  Matrix sigmaVectors;
  Matrix interactionMatrix;
  Eigen::DiagonalMatrix<double, -1, -1> preConditioner;
  Eigen::DiagonalMatrix<double, -1, -1> diagonal;
  Vector eigvals;
  Matrix subspaceEigvecs;
  Matrix RitzVectors;
  Matrix residuals;
  Matrix newDirection;

  int numEigpairs;
  int maxSubspaceDim;
  int fullDim;
  int subspaceDim;

 public:
  int getSubspaceDim() const {
    return subspaceDim;
  }

 private:
  bool guessProvided = false;

 public:
  Davidson(int numberOfEigenPairs, int maxSubspaceDimension, int fullDimension)
    : numEigpairs(numberOfEigenPairs), maxSubspaceDim(maxSubspaceDimension), fullDim(fullDimension) {
  }

  int getIterations() const {
    return iterations;
  }

  const Vector& getEigvals() const {
    return eigvals;
  }

  const Matrix& getEigvecs() const {
    return RitzVectors;
  }

  const Matrix& getProjector() const {
    return projector;
  }

  double getError() const {
    return error;
  }

  void setMaxIterations(int maxIt) {
    maxIterations = maxIt;
  }

  void setDiagonalPreconditioner(Eigen::DiagonalMatrix<double, -1, -1>&& preCond) {
    diagonal = preCond;
    useDiagonalPreconditioner = true;
  }

  void setGuess(Eigen::MatrixXd&& guess) {
    projector = guess;
    guessProvided = true;
    if (projector.cols() < numEigpairs) {
      throw std::runtime_error("Provide at least as many trial vectors, as eigenpairs are requested!");
    }
  }

  void setThresh(double thresh) {
    Davidson::thresh = thresh;
  }

  /**
   * @brief Main routine
   * @tparam SigmaVectorEvaluator Class that contains the sigma-vector evaluator.
   *         --> Must have the following function implemented: `SigmaVectorEvaluator.evaluate(Eigen::MatrixXd) ->
   * Eigen::MatrixXd`
   * @param evaluator
   * @param verbose
   */
  template<class SigmaVectorEvaluator>
  auto compute(SigmaVectorEvaluator& evaluator, bool verbose = true) -> void {
    if (verbose) {
      std::cout << "----------------------------\n";
      std::cout << " Davidson\n";
      std::cout << "----------------------------\n";
      std::cout << " Iteration   Error          \n";
      std::cout << "----------------------------\n";
    }

    iterations = 0;

    residuals.resize(fullDim, numEigpairs);
    RitzVectors.resize(fullDim, numEigpairs);
    newDirection.resize(fullDim, numEigpairs);

    // int sizeOfGuess=0;

    if (!guessProvided) {
      projector = Eigen::MatrixXd::Identity(fullDim, numEigpairs);
      subspaceDim = numEigpairs;
    }
    else {
      subspaceDim = projector.cols();
    }

    Eigen::SelfAdjointEigenSolver<Matrix> subspaceDiagonalizer;

    Eigen::DiagonalMatrix<double, -1, -1> Id;
    Id.setIdentity(fullDim);

    while (true) {
      sigmaVectors = evaluator.evaluate(projector);

      interactionMatrix = projector.transpose() * sigmaVectors;

      subspaceDiagonalizer.compute(interactionMatrix);

      eigvals = subspaceDiagonalizer.eigenvalues().block(0, 0, numEigpairs, 1);
      subspaceEigvecs = subspaceDiagonalizer.eigenvectors().block(0, 0, subspaceDim, numEigpairs);

      for (auto i = 0; i < numEigpairs; ++i) {
        RitzVectors.col(i) = projector * subspaceEigvecs.col(i);
        residuals.col(i) = eigvals(i) * RitzVectors.col(i) - sigmaVectors * subspaceEigvecs.col(i);
      }

      iterations++;

      error = residuals.norm();

      if (verbose) {
        std::cout << std::right << std::setw(10) << iterations << std::string(3, ' ') << std::setw(15)
                  << std::scientific << std::setprecision(5) << error << std::endl;
      }

      if (error < thresh) {
        break;
      }
      if (iterations == maxIterations) {
        std::cout << "Davidson diagonalization failed after " << iterations << " iterations." << std::endl;
        break;
      }

      /* Insert preconditioner here */
      //
      if (useDiagonalPreconditioner) {
        for (auto i = 0; i < numEigpairs; ++i) {
          preConditioner = (eigvals(i) * Id.diagonal() - diagonal.diagonal()).cwiseInverse().asDiagonal();
          newDirection.col(i) = preConditioner * residuals.col(i);
        }
      }
      else {
        newDirection = residuals;
      }

      initiateModifiedGramSchmidt();
    }

    if (verbose) {
      std::cout << "----------------------------\n";
    }
  }

  const Matrix& getSubspaceEigvecs() const {
    return subspaceEigvecs;
  }

  auto initiateModifiedGramSchmidt() -> void {
    if (subspaceDim <= maxSubspaceDim - numEigpairs) {
      int oldSubspaceDim = subspaceDim;

      subspaceDim = subspaceDim + numEigpairs;

      projector.conservativeResize(fullDim, subspaceDim);

      for (auto i = 0; i < numEigpairs; ++i) {
        projector.col(oldSubspaceDim + i) = newDirection.col(i);
      }
    }
    else {
      subspaceDim = 2 * numEigpairs;

      projector.resize(fullDim, subspaceDim);

      for (auto i = 0; i < numEigpairs; ++i) {
        projector.col(i) = RitzVectors.col(i);
      }
      for (auto i = numEigpairs; i < subspaceDim; ++i) {
        projector.col(i) = newDirection.col(i - numEigpairs);
      }
    }

    projector = modifiedGramSchmidt(projector);
  }

  /**
   * @brief Modified Gram-Schmidt method that can handle rectangular matrices.
   * @param vectorSet
   * @return
   */
  static auto modifiedGramSchmidt(const Eigen::MatrixXd& vectorSet) -> Eigen::MatrixXd {
    // Calculate modified Gram-Schmidt QR decomposition
    Matrix Q(vectorSet.rows(), vectorSet.cols());

    for (int j = 0; j < vectorSet.cols(); ++j) {
      Q.col(j) = vectorSet.col(j);
      for (int i = 0; i < j; ++i) {
        Q.col(j) -= Q.col(i) * (Q.col(i).adjoint() * Q.col(j));
      }
      Q.col(j).normalize();
    }

    return Q;
  }
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_DAVIDSON_H
