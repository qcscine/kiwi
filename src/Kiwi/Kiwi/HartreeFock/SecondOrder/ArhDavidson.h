/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_ARHDAVIDSON_H
#define KIWI_ARHDAVIDSON_H

#include <Kiwi/HartreeFock/SecondOrder/SigmaVectorEvaluator.h>
#include <Eigen/Eigenvalues>
#include <cstddef>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {
namespace Arh {

/**
 * https://epubs.siam.org/doi/abs/10.1137/0915004
 * SIAM J. Sci. Comput., 15(1), 62–76. (15 pages)
 */
class Davidson {
  using Matrix = Eigen::MatrixXd;
  using Vector = Eigen::VectorXd;

 private:
  double thresh = 1e-10;

  std::shared_ptr<SigmaVectorEvaluator> sigmaVectorEvaluator;
  int maxIterations;
  int fullDim;
  int subspaceDim;
  int iterations;
  double error;
  bool successful = false;

  bool precondition = true;

  Eigen::DiagonalMatrix<double, -1, -1> diagonal;
  double eigval;
  Vector subspaceEigvec;
  Vector ritzVector;
  Vector residual;
  Vector newDirection;

 public:
  Davidson(std::shared_ptr<SigmaVectorEvaluator> evaluator)
    : sigmaVectorEvaluator(std::move(evaluator)),
      maxIterations(sigmaVectorEvaluator->maxRedSpaceDim - 2),
      fullDim(sigmaVectorEvaluator->fullDim) {
  }
  void setPrecondition(bool precond) {
    precondition = precond;
  }
  bool isSuccessful() const {
    return successful;
  }
  void setThresh(double threshold) {
    thresh = threshold;
  }

  int getIterations() const {
    return iterations;
  }

  const double& getEigval() const {
    return eigval;
  }

  const Vector& getEigvec() const {
    return ritzVector;
  }

  double getError() const {
    return error;
  }

  void setMaxIterations(int maxIt) {
    maxIterations = maxIt;
  }

  auto restartFromLastSolution(bool verbose = true) -> void {
    sigmaVectorEvaluator->addTrialVectorAndRestart(ritzVector);

    compute(verbose);
  }

  /**
   * Main routine.
   */
  auto compute(bool verbose = true) -> void {
    if (verbose) {
      std::cout << "------------------------------\n";
      std::cout << " Davidson\n";
      std::cout << "------------------------------\n";
      std::cout << " Iteration   Error          \n";
      std::cout << "------------------------------\n";
    }

    iterations = 0;

    residual.resize(fullDim - 1);
    ritzVector.resize(fullDim);

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> precond;

    precondition = true;

    if (precondition) {
      sigmaVectorEvaluator->evaluatePreconditioner();
      precond.setIdentity(sigmaVectorEvaluator->fullDim);
    }

    iterations += sigmaVectorEvaluator->trialVecDim;

    while (true) {
      subspaceDim = sigmaVectorEvaluator->redSpaceDim;

      ReducedHessianOptimizer reducedHessianOptimizer(sigmaVectorEvaluator->redSpaceTRAH);
      reducedHessianOptimizer.evaluate();

      eigval = reducedHessianOptimizer.getEigenvalue();
      subspaceEigvec = reducedHessianOptimizer.getEigenvector();

      ritzVector = sigmaVectorEvaluator->trialVectors * subspaceEigvec;

      residual = (sigmaVectorEvaluator->sigmaVectors * subspaceEigvec - eigval * ritzVector);

      iterations++;

      error = residual.norm();

      if (verbose) {
        std::cout << std::right << std::setw(10) << iterations << std::string(3, ' ') << std::setw(15)
                  << std::scientific << std::setprecision(5) << error << std::endl;
      }

      if (error < thresh) {
        successful = true;
        break;
      }
      if (iterations == maxIterations) {
        successful = false;
        std::cout << "ModifiedDavidson failed after " << iterations << " iterations." << std::endl;
        break;
      }

      if (precondition) {
        // precond = sigmaVectorEvaluator->getPreconditioner(true, eigval);
        precond = sigmaVectorEvaluator->getPreconditioner(false);
      }

      if (precondition) {
        newDirection = residual;
        newDirection = (precond * newDirection).eval();
      }
      else {
        newDirection = residual;
      }

      if (sigmaVectorEvaluator->redSpaceDim == sigmaVectorEvaluator->maxRedSpaceDim) {
        sigmaVectorEvaluator->addTrialVectorAndRestart(newDirection);
      }
      else {
        sigmaVectorEvaluator->addTrialVector(newDirection);
      }
    }

    if (verbose) {
      std::cout << "----------------------------\n";
      std::cout << "Lowest eigenvalue\n";
      std::cout << std::fixed << std::setprecision(10) << eigval << std::endl;
      std::cout << "----------------------------\n";
    }
  }

 private:
  class ReducedHessianOptimizer {
    const Matrix& reducedSpaceAugmentedHessian;
    Vector eigvec;
    Vector eigvals;
    Eigen::Index min;
    int dim;
    Eigen::SelfAdjointEigenSolver<Matrix> subspaceDiagonalizer;
    Matrix tmpHess;

   public:
    ReducedHessianOptimizer(const Matrix& augHess)
      : reducedSpaceAugmentedHessian(augHess), dim(reducedSpaceAugmentedHessian.rows()) {
    }

    auto evaluate() -> void {
      tmpHess = reducedSpaceAugmentedHessian.block(1, 1, dim - 1, dim - 1);

      subspaceDiagonalizer.compute(tmpHess);

      eigvals = subspaceDiagonalizer.eigenvalues().real();

      eigvals.minCoeff(&min);
      eigvec = subspaceDiagonalizer.eigenvectors().real().col(min);
    }

    auto getEigenvector() -> Vector {
      return eigvec;
    }

    auto getEigenvalue() -> double {
      return eigvals(min);
    }
  };
};

} // namespace Arh
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_ARHDAVIDSON_H
