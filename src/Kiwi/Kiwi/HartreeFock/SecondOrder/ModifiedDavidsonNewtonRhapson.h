/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_MODIFIEDDAVIDSONNEWTONRHAPSON_H
#define KIWI_MODIFIEDDAVIDSONNEWTONRHAPSON_H

#include <Kiwi/HartreeFock/SecondOrder/SigmaVectorEvaluator.h>
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
class ModifiedDavidsonNewtonRhapson {
  using Matrix = Eigen::MatrixXd;
  using Vector = Eigen::VectorXd;

 private:
  double thresh = 1e-10;

  std::shared_ptr<SigmaVectorEvaluator> sigmaVectorEvaluator;
  int maxIterations;
  int subspaceDim;
  int iterations;
  double error;
  bool successful = false;

  bool precondition = true;

  Eigen::DiagonalMatrix<double, -1, -1> diagonal;
  double eigval;
  Vector subspaceEigvec;
  Vector diminishedRitzVector;
  Vector residual;
  Vector newDirection;

 public:
  ModifiedDavidsonNewtonRhapson(std::shared_ptr<SigmaVectorEvaluator> evaluator)
    : sigmaVectorEvaluator(std::move(evaluator)), maxIterations(sigmaVectorEvaluator->maxRedSpaceDim - 2) {
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
    return diminishedRitzVector;
  }

  double getError() const {
    return error;
  }

  void setMaxIterations(int maxIt) {
    maxIterations = maxIt;
  }

  auto restartFromLastSolution(bool verbose = true) -> void {
    sigmaVectorEvaluator->addTrialVectorAndRestart(diminishedRitzVector);
    compute(verbose);
  }

  /**
   * Main routine.
   */
  auto compute(bool verbose = true) -> void {
    if (verbose) {
      std::cout << "----------------------------\n";
      std::cout << "Davidson-Newton-Rhapson\n";
      std::cout << "-----------------------------\n";
      std::cout << " Iteration   Error          \n";
      std::cout << "-----------------------------\n";
    }

    iterations = 0;

    residual.resize(sigmaVectorEvaluator->fullDim);
    diminishedRitzVector.resize(sigmaVectorEvaluator->fullDim);

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> precond;

    precondition = true;

    if (precondition) {
      sigmaVectorEvaluator->evaluatePreconditioner();
      precond = sigmaVectorEvaluator->getPreconditioner(false);
    }

    while (true) {
      subspaceDim = sigmaVectorEvaluator->redSpaceDim;

      AlphaZeroOptimizer alphaZeroOptimizer(sigmaVectorEvaluator->redSpaceTRAH);
      eigval = 0;
      alphaZeroOptimizer.evaluate();
      subspaceEigvec.resize(subspaceDim);
      subspaceEigvec.col(0)(0) = 0;
      subspaceEigvec.segment(1, subspaceDim - 1) = alphaZeroOptimizer.getSolution();

      diminishedRitzVector = sigmaVectorEvaluator->trialVectors * subspaceEigvec.segment(1, subspaceDim - 1);

      residual = sigmaVectorEvaluator->getGradient() +
                 sigmaVectorEvaluator->sigmaVectors * subspaceEigvec.segment(1, subspaceDim - 1);

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
        std::cout << "ModifiedDavidsonNewtonRhapson failed after " << iterations << " iterations." << std::endl;
        break;
      }

      if (precondition) {
        newDirection = precond * residual;
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
  }

 private:
  class AlphaZeroOptimizer {
    const Matrix& reducedSpaceAugmentedHessian;
    Vector solution;
    int dim;

   public:
    AlphaZeroOptimizer(const Matrix& augHess)
      : reducedSpaceAugmentedHessian(augHess), dim(reducedSpaceAugmentedHessian.rows()) {
    }

    auto evaluate() -> void {
      solution = reducedSpaceAugmentedHessian.block(1, 1, dim - 1, dim - 1)
                     .fullPivLu()
                     .solve(-reducedSpaceAugmentedHessian.block(1, 0, dim - 1, 1));
    }

    auto getSolution() -> Vector {
      return solution;
    }
  };
};

} // namespace Arh
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_MODIFIEDDAVIDSONNEWTONRHAPSON_H
