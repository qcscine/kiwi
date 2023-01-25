/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_MODIFIEDDAVIDSON_H
#define KIWI_MODIFIEDDAVIDSON_H

#include <Kiwi/HartreeFock/SecondOrder/SigmaVectorEvaluator.h>
#include <Kiwi/KiwiOpt/Bisection.h>
#include <Eigen/Eigenvalues>
#include <chrono>
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
class ModifiedDavidson {
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

  double trustRadius;
  double alpha;

  bool precondition = true;

  Eigen::DiagonalMatrix<double, -1, -1> diagonal;
  double eigval;
  Vector subspaceEigvec;
  Vector diminishedRitzVector;
  Vector residual;
  Vector newDirection;

  bool zeroLevelShift_ = false;

 public:
  ModifiedDavidson(std::shared_ptr<SigmaVectorEvaluator> evaluator)
    : sigmaVectorEvaluator(std::move(evaluator)),
      maxIterations(sigmaVectorEvaluator->maxRedSpaceDim - 2),
      fullDim(sigmaVectorEvaluator->fullDim) {
  }
  void setPrecondition(bool precond) {
    precondition = precond;
  }
  double getAlpha() const {
    return alpha;
  }
  bool isSuccessful() const {
    return successful;
  }
  void setThresh(double threshold) {
    thresh = threshold;
  }
  bool zeroLevelShift() const {
    return zeroLevelShift_;
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

  void setTrustRadius(double trRad) {
    trustRadius = trRad;
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
      std::cout << "--------------------------------------------\n";
      std::cout << " ModifiedDavidson\n";
      std::cout << "--------------------------------------------\n";
      std::cout << " Iteration   Error          Alpha           \n";
      std::cout << "--------------------------------------------\n";
    }

    iterations = 0;

    residual.resize(fullDim - 1);
    diminishedRitzVector.resize(fullDim);

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> precond;

    precondition = true;

    if (precondition) {
      sigmaVectorEvaluator->evaluatePreconditioner();
      precond.setIdentity(sigmaVectorEvaluator->fullDim);
    }

    double bisectionThresh = 1.;

    iterations += sigmaVectorEvaluator->trialVecDim;

    while (true) {
      subspaceDim = sigmaVectorEvaluator->redSpaceDim;
      double alphaMin = 1;
      double alphaMax = 1000;

      bool doBisection = false;

      // Get alpha:
      {
        AlphaOptimization alphaOptimization(sigmaVectorEvaluator->redSpaceTRAH, trustRadius);

        double normMinusTrustRadiusMin = alphaOptimization.evaluate(alphaMin);

        doBisection = normMinusTrustRadiusMin >= 0;
        // isAlphaZero = (normMinusTrustRadiusMin < 0);
      }

      if (doBisection) {
        // if (!isAlphaZero) {
        zeroLevelShift_ = false;
        AlphaOptimization alphaOptimization(sigmaVectorEvaluator->redSpaceTRAH, trustRadius);

        // TODO: check if [alphaMin,alphaMax] is a valid interval. If not do bracketAlpha
        bracketAlpha(alphaMin, alphaMax);

        Optimization::Bisection bisection;
        bisection.setMaxIterations(100);
        bisection.setLowerBound(alphaMin);
        bisection.setUpperBound(alphaMax);

        bisectionThresh = error * 0.1;

        if (bisectionThresh > 1e-4) {
          bisectionThresh = 1e-4;
        }

        bisection.setTolerance(bisectionThresh);

        alpha = bisection.template compute(alphaOptimization, false);

        if (!bisection.hasBoundaryError()) {
          alphaOptimization.evaluate(alpha);

          eigval = alphaOptimization.getEigenvalue();
          subspaceEigvec = alphaOptimization.getEigenvector();
          subspaceEigvec.segment(1, subspaceDim - 1) /= alpha;
        }
        else {
          throw std::runtime_error("Boundary error in bisection!");
        }
      }
      else {
        alpha = alphaMin;

        AlphaOptimization alphaOptimization(sigmaVectorEvaluator->redSpaceTRAH, trustRadius);
        alphaOptimization.evaluate(alpha);

        eigval = alphaOptimization.getEigenvalue();
        subspaceEigvec = alphaOptimization.getEigenvector();
        subspaceEigvec.segment(1, subspaceDim - 1) /= alpha;
      }

      diminishedRitzVector = sigmaVectorEvaluator->trialVectors * subspaceEigvec.segment(1, subspaceDim - 1);

      residual =
          sigmaVectorEvaluator->getGradient() +
          (sigmaVectorEvaluator->sigmaVectors * subspaceEigvec.segment(1, subspaceDim - 1) - eigval * diminishedRitzVector);

      iterations++;

      error = residual.norm();

      if (verbose) {
        std::cout << std::right << std::setw(10) << iterations << std::string(3, ' ') << std::setw(15)
                  << std::scientific << std::setprecision(5) << error << std::setw(15) << std::scientific
                  << std::setprecision(5) << alpha << std::endl;
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
        precond = sigmaVectorEvaluator->getPreconditioner(true, eigval);
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
  }

 private:
  void bracketAlpha(const double& alphaMin, double& alphaMax) {
    double upperBound = 1e3;

    double factor = 10;
    alphaMax = alphaMin + factor;

    int Nmax = 500;

    AlphaOptimization alphaOptimization(sigmaVectorEvaluator->redSpaceTRAH, trustRadius);

    auto fMin = alphaOptimization.evaluate(alphaMin);
    auto fMax = alphaOptimization.evaluate(alphaMax);
    if (fMin * fMax < 0) {
      return;
    }

    for (int i = 1; i < Nmax; ++i) {
      alphaMax += factor;
      fMax = alphaOptimization.evaluate(alphaMax);
      if (fMin * fMax < 0) {
        return;
      }
      if (alphaMax > upperBound) {
        if (fMin < 0 && fMax < 0) {
          zeroLevelShift_ = true;
          return;
        }
        upperBound *= 10;
      }
      factor *= 10;
    }
  }

  class AlphaOptimization {
    const Matrix& reducedSpaceAugmentedHessian;
    Vector eigvec;
    Vector eigvals;
    Eigen::Index min;
    int dim;
    Eigen::EigenSolver<Matrix> subspaceDiagonalizer;
    double trustRadius;
    Matrix tmpHess;
    double alpha;
    double norm;

   public:
    AlphaOptimization(const Matrix& augHess, double tr)
      : reducedSpaceAugmentedHessian(augHess), dim(reducedSpaceAugmentedHessian.rows()), trustRadius(tr) {
    }

    auto evaluate(double scaling) -> double {
      alpha = scaling;
      tmpHess = reducedSpaceAugmentedHessian;

      tmpHess(0, 1) *= alpha;
      tmpHess(1, 0) *= alpha;

      subspaceDiagonalizer.compute(tmpHess);

      eigvals = subspaceDiagonalizer.eigenvalues().real();

      eigvals.minCoeff(&min);
      eigvec = subspaceDiagonalizer.eigenvectors().real().col(min);
      eigvec = eigvec / eigvec(0);

      norm = std::sqrt(eigvec.block(1, 0, dim - 1, 1).squaredNorm() / (alpha * alpha));
      return norm - trustRadius;
    }

    auto getEigenvector() -> Vector {
      return eigvec;
    }

    auto getEigenvalue() -> double {
      return eigvals(min);
    }

    auto getNorm() -> double {
      return norm;
    }
  };
};

} // namespace Arh
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_MODIFIEDDAVIDSON_H
