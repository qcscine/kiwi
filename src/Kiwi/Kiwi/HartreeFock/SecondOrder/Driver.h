/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_DRIVER_H
#define KIWI_DRIVER_H

#include <Kiwi/HartreeFock/SecondOrder/ArhBuilder.h>
#include <Kiwi/HartreeFock/SecondOrder/ModifiedDavidson.h>
#include <Kiwi/HartreeFock/SecondOrder/ModifiedDavidsonNewtonRhapson.h>
#include <Kiwi/KiwiOpt/PreconditionedDirectConjugateGradient.h>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {
namespace Arh {

class Driver {
  std::shared_ptr<SigmaVectorEvaluator> arhSigmaVectorEvaluator;
  // Shift
  double mu;
  // Solution to the (shifted) Newton equations: a^-1 * solution
  Eigen::VectorXd solution;

  bool success = false;

 public:
  bool wasSuccessful() const {
    return success;
  }

 private:
  int iterations;
  double error;

 public:
  double getError() const {
    return error;
  }

 private:
  double thresh;

 public:
  double getMu() const {
    return mu;
  }

  int getIterations() const {
    return iterations;
  }
  const Eigen::VectorXd& getSolution() const {
    return solution;
  }

 private:
  double trustRadius = 0.5;

 public:
  void setTrustRadius(double trustRadius) {
    Driver::trustRadius = trustRadius;
  }

  Driver(std::shared_ptr<SigmaVectorEvaluator> sigmaVectorEvaluator, double threshold)
    : arhSigmaVectorEvaluator(std::move(sigmaVectorEvaluator)), thresh(threshold) {
  }

  auto evaluateTRAH(int maxIt = 50) -> void {
    Arh::ModifiedDavidson davidson(arhSigmaVectorEvaluator);
    davidson.setTrustRadius(trustRadius);
    davidson.setMaxIterations(maxIt);
    davidson.setThresh(thresh);

    davidson.compute(false);
    solution = davidson.getEigvec();
    iterations = davidson.getIterations();
    error = davidson.getError();
    mu = davidson.getEigval();
  }

  auto evaluateNR(int maxIt = 50) -> void {
    Arh::ModifiedDavidsonNewtonRhapson davidsonNewtonRhapson(arhSigmaVectorEvaluator);

    davidsonNewtonRhapson.setMaxIterations(maxIt);
    davidsonNewtonRhapson.setThresh(thresh);
    davidsonNewtonRhapson.compute(false);

    success = davidsonNewtonRhapson.isSuccessful();
    solution = davidsonNewtonRhapson.getEigvec();
    iterations = davidsonNewtonRhapson.getIterations();
    error = davidsonNewtonRhapson.getError();
    mu = 0;
  }

  auto evaluatePCG(int maxIt = 50) -> void {
    Eigen::VectorXd guess = arhSigmaVectorEvaluator->getGradient();
    Eigen::VectorXd image = -arhSigmaVectorEvaluator->getGradient();

    Optimization::PreconditionedDirectConjugateGradient<double> pcg(std::move(guess));

    pcg.updateImage(std::move(image));

    pcg.setMaxIterations(maxIt);
    pcg.setThresh(thresh);
    arhSigmaVectorEvaluator->evaluatePreconditioner();
    pcg.optimize(*arhSigmaVectorEvaluator, arhSigmaVectorEvaluator->getPreconditioner());

    success = pcg.getError() < thresh;
    solution = pcg.getSolution();
    iterations = pcg.getIterations();
    error = pcg.getError();
    mu = 0;
  }
};

} // namespace Arh
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_DRIVER_H
