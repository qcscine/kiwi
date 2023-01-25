/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_BISECTION_H
#define KIWI_BISECTION_H

#include "Kiwi/KiwiUtils/Data.h"
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {
namespace Optimization {

/**
 * @class Bisection @file Bisction
 * @brief Performs a bisection of a function in an interval [a,b].
 */
class Bisection {
 public:
  double getLowerBound() const {
    return lowerBound;
  }

  void setLowerBound(double lb) {
    lowerBound = lb;
  }

  double getUpperBound() const {
    return upperBound;
  }

  void setUpperBound(double ub) {
    upperBound = ub;
  }

  double getTolerance() const {
    return tolerance;
  }

  void setTolerance(double tol) {
    tolerance = tol;
  }

  int getMaxIterations() const {
    return maxIterations;
  }

  int getIterations() const {
    return iteration;
  }

  void setMaxIterations(int maxIt) {
    maxIterations = maxIt;
  }

  bool wasSuccessful() const {
    return wasSuccessful_;
  }

 private:
  double lowerBound = 1;
  double midPoint;
  double upperBound = 1000;

  double valueLB;
  double valueMidPoint;
  double valueUB;

  double tolerance = 1e-10;
  int iteration = 0;
  int maxIterations = 500;
  bool wasSuccessful_ = false;
  bool hasBoundaryError_ = false;

 public:
  bool hasBoundaryError() const {
    return hasBoundaryError_;
  }

 public:
  Bisection() = default;

  static auto sgn(double val) -> int {
    return (double(0) < val) - (val < double(0));
  }

  /**
   * @brief Bisection algorithm.
   * @tparam FunctionEvaluator Class that contains the function of which we want to find the zero.
   *         --> The FunctionEvaluator must have a method `FunctionEvaluator.evaluate(double) -> double`.
   * @param evaluator
   * @param verbose
   * @return solution of the bisection algorithm.
   */
  template<class FunctionEvaluator>
  auto compute(FunctionEvaluator& evaluator, bool verbose = true) -> double {
    if (verbose) {
      std::cout << "----------------------------\n";
      std::cout << " Bisection\n";
      std::cout << "----------------------------\n";
      std::cout << " Iteration   Error          \n";
      std::cout << "----------------------------\n";
    }

    if (upperBound < lowerBound) {
      throw std::runtime_error("Error in Bisection. Upper bound is lower than lower bound.");
    }

    valueLB = evaluator.evaluate(lowerBound);
    valueUB = evaluator.evaluate(upperBound);

    if (sgn(valueUB) == sgn(valueLB)) {
      std::cout << "Lower bound = " << valueLB << std::endl;
      std::cout << "Upper bound = " << valueUB << std::endl;
      std::cout << "Error in Bisection. Values at upper bound and lower bound have the same sign." << std::endl;
      hasBoundaryError_ = true;
      return 0;
    }

    iteration = 0;

    for (; iteration < maxIterations; ++iteration) {
      midPoint = (lowerBound + upperBound) / 2;

      valueMidPoint = evaluator.evaluate(midPoint);

      if (verbose) {
        std::cout << std::right << std::setw(10) << iteration << std::string(3, ' ') << std::setw(15) << std::scientific
                  << std::setprecision(5) << valueMidPoint << std::endl;
      }

      if (std::abs(valueMidPoint) < tolerance || upperBound - lowerBound / 2 < tolerance) {
        wasSuccessful_ = true;
        break;
      }

      if (sgn(valueMidPoint) == sgn(valueLB)) {
        lowerBound = midPoint;
      }
      else {
        upperBound = midPoint;
      }
    }

    if (verbose) {
      std::cout << "----------------------------\n";
    }

    return midPoint;
  }
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_BISECTION_H
