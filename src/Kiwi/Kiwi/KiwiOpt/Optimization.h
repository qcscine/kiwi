/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_OPTIMIZATION_H
#define KIWI_OPTIMIZATION_H

#include <Kiwi/KiwiOpt/MoreThuente.h>
#include <Kiwi/KiwiOpt/Optimizer.h>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {
namespace Optimization {

/**
 * @brief Main routine of the optimization.
 */
class Optimization {
 private:
  double thresh_ = 1e-16;

  std::shared_ptr<Optimizer> optimizer;

  int maxIterations_ = 100;
  int iterations_;

  double error_;

  double alpha_0 = 1.0;

  bool successful_ = false;

 public:
  Optimization(std::shared_ptr<Optimizer> opt) : optimizer(std::move(opt)) {
  }

  bool wasSuccessful() const {
    return successful_;
  }

  void setMaxIterations(int maxIterations) {
    maxIterations_ = maxIterations;
  }

  int getIterations() const {
    return iterations_;
  }

  double getError() const {
    return error_;
  }

  void setThresh(double thresh) {
    thresh_ = thresh;
  }

  auto optimize() -> void {
    successful_ = false;

    iterations_ = 1;
    int total_micro_iterations = 0;

    MoreThuente<Optimizer> moreThuente(*optimizer);

    std::cout << std::left << std::setw(10) << "Macro-it" << std::setw(10) << "Micro-it" << std::setw(20) << "Value"
              << std::setw(20) << "||g||" << std::endl;
    std::cout << std::string(2 * 16 + 2 * 20 + 2, '-') << std::endl;

    double currentThresh = 0.01;

    while (true) {
      if (iterations_ == maxIterations_) {
        std::cout << "GradientDescent did not convergege after " << iterations_ << " iterations.\n";
        break;
      }

      optimizer->evaluateDirection();

      moreThuente.setTol(currentThresh);
      moreThuente(alpha_0, optimizer->getZeroValue(), optimizer->getZeroDerivative());
      total_micro_iterations += moreThuente.getIterations();

      ++iterations_;

      optimizer->applyUpdate();

      error_ = optimizer->getError();

      if (error_ < currentThresh) {
        currentThresh = 0.1 * error_;
      }

      std::cout << std::right << std::setw(10) << iterations_ << std::setw(10) << moreThuente.getIterations();
      std::cout << std::setw(20) << std::fixed << std::setprecision(10) << optimizer->getLineSearchValue();
      std::cout << std::setw(20) << std::scientific << std::setprecision(10) << optimizer->getError();
      std::cout << std::endl;

      if (error_ < thresh_) {
        break;
      }
    }
    std::cout << std::endl
              << "Total number of gradient and function evaluations = " << total_micro_iterations << std::endl;
  }
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_OPTIMIZATION_H
