/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_INTERFACE_H
#define KIWI_INTERFACE_H

#include <Kiwi/KiwiOpt/ObjectiveFunction.h>
#include <Kiwi/KiwiOpt/Optimizer.h>
#include <Eigen/Dense>
#include <deque>
#include <memory>

namespace Scine {
namespace Kiwi {
namespace Optimization {

class LineSearchInterface : public Kiwi::Optimization::Optimizer {
  std::shared_ptr<ObjectiveFunction> objectiveFunction;

  Eigen::VectorXd x_0;
  double value_0;
  Eigen::VectorXd gradient_0;

 private:
  Eigen::VectorXd direction;

  double alpha_n;
  Eigen::VectorXd x_alpha;
  Eigen::VectorXd gradient_alpha;

 public:
  virtual ~LineSearchInterface() = default;
  LineSearchInterface(std::shared_ptr<ObjectiveFunction> function, Eigen::VectorXd guess)
    : objectiveFunction(std::move(function)), x_0(std::move(guess)) {
    // gradient descent:
    objectiveFunction->evaluate(x_0);
    value_0 = objectiveFunction->getValue();
    gradient_0 = objectiveFunction->getGradient();
    direction = -gradient_0;
    phi_0_value = value_0;
    phi_0_derivative = gradient_0.dot(direction);
  }

  auto evaluate(double alpha) -> void override {
    alpha_n = alpha;
    x_alpha = x_0 + alpha_n * direction;
    objectiveFunction->evaluate(x_alpha);
    phi_alpha_value = objectiveFunction->getValue();
    gradient_alpha = objectiveFunction->getGradient();
    phi_alpha_derivative = gradient_alpha.dot(direction);
  }

  auto evaluateDirection() -> void override {
    // gradient descent:
    direction = -gradient_0;
    phi_0_value = value_0;
    phi_0_derivative = gradient_0.dot(direction);
    phi_alpha_value = phi_0_value;
    phi_alpha_derivative = phi_0_derivative;
  }

  auto setDirection(const Eigen::VectorXd& new_direction) -> void {
    direction = new_direction;
    phi_0_value = value_0;
    phi_0_derivative = gradient_0.dot(direction);
    phi_alpha_value = phi_0_value;
    phi_alpha_derivative = phi_0_derivative;
  }

  auto applyUpdate() -> void override {
    x_0 = x_alpha;
    value_0 = phi_alpha_value;
    gradient_0 = gradient_alpha;
  }

  auto getError() -> double override {
    return gradient_0.norm();
  }

  auto printIteration() -> void override {
  }

  auto getResult() const -> Eigen::VectorXd {
    return x_alpha;
  }
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_INTERFACE_H
