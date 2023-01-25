/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_OBJECTIVEFUNCTION_H
#define KIWI_OBJECTIVEFUNCTION_H

#include <Eigen/Dense>

namespace Scine {
namespace Kiwi {
namespace Optimization {

class ObjectiveFunction {
 protected:
  Eigen::VectorXd gradient;
  Eigen::VectorXd parameters;
  Eigen::MatrixXd Hessian;

  double value;

 public:
  virtual auto evaluate(const Eigen::VectorXd& x) -> void = 0;

  auto getGradient() const -> const Eigen::VectorXd& {
    return gradient;
  }

  auto getParameters() const -> const Eigen::VectorXd& {
    return parameters;
  }

  auto getValue() const -> double {
    return value;
  }

  auto getHessian() const -> const Eigen::MatrixXd& {
    return Hessian;
  }
};

} // namespace Optimization
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_OBJECTIVEFUNCTION_H