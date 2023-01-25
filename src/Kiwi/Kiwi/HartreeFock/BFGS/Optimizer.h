/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_FIRSTORDEROPTIMIZER_H
#define KIWI_FIRSTORDEROPTIMIZER_H

#include <Kiwi/HartreeFock/BFGS/Interface.h>
#include <Kiwi/KiwiOpt/Optimizer.h>
#include <iomanip>
#include <iostream>
#include <utility>

namespace Scine {
namespace Kiwi {
namespace BFGS {

template<SymmetryType Symmetry>
class Optimizer : public Optimization::Optimizer {
 private:
  std::shared_ptr<BFGS::Interface<Symmetry>> interface;
  const std::vector<int> dim;
  // Number of different types and spins
  const int numElements;

  std::vector<int> offset;

  Eigen::VectorXd gradient_0;

  Eigen::VectorXd gradient_alpha;

  Eigen::VectorXd direction;

  double gradientNorm;

  int fullDim;

  // Preconditioner
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> precond;
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> eye;

 public:
  // Dimension of all bases combined

  Optimizer(std::shared_ptr<BFGS::Interface<Symmetry>> firstOrderInterface)
    : interface(std::move(firstOrderInterface)),
      dim(interface->getDim()),
      numElements(dim.size()),
      offset(interface->getOffset()),
      fullDim(std::accumulate(dim.begin(), dim.end(), 0)) {
    gradient_0.resize(fullDim);
    gradient_alpha.resize(fullDim);
    direction.resize(fullDim);

    interface->evaluateGradient();
    for (int i = 0; i < numElements; ++i) {
      gradient_0.segment(offset[i], dim[i]) = interface->getGradient(i);
    }
    interface->evaluateDirection();
    populateDirection();

    phi_0_value = interface->getEnergyAlpha();
    phi_0_derivative = gradient_0.dot(direction);
    std::cout << "phi(0)' = " << phi_0_derivative << std::endl;

    gradientNorm = gradient_0.lpNorm<2>();
  }

  [[nodiscard]] double getError() const {
    return gradientNorm;
  }

  auto evaluateDirection() -> void final {
    interface->evaluateDirection();

    populateDirection();

    phi_0_derivative = gradient_0.dot(direction);
  }

  auto getError() -> double final {
    return gradientNorm;
  }

  auto evaluate(double alpha) -> void final {
    current_alpha = alpha;

    interface->evaluate(alpha);
    populateGradient();

    phi_alpha_value = interface->getEnergyAlpha();
    phi_alpha_derivative = gradient_alpha.dot(direction);
  }

  auto applyUpdate() -> void final {
    // The last point of the last iteration will be the new zero-point
    // Optimization::Optimizer<double>::
    // Optimization::Optimizer<double>::

    interface->applyUpdate(true, current_alpha);
    // evaluate(0);
    // interface->applyUpdate(false, 0);

    phi_0_value = interface->getEnergyAlpha();

    for (int i = 0; i < numElements; ++i) {
      gradient_0.segment(offset[i], dim[i]) = interface->getGradient(i);
    }

    gradientNorm = gradient_0.lpNorm<2>();
  }

  auto printIteration() -> void {
    std::cout << "Hello there\n";
  }

 private:
  auto populateGradient() -> void {
    for (int i = 0; i < numElements; ++i) {
      gradient_alpha.segment(offset[i], dim[i]) = interface->getGradient(i);
    }
  }

  auto populateDirection() -> void {
    for (int i = 0; i < numElements; ++i) {
      direction.segment(offset[i], dim[i]) = interface->getDirection(i);
    }
  }
};

} // namespace BFGS
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_FIRSTORDEROPTIMIZER_H
