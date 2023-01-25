/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* internal */
#include <Kiwi/HartreeFock/FockMatrixBuilder.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/HartreeFock/Projector/IdentityProjector.h>
#include <Kiwi/KiwiUtils/Data.h>

//#include <Kiwi/HartreeFock/FockMatrix/BoHelper.h>
//#include <Kiwi/HartreeFock/FockMatrix/PreBoHelper.h>
/* external */
#include <iomanip>
#include <iostream>

namespace Scine {
namespace Kiwi {

auto IdentityProjector::init(bool verbose) -> void {
  // builder_ = std::make_shared<std::map<Utils::ElementType, FockMatrixBuilder>>();
  boHelper_ = std::make_unique<std::map<Utils::ElementType, FockMatrix::BoHelper>>();

  if (molecule_->size() > 1) {
    preBoHelper_ = std::make_unique<FockMatrix::PreBoHelper>(data_, molecule_, data_->hartreeFockData);
    preBoHelper_->updateL(data_->D, data_->integralDirect);
  }

  for (auto const& elem : *molecule_) {
    boHelper_->insert({elem.first, FockMatrix::BoHelper(elem.first, data_, molecule_, data_->hartreeFockData)});
  }
  for (auto& helper : *boHelper_) {
    helper.second.update(data_->D, false, 0);
  }

  auto energy = HartreeFockUtils::getEnergy(molecule_, data_);

  if (verbose) {
    std::cout << std::endl
              << "Guess energy = " << std::setprecision(10) << std::fixed << energy << " Ha" << std::endl
              << std::endl;
  }

  energy_ = energy;
}

auto IdentityProjector::updateFockMatrices(bool incremental, std::size_t numFockRebuilds) -> void {
  if (molecule_->size() > 1) {
    preBoHelper_->updateL(data_->D, data_->integralDirect);
  }
  for (auto& helper : *boHelper_) {
    helper.second.update(incrementalDensityMatrix_, incremental, numFockRebuilds);
  }
}

auto IdentityProjector::evaluateEnergy() -> void {
  energy_ = HartreeFockUtils::getEnergy(molecule_, data_);
}
auto IdentityProjector::getErrorMap() -> std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble> {
  errorMap_ = HartreeFockUtils::getDensityGradientErrorMapOrtho(molecule_, data_);
  return errorMap_;
}

auto IdentityProjector::getFockMatrixFormationTime() -> std::chrono::duration<double, std::milli> {
  fockMatrixFormationTime_ = std::chrono::duration<double, std::milli>::zero();
  ;

  if (molecule_->size() > 1) {
    fockMatrixFormationTime_ += preBoHelper_->getElapsedTime();
  }

  for (auto& b : *boHelper_) {
    fockMatrixFormationTime_ += b.second.getElapsedTime();
  }
  return fockMatrixFormationTime_;
}

} // namespace Kiwi
} // namespace Scine
