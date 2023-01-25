/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_REFERENCECALCULATIONTASK_H
#define KIWI_REFERENCECALCULATIONTASK_H

#include "Task.h"
#include <Kiwi/HartreeFock/HartreeFockSettings.h>

namespace Scine {
namespace Kiwi {

class ReferenceCalculationTask : public Task {
 private:
  HartreeFockSettings settings_;
  TRAHSettings trahSettings_;

  bool doScf_ = true;

 public:
  // ReferenceCalculationTask() = default;

  ReferenceCalculationTask(YAML::Node& input);

  // ReferenceCalculationTask(const ReferenceCalculationTask&) = delete;
  // ReferenceCalculationTask& operator=(const ReferenceCalculationTask&) = delete;
  //~ReferenceCalculationTask() override = default;

  auto name() const -> const std::string final;

  auto run() -> void final;

  auto settings(YAML::Node& input) -> void;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_REFERENCECALCULATIONTASK_H
