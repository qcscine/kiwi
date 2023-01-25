/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_TASKFACTORY_H
#define KIWI_TASKFACTORY_H

#include "Ao2MoTask.h"
#include "NaturalOrbitalsTask.h"
#include "ParticleDensityTask.h"
#include "ReferenceCalculationTask.h"
#include "Task.h"
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace YAML {
class Node;
}

namespace Scine {
namespace Kiwi {

/**
 * @brief A factory generating Tasks by name.
 */
class TaskFactory {
 public:
  /// @brief Has only static functions.
  TaskFactory() = delete;
  /**
   * @brief Contstructs a Task with a given set of input and output systems.
   *
   * @param name   The name of the requested task.
   * @return std::unique_ptr<Task> The requested task.
   */
  static std::unique_ptr<Task> produce(std::string name, YAML::Node& input) {
    std::unique_ptr<Task> task;
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    if (name == "HF" || name == "HARTREE FOCK" || name == "HARTREE-FOCK") {
      task = std::make_unique<ReferenceCalculationTask>(input);
    }
    else if (name == "READ" || name == "READ MOS") {
      throw std::runtime_error("MO reading not implemented, yet.");
    }
    else if (name == "AO2MO" || name == "AO 2 MO" || name == "AO TO MO" || name == "AOTOMO") {
      task = std::make_unique<Ao2MoTask>(input);
    }
    else if (name == "PARTICLE DENSITY" || name == "DENSITY") {
      task = std::make_unique<ParticleDensityTask>(input);
    }
    else if (name == "NATURAL ORBITALS" || name == "NATURALORBITALS") {
      task = std::make_unique<NaturalOrbitalsTask>(input);
    }
    else {
      throw std::runtime_error("The requested task '" + name + "' is not available.\n");
    }
    return task;
  }
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_TASKFACTORY_H
