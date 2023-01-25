/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_CALCULATIONMANAGER_H
#define KIWI_CALCULATIONMANAGER_H

#include <map>
#include <memory>
#include <ostream>
#include <vector>

namespace YAML {
class Node;
}

namespace Scine {
namespace Kiwi {

class Molecule;
class Task;

class CalculationManager {
 public:
  CalculationManager(YAML::Node& input, const std::string& inputFileName);

  CalculationManager(const CalculationManager&) = delete;

  CalculationManager& operator=(const CalculationManager&) = delete;

  ~CalculationManager();

  auto initializeTasks() -> void;

  auto executeTasks() -> void;

 private:
  std::vector<std::unique_ptr<Task>> tasks_;

  std::shared_ptr<Molecule> molecule_;

  YAML::Node& input_;

  const std::string inputFileName_;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_CALCULATIONMANAGER_H
