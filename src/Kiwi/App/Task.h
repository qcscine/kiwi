/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_TASK_H
#define KIWI_TASK_H

#include <memory>

namespace YAML {
class Node;
}

namespace Scine {
namespace Kiwi {

class Data;
class Molecule;

/**
 * @class Task
 * @brief Abstract base class for all Tasks
 */
class Task {
 protected:
  std::shared_ptr<Molecule> molecule;
  std::shared_ptr<Data> data;
  bool dataIsSet = false;

 public:
  explicit Task() = default;

  virtual auto name() const -> const std::string = 0;

  virtual auto run() -> void = 0;

  virtual auto setData(const std::shared_ptr<Data>& dat) -> void final {
    data = dat;
    dataIsSet = true;
  };

  virtual auto setMolecule(const std::shared_ptr<Molecule>& mol) -> void final {
    molecule = mol;
  };

  virtual std::shared_ptr<Molecule> getMolecule() const final {
    return molecule;
  };

  virtual std::shared_ptr<Data> getData() const final {
    return data;
  };
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_TASK_H
