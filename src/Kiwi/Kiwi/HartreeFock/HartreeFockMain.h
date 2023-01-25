/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_HARTREEFOCKMAIN_H
#define KIWI_HARTREEFOCKMAIN_H

#include <Utils/Geometry/ElementTypes.h>
#include <memory>

namespace Scine {
namespace Kiwi {

class Molecule;
class Data;
class TRAHSettings;
class HartreeFockSettings;

class HartreeFockMain {
 public:
  HartreeFockMain(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data, HartreeFockSettings& settings,
                  bool hasGuess = false, bool verbose = true);

  auto setTrahSettings(const std::shared_ptr<TRAHSettings>& trahSettings) -> void {
    trahSettings_ = trahSettings;
  };

  auto scf() -> void;

  static auto makeGuess(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data, HartreeFockSettings& settings,
                        bool verbose = true) -> void;

  auto singleIteration() -> void;

  auto writeOrbitals() -> void;

 private:
  std::shared_ptr<Molecule> molecule_;
  std::shared_ptr<Data> data_;
  HartreeFockSettings& settings_;
  std::shared_ptr<TRAHSettings> trahSettings_;
  bool verbose_ = true;

 public:
  void setVerbose(bool verbose);
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_HARTREEFOCKMAIN_H
