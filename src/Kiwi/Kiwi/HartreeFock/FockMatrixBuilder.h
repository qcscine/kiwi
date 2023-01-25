/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_FOCKMATRIXBUILDER_H
#define KIWI_FOCKMATRIXBUILDER_H

#include <Kiwi/KiwiUtils/Data.h>
#include <chrono>

namespace Scine {
namespace Kiwi {

class FockMatrixBuilder {
 public:
  FockMatrixBuilder() = default;

  FockMatrixBuilder(Utils::ElementType type, std::shared_ptr<Data> data, std::shared_ptr<Molecule> molecule,
                    std::shared_ptr<HartreeFockData> hartreeFockData);

  // auto update(std::map<Utils::ElementType, Utils::DensityMatrix>& D, bool incremental = false, std::size_t
  // numRebuilds = 0)
  //    -> void;

  auto updateGeneralized(std::map<Utils::ElementType, Utils::DensityMatrix>& D, bool verbose = false) -> void;

  auto init() -> void;

  auto initGeneralized() -> void;

  auto getElapsedTime() -> std::chrono::duration<double, std::milli> {
    return elapsedTime_ / fockMatrixFunctionCalls_;
  }

 private:
  Utils::ElementType type_;

  std::shared_ptr<Data> data_;

  std::shared_ptr<Molecule> molecule_;

  std::shared_ptr<HartreeFockData> hartreeFockData_;

 public:
  const std::shared_ptr<HartreeFockData>& getHartreeFockData() const;

 private:
  template<SpinSymmetry>
  auto buildJKmatrices(std::map<Utils::ElementType, Utils::DensityMatrix>& D) -> void;

  template<SpinSymmetry>
  auto buildFockmatrices(std::map<Utils::ElementType, Utils::DensityMatrix>& D) -> void;

  auto buildLmatrices(std::map<Utils::ElementType, Utils::DensityMatrix>& D) -> void;

  // auto buildJKmatricesDirect(std::map<Utils::ElementType, Utils::DensityMatrix>& D) -> void;

  // auto buildLmatricesDirect(std::map<Utils::ElementType, Utils::DensityMatrix>& D) -> void;

  std::chrono::duration<double, std::milli> elapsedTime_;

  int fockMatrixFunctionCalls_ = 0;

  bool formIncremental_;

  std::size_t numRebuilds_;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_FOCKMATRIXBUILDER_H
