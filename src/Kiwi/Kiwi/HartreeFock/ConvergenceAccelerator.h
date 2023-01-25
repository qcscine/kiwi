/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_CONVERGENCEACCELERATOR_H
#define KIWI_CONVERGENCEACCELERATOR_H

#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Scf/ConvergenceAccelerators/ConvergenceAcceleratorFactory.h>
#include <Utils/Scf/ConvergenceAccelerators/Ediis.h>
#include <Utils/Scf/ConvergenceAccelerators/FockDiis.h>
#include <Utils/Scf/ConvergenceAccelerators/FockSimple.h>
#include <iostream>
#include <memory>

namespace Scine {
namespace Kiwi {

enum class AcceleratorType { Diis, EDiis, NucleraElectronicDiis };

/**
 * @class ConvergenceAccelerator
 * @brief This class is an interface and handles the energy- and regular-DIIS for electrons and nuclei.
 *        I don't particularly like this interface, but it works.
 */
class ConvergenceAccelerator {
 private:
  std::map<Utils::ElementType, std::unique_ptr<Utils::Ediis>> _ediisVector;
  std::map<Utils::ElementType, std::size_t> _eDiisCounter;
  std::map<Utils::ElementType, std::unique_ptr<Utils::FockDiis>> _fockDiisVector;
  std::map<Utils::ElementType, std::size_t> _fockDiisCounter;
  std::map<Utils::ElementType, std::unique_ptr<Utils::FockSimple>> _fockSimpleVector;
  std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble> _errorMap;
  std::map<Utils::ElementType, Utils::SpinAdaptedMatrix> _oldFockMatrix;
  std::map<Utils::ElementType, bool> _fockDiisStarted;
  std::map<Utils::ElementType, bool> _ediisStarted;

  double _electronicEDiisThresh = 0.2;
  double _electronicFockDiisThresh = 0.02;
  double _nuclearEDiisThresh = 0.2;
  double _nuclearFockDiisThresh = 0.02;
  int _subspaceSize = 15;

  bool first = true;

  std::size_t _restartAt = 15;

  std::shared_ptr<Molecule> _molecule;
  std::shared_ptr<Data> _data;
  Utils::scf_mixer_t _mixer;
  bool verbose_;

  inline auto print(Utils::ElementType type, const std::string& which) {
    if (verbose_) {
      std::cout << "                      *** Initiating " << which << " for type " << Utils::ElementInfo::symbol(type)
                << " ***" << std::endl;
    }
  }

 public:
  inline void setElectronicEDiisThresh(double electronicEDiisThresh) {
    _electronicEDiisThresh = electronicEDiisThresh;
  }
  inline void setElectronicFockDiisThresh(double electronicFockDiisThresh) {
    _electronicFockDiisThresh = electronicFockDiisThresh;
  }
  inline void setNuclearEDiisThresh(double nuclearEDiisThresh) {
    _nuclearEDiisThresh = nuclearEDiisThresh;
  }
  inline void setNuclearFockDiisThresh(double nuclearFockDiisThresh) {
    _nuclearFockDiisThresh = nuclearFockDiisThresh;
  }
  inline void setSubspaceSize(int subspaceSize) {
    _subspaceSize = subspaceSize;
  }
  inline void setMaxIterRestart(int maxIterRestart) {
    _restartAt = maxIterRestart;
  }

  ConvergenceAccelerator(std::shared_ptr<Molecule> molecule, std::shared_ptr<Data> data, Utils::scf_mixer_t mixer,
                         bool verbose = true);

  auto update(const std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble>& errorMap, std::size_t restartAt = 15)
      -> void;

  auto update(const std::map<Utils::ElementType, HartreeFockUtils::SpinAdaptedDouble>& errorMap,
              Utils::ElementType type, std::size_t restartAt = 15) -> void;

  auto init() -> void;

  auto reset() -> void;

 private:
  auto _check(Utils::ElementType type, double FockDiisThresh, double EDiisThresh, bool isRestricted) -> Utils::scf_mixer_t;

  inline auto _checkType(Utils::ElementType type, bool isRestricted, double thresh) -> bool {
    if (!isRestricted) {
      if (_errorMap[type].alpha < thresh && _errorMap[type].beta < thresh) {
        return true;
      }
    }
    else {
      if (_errorMap[type].restricted < thresh) {
        return true;
      }
    }
    return false;
  }
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_CONVERGENCEACCELERATOR_H
