/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_HARTREEFOCKSETTINGS_H
#define KIWI_HARTREEFOCKSETTINGS_H

#include <Utils/Scf/ConvergenceAccelerators/ConvergenceAcceleratorFactory.h>

namespace Scine {
namespace Kiwi {

enum class NuclearElectronicSCF { Serial, Alternating, TRAH, BFGS };

enum class InitialGuessSCF { Core, Hueckel, BO, SAD, SadHueckel, SADNO, SND, Read };

struct HartreeFockSettings {
  std::size_t maxIterations = 500;
  std::size_t maxIterationsNested = 20;

  double energyThreshold = 1e-9;

  double nestedEnergyThreshold(double n = 1) const {
    return energyThreshold * n;
  }

  double fockMatrixNormThresh = 1e-5;

  Utils::scf_mixer_t accelerator = Utils::scf_mixer_t::fock_diis;

  double electronicDiisThresh = 1.0;

  double nuclearDiisThresh = 1.0;

  double electronicEDiisThresh = 20.0;

  double nuclearEDiisThresh = 20.0;

  double loewdinThresh = 1e-9;

  InitialGuessSCF guess = InitialGuessSCF::SAD;
  InitialGuessSCF nuclearGuess = InitialGuessSCF::SND;

  NuclearElectronicSCF neScfType = NuclearElectronicSCF::TRAH;

  bool verbose = true;

  bool stabilityAnalysis = false;

  int numberOfPerturbations = 0;

  bool useOrbitalSteering = false;

  int numberOfMixes = 10;

  double maxAngle = 0.2;

  bool integralDirect = true;

  bool disableIncrementalFock = true;

  std::size_t resetIncremental = 10;
};

struct TRAHSettings {
  enum class TrahOptimizer { ARH, Newton };

  double initialTrustRadius = 0.5;

  double maxTrustRadius = 1.0;

  int maxDavidsonIterations = 32;

  int maxDavidsonSubspaceDimension = 70;

  double microThreshGradScaling = 0.1;

  double microMinThresh = 0.01;

  double localRegionThresh = 0.001;

  double secondOrderThresh = 0.001;

  double minimumAlpha = 1;

  bool useSecondStartVector = true;

  bool useWhiteNoiseSecondStartVec = false;

  bool useThirdStartVector = true;

  bool useWhiteNoiseThirdStartVec = false;

  bool auto2ndOrder = false;

  TrahOptimizer opt = TrahOptimizer::ARH;

  int maxArhDim = 25;
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_HARTREEFOCKSETTINGS_H
