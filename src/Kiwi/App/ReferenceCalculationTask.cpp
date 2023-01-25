/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ReferenceCalculationTask.h"
#include "Keywords.h"
#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Kiwi/KiwiUtils/Loewdin.h>
#include <Utils/IO/Yaml.h>
#include <yaml-cpp/yaml.h>

Scine::Kiwi::ReferenceCalculationTask::ReferenceCalculationTask(YAML::Node& input) {
  settings(input);
}

auto Scine::Kiwi::ReferenceCalculationTask::name() const -> const std::string {
  return "Hartree-Fock";
}

auto Scine::Kiwi::ReferenceCalculationTask::settings(YAML::Node& input) -> void {
  Scine::Utils::checkYamlKeyRecognition(input, Scine::Kiwi::Keywords::referenceCalculationSettings);

  if (input["print keywords"]) {
    auto printKeywords = input["print keywords"].as<bool>();
    if (printKeywords) {
      printAllowedKeywords(Scine::Kiwi::Keywords::referenceCalculationSettings);
    }
  }
  if (input["scf"]) {
    doScf_ = input["scf"].as<bool>();
  }
  if (input["max iter"]) {
    settings_.maxIterations = input["max iter"].as<std::size_t>();
  }
  if (input["max iter nested"]) {
    settings_.maxIterationsNested = input["max iter nested"].as<std::size_t>();
  }
  if (input["thresh"]) {
    settings_.energyThreshold = input["thresh"].as<double>();
  }
  if (input["reset incremental"]) {
    settings_.resetIncremental = input["reset incremental"].as<std::size_t>();
  }
  if (input["disable incremental"]) {
    settings_.disableIncrementalFock = input["disable incremental"].as<bool>();
  }
  if (input["integral direct"]) {
    settings_.integralDirect = input["integral direct"].as<bool>();
  }
  if (input["accelerator"]) {
    auto which = input["accelerator"].as<std::string>();
    if (which == "ediis diis") {
      settings_.accelerator = Utils::scf_mixer_t::ediis_diis;
    }
    else if (which == "ediis") {
      settings_.accelerator = Utils::scf_mixer_t::ediis;
    }
    else if (which == "diis") {
      settings_.accelerator = Utils::scf_mixer_t::fock_diis;
    }
    else if (which == "none") {
      settings_.accelerator = Utils::scf_mixer_t::none;
    }
    else {
      throw std::runtime_error("Don't know this accelerator. Available are: `diis`, `ediis`, `ediis diis`, `none`.");
    }
  }
  if (input["thresh fock"]) {
    settings_.fockMatrixNormThresh = input["thresh fock"].as<double>();
  }
  if (input["mix angle"]) {
    settings_.maxAngle = input["mix angle"].as<double>();
  }
  if (input["mixes"]) {
    settings_.numberOfMixes = input["mixes"].as<int>();
  }
  if (input["perturbations"]) {
    settings_.numberOfPerturbations = input["perturbations"].as<int>();
  }
  if (input["loewdin thresh"]) {
    settings_.loewdinThresh = input["loewdin thresh"].as<double>();
  }
  if (input["diis thresh"]) {
    settings_.electronicDiisThresh = input["diis thresh"].as<double>();
  }
  if (input["nuclear diis thresh"]) {
    settings_.nuclearDiisThresh = input["nuclear diis thresh"].as<double>();
  }
  if (input["ediis thresh"]) {
    settings_.electronicEDiisThresh = input["ediis thresh"].as<double>();
  }
  if (input["nuclear ediis thresh"]) {
    settings_.nuclearEDiisThresh = input["nuclear ediis thresh"].as<double>();
  }
  if (input["nuclear guess"]) {
    if (input["nuclear guess"].as<std::string>() == "core") {
      settings_.nuclearGuess = InitialGuessSCF::Core;
    }
    else if (input["nuclear guess"].as<std::string>() == "snd") {
      settings_.nuclearGuess = InitialGuessSCF::SND;
    }
    else {
      throw std::runtime_error("Valid options for nuclear guess: core, snd.");
    }
  }
  if (input["guess"]) {
    if (input["guess"].as<std::string>() == "read") {
      settings_.guess = InitialGuessSCF::Read;
    }
    if (input["guess"].as<std::string>() == "core") {
      settings_.guess = InitialGuessSCF::Core;
    }
    if (input["guess"].as<std::string>() == "sad") {
      settings_.guess = InitialGuessSCF::SAD;
    }
    if (input["guess"].as<std::string>() == "sadno") {
      settings_.guess = InitialGuessSCF::SADNO;
    }
    if (input["guess"].as<std::string>() == "bo") {
      settings_.guess = InitialGuessSCF::BO;
    }
    else if (input["guess"].as<std::string>() == "sad hueckel") {
      settings_.guess = InitialGuessSCF::SadHueckel;
    }
    else if (input["guess"].as<std::string>() == "hueckel") {
      settings_.guess = InitialGuessSCF::Hueckel;
    }
  }
  if (input["verbose"]) {
    settings_.verbose = input["verbose"].as<bool>();
  }
  if (input["stability analysis"]) {
    settings_.stabilityAnalysis = input["stability analysis"].as<bool>();
  }
  if (input["scf type"]) {
    if (input["scf type"].as<std::string>() == "serial") {
      settings_.neScfType = NuclearElectronicSCF::Serial;
    }
    else if (input["scf type"].as<std::string>() == "bfgs") {
      settings_.neScfType = NuclearElectronicSCF::BFGS;
    }
    else if (input["scf type"].as<std::string>() == "alternating") {
      settings_.neScfType = NuclearElectronicSCF::Alternating;
    }
    else if (input["scf type"].as<std::string>() == "trah") {
      settings_.neScfType = NuclearElectronicSCF::TRAH;
    }
    else {
      throw std::runtime_error("Do not know this scf type.");
    }
  }
  if (auto trah = input["trah settings"]) {
    Scine::Utils::checkYamlKeyRecognition(trah, Scine::Kiwi::Keywords::trahSettings);

    if (input["print keywords"]) {
      auto printKeywords = input["print keywords"].as<bool>();
      if (printKeywords) {
        printAllowedKeywords(Scine::Kiwi::Keywords::trahSettings);
      }
    }
    if (trah["optimizer"]) {
      if (trah["optimizer"].as<std::string>() == "arh") {
        trahSettings_.opt = TRAHSettings::TrahOptimizer::ARH;
      }
      else if (trah["optimizer"].as<std::string>() == "newton") {
        trahSettings_.opt = TRAHSettings::TrahOptimizer::Newton;
      }
      else {
        throw std::runtime_error("Don't know this TRAH optimizer.");
      }
    }
    if (trah["initial trust radius"]) {
      trahSettings_.initialTrustRadius = trah["initial trust radius"].as<double>();
    }
    if (trah["max trust radius"]) {
      trahSettings_.maxTrustRadius = trah["max trust radius"].as<double>();
    }
    if (trah["max davidson iterations"]) {
      trahSettings_.maxDavidsonIterations = trah["max davidson iterations"].as<int>();
    }
    if (trah["max davidson subspace dim"]) {
      trahSettings_.maxDavidsonSubspaceDimension = trah["max davidson subspace dim"].as<int>();
    }
    if (trah["grad scaling"]) {
      trahSettings_.microThreshGradScaling = trah["grad scaling"].as<double>();
    }
    if (trah["max arh dim"]) {
      trahSettings_.maxArhDim = trah["max arh dim"].as<int>();
    }
    if (trah["local thresh"]) {
      trahSettings_.localRegionThresh = trah["local thresh"].as<double>();
    }
    if (trah["2nd order thresh"]) {
      trahSettings_.secondOrderThresh = trah["2nd order thresh"].as<double>();
    }
    if (trah["min thresh"]) {
      trahSettings_.microMinThresh = trah["min thresh"].as<double>();
    }
    if (trah["auto 2nd order"]) {
      trahSettings_.auto2ndOrder = trah["auto 2nd order"].as<bool>();
    }
    if (trah["2nd start vector"]) {
      trahSettings_.useSecondStartVector = trah["2nd start vector"].as<bool>();
    }
    if (trah["2nd start vector noise"]) {
      trahSettings_.useWhiteNoiseSecondStartVec = trah["2nd start vector noise"].as<bool>();
    }
    if (trah["3rd start vector"]) {
      trahSettings_.useThirdStartVector = trah["3rd start vector"].as<bool>();
    }
    if (trah["3rd start vector noise"]) {
      trahSettings_.useWhiteNoiseThirdStartVec = trah["3rd start vector noise"].as<bool>();
    }
  }
}

auto Scine::Kiwi::ReferenceCalculationTask::run() -> void {
  std::cout << "\n-----------------------------------------------------------\n";
  std::cout << "  Integral evaluation\n";
  std::cout << "-----------------------------------------------------------\n\n";

  bool hasGuess = true;

  // omp_set_num_threads(4);

  if (!this->dataIsSet) {
    hasGuess = false;
    this->data = std::make_shared<Data>(this->molecule);
    this->data->oneBodyIntegrals();

    this->data->integralDirect = settings_.integralDirect;
    this->data->twoBodyIntegrals();
    this->data->makeExchange();

    std::cout << "\n-----------------------------------------------------------\n";
    std::cout << "  Iterative canonical Loewdin orthonormalization\n";
    std::cout << "-----------------------------------------------------------\n\n";

    for (auto const& elem : *this->molecule) {
      std::cout << "Particle type: " << Utils::ElementInfo::symbol(elem.first) << std::endl;
      this->data->X[elem.first] = Kiwi::Loewdin(this->data->S[elem.first], settings_.loewdinThresh, true);
      this->molecule->at(elem.first).LMO = this->data->X[elem.first].cols();
      std::cout << std::endl;
    }
  }

  // omp_set_num_threads(4);

  std::cout << "\n-----------------------------------------------------------\n";
  std::cout << "  Self-consistent field procedure\n";
  std::cout << "-----------------------------------------------------------\n\n";

  HartreeFockMain hartreeFockMain(this->molecule, this->data, settings_, hasGuess);

  std::shared_ptr<TRAHSettings> ptrTrahSettings = std::make_shared<TRAHSettings>(trahSettings_);

  hartreeFockMain.setTrahSettings(ptrTrahSettings);

  if (doScf_) {
    hartreeFockMain.scf();
  }

  hartreeFockMain.writeOrbitals();
}
