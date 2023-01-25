/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ParticleDensityTask.h"
#include "Keywords.h"
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Kiwi/KiwiUtils/ParticleDensity/Evaluator.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/Yaml.h>
#include <yaml-cpp/yaml.h>
#include <iomanip> // std::setprecision
#include <iostream>

using namespace Scine;
using namespace Kiwi;

ParticleDensityTask::ParticleDensityTask(YAML::Node& input) {
  settings(input);
}

auto ParticleDensityTask::name() const -> const std::string {
  return "ParticleDensity";
}

auto ParticleDensityTask::settings(YAML::Node& input) -> void {
  Scine::Utils::checkYamlKeyRecognition(input, Scine::Kiwi::Keywords::particleDensitySettings);

  std::vector<std::string> requiredKeywords = {"format", "mode", "particle type", "type"};

  for (const auto& word : requiredKeywords) {
    if (!input[word]) {
      throw std::runtime_error("Keyword " + word + " does not exist. Terminating program!");
    }
    if (input[word].IsNull()) {
      throw std::runtime_error("Keyword " + word + " is empty. Terminating program!");
    }
  }

  if (input["print keywords"]) {
    auto printKeywords = input["print keywords"].as<bool>();
    if (printKeywords) {
      printAllowedKeywords(Scine::Kiwi::Keywords::particleDensitySettings);
    }
  }
  if (input["particle type"]) {
    settings_.type = Utils::ElementInfo::elementTypeForSymbol(input["particle type"].as<std::string>());
  }
  if (input["rdm file"]) {
    rdmFileNameRestricted_ = input["rdm file"].as<std::string>();
    useExternalRdm_ = true;
  }
  if (input["rdm file alpha"]) {
    rdmFileNameAlpha_ = input["rdm file alpha"].as<std::string>();
    useExternalRdm_ = true;
  }
  if (input["rdm file beta"]) {
    rdmFileNameBeta_ = input["rdm file beta"].as<std::string>();
    useExternalRdm_ = true;
  }
  if (input["xmin"]) {
    settings_.xmin = Utils::toBohr(Utils::Angstrom(input["xmin"].as<double>()));
  }
  if (input["ymin"]) {
    settings_.ymin = Utils::toBohr(Utils::Angstrom(input["ymin"].as<double>()));
  }
  if (input["zmin"]) {
    settings_.zmin = Utils::toBohr(Utils::Angstrom(input["zmin"].as<double>()));
  }
  if (input["xmax"]) {
    settings_.xmax = Utils::toBohr(Utils::Angstrom(input["xmax"].as<double>()));
  }
  if (input["ymax"]) {
    settings_.ymax = Utils::toBohr(Utils::Angstrom(input["ymax"].as<double>()));
  }
  if (input["zmax"]) {
    settings_.zmax = Utils::toBohr(Utils::Angstrom(input["zmax"].as<double>()));
  }
  if (input["NX"]) {
    settings_.N3DX = input["NX"].as<int>();
  }
  if (input["NY"]) {
    settings_.N3DY = input["NY"].as<int>();
  }
  if (input["NZ"]) {
    settings_.N3DZ = input["NZ"].as<int>();
  }
  if (input["N"]) {
    settings_.N1D = input["N"].as<int>();
  }
  if (input["mo index"]) {
    settings_.moIndex = input["mo index"].as<int>();
  }
  if (input["ms index"]) {
    settings_.msIndex = input["ms index"].as<int>();
  }
  if (input["mode"]) {
    auto modeStr = input["mode"].as<std::string>();
    std::transform(modeStr.begin(), modeStr.end(), modeStr.begin(), ::toupper);
    if (modeStr == "DENSITY" || modeStr == "DENS") {
      settings_.densityMode = ParticleDensity::DensityMode::PartDens;
    }
    else if (modeStr == "MO" || modeStr == "MOS") {
      settings_.densityMode = ParticleDensity::DensityMode::MO;
    }
    else {
      throw std::runtime_error("... don't know this mode.");
    }
  }
  if (input["format"]) {
    auto formatStr = input["format"].as<std::string>();
    std::transform(formatStr.begin(), formatStr.end(), formatStr.begin(), ::toupper);
    if (formatStr == "GRID") {
      settings_.format = ParticleDensity::Format::Grid;
    }
    else if (formatStr == "CUBE") {
      settings_.format = ParticleDensity::Format::Cube;
    }
    else {
      throw std::runtime_error("... don't know this format.");
    }
  }
  if (input["N"]) {
    settings_.N3DX = input["N"].as<double>();
    settings_.N3DY = input["N"].as<double>();
    settings_.N3DZ = input["N"].as<double>();
  }
}

auto ParticleDensityTask::run() -> void {
  std::cout << "\n";

  settings_.filename = molecule->getInputFileName();
  settings_.filename += ".particle_denstiy";
  if (!useExternalRdm_) {
    settings_.filename += ".HF.";
  }
  if (useExternalRdm_) {
    settings_.filename += ".RDM.";
  }
  settings_.filename += Utils::ElementInfo::symbol(settings_.type);
  settings_.filename += ".";
  switch (settings_.format) {
    case ParticleDensity::Format::Grid:
      settings_.filename += "grid";
      break;
    case ParticleDensity::Format::Cube:
      settings_.filename += "cube";
      break;
  }

  ParticleDensity::Evaluator evaluator(settings_, this->data, this->molecule);
  Utils::SpinAdaptedMatrix D;
  if (useExternalRdm_) {
    if (this->molecule->at(settings_.type).isRestricted) {
      Eigen::MatrixXd RDM = csvFileReader(rdmFileNameRestricted_);
      D.setRestrictedMatrix(RDM);
    }
    else {
      Eigen::MatrixXd RDMa = csvFileReader(rdmFileNameAlpha_);
      std::cout << "RDM-alpha\n" << RDMa << std::endl << std::endl;
      D.setAlphaMatrix(RDMa);
      if (this->molecule->at(settings_.type).msVector[1] > 0) {
        std::ifstream f(rdmFileNameBeta_.c_str());
        if (!f.good()) {
          throw std::runtime_error("Problem with RDM file for beta particles.");
        }
        Eigen::MatrixXd RDMb = csvFileReader(rdmFileNameBeta_);
        D.setBetaMatrix(RDMb);
      }
    }
  }
  else {
    if (this->molecule->at(settings_.type).isRestricted) {
      Eigen::MatrixXd RDM =
          Eigen::MatrixXd::Zero(this->molecule->at(settings_.type).LAO, this->molecule->at(settings_.type).LAO);
      for (auto i = 0; i < int(this->molecule->at(settings_.type).N / 2); ++i) {
        RDM(i, i) = 2;
      }
      D.setRestrictedMatrix(RDM);
    }
    else {
      Eigen::MatrixXd RDMa =
          Eigen::MatrixXd::Zero(this->molecule->at(settings_.type).LAO, this->molecule->at(settings_.type).LAO);
      for (auto i = 0; i < int(this->molecule->at(settings_.type).msVector[0]); ++i) {
        RDMa(i, i) = 1;
      }
      D.setAlphaMatrix(RDMa);
      if (this->molecule->at(settings_.type).msVector[0] > 0) {
        Eigen::MatrixXd RDMb =
            Eigen::MatrixXd::Zero(this->molecule->at(settings_.type).LAO, this->molecule->at(settings_.type).LAO);
        for (auto i = 0; i < int(this->molecule->at(settings_.type).msVector[1]); ++i) {
          RDMb(i, i) = 1;
        }
        D.setBetaMatrix(RDMb);
      }
    }
  }
  evaluator.setDensityMatrix(D);

  switch (settings_.format) {
    case ParticleDensity::Format::Grid: {
      evaluator.constructGrid();
      evaluator.runCalculationGrid();
      evaluator.writeGrid();
    } break;
    case ParticleDensity::Format::Cube: {
      evaluator.constructCube();
      evaluator.runCalculationCube();
      evaluator.writeCube();
    } break;
  }
  std::cout << "\n";
}
