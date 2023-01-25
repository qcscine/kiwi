/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NaturalOrbitalsTask.h"
#include "Keywords.h"
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Kiwi/KiwiUtils/NaturalOrbitals.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/Yaml.h>
#include <yaml-cpp/yaml.h>
#include <iomanip> // std::setprecision
#include <iostream>

using namespace Scine;
using namespace Kiwi;

NaturalOrbitalsTask::NaturalOrbitalsTask(YAML::Node& input) {
  settings(input);
}

auto NaturalOrbitalsTask::name() const -> const std::string {
  return "NaturalOrbitals";
}

auto NaturalOrbitalsTask::settings(YAML::Node& input) -> void {
  Scine::Utils::checkYamlKeyRecognition(input, Scine::Kiwi::Keywords::naturalOrbitalSettings);

  std::vector<std::string> requiredKeywords = {"type", "rdm files"};

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
      printAllowedKeywords(Scine::Kiwi::Keywords::naturalOrbitalSettings);
    }
  }
  for (auto const& it : input["rdm files"]) {
    auto tp = Utils::ElementInfo::elementTypeForSymbol(it.first.as<std::string>());
    elementType_LineMap_[tp] = it.second.as<std::string>();
  }
}

auto NaturalOrbitalsTask::run() -> void {
  NaturalOrbitals naturalOrbitals(data, molecule, elementType_LineMap_);

  naturalOrbitals.generateNaturalOrbitals();

  std::cout << "\n"
            << "Natural orbitals were constructed from external RDMs.\n";
  std::cout << "Attention! External natural orbitals can only be used for:\n";
  std::cout << "  - Obtaining an FCIDUMP file in terms of the natural orbitals.\n";
  std::cout << "  - Obtaining particle densities with an external RDM.\n";
  std::cout << "\n\n";
}
