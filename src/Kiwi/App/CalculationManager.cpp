/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CalculationManager.h"
#include "TaskFactory.h"
#include "boost/exception/diagnostic_information.hpp"
#include "boost/filesystem.hpp"
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Kiwi/KiwiUtils/Molecule.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/Yaml.h>
#include <yaml-cpp/yaml.h>
#include <iostream>

using namespace Scine;

using namespace Kiwi;

CalculationManager::CalculationManager(YAML::Node& input, const std::string& inputFileName)
  : input_(input), inputFileName_(inputFileName) {
  // Check for invalid top level input_ sections
  std::vector<std::string> requiredKeywords{"tasks", "xyz file", "basis"};
  std::vector<std::string> allowedKeywords{"tasks", "xyz file",      "molecule",
                                           "basis", "nuclear basis", "xyz nuclear centers"};
  Scine::Utils::checkYamlKeyRecognition(input_, allowedKeywords);
  for (const auto& word : requiredKeywords) {
    if (!input_[word]) {
      throw std::runtime_error("Keyword " + word + " does not exist. Terminating program!");
    }
    if (input_[word].IsNull()) {
      throw std::runtime_error("Keyword " + word + " is empty. Terminating program!");
    }
  }

  // Load xyz file
  // `structureBoolVectorPair` stores which nuclei within the given structure should be treated quantum mechanically.
  std::pair<Utils::AtomCollection, std::vector<bool>> structureBoolVectorPair;

  MoleculeSettings molSettings;

  molSettings.multiplicity = 1;
  molSettings.charge = 0;
  molSettings.useHighSpinApproximation = true;
  molSettings.isRestricted = true;
  molSettings.onlyElectrons = true;

  bool uncontract = false;

  const std::string xyzFileName = input_["xyz file"].as<std::string>();

  if (!boost::filesystem::exists(xyzFileName)) {
    throw std::runtime_error("Specified xyz file file does not exist!");
  }
  std::ifstream xyfFile(xyzFileName);
  if (xyfFile.is_open()) {
    structureBoolVectorPair = Utils::XyzStreamHandler::readNuclearElectronic(xyfFile);
    xyfFile.close();
  }
  else {
    throw std::runtime_error("Unable to open xyz file.");
  }

  std::cout << "--------------------------------------------------------------------------------" << std::endl;
  std::cout << "                              Geometry (Angstrom)" << std::endl;
  std::cout << "--------------------------------------------------------------------------------" << std::endl;

  std::cout << std::setprecision(12) << std::fixed;

  for (auto i = 0UL; i < structureBoolVectorPair.second.size(); ++i) {
    auto atom = structureBoolVectorPair.first.at(i);
    std::cout << " " << std::left << std::setw(8) << Utils::ElementInfo::symbol(atom.getElementType()) << std::setw(20)
              << std::right << Utils::toAngstrom(Utils::Bohr(atom.getPosition()[0])) << std::setw(20)
              << Utils::toAngstrom(Utils::Bohr(atom.getPosition()[1])) << std::setw(20)
              << Utils::toAngstrom(Utils::Bohr(atom.getPosition()[2]));
    if (structureBoolVectorPair.second.at(i)) {
      std::cout << std::setw(8) << "Q";
    }
    std::cout << std::endl;
  }

  std::cout << "--------------------------------------------------------------------------------" << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;

  // Check if the molecule_ contains quantum nuclei:
  if (std::find(structureBoolVectorPair.second.begin(), structureBoolVectorPair.second.end(), true) !=
      structureBoolVectorPair.second.end()) {
    molSettings.onlyElectrons = false;
  }

  // Get molecule_ input_
  if (input_["molecule"]) {
    const auto& moleculeInput = input_["molecule"];
    if (moleculeInput.size() == 0) {
      throw std::runtime_error("The molecule input is empty.");
    }
    // Check for invalid task input_ sections
    std::vector<std::string> allowedKeywordsMol = {"mult",      "multiplicity",   "charge",    "restricted",
                                                   "high spin", "pure spherical", "uncontract"};
    Scine::Utils::checkYamlKeyRecognition(moleculeInput, allowedKeywordsMol);
    if (moleculeInput["mult"]) {
      molSettings.multiplicity = moleculeInput["mult"].as<std::size_t>();
    }
    if (moleculeInput["multiplicity"]) {
      molSettings.multiplicity = moleculeInput["multiplicity"].as<std::size_t>();
    }
    if (moleculeInput["charge"]) {
      molSettings.charge = moleculeInput["charge"].as<int>();
    }
    if (moleculeInput["high spin"]) {
      molSettings.useHighSpinApproximation = moleculeInput["high spin"].as<bool>();
    }
    if (moleculeInput["restricted"]) {
      molSettings.isRestricted = moleculeInput["restricted"].as<bool>();
    }
    if (moleculeInput["pure spherical"]) {
      molSettings.usePureSpherical = moleculeInput["pure spherical"].as<bool>();
    }
    if (moleculeInput["uncontract"]) {
      uncontract = moleculeInput["uncontract"].as<bool>();
    }
  }
  std::unordered_map<Utils::ElementType, std::string> basisSetNameMap;
  std::unordered_map<Utils::ElementType, std::string> electronBasisSetNameMap;
  // Default basis set input_:
  if (input_["basis"].Type() == YAML::NodeType::Scalar) {
    auto name = input_["basis"].as<std::string>();
    // This is the BO-branch:
    if (molSettings.onlyElectrons) {
      molecule_ = std::make_unique<Kiwi::Molecule>(name, structureBoolVectorPair, molSettings);
    }
    else {
      basisSetNameMap[Utils::ElementType::E] = name;
    }
  }
  else if (input_["basis"].Type() == YAML::NodeType::Map) {
    for (auto const& it : input_["basis"]) {
      auto tp = Utils::ElementInfo::elementTypeForSymbol(it.first.as<std::string>());
      electronBasisSetNameMap[tp] = it.second.as<std::string>();
    }
    // This is a BO-branch
    if (molSettings.onlyElectrons) {
      molecule_ = std::make_unique<Kiwi::Molecule>(electronBasisSetNameMap, structureBoolVectorPair, molSettings);
    }
  }
  else {
    throw std::runtime_error("Format not recognized.");
  }
  // Pre-BO branches:
  if (!molSettings.onlyElectrons) {
    if (!input_["nuclear basis"]) {
      throw std::runtime_error("Please provide nuclear basis.");
    }
    for (auto const& it : input_["nuclear basis"]) {
      auto tp = Utils::ElementInfo::elementTypeForSymbol(it.first.as<std::string>());
      basisSetNameMap[tp] = it.second.as<std::string>();
    }
    if (!electronBasisSetNameMap.empty()) {
      molecule_ =
          std::make_unique<Kiwi::Molecule>(electronBasisSetNameMap, basisSetNameMap, structureBoolVectorPair, molSettings);
    }
    else {
      molecule_ = std::make_unique<Kiwi::Molecule>(basisSetNameMap, structureBoolVectorPair, molSettings);
    }
    if (input_["xyz nuclear centers"]) {
      const std::string xyzAdditionalNuclearCenters = input_["xyz nuclear centers"].as<std::string>();
      std::pair<Utils::AtomCollection, std::vector<bool>> structureBoolVectorPairAdditionalCenters;

      if (!boost::filesystem::exists(xyzAdditionalNuclearCenters)) {
        throw std::runtime_error("Specified xyz file with name '" + xyzAdditionalNuclearCenters + "'  does not exist!");
      }
      std::ifstream xyzAdditionalNuclearCentersFile(xyzAdditionalNuclearCenters);
      if (xyzAdditionalNuclearCentersFile.is_open()) {
        structureBoolVectorPairAdditionalCenters =
            Utils::XyzStreamHandler::readNuclearElectronic(xyzAdditionalNuclearCentersFile);
        xyzAdditionalNuclearCentersFile.close();
      }
      else {
        throw std::runtime_error("Unable to open xyz file.");
      }

      molecule_->addAdditionalNuclearCenters(structureBoolVectorPairAdditionalCenters.first, basisSetNameMap);
    }
  }

  // Un-contraction routine
  if (uncontract) {
    molecule_->uncontractElectronicBasis();
  }

  std::string baseName = inputFileName_;
  std::string extension = ".yml";

  std::size_t ind = baseName.find(extension); // Find the starting position of substring in the string
  if (ind != std::string::npos) {
    baseName.erase(ind, extension.length()); // erase function takes two parameter, the starting index in the string
                                             // from where you want to erase characters and total no of characters you
                                             // want to erase.
  }

  molecule_->setInputFileName(baseName);

  //
  // The molecule_ object has been created. We can proceed with the tasks now.
  //
}

void CalculationManager::executeTasks() {
  std::cout << ".==============================================================================." << std::endl;
  std::cout << "|                             Running Tasks                                    |" << std::endl;
  std::cout << "'=============================================================================='" << std::endl;

  auto& clock = Clock::getInstance();

  for (size_t i = 0; i < tasks_.size(); i++) {
    clock.time("task" + std::to_string(i));

    if (i != 0) {
      tasks_[i]->setData(tasks_[i - 1]->getData());
    }
    auto name = tasks_[i]->name();
    printf("  .===============%s==.\n", std::string(name.size(), '=').c_str());
    printf("  |  Task %3d     %s  |\n", static_cast<int>(i + 1), name.c_str());
    printf("  '===============%s=='\n", std::string(name.size(), '=').c_str());
    tasks_[i]->run();
    std::cout << "\nTask " + std::to_string(i + 1) + " took ";
    clock.time("task" + std::to_string(i));
    std::cout << "\n\n";
  }
}

auto CalculationManager::initializeTasks() -> void {
  auto tasksInput = input_["tasks"];

  tasks_.reserve(tasksInput.size());

  for (auto&& current : tasksInput) {
    std::string type = current["type"].as<std::string>();

    std::unique_ptr<Task> task;
    // Generate task
    task = TaskFactory::produce(type, current);

    task->setMolecule(molecule_);

    tasks_.emplace_back(std::move(task));
  }
}

CalculationManager::~CalculationManager() = default;
