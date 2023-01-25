/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
/* internal */
#include "CalculationManager.h"
#include <Kiwi/KiwiUtils/GeneralUtility.h>

/* Boost program arguments */
#include "boost/program_options.hpp"
/* External dependencies */
#include <yaml-cpp/yaml.h>
#include <boost/filesystem/operations.hpp>
#include <iostream>

#define VERSION_NUMBER "1.0.0"

using namespace Scine;

int main(int argc, char* argv[]) {
  auto& clock = Kiwi::Clock::getInstance();
  clock.time("main");

  namespace po = boost::program_options;

  // Arguments
  po::options_description optionsDescription("Recognized options");
  optionsDescription.add_options()("help,h", "Produce this help message")("config,c", po::value<std::string>(),
                                                                          "YAML input file to read");

  po::variables_map optionsMap;
  po::positional_options_description positionalDescription;
  positionalDescription.add("config", 1);
  po::store(po::command_line_parser(argc, argv)
                .options(optionsDescription)
                .positional(positionalDescription)
                .style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise)
                .run(),
            optionsMap);
  po::notify(optionsMap);

  // Handle help
  if (optionsMap.count("help") > 0 || optionsMap.count("config") == 0) {
    std::cout << optionsDescription << std::endl;
    return 1;
  }
  // Load input file
  const std::string filename = optionsMap["config"].as<std::string>();
  if (!boost::filesystem::exists(filename)) {
    std::cout << "Specified config file does not exist!\n";
    return 1;
  }
  auto input = YAML::LoadFile(filename);
  if (input.size() == 0) {
    std::cout << "Input file is empty." << std::endl;
    return 1;
  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << R"(.==============================================================================.)" << std::endl;
  std::cout << R"(|          .-------------------------------------------------------.           |)" << std::endl;
  std::cout << R"(|          |  ___  ____       _____     _____  _____      _____    |           |)" << std::endl;
  std::cout << R"(|          | |_  ||_  _|     |_   _|   |_   _||_   _|    |_   _|   |           |)" << std::endl;
  std::cout << R"(|          |   | |_/ /         | |       | | /\ | |        | |     |           |)" << std::endl;
  std::cout << R"(|          |   |  __'.         | |       | |/  \| |        | |     |           |)" << std::endl;
  std::cout << R"(|          |  _| |  \ \_      _| |_      |   /\   |       _| |_    |           |)" << std::endl;
  std::cout << R"(|          | |____||____|    |_____|     |__/  \__|      |_____|   |           |)" << std::endl;
  std::cout << R"(|          |                                                       |           |)" << std::endl;
  std::cout << R"(|          '-------------------------------------------------------'           |)" << std::endl;
  std::cout << R"('==============================================================================')" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "Version " << VERSION_NUMBER << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "With contributions from:" << std::endl;
  std::cout << "\tRobin Feldmann" << std::endl;
  std::cout << "\tFrancesco Bosia" << std::endl;
  std::cout << "\tMarkus Reiher" << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "--------------------------------------------------------------------------------" << std::endl;
  std::cout << "                              Input file" << std::endl;
  std::cout << "--------------------------------------------------------------------------------" << std::endl;
  std::cout << input << std::endl;
  std::cout << "--------------------------------------------------------------------------------" << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;

  Kiwi::CalculationManager manager(input, filename);

  manager.initializeTasks();

  manager.executeTasks();

  std::cout << "               *** Successfully terminated program. ***" << std::endl;
  std::cout << "Total run time:\t";
  clock.time("main");
  return 0;
}
