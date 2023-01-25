/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/KiwiUtils/Loewdin.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>

namespace Scine {

using namespace testing;

class PreBoRoothaanHallTest : public Test {};

TEST_F(PreBoRoothaanHallTest, Test_H2_RHF_HighSpin_Alternating) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Alternating;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724749729, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_H2_RHF_LowSpin_Alternating) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.useHighSpinApproximation = false;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Alternating;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_H2_UHF_HighSpin_Alternating) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.useHighSpinApproximation = true;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Alternating;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_H2_UHF_LowSpin_Alternating) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.useHighSpinApproximation = false;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Alternating;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_HCN_RHF_Alternating) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0 Q\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = true;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Alternating;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from blueberry
  EXPECT_THAT(energy, DoubleNear(-92.75837640221224, 5e-5));
}

TEST_F(PreBoRoothaanHallTest, Test_HCN_UHF_Alternating) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0 Q\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = true;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Alternating;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from blueberry
  EXPECT_THAT(energy, DoubleNear(-92.75837640221224, 5e-5));
}

TEST_F(PreBoRoothaanHallTest, Test_H_Alternating) {
  std::stringstream xyzInput("1\n\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  moleculeSettings.multiplicity = 2;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = true;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Alternating;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-0.490849245216, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_H2_RHF_HighSpin_Serial) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_H2_RHF_LowSpin_Serial) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.useHighSpinApproximation = false;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_H2_UHF_HighSpin_Serial) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.useHighSpinApproximation = true;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_H2_UHF_LowSpin_Serial) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.useHighSpinApproximation = false;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoRoothaanHallTest, Test_HCN_RHF_Serial) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0 Q\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = true;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from blueberry
  EXPECT_THAT(energy, DoubleNear(-92.75837640221224, 5e-5));
}

TEST_F(PreBoRoothaanHallTest, Test_HCN_UHF_Serial) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0 Q\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = true;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from blueberry
  EXPECT_THAT(energy, DoubleNear(-92.75837640221224, 5e-5));
}

TEST_F(PreBoRoothaanHallTest, Test_H_Serial) {
  std::stringstream xyzInput("1\n\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  moleculeSettings.multiplicity = 2;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = true;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::BO;
  settings.nuclearGuess = Kiwi::InitialGuessSCF::SND;
  settings.maxIterations = 500;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-0.490849245216, 1e-7));
}

} // namespace Scine
