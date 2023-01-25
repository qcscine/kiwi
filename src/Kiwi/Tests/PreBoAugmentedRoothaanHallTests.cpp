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

class PreBoArhTest : public Test {};

//
// TODO Symmetry test of the Hessian and eigenvalue test.
//

TEST_F(PreBoArhTest, Test_H2_RHF_HighSpin_Arh) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoArhTest, Test_H2_RHF_LowSpin_Arh) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoArhTest, Test_H2_UHF_HighSpin_Arh) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoArhTest, Test_H2_UHF_LowSpin_Arh) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoArhTest, Test_H_Arh) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-0.490849245216, 1e-7));
}

TEST_F(PreBoArhTest, Test_H2_RHF_HighSpin_Newton) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoArhTest, Test_H2_RHF_LowSpin_Newton) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoArhTest, Test_H2_UHF_HighSpin_Newton) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoArhTest, Test_H2_UHF_LowSpin_Newton) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(PreBoArhTest, Test_H_Newton) {
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
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-0.490849245216, 1e-7));
}

} // namespace Scine
