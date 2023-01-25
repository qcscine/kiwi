/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/FockMatrixBuilder.h>
#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <Kiwi/KiwiUtils/Loewdin.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>

namespace Scine {

using namespace testing;

class BoDirectHFTest : public Test {};

TEST_F(BoDirectHFTest, Test_HCN_RHF) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::Hueckel;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-92.797276770625, 1e-8));
}

TEST_F(BoDirectHFTest, Test_H2_UHF_Triplet) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.0     0.0 0.0\n"
                             "H    0.7   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  moleculeSettings.multiplicity = 3;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::Core;
  settings.electronicDiisThresh = 0.2;
  // settings.useOrbitalSteering = true;
  // settings.numberOfPerturbations = 1;
  // settings.guess = Kiwi::InitialGuessSCF::SAD;
  settings.accelerator = Utils::scf_mixer_t::ediis_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-0.7420515415, 5e-7));
}

TEST_F(BoDirectHFTest, Test_H2plus_UHF) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.0     0.0 0.0\n"
                             "H    0.7   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  moleculeSettings.multiplicity = 2;
  moleculeSettings.charge = 1;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.electronicDiisThresh = 0.2;
  settings.guess = Kiwi::InitialGuessSCF::Core;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-0.551832212171, 5e-7));
}

TEST_F(BoDirectHFTest, Test_HCN_UHF_Doublet) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  moleculeSettings.multiplicity = 2;
  moleculeSettings.charge = 1;
  auto mol = std::make_shared<Kiwi::Molecule>("cc-pvdz", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.electronicDiisThresh = 0.2;
  settings.guess = Kiwi::InitialGuessSCF::Hueckel;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  // trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca
  EXPECT_THAT(energy, DoubleNear(-92.4251059200, 1e-7));
}

TEST_F(BoDirectHFTest, Test_Benzene_RHF_SAD) {
  std::stringstream xyzInput("12\n"
                             "benzene example\n"
                             "  C        0.00000        1.40272        0.00000\n"
                             "  H        0.00000        2.49029        0.00000\n"
                             "  C       -1.21479        0.70136        0.00000\n"
                             "  H       -2.15666        1.24515        0.00000\n"
                             "  C       -1.21479       -0.70136        0.00000\n"
                             "  H       -2.15666       -1.24515        0.00000\n"
                             "  C        0.00000       -1.40272        0.00000\n"
                             "  H        0.00000       -2.49029        0.00000\n"
                             "  C        1.21479       -0.70136        0.00000\n"
                             "  H        2.15666       -1.24515        0.00000\n"
                             "  C        1.21479        0.70136        0.00000\n"
                             "  H        2.15666        1.24515        0.00000");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>("6-31g", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.electronicDiisThresh = 1.0;
  settings.guess = Kiwi::InitialGuessSCF::SAD;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-230.622311571224, 1e-7));
}

TEST_F(BoDirectHFTest, Test_HCN_RHF_Newton) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::Hueckel;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-92.797276770625, 1e-8));
}

TEST_F(BoDirectHFTest, Test_H2_UHF_Triplet_Newton) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.0     0.0 0.0\n"
                             "H    0.7   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  moleculeSettings.multiplicity = 3;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::Core;
  settings.electronicDiisThresh = 0.2;
  // settings.useOrbitalSteering = true;
  // settings.numberOfPerturbations = 1;
  // settings.guess = Kiwi::InitialGuessSCF::SAD;
  settings.accelerator = Utils::scf_mixer_t::ediis_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-0.7420515415, 5e-7));
}

TEST_F(BoDirectHFTest, Test_H2plus_UHF_Newton) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.0     0.0 0.0\n"
                             "H    0.7   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  moleculeSettings.multiplicity = 2;
  moleculeSettings.charge = 1;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.electronicDiisThresh = 0.2;
  settings.guess = Kiwi::InitialGuessSCF::Core;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-0.551832212171, 5e-7));
}

TEST_F(BoDirectHFTest, Test_HCN_UHF_Doublet_Newton) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  moleculeSettings.multiplicity = 2;
  moleculeSettings.charge = 1;
  auto mol = std::make_shared<Kiwi::Molecule>("cc-pvdz", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.electronicDiisThresh = 0.2;
  settings.guess = Kiwi::InitialGuessSCF::Hueckel;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca
  EXPECT_THAT(energy, DoubleNear(-92.4251059200, 1e-7));
}

TEST_F(BoDirectHFTest, Test_Benzene_RHF_SAD_Newton) {
  std::stringstream xyzInput("12\n"
                             "benzene example\n"
                             "  C        0.00000        1.40272        0.00000\n"
                             "  H        0.00000        2.49029        0.00000\n"
                             "  C       -1.21479        0.70136        0.00000\n"
                             "  H       -2.15666        1.24515        0.00000\n"
                             "  C       -1.21479       -0.70136        0.00000\n"
                             "  H       -2.15666       -1.24515        0.00000\n"
                             "  C        0.00000       -1.40272        0.00000\n"
                             "  H        0.00000       -2.49029        0.00000\n"
                             "  C        1.21479       -0.70136        0.00000\n"
                             "  H        2.15666       -1.24515        0.00000\n"
                             "  C        1.21479        0.70136        0.00000\n"
                             "  H        2.15666        1.24515        0.00000");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>("6-31g", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.electronicDiisThresh = 1.0;
  settings.guess = Kiwi::InitialGuessSCF::SAD;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;

  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  trahSettings->opt = Kiwi::TRAHSettings::TrahOptimizer::Newton;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-230.622311571224, 1e-7));
}

TEST_F(BoDirectHFTest, Test_HCN_RHF_Incremental_RH) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::Hueckel;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  settings.disableIncrementalFock = false;
  settings.resetIncremental = 5;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-92.797276770625, 1e-8));
}

TEST_F(BoDirectHFTest, Test_HCN_UHF_Incremental_RH) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->oneBodyIntegrals();
  data->twoBodyIntegrals();
  data->makeExchange();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }

  Kiwi::HartreeFockSettings settings;
  settings.guess = Kiwi::InitialGuessSCF::Hueckel;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  settings.disableIncrementalFock = false;
  settings.resetIncremental = 5;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-92.797276770625, 1e-8));
}

} // namespace Scine
