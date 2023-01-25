/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/HartreeFock/InitialGuess/Read.h>
#include <Kiwi/KiwiUtils/AO2MO.h>
#include <Kiwi/KiwiUtils/Loewdin.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>

namespace Scine {

using namespace testing;

class IOTest : public Test {};

// TODO: Test density evaluation with and without natural orbitals, and with and without external RDM.

TEST_F(IOTest, Test_MolecularOrbitalsIO_HCN_RHF) {
  std::stringstream xyzInput("3\n\n"
                             "H    0.0     0.0 0.0\n"
                             "C    1.0378  0.0 0.0\n"
                             "N    2.204   0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-svp", atoms, moleculeSettings);
  mol->setInputFileName("HCN");

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
  settings.guess = Kiwi::InitialGuessSCF::SAD;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;
  {
    Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

    hartreeFockMain.scf();

    hartreeFockMain.writeOrbitals();
  }

  settings.guess = Kiwi::InitialGuessSCF::Read;
  {
    Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

    hartreeFockMain.scf();
  }

  Kiwi::Read densityMatrixIo(data, mol);

  densityMatrixIo.cleanUp();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-92.797276770625, 1e-8));
}

TEST_F(IOTest, Test_MolecularOrbitalsIO_H2_RHF_HighSpin) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0 Q\n"
                             "H    0.0 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-svp";
  basisSetsNames[Utils::ElementType::H] = "pb4-d";

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);
  mol->setInputFileName("H2");

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
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  settings.maxIterations = 500;
  {
    Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
    auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
    trahSettings->initialTrustRadius = 0.1;
    trahSettings->maxTrustRadius = 0.2;
    hartreeFockMain.setTrahSettings(trahSettings);

    hartreeFockMain.scf();

    hartreeFockMain.writeOrbitals();
  }

  settings.guess = Kiwi::InitialGuessSCF::Read;
  {
    Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
    auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
    trahSettings->initialTrustRadius = 0.1;
    trahSettings->maxTrustRadius = 0.2;
    hartreeFockMain.setTrahSettings(trahSettings);

    hartreeFockMain.scf();
  }

  Kiwi::Read densityMatrixIo(data, mol);
  densityMatrixIo.cleanUp();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(IOTest, Test_MolecularOrbitalsIO_H2_RHF_LowSpin) {
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
  mol->setInputFileName("H2");

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
  settings.neScfType = Kiwi::NuclearElectronicSCF::Serial;
  settings.maxIterations = 500;
  {
    Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
    auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
    trahSettings->initialTrustRadius = 0.1;
    trahSettings->maxTrustRadius = 0.2;
    hartreeFockMain.setTrahSettings(trahSettings);

    hartreeFockMain.scf();

    hartreeFockMain.writeOrbitals();
  }

  settings.guess = Kiwi::InitialGuessSCF::Read;
  {
    Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
    auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
    trahSettings->initialTrustRadius = 0.1;
    trahSettings->maxTrustRadius = 0.2;
    hartreeFockMain.setTrahSettings(trahSettings);

    hartreeFockMain.scf();
  }

  Kiwi::Read densityMatrixIo(data, mol);
  densityMatrixIo.cleanUp();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  EXPECT_THAT(energy, DoubleNear(-1.0724750406, 1e-7));
}

TEST_F(IOTest, Test_FCIDUMP) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0\n"
                             "H    0.0 0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);
  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "def2-tzvp";

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);
  mol->setInputFileName("H2");

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
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  trahSettings->initialTrustRadius = 0.1;
  trahSettings->maxTrustRadius = 0.2;
  hartreeFockMain.setTrahSettings(trahSettings);

  hartreeFockMain.scf();

  Kiwi::AO2MO ao2mo(data, true);

  ao2mo.perform();

  ao2mo.write(1e-16);

  const auto& H_ref = ao2mo.getOneBody()->at(Utils::ElementType::E);
  const auto& V_ref = ao2mo.getTwoBody()->at({Utils::ElementType::E, Utils::ElementType::E});

  Kiwi::BoFciDumpData boFciDumpData;
  int L = mol->at(Utils::ElementType::E).LMO;

  std::string fname = "FCIDUMP.H2";
  std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Kiwi::BoFciDumpData> Core_Eri_data = Kiwi::BoFciDumper::read(fname);

  auto const& H = get<0>(Core_Eri_data);
  auto const& V = get<1>(Core_Eri_data);

  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
      EXPECT_THAT(H(i, j), DoubleNear(H_ref(i, j), 1e-10));
    }
  }
  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
      for (int k = 0; k < L; ++k) {
        for (int l = 0; l < L; ++l) {
          EXPECT_THAT(V(i * L + j, k * L + l), DoubleNear(V_ref(i * L + j, k * L + l), 1e-10));
        }
      }
    }
  }

  std::remove("FCIDUMP.H2");
}

} // namespace Scine
