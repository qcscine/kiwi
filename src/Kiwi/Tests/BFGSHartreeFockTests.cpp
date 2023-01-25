/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/BFGS/Optimizer.h>
#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/KiwiOpt/Optimization.h>
#include <Kiwi/KiwiUtils/Loewdin.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>

namespace Scine {

using namespace testing;

class BFGSTest : public Test {};

TEST_F(BFGSTest, Test_Evaluate_Function) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0\n"
                             "H    0.0 0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = true;
  moleculeSettings.multiplicity = 1;
  moleculeSettings.charge = 0;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-tzvp", atoms, moleculeSettings);

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
  settings.disableIncrementalFock = false;
  settings.neScfType = Kiwi::NuclearElectronicSCF::TRAH;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  // hartreeFockMain.scf();

  auto interface = std::make_shared<Kiwi::BFGS::Interface<Kiwi::SymmetryType::None>>(mol, data);

  auto opt = std::make_shared<Kiwi::BFGS::Optimizer<Kiwi::SymmetryType::None>>(interface);

  auto energy_before = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  auto alpha = double(1);

  opt->evaluate(alpha);

  auto energy_after = opt->getLineSearchValue();

  ASSERT_TRUE(energy_after < energy_before);
}

TEST_F(BFGSTest, Test_Apply_Update) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0\n"
                             "H    0.0 0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = true;
  moleculeSettings.multiplicity = 1;
  moleculeSettings.charge = 0;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-tzvp", atoms, moleculeSettings);

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
  settings.disableIncrementalFock = false;
  settings.neScfType = Kiwi::NuclearElectronicSCF::TRAH;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  hartreeFockMain.singleIteration();

  auto interface = std::make_shared<Kiwi::BFGS::Interface<Kiwi::SymmetryType::None>>(mol, data);

  auto opt = std::make_shared<Kiwi::BFGS::Optimizer<Kiwi::SymmetryType::None>>(interface);

  auto alpha = double(1);

  auto gradientOld = interface->getGradient(0);
  auto directionOld = interface->getGradient(0);

  opt->evaluate(alpha);

  opt->applyUpdate();

  auto gradientNew = interface->getGradient(0);
  auto directionNew = interface->getGradient(0);

  auto gradient_diff = (gradientOld - gradientNew).norm();
  std::cout << "gradient diff = " << gradient_diff << std::endl;
  auto direction_diff = (directionOld - directionNew).norm();
  std::cout << "direction diff = " << direction_diff << std::endl;

  ASSERT_TRUE(gradient_diff > 1e-12);
  ASSERT_TRUE(direction_diff > 1e-12);
}

TEST_F(BFGSTest, Test_H2_RHF) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.7 0.0 0.0\n"
                             "H    0.0 0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = true;
  moleculeSettings.multiplicity = 1;
  moleculeSettings.charge = 0;
  auto mol = std::make_shared<Kiwi::Molecule>("def2-tzvp", atoms, moleculeSettings);

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
  settings.guess = Kiwi::InitialGuessSCF::Hueckel;
  settings.disableIncrementalFock = false;
  settings.neScfType = Kiwi::NuclearElectronicSCF::TRAH;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);
  hartreeFockMain.singleIteration();

  auto interface = std::make_shared<Kiwi::BFGS::Interface<Kiwi::SymmetryType::None>>(mol, data);

  auto opt = std::make_shared<Kiwi::BFGS::Optimizer<Kiwi::SymmetryType::None>>(interface);

  Kiwi::Optimization::Optimization optimization(opt);
  optimization.setThresh(settings.fockMatrixNormThresh);

  optimization.optimize();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-1.131600453888, 1e-8));
}

TEST_F(BFGSTest, Test_Benzene_RHF) {
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
  moleculeSettings.isRestricted = true;
  moleculeSettings.multiplicity = 1;
  moleculeSettings.charge = 0;
  auto mol = std::make_shared<Kiwi::Molecule>("6-31g", atoms, moleculeSettings);

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
  settings.electronicDiisThresh = 1.0;
  settings.guess = Kiwi::InitialGuessSCF::Hueckel;
  settings.accelerator = Utils::scf_mixer_t::fock_diis;
  settings.neScfType = Kiwi::NuclearElectronicSCF::BFGS;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings, false, true);

  hartreeFockMain.scf();

  auto energy = Kiwi::HartreeFockUtils::getEnergy(mol, data);

  // Reference energy from orca:
  EXPECT_THAT(energy, DoubleNear(-230.622311571224, 1e-7));
}

} // namespace Scine
