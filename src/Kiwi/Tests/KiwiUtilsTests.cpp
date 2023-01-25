/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/HartreeFockMain.h>
#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/KiwiOpt/Bisection.h>
#include <Kiwi/KiwiOpt/Davidson.h>
#include <Kiwi/KiwiOpt/DirectConjugateGradient.h>
#include <Kiwi/KiwiOpt/MoreThuente.h>
#include <Kiwi/KiwiOpt/Optimization.h>
#include <Kiwi/KiwiOpt/PreconditionedDirectConjugateGradient.h>
#include <Kiwi/KiwiUtils/AO2MO.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <Kiwi/KiwiUtils/Loewdin.h>
#include <Kiwi/KiwiUtils/Molecule.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <utility>

namespace Scine {

using namespace testing;

class MoleculeTest : public Test {};

class CalculationDataTest : public Test {};

class OptimizationTest : public Test {};

class LoewdinTest : public Test {};

TEST_F(LoewdinTest, SingleStepLoewdinTest) {
  std::stringstream xyzInput("2\n\n"
                             "H    0.025 0.0 0.0 Q\n"
                             "H    -0.025 0.0 0.0 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  std::unordered_map<Utils::ElementType, std::string> basisSetsNames;
  basisSetsNames[Utils::ElementType::E] = "cc-pvtz";
  basisSetsNames[Utils::ElementType::H] = "def2-svp";

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.isRestricted = false;
  auto mol_unrest = std::make_shared<Kiwi::Molecule>(basisSetsNames, atoms, moleculeSettings);
  auto data_unrest = std::make_shared<Kiwi::Data>(mol_unrest);

  data_unrest->oneBodyIntegrals();

  for (auto const& elem : *mol_unrest) {
    data_unrest->X[elem.first] = Kiwi::SingleStepLoewdin(data_unrest->S[elem.first], 1e-6, true);
    const auto& X = data_unrest->X[elem.first];
    const auto& S = data_unrest->S[elem.first];
    Eigen::MatrixXd testMat = X.transpose() * S.selfadjointView<Eigen::Lower>() * X;
    EXPECT_EQ(testMat.cols(), testMat.rows());
    // Generate testMat matrix
    for (auto i = 0L; i < testMat.cols(); i++) {
      // Diminish by one -> Matrix should be zero.
      testMat(i, i) -= 1.0;
    }

    auto loewdinError = testMat.norm(); // Search error in transformation of the overlap
    EXPECT_THAT(loewdinError, DoubleNear(0.0, 1e-9));
  }
}

TEST_F(MoleculeTest, constructorTestDefault) {
  std::stringstream xyzInput("3\n\n"
                             "C12    -4.0410150   -1.2118929   -0.0394793 Q\n"
                             "H    -4.1063106   -1.6202046   -1.0499636 Q\n"
                             "H    -0.8426606   -1.1792649    0.1904767 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  std::unordered_map<Utils::ElementType, std::string> nameMap = {
      {Utils::ElementType::E, "sto-3g"}, {Utils::ElementType::C12, "sto-3g"}, {Utils::ElementType::H, "sto-3g"}};

  Kiwi::MoleculeSettings moleculeSettings;
  Kiwi::Molecule mol(nameMap, atoms, moleculeSettings);

  EXPECT_EQ(mol.at(Utils::ElementType::E).N, 8);
  EXPECT_EQ(mol.at(Utils::ElementType::E).msVector[0], 4);
  EXPECT_EQ(mol.at(Utils::ElementType::E).msVector[1], 4);

  EXPECT_EQ(mol.at(Utils::ElementType::C12).N, 1);
  EXPECT_EQ(mol.at(Utils::ElementType::C12).msVector[0], 1);

  EXPECT_EQ(mol.at(Utils::ElementType::H).N, 2);
  EXPECT_EQ(mol.at(Utils::ElementType::H).msVector[0], 2);
  EXPECT_EQ(mol.at(Utils::ElementType::H).msVector[1], 0);
}

TEST_F(MoleculeTest, constructorTestOpenShell) {
  std::stringstream xyzInput("3\n\n"
                             "C12    -4.0410150   -1.2118929   -0.0394793 Q\n"
                             "H    -4.1063106   -1.6202046   -1.0499636 Q\n"
                             "H    -0.8426606   -1.1792649    0.1904767 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  std::unordered_map<Utils::ElementType, std::string> nameMap = {
      {Utils::ElementType::E, "sto-3g"}, {Utils::ElementType::C12, "sto-3g"}, {Utils::ElementType::H, "sto-3g"}};

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.multiplicity = 3;
  moleculeSettings.isRestricted = false;
  Kiwi::Molecule mol(nameMap, atoms, moleculeSettings);

  EXPECT_EQ(mol.at(Utils::ElementType::E).N, 8);
  EXPECT_EQ(mol.at(Utils::ElementType::E).msVector[0], 5);
  EXPECT_EQ(mol.at(Utils::ElementType::E).msVector[1], 3);

  EXPECT_EQ(mol.at(Utils::ElementType::C12).N, 1);
  EXPECT_EQ(mol.at(Utils::ElementType::C12).msVector[0], 1);

  EXPECT_EQ(mol.at(Utils::ElementType::H).N, 2);
  EXPECT_EQ(mol.at(Utils::ElementType::H).msVector[0], 2);
  EXPECT_EQ(mol.at(Utils::ElementType::H).msVector[1], 0);
}

TEST_F(MoleculeTest, constructorTestNuclearClosedShell) {
  std::stringstream xyzInput("3\n\n"
                             "C12    -4.0410150   -1.2118929   -0.0394793 Q\n"
                             "H    -4.1063106   -1.6202046   -1.0499636 Q\n"
                             "H    -0.8426606   -1.1792649    0.1904767 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  std::unordered_map<Utils::ElementType, std::string> nameMap = {
      {Utils::ElementType::E, "sto-3g"}, {Utils::ElementType::C12, "sto-3g"}, {Utils::ElementType::H, "sto-3g"}};

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.useHighSpinApproximation = false;
  Kiwi::Molecule mol(nameMap, atoms, moleculeSettings);

  EXPECT_EQ(mol.at(Utils::ElementType::E).N, 8);
  EXPECT_EQ(mol.at(Utils::ElementType::E).msVector[0], 4);
  EXPECT_EQ(mol.at(Utils::ElementType::E).msVector[1], 4);

  EXPECT_EQ(mol.at(Utils::ElementType::C12).N, 1);
  EXPECT_EQ(mol.at(Utils::ElementType::C12).msVector[0], 1);

  EXPECT_EQ(mol.at(Utils::ElementType::H).N, 2);
  EXPECT_EQ(mol.at(Utils::ElementType::H).msVector[0], 1);
  EXPECT_EQ(mol.at(Utils::ElementType::H).msVector[1], 1);
}

TEST_F(MoleculeTest, constructorTestCharged) {
  std::stringstream xyzInput("3\n\n"
                             "C12    -4.0410150   -1.2118929   -0.0394793 Q\n"
                             "H    -4.1063106   -1.6202046   -1.0499636 Q\n"
                             "H    -0.8426606   -1.1792649    0.1904767 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  std::unordered_map<Utils::ElementType, std::string> nameMap = {
      {Utils::ElementType::E, "sto-3g"}, {Utils::ElementType::C12, "sto-3g"}, {Utils::ElementType::H, "sto-3g"}};

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.multiplicity = 2;
  moleculeSettings.charge = -1;
  moleculeSettings.isRestricted = false;
  Kiwi::Molecule mol(nameMap, atoms, moleculeSettings);

  EXPECT_EQ(mol.at(Utils::ElementType::E).N, 9);
  EXPECT_EQ(mol.at(Utils::ElementType::E).msVector[0], 5);
  EXPECT_EQ(mol.at(Utils::ElementType::E).msVector[1], 4);
}

TEST_F(CalculationDataTest, TestsMap) {
  auto ptr_data = std::make_shared<Kiwi::Data>();
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::E, Utils::ElementType::E)] = Kiwi::Matrix::Zero(1, 1);
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::E, Utils::ElementType::H)] = Kiwi::Matrix::Zero(1, 2);
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::H, Utils::ElementType::E)] = Kiwi::Matrix::Zero(1, 3);
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::H, Utils::ElementType::H)] = Kiwi::Matrix::Zero(1, 4);
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::C12, Utils::ElementType::C12)] = Kiwi::Matrix::Zero(1, 5);
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::E, Utils::ElementType::C12)] = Kiwi::Matrix::Zero(1, 6);
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::C12, Utils::ElementType::E)] = Kiwi::Matrix::Zero(1, 7);
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::H, Utils::ElementType::C12)] = Kiwi::Matrix::Zero(1, 8);
  ptr_data->Coulomb[Kiwi::Data::getUnique(Utils::ElementType::C12, Utils::ElementType::H)] = Kiwi::Matrix::Zero(1, 9);

  std::vector<Utils::ElementType> list = {Utils::ElementType::E, Utils::ElementType::H, Utils::ElementType::C12};

  for (auto const& first : list) {
    for (auto const& second : list) {
      int left = ptr_data->Coulomb[Kiwi::Data::getUnique(first, second)].cols();
      int right = ptr_data->Coulomb[Kiwi::Data::getUnique(second, first)].cols();
      EXPECT_EQ(left, right);
    }
  }
}

TEST_F(CalculationDataTest, TestsConstructor) {
  std::stringstream xyzInput("3\n\n"
                             "C12    -4.0410150   -1.2118929   -0.0394793 Q\n"
                             "H    -4.1063106   -1.6202046   -1.0499636 Q\n"
                             "H    -0.8426606   -1.1792649    0.1904767 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  std::unordered_map<Utils::ElementType, std::string> nameMap = {
      {Utils::ElementType::E, "sto-3g"}, {Utils::ElementType::C12, "sto-3g"}, {Utils::ElementType::H, "sto-3g"}};

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.multiplicity = 2;
  moleculeSettings.charge = -1;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>(nameMap, atoms, moleculeSettings);

  auto data = Kiwi::Data(mol);

  EXPECT_EQ(data.uniquePairs.size(), 6);
}

TEST_F(CalculationDataTest, CanCalculateIntegrals) {
  std::stringstream xyzInput("3\n\n"
                             "C12    -4.0410150   -1.2118929   -0.0394793 Q\n"
                             "H    -4.1063106   -1.6202046   -1.0499636 Q\n"
                             "H    -0.8426606   -1.1792649    0.1904767 Q");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  std::unordered_map<Utils::ElementType, std::string> nameMap = {
      {Utils::ElementType::E, "def2-svp"}, {Utils::ElementType::C12, "sto-3g"}, {Utils::ElementType::H, "sto-3g"}};

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.multiplicity = 2;
  moleculeSettings.charge = -1;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>(nameMap, atoms, moleculeSettings);

  auto data = Kiwi::Data(mol);

  data.oneBodyIntegrals();
  data.twoBodyIntegrals();
  data.makeExchange();
}

TEST_F(CalculationDataTest, AO2MOTest) {
  std::map<std::array<int, 4>, double> pyScfMOEriMap2;

  /** pyscf source code:
   *
   * from pyscf import gto, scf, ao2mo
   * import numpy as np
   * HH = 'H 0 0 0; H 1.2 0 0'
   *
   * mol = gto.M(atom=HH, basis='6-31g*')
   * mf = scf.RHF(mol)
   * mf.init_guess = 'hcore'
   * mf.diis = False
   *
   * mf.kernel()
   *
   * print(mf.get_fock())
   * print(mf.energy_elec())
   * print(mf.make_rdm1())
   * mo_ints = ao2mo.kernel(mol, mf.mo_coeff, compact=False)
   * print(mo_ints)
   *
   * dim=mf.mo_coeff.shape[0]
   *
   * for i in range(0,dim):
   *     for j in range(0,dim):
   *         for k in range(0,dim):
   *             for l in range(0,dim):
   *                 index1 = dim*i+j
   *                 index2 = dim*k+l
   *                 if (mo_ints[index1,index2] != 0):
   *                     print("pyScfMOEriMap2.insert({{{", i, ",", j, ",", k,
   *                         ",", l, "}}, ", mo_ints[index1,index2],
   *                         "});")
   *
   * for i in range(0,dim):
   *     for j in range(0,dim):
   *         print("%.16f ," % mf.mo_coeff[i,j] )
   *
   * np.set_printoptions(precision=16)
   * print(mol.atom_coords())
   *
   */

  pyScfMOEriMap2.insert({{{0, 0, 0, 0}}, 0.5241516710740569});
  pyScfMOEriMap2.insert({{{0, 0, 0, 1}}, 5.129044490757549e-17});
  pyScfMOEriMap2.insert({{{0, 0, 0, 2}}, 0.11421947189024337});
  pyScfMOEriMap2.insert({{{0, 0, 0, 3}}, -6.152426396084604e-16});
  pyScfMOEriMap2.insert({{{0, 0, 1, 0}}, 7.190029527983987e-17});
  pyScfMOEriMap2.insert({{{0, 0, 1, 1}}, 0.43928817323738495});
  pyScfMOEriMap2.insert({{{0, 0, 1, 2}}, -5.781139048104144e-16});
  pyScfMOEriMap2.insert({{{0, 0, 1, 3}}, -0.12802272602631076});
  pyScfMOEriMap2.insert({{{0, 0, 2, 0}}, 0.11421947189024341});
  pyScfMOEriMap2.insert({{{0, 0, 2, 1}}, -5.324464406350075e-16});
  pyScfMOEriMap2.insert({{{0, 0, 2, 2}}, 0.4867453413927506});
  pyScfMOEriMap2.insert({{{0, 0, 2, 3}}, -2.125548325963349e-16});
  pyScfMOEriMap2.insert({{{0, 0, 3, 0}}, -5.522988836325202e-16});
  pyScfMOEriMap2.insert({{{0, 0, 3, 1}}, -0.12802272602631073});
  pyScfMOEriMap2.insert({{{0, 0, 3, 2}}, -2.728481138705354e-16});
  pyScfMOEriMap2.insert({{{0, 0, 3, 3}}, 0.5188746251919503});
  pyScfMOEriMap2.insert({{{0, 1, 0, 0}}, 5.750894739670236e-17});
  pyScfMOEriMap2.insert({{{0, 1, 0, 1}}, 0.12298082684985662});
  pyScfMOEriMap2.insert({{{0, 1, 0, 2}}, -3.0301438461884693e-16});
  pyScfMOEriMap2.insert({{{0, 1, 0, 3}}, -0.07747835724448596});
  pyScfMOEriMap2.insert({{{0, 1, 1, 0}}, 0.12298082684985663});
  pyScfMOEriMap2.insert({{{0, 1, 1, 1}}, 2.2007188932312939e-16});
  pyScfMOEriMap2.insert({{{0, 1, 1, 2}}, 0.03266679248305437});
  pyScfMOEriMap2.insert({{{0, 1, 1, 3}}, -2.394064809792865e-16});
  pyScfMOEriMap2.insert({{{0, 1, 2, 0}}, -3.126936559011105e-16});
  pyScfMOEriMap2.insert({{{0, 1, 2, 1}}, 0.03266679248305437});
  pyScfMOEriMap2.insert({{{0, 1, 2, 2}}, -1.4359958592042849e-15});
  pyScfMOEriMap2.insert({{{0, 1, 2, 3}}, -0.13878148169674515});
  pyScfMOEriMap2.insert({{{0, 1, 3, 0}}, -0.07747835724448596});
  pyScfMOEriMap2.insert({{{0, 1, 3, 1}}, -2.247294912558349e-16});
  pyScfMOEriMap2.insert({{{0, 1, 3, 2}}, -0.13878148169674512});
  pyScfMOEriMap2.insert({{{0, 1, 3, 3}}, 1.8395876202956174e-15});
  pyScfMOEriMap2.insert({{{0, 2, 0, 0}}, 0.11421947189024342});
  pyScfMOEriMap2.insert({{{0, 2, 0, 1}}, -3.750476287542659e-16});
  pyScfMOEriMap2.insert({{{0, 2, 0, 2}}, 0.08454808151948119});
  pyScfMOEriMap2.insert({{{0, 2, 0, 3}}, -3.394102315216874e-17});
  pyScfMOEriMap2.insert({{{0, 2, 1, 0}}, -3.580955298938201e-16});
  pyScfMOEriMap2.insert({{{0, 2, 1, 1}}, 0.07663203050639272});
  pyScfMOEriMap2.insert({{{0, 2, 1, 2}}, -6.855109115151029e-16});
  pyScfMOEriMap2.insert({{{0, 2, 1, 3}}, -0.08038181988424592});
  pyScfMOEriMap2.insert({{{0, 2, 2, 0}}, 0.08454808151948118});
  pyScfMOEriMap2.insert({{{0, 2, 2, 1}}, -6.882532393661093e-16});
  pyScfMOEriMap2.insert({{{0, 2, 2, 2}}, 0.11755793467779554});
  pyScfMOEriMap2.insert({{{0, 2, 2, 3}}, 7.809884510997569e-16});
  pyScfMOEriMap2.insert({{{0, 2, 3, 0}}, -2.2663145039027866e-17});
  pyScfMOEriMap2.insert({{{0, 2, 3, 1}}, -0.08038181988424592});
  pyScfMOEriMap2.insert({{{0, 2, 3, 2}}, 7.70056540469179e-16});
  pyScfMOEriMap2.insert({{{0, 2, 3, 3}}, 0.14972481025996567});
  pyScfMOEriMap2.insert({{{0, 3, 0, 0}}, -6.10214475713126e-16});
  pyScfMOEriMap2.insert({{{0, 3, 0, 1}}, -0.07747835724448594});
  pyScfMOEriMap2.insert({{{0, 3, 0, 2}}, -3.6319484411394934e-17});
  pyScfMOEriMap2.insert({{{0, 3, 0, 3}}, 0.09471344004943412});
  pyScfMOEriMap2.insert({{{0, 3, 1, 0}}, -0.07747835724448593});
  pyScfMOEriMap2.insert({{{0, 3, 1, 1}}, -4.4352813927667497e-16});
  pyScfMOEriMap2.insert({{{0, 3, 1, 2}}, -0.06118960215791901});
  pyScfMOEriMap2.insert({{{0, 3, 1, 3}}, 1.0314674695705134e-15});
  pyScfMOEriMap2.insert({{{0, 3, 2, 0}}, -2.3006429080291718e-17});
  pyScfMOEriMap2.insert({{{0, 3, 2, 1}}, -0.061189602157919});
  pyScfMOEriMap2.insert({{{0, 3, 2, 2}}, 6.62785125997328e-16});
  pyScfMOEriMap2.insert({{{0, 3, 2, 3}}, 0.12014457557547449});
  pyScfMOEriMap2.insert({{{0, 3, 3, 0}}, 0.0947134400494341});
  pyScfMOEriMap2.insert({{{0, 3, 3, 1}}, 1.031577044033484e-15});
  pyScfMOEriMap2.insert({{{0, 3, 3, 2}}, 0.12014457557547445});
  pyScfMOEriMap2.insert({{{0, 3, 3, 3}}, -1.9239412992225977e-15});
  pyScfMOEriMap2.insert({{{1, 0, 0, 0}}, 7.305738330580873e-17});
  pyScfMOEriMap2.insert({{{1, 0, 0, 1}}, 0.12298082684985664});
  pyScfMOEriMap2.insert({{{1, 0, 0, 2}}, -3.50873330935985e-16});
  pyScfMOEriMap2.insert({{{1, 0, 0, 3}}, -0.07747835724448594});
  pyScfMOEriMap2.insert({{{1, 0, 1, 0}}, 0.12298082684985667});
  pyScfMOEriMap2.insert({{{1, 0, 1, 1}}, 1.6424952818057248e-16});
  pyScfMOEriMap2.insert({{{1, 0, 1, 2}}, 0.03266679248305439});
  pyScfMOEriMap2.insert({{{1, 0, 1, 3}}, -3.025177058810046e-16});
  pyScfMOEriMap2.insert({{{1, 0, 2, 0}}, -3.610096117297215e-16});
  pyScfMOEriMap2.insert({{{1, 0, 2, 1}}, 0.032666792483054374});
  pyScfMOEriMap2.insert({{{1, 0, 2, 2}}, -1.438849103245135e-15});
  pyScfMOEriMap2.insert({{{1, 0, 2, 3}}, -0.13878148169674503});
  pyScfMOEriMap2.insert({{{1, 0, 3, 0}}, -0.07747835724448594});
  pyScfMOEriMap2.insert({{{1, 0, 3, 1}}, -2.873351413466495e-16});
  pyScfMOEriMap2.insert({{{1, 0, 3, 2}}, -0.138781481696745});
  pyScfMOEriMap2.insert({{{1, 0, 3, 3}}, 1.6900708114783353e-15});
  pyScfMOEriMap2.insert({{{1, 1, 0, 0}}, 0.43928817323738506});
  pyScfMOEriMap2.insert({{{1, 1, 0, 1}}, 2.721134994789749e-16});
  pyScfMOEriMap2.insert({{{1, 1, 0, 2}}, 0.07663203050639279});
  pyScfMOEriMap2.insert({{{1, 1, 0, 3}}, -4.1519675163054267e-16});
  pyScfMOEriMap2.insert({{{1, 1, 1, 0}}, 1.9956017248831904e-16});
  pyScfMOEriMap2.insert({{{1, 1, 1, 1}}, 0.4096388175479753});
  pyScfMOEriMap2.insert({{{1, 1, 1, 2}}, -1.1295603843903856e-16});
  pyScfMOEriMap2.insert({{{1, 1, 1, 3}}, -0.08142205002071352});
  pyScfMOEriMap2.insert({{{1, 1, 2, 0}}, 0.07663203050639275});
  pyScfMOEriMap2.insert({{{1, 1, 2, 1}}, -1.2564053429859476e-16});
  pyScfMOEriMap2.insert({{{1, 1, 2, 2}}, 0.41635882599304774});
  pyScfMOEriMap2.insert({{{1, 1, 2, 3}}, -1.4745406145892995e-16});
  pyScfMOEriMap2.insert({{{1, 1, 3, 0}}, -4.704444476927934e-16});
  pyScfMOEriMap2.insert({{{1, 1, 3, 1}}, -0.08142205002071347});
  pyScfMOEriMap2.insert({{{1, 1, 3, 2}}, -1.0408824007153873e-16});
  pyScfMOEriMap2.insert({{{1, 1, 3, 3}}, 0.4455482953329862});
  pyScfMOEriMap2.insert({{{1, 2, 0, 0}}, -5.944279886205798e-16});
  pyScfMOEriMap2.insert({{{1, 2, 0, 1}}, 0.03266679248305431});
  pyScfMOEriMap2.insert({{{1, 2, 0, 2}}, -6.817668979694943e-16});
  pyScfMOEriMap2.insert({{{1, 2, 0, 3}}, -0.06118960215791902});
  pyScfMOEriMap2.insert({{{1, 2, 1, 0}}, 0.03266679248305431});
  pyScfMOEriMap2.insert({{{1, 2, 1, 1}}, -3.0691054227465586e-16});
  pyScfMOEriMap2.insert({{{1, 2, 1, 2}}, 0.04696781757802746});
  pyScfMOEriMap2.insert({{{1, 2, 1, 3}}, 1.9063297697549776e-16});
  pyScfMOEriMap2.insert({{{1, 2, 2, 0}}, -6.740039616587148e-16});
  pyScfMOEriMap2.insert({{{1, 2, 2, 1}}, 0.046967817578027454});
  pyScfMOEriMap2.insert({{{1, 2, 2, 2}}, -9.722012152976426e-16});
  pyScfMOEriMap2.insert({{{1, 2, 2, 3}}, -0.0642495577079315});
  pyScfMOEriMap2.insert({{{1, 2, 3, 0}}, -0.061189602157919035});
  pyScfMOEriMap2.insert({{{1, 2, 3, 1}}, 1.8654454227894157e-16});
  pyScfMOEriMap2.insert({{{1, 2, 3, 2}}, -0.0642495577079315});
  pyScfMOEriMap2.insert({{{1, 2, 3, 3}}, 4.718372790022788e-16});
  pyScfMOEriMap2.insert({{{1, 3, 0, 0}}, -0.12802272602631073});
  pyScfMOEriMap2.insert({{{1, 3, 0, 1}}, -2.7063228870451843e-16});
  pyScfMOEriMap2.insert({{{1, 3, 0, 2}}, -0.08038181988424582});
  pyScfMOEriMap2.insert({{{1, 3, 0, 3}}, 8.51321346936912e-16});
  pyScfMOEriMap2.insert({{{1, 3, 1, 0}}, -2.656657878323823e-16});
  pyScfMOEriMap2.insert({{{1, 3, 1, 1}}, -0.08142205002071341});
  pyScfMOEriMap2.insert({{{1, 3, 1, 2}}, 2.632036848200921e-16});
  pyScfMOEriMap2.insert({{{1, 3, 1, 3}}, 0.08637200526415965});
  pyScfMOEriMap2.insert({{{1, 3, 2, 0}}, -0.08038181988424582});
  pyScfMOEriMap2.insert({{{1, 3, 2, 1}}, 2.704096091540828e-16});
  pyScfMOEriMap2.insert({{{1, 3, 2, 2}}, -0.12767413616212886});
  pyScfMOEriMap2.insert({{{1, 3, 2, 3}}, 4.24425319675475e-16});
  pyScfMOEriMap2.insert({{{1, 3, 3, 0}}, 8.587102059679663e-16});
  pyScfMOEriMap2.insert({{{1, 3, 3, 1}}, 0.08637200526415965});
  pyScfMOEriMap2.insert({{{1, 3, 3, 2}}, 4.926430294172555e-16});
  pyScfMOEriMap2.insert({{{1, 3, 3, 3}}, -0.15246130026914106});
  pyScfMOEriMap2.insert({{{2, 0, 0, 0}}, 0.11421947189024338});
  pyScfMOEriMap2.insert({{{2, 0, 0, 1}}, -3.7489202379861636e-16});
  pyScfMOEriMap2.insert({{{2, 0, 0, 2}}, 0.08454808151948118});
  pyScfMOEriMap2.insert({{{2, 0, 0, 3}}, 1.1492377106504773e-17});
  pyScfMOEriMap2.insert({{{2, 0, 1, 0}}, -3.778846614308411e-16});
  pyScfMOEriMap2.insert({{{2, 0, 1, 1}}, 0.07663203050639276});
  pyScfMOEriMap2.insert({{{2, 0, 1, 2}}, -7.063761472763791e-16});
  pyScfMOEriMap2.insert({{{2, 0, 1, 3}}, -0.08038181988424573});
  pyScfMOEriMap2.insert({{{2, 0, 2, 0}}, 0.08454808151948116});
  pyScfMOEriMap2.insert({{{2, 0, 2, 1}}, -7.086809086887071e-16});
  pyScfMOEriMap2.insert({{{2, 0, 2, 2}}, 0.11755793467779584});
  pyScfMOEriMap2.insert({{{2, 0, 2, 3}}, 6.864051346292238e-16});
  pyScfMOEriMap2.insert({{{2, 0, 3, 0}}, 2.746646861375418e-18});
  pyScfMOEriMap2.insert({{{2, 0, 3, 1}}, -0.08038181988424575});
  pyScfMOEriMap2.insert({{{2, 0, 3, 2}}, 6.76997675677535e-16});
  pyScfMOEriMap2.insert({{{2, 0, 3, 3}}, 0.14972481025996606});
  pyScfMOEriMap2.insert({{{2, 1, 0, 0}}, -5.309652400645142e-16});
  pyScfMOEriMap2.insert({{{2, 1, 0, 1}}, 0.032666792483054326});
  pyScfMOEriMap2.insert({{{2, 1, 0, 2}}, -7.262618940129003e-16});
  pyScfMOEriMap2.insert({{{2, 1, 0, 3}}, -0.061189602157919035});
  pyScfMOEriMap2.insert({{{2, 1, 1, 0}}, 0.032666792483054326});
  pyScfMOEriMap2.insert({{{2, 1, 1, 1}}, -3.5926345646525915e-16});
  pyScfMOEriMap2.insert({{{2, 1, 1, 2}}, 0.04696781757802742});
  pyScfMOEriMap2.insert({{{2, 1, 1, 3}}, 5.466336608275374e-17});
  pyScfMOEriMap2.insert({{{2, 1, 2, 0}}, -7.364540067938799e-16});
  pyScfMOEriMap2.insert({{{2, 1, 2, 1}}, 0.046967817578027427});
  pyScfMOEriMap2.insert({{{2, 1, 2, 2}}, -1.1005707525199263e-15});
  pyScfMOEriMap2.insert({{{2, 1, 2, 3}}, -0.06424955770793161});
  pyScfMOEriMap2.insert({{{2, 1, 3, 0}}, -0.061189602157919035});
  pyScfMOEriMap2.insert({{{2, 1, 3, 1}}, 5.4122039521018457e-17});
  pyScfMOEriMap2.insert({{{2, 1, 3, 2}}, -0.06424955770793161});
  pyScfMOEriMap2.insert({{{2, 1, 3, 3}}, 1.3839802129043326e-16});
  pyScfMOEriMap2.insert({{{2, 2, 0, 0}}, 0.48674534139275055});
  pyScfMOEriMap2.insert({{{2, 2, 0, 1}}, -1.3529106966456838e-15});
  pyScfMOEriMap2.insert({{{2, 2, 0, 2}}, 0.11755793467779589});
  pyScfMOEriMap2.insert({{{2, 2, 0, 3}}, 6.235722586757672e-16});
  pyScfMOEriMap2.insert({{{2, 2, 1, 0}}, -1.4155705750459783e-15});
  pyScfMOEriMap2.insert({{{2, 2, 1, 1}}, 0.416358825993048});
  pyScfMOEriMap2.insert({{{2, 2, 1, 2}}, -1.1332955448776432e-15});
  pyScfMOEriMap2.insert({{{2, 2, 1, 3}}, -0.12767413616212892});
  pyScfMOEriMap2.insert({{{2, 2, 2, 0}}, 0.11755793467779581});
  pyScfMOEriMap2.insert({{{2, 2, 2, 1}}, -1.2024124504104474e-15});
  pyScfMOEriMap2.insert({{{2, 2, 2, 2}}, 0.48240812137449063});
  pyScfMOEriMap2.insert({{{2, 2, 2, 3}}, 1.771683264600496e-15});
  pyScfMOEriMap2.insert({{{2, 2, 3, 0}}, 6.398142519332816e-16});
  pyScfMOEriMap2.insert({{{2, 2, 3, 1}}, -0.12767413616212897});
  pyScfMOEriMap2.insert({{{2, 2, 3, 2}}, 1.8356257422537076e-15});
  pyScfMOEriMap2.insert({{{2, 2, 3, 3}}, 0.5168722520968734});
  pyScfMOEriMap2.insert({{{2, 3, 0, 0}}, -1.4118179407985874e-16});
  pyScfMOEriMap2.insert({{{2, 3, 0, 1}}, -0.138781481696745});
  pyScfMOEriMap2.insert({{{2, 3, 0, 2}}, 6.881415829148346e-16});
  pyScfMOEriMap2.insert({{{2, 3, 0, 3}}, 0.12014457557547478});
  pyScfMOEriMap2.insert({{{2, 3, 1, 0}}, -0.13878148169674503});
  pyScfMOEriMap2.insert({{{2, 3, 1, 1}}, -3.3280072395784165e-16});
  pyScfMOEriMap2.insert({{{2, 3, 1, 2}}, -0.06424955770793182});
  pyScfMOEriMap2.insert({{{2, 3, 1, 3}}, 3.7761030903612437e-16});
  pyScfMOEriMap2.insert({{{2, 3, 2, 0}}, 6.878765017297762e-16});
  pyScfMOEriMap2.insert({{{2, 3, 2, 1}}, -0.0642495577079318});
  pyScfMOEriMap2.insert({{{2, 3, 2, 2}}, 1.856374402575028e-15});
  pyScfMOEriMap2.insert({{{2, 3, 2, 3}}, 0.19560058197078858});
  pyScfMOEriMap2.insert({{{2, 3, 3, 0}}, 0.12014457557547477});
  pyScfMOEriMap2.insert({{{2, 3, 3, 1}}, 3.4793285115174377e-16});
  pyScfMOEriMap2.insert({{{2, 3, 3, 2}}, 0.19560058197078864});
  pyScfMOEriMap2.insert({{{2, 3, 3, 3}}, -2.4332677506218387e-15});
  pyScfMOEriMap2.insert({{{3, 0, 0, 0}}, -5.686451876643675e-16});
  pyScfMOEriMap2.insert({{{3, 0, 0, 1}}, -0.07747835724448592});
  pyScfMOEriMap2.insert({{{3, 0, 0, 2}}, 4.766055344570432e-17});
  pyScfMOEriMap2.insert({{{3, 0, 0, 3}}, 0.09471344004943409});
  pyScfMOEriMap2.insert({{{3, 0, 1, 0}}, -0.07747835724448592});
  pyScfMOEriMap2.insert({{{3, 0, 1, 1}}, -4.296503514688605e-16});
  pyScfMOEriMap2.insert({{{3, 0, 1, 2}}, -0.06118960215791895});
  pyScfMOEriMap2.insert({{{3, 0, 1, 3}}, 9.412618488197194e-16});
  pyScfMOEriMap2.insert({{{3, 0, 2, 0}}, 1.6506586529194315e-17});
  pyScfMOEriMap2.insert({{{3, 0, 2, 1}}, -0.061189602157918944});
  pyScfMOEriMap2.insert({{{3, 0, 2, 2}}, 6.735258453335178e-16});
  pyScfMOEriMap2.insert({{{3, 0, 2, 3}}, 0.1201445755754747});
  pyScfMOEriMap2.insert({{{3, 0, 3, 0}}, 0.0947134400494341});
  pyScfMOEriMap2.insert({{{3, 0, 3, 1}}, 9.206325440487932e-16});
  pyScfMOEriMap2.insert({{{3, 0, 3, 2}}, 0.12014457557547471});
  pyScfMOEriMap2.insert({{{3, 0, 3, 3}}, -2.1232692040182922e-15});
  pyScfMOEriMap2.insert({{{3, 1, 0, 0}}, -0.12802272602631076});
  pyScfMOEriMap2.insert({{{3, 1, 0, 1}}, -2.426433056554152e-16});
  pyScfMOEriMap2.insert({{{3, 1, 0, 2}}, -0.08038181988424582});
  pyScfMOEriMap2.insert({{{3, 1, 0, 3}}, 8.847769778053861e-16});
  pyScfMOEriMap2.insert({{{3, 1, 1, 0}}, -2.416359606674466e-16});
  pyScfMOEriMap2.insert({{{3, 1, 1, 1}}, -0.08142205002071345});
  pyScfMOEriMap2.insert({{{3, 1, 1, 2}}, 1.8972909591796262e-16});
  pyScfMOEriMap2.insert({{{3, 1, 1, 3}}, 0.08637200526415967});
  pyScfMOEriMap2.insert({{{3, 1, 2, 0}}, -0.08038181988424585});
  pyScfMOEriMap2.insert({{{3, 1, 2, 1}}, 2.0753165039088153e-16});
  pyScfMOEriMap2.insert({{{3, 1, 2, 2}}, -0.12767413616212914});
  pyScfMOEriMap2.insert({{{3, 1, 2, 3}}, 4.833641574176644e-16});
  pyScfMOEriMap2.insert({{{3, 1, 3, 0}}, 8.919915009748213e-16});
  pyScfMOEriMap2.insert({{{3, 1, 3, 1}}, 0.08637200526415967});
  pyScfMOEriMap2.insert({{{3, 1, 3, 2}}, 5.034323028029999e-16});
  pyScfMOEriMap2.insert({{{3, 1, 3, 3}}, -0.1524613002691411});
  pyScfMOEriMap2.insert({{{3, 2, 0, 0}}, -1.757481128500251e-16});
  pyScfMOEriMap2.insert({{{3, 2, 0, 1}}, -0.138781481696745});
  pyScfMOEriMap2.insert({{{3, 2, 0, 2}}, 6.936538547983168e-16});
  pyScfMOEriMap2.insert({{{3, 2, 0, 3}}, 0.12014457557547467});
  pyScfMOEriMap2.insert({{{3, 2, 1, 0}}, -0.138781481696745});
  pyScfMOEriMap2.insert({{{3, 2, 1, 1}}, -2.3511159195836706e-16});
  pyScfMOEriMap2.insert({{{3, 2, 1, 2}}, -0.06424955770793159});
  pyScfMOEriMap2.insert({{{3, 2, 1, 3}}, 5.574438916532033e-16});
  pyScfMOEriMap2.insert({{{3, 2, 2, 0}}, 6.860284742582298e-16});
  pyScfMOEriMap2.insert({{{3, 2, 2, 1}}, -0.06424955770793161});
  pyScfMOEriMap2.insert({{{3, 2, 2, 2}}, 1.753061247655299e-15});
  pyScfMOEriMap2.insert({{{3, 2, 2, 3}}, 0.19560058197078883});
  pyScfMOEriMap2.insert({{{3, 2, 3, 0}}, 0.12014457557547467});
  pyScfMOEriMap2.insert({{{3, 2, 3, 1}}, 5.503667127713373e-16});
  pyScfMOEriMap2.insert({{{3, 2, 3, 2}}, 0.1956005819707888});
  pyScfMOEriMap2.insert({{{3, 2, 3, 3}}, -2.188711754127681e-15});
  pyScfMOEriMap2.insert({{{3, 3, 0, 0}}, 0.5188746251919502});
  pyScfMOEriMap2.insert({{{3, 3, 0, 1}}, 1.822686578303778e-15});
  pyScfMOEriMap2.insert({{{3, 3, 0, 2}}, 0.14972481025996576});
  pyScfMOEriMap2.insert({{{3, 3, 0, 3}}, -2.3001008422707206e-15});
  pyScfMOEriMap2.insert({{{3, 3, 1, 0}}, 1.8358916752164607e-15});
  pyScfMOEriMap2.insert({{{3, 3, 1, 1}}, 0.44554829533298634});
  pyScfMOEriMap2.insert({{{3, 3, 1, 2}}, 4.527695169115665e-16});
  pyScfMOEriMap2.insert({{{3, 3, 1, 3}}, -0.15246130026914065});
  pyScfMOEriMap2.insert({{{3, 3, 2, 0}}, 0.1497248102599658});
  pyScfMOEriMap2.insert({{{3, 3, 2, 1}}, 4.136990792896897e-16});
  pyScfMOEriMap2.insert({{{3, 3, 2, 2}}, 0.516872252096874});
  pyScfMOEriMap2.insert({{{3, 3, 2, 3}}, -2.2377483288302723e-15});
  pyScfMOEriMap2.insert({{{3, 3, 3, 0}}, -2.3030047574572572e-15});
  pyScfMOEriMap2.insert({{{3, 3, 3, 1}}, -0.15246130026914062});
  pyScfMOEriMap2.insert({{{3, 3, 3, 2}}, -2.2164074589377752e-15});
  pyScfMOEriMap2.insert({{{3, 3, 3, 3}}, 0.5701251831541925});

  std::stringstream xyzInput("2\n\n"
                             "H 0 0 0\n"
                             "H 1.2 0 0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  auto mol = std::make_shared<Kiwi::Molecule>("6-31g*", atoms, moleculeSettings);

  auto data = std::make_shared<Kiwi::Data>(mol);

  data->integralDirect = false;
  data->oneBodyIntegrals();
  data->twoBodyIntegrals();

  for (auto const& elem : *mol) {
    data->X[elem.first] = Kiwi::Loewdin(data->S[elem.first], 1e-10, true);
    mol->at(elem.first).LMO = data->X[elem.first].cols();
  }
  auto dim = mol->at(Utils::ElementType::E).LMO;

  Kiwi::HartreeFockSettings settings;
  Kiwi::HartreeFockMain hartreeFockMain(mol, data, settings);

  Eigen::MatrixXd moMat(dim, dim);
  moMat << 0.2625636271904276, 0.1924218495385382, 0.9148555725536245, -0.9374180032850139, 0.3703828826631277,
      1.0056062593622421, -0.6986005248033106, 1.1369107558011013, 0.2625636271904279, -0.1924218495385377,
      0.9148555725536341, 0.9374180032850031, 0.3703828826631277, -1.0056062593622421, -0.6986005248033225,
      -1.1369107558010931;

  Utils::MolecularOrbitals mos = Utils::MolecularOrbitals::createFromRestrictedCoefficients(moMat);

  data->C.at(Utils::ElementType::E) = mos;

  Kiwi::AO2MO ao2mo(data, true);
  ao2mo.perform();
  auto integrals = ao2mo.getTwoBody()->at({Utils::ElementType::E, Utils::ElementType::E});

  for (auto mo : pyScfMOEriMap2) {
    std::array<int, 4> index = mo.first;
    EXPECT_THAT(integrals(index[0] * dim + index[1], index[2] * dim + index[3]), DoubleNear(pyScfMOEriMap2[mo.first], 1e-8));
  }
}

TEST_F(CalculationDataTest, CanPerformLoewdinOrtho) {
  std::stringstream xyzInput("4\n\n"
                             "H    0.6 0.0 0.0\n"
                             "H    0.4 0.0 0.0\n"
                             "H    0.2 0.0 0.0\n"
                             "H    0.0 0.0 0.0");

  auto atoms = Utils::XyzStreamHandler::readNuclearElectronic(xyzInput);

  Kiwi::MoleculeSettings moleculeSettings;
  moleculeSettings.multiplicity = 2;
  moleculeSettings.charge = -1;
  moleculeSettings.isRestricted = false;
  auto mol = std::make_shared<Kiwi::Molecule>("6-31g", atoms, moleculeSettings);

  auto data = Kiwi::Data(mol);

  data.oneBodyIntegrals();
  data.twoBodyIntegrals();

  data.X[Utils::ElementType::E] = Kiwi::Loewdin(data.S[Utils::ElementType::E], 1e-12, true);

  auto const& X = data.X[Utils::ElementType::E];
  auto const& S = data.S[Utils::ElementType::E];

  Eigen::MatrixXd test = X.transpose() * S * X;
  // Generate TEST matrix
  for (auto i = 0; i < test.rows(); i++) {
    // Diminish by one -> Matrix should be zero.
    test(i, i) -= 1.0;
  }
  auto loewdinError = test.norm();
  ASSERT_TRUE(loewdinError < 1e-12);
}

} // namespace Scine