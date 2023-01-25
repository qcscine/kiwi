/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Kiwi/HartreeFock/HartreeFockSettings.h>
#include <Kiwi/HartreeFock/HartreeFockUtils.h>
#include <Kiwi/HartreeFock/InitialGuess/InitialGuess.h>
#include <Kiwi/HartreeFock/InitialGuess/SAD.h>
#include <Kiwi/HartreeFock/InitialGuess/SND.h>
#include <Kiwi/HartreeFock/SecondOrder.h>
#include <Kiwi/HartreeFock/SerialScf.h>
#include <Kiwi/KiwiUtils/Data.h>
#include <Kiwi/KiwiUtils/GeneralUtility.h>
#include <LibintIntegrals/LibintIntegrals.h>
#include <Utils/DataStructures/IntegralSpecifier.h>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#define UNUSED(expr) \
  do {               \
    (void)(expr);    \
  } while (0)

namespace Scine {
namespace Kiwi {

template<>
auto InitialGuess::perform<InitialGuessSCF::Core>(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data,
                                                  const Utils::ElementType type, Utils::DensityMatrix& D,
                                                  Utils::MolecularOrbitals& C, bool verbose) -> void {
  if (verbose) {
    std::cout << "Performing Core guess for type " << Utils::ElementInfo::symbol(type) << "...\n";
  }

  Matrix tmp = CoreGuessCoefficients(data->H[type], data->X[type]);

  if (mol->at(type).isRestricted) {
    C = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Matrix>(std::move(tmp));
  }
  else {
    Matrix tmpcp = tmp;
    if (type != Utils::ElementType::E && mol->at(type).msVector[1] > 0) {
      double k = 0.9;
      double scaling = 1 / std::sqrt(1 + k * k);
      int HOMO = mol->at(type).msVector[0] - 1;
      int LUMO = mol->at(type).msVector[0];
      for (auto i = 0UL; i < mol->at(type).LMO; ++i) {
        auto oldHOMO = tmp(i, HOMO);
        auto oldLUMO = tmp(i, LUMO);
        tmp(i, HOMO) = scaling * (oldHOMO + k * oldLUMO);
        tmp(i, LUMO) = scaling * (-k * oldHOMO + oldLUMO);
      }
    }
    C = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Matrix>(std::move(tmp), std::move(tmpcp));
  }
  D = HartreeFockUtils::makeDensity(type, mol, C);
}

template<>
auto InitialGuess::perform<InitialGuessSCF::Hueckel>(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data,
                                                     const Utils::ElementType type, Utils::DensityMatrix& D,
                                                     Utils::MolecularOrbitals& C, bool verbose) -> void {
  if (type != Utils::ElementType::E) {
    throw std::runtime_error("Hueckel guess was called for nuclei!");
  }

  if (verbose) {
    std::cout << "Performing Hueckel guess for type " << Utils::ElementInfo::symbol(type) << "...\n";
  }

  Eigen::MatrixXd tmp;
  tmp = HueckelSTO3gGuessCoefficients(mol, data);
  if (mol->isRestricted()) {
    C = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Matrix>(std::move(tmp));
  }
  else {
    auto tmpcp = tmp;
    C = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Matrix>(std::move(tmp), std::move(tmpcp));
  }
  D = HartreeFockUtils::makeDensity(type, mol, C);
}

template<>
auto InitialGuess::perform<InitialGuessSCF::SND>(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data,
                                                 const Utils::ElementType type, Utils::DensityMatrix& D,
                                                 Utils::MolecularOrbitals& C, bool verbose) -> void {
  UNUSED(data);

  if (type == Utils::ElementType::E) {
    throw std::runtime_error("SND guess was called for electrons!");
  }

  if (verbose) {
    std::cout << "Performing SND guess for type " << Utils::ElementInfo::symbol(type) << "...\n";
    auto& clock = Clock::getInstance();
    clock.time("SND");
  }

  //
  // source: me
  //
  auto snd = Kiwi::SNDGuess(mol->at(type).positions, mol->at(type).basisSet, mol->at(Utils::ElementType::E).basisSet,
                            mol->at(type).msVector[0], mol->at(type).msVector[1]);
  snd.runScf();

  C = snd.getMOs();
  // This is necessary such that the density matrix was built once properly.
  D = HartreeFockUtils::makeDensity(type, mol, C);
  D.getAlphaMatrix() = snd.getDensity().alphaMatrix();
  D.getBetaMatrix() = snd.getDensity().betaMatrix();
  D.getRestrictedMatrix() = snd.getDensity().restrictedMatrix();

  if (verbose) {
    std::cout << "SND guess took ";
    auto& clock = Clock::getInstance();
    clock.time("SND");
  }
}

template<>
auto InitialGuess::perform<InitialGuessSCF::SAD>(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data,
                                                 const Utils::ElementType type, Utils::DensityMatrix& D,
                                                 Utils::MolecularOrbitals& C, bool verbose) -> void {
  UNUSED(data);

  if (type != Utils::ElementType::E) {
    throw std::runtime_error("SAD guess was called for nuclei!");
  }

  if (verbose) {
    std::cout << "Performing SAD guess for type " << Utils::ElementInfo::symbol(type) << "...\n";
    auto& clock = Clock::getInstance();
    clock.time("SAD");
  }

  //
  // https://pubs.rsc.org/en/content/articlehtml/2009/cp/b901987a
  //

  auto sad = SADGuess(mol->getGeometry(), mol->isRestricted(), mol->at(Utils::ElementType::E).basisSet);
  sad.computeInitialDensityMatrices();
  sad.runScf();

  C = sad.getMOs();
  D = HartreeFockUtils::makeDensity(type, mol, C);
  // This is necessary such that the density matrix was built once properly.
  D.getRestrictedMatrix() = sad.getDensity().restrictedMatrix();
  D.getAlphaMatrix() = sad.getDensity().alphaMatrix();
  D.getBetaMatrix() = sad.getDensity().betaMatrix();

  if (verbose) {
    std::cout << "SAD guess took ";
    auto& clock = Clock::getInstance();
    clock.time("SAD");
  }
}

template<>
auto InitialGuess::perform<InitialGuessSCF::SADNO>(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data,
                                                   const Utils::ElementType type, Utils::DensityMatrix& D,
                                                   Utils::MolecularOrbitals& C, bool verbose) -> void {
  if (type != Utils::ElementType::E) {
    throw std::runtime_error("SADNO guess was called for nuclei!");
  }

  if (verbose) {
    std::cout << "Performing SADNO guess for type " << Utils::ElementInfo::symbol(type) << "...\n";
    auto& clock = Clock::getInstance();
    clock.time("SADNO");
  }

  //
  // https://pubs.rsc.org/en/content/articlehtml/2009/cp/b901987a
  // https://pubs.acs.org/doi/10.1021/acs.jctc.8b01089
  //

  auto sad = SADGuess(mol->getGeometry(), mol->isRestricted(), mol->at(Utils::ElementType::E).basisSet);
  sad.computeInitialDensityMatrices();
  sad.runScf();

  D = sad.getDensity();

  const auto& X = data->X.at(Utils::ElementType::E);

  if (mol->isRestricted()) {
    Eigen::MatrixXd D_OAO = X.transpose() * D.restrictedMatrix() * X;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(D_OAO);
    Eigen::MatrixXd tmp = X * es.eigenvectors().rowwise().reverse();
    C = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Matrix>(std::move(tmp));
  }
  else {
    Eigen::MatrixXd D_OAO = X.transpose() * D.alphaMatrix() * X;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(D_OAO);
    Eigen::MatrixXd tmp = X * es.eigenvectors().rowwise().reverse();
    auto tmpcp = tmp;
    C = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Matrix>(std::move(tmp), std::move(tmpcp));
  }

  D = HartreeFockUtils::makeDensity(type, mol, C);

  if (verbose) {
    std::cout << "SADNO guess took ";
    auto& clock = Clock::getInstance();
    clock.time("SADNO");
  }
}

template<>
auto InitialGuess::perform<InitialGuessSCF::SadHueckel>(const std::shared_ptr<Molecule>& mol,
                                                        const std::shared_ptr<Data>& data, const Utils::ElementType type,
                                                        Utils::DensityMatrix& D, Utils::MolecularOrbitals& C, bool verbose)
    -> void {
  if (type != Utils::ElementType::E) {
    throw std::runtime_error("SadHueckel guess was called for nuclei!");
  }

  if (verbose) {
    std::cout << "Performing SAD-Hueckel guess for type " << Utils::ElementInfo::symbol(type) << "...\n";
  }

  //
  // https://pubs.rsc.org/en/content/articlehtml/2009/cp/b901987a
  // https://pubs.acs.org/doi/10.1021/acs.jctc.8b01089
  // https://www.sciencedirect.com/science/article/pii/S0009261412001996
  //

  auto sad = SADGuess(mol->getGeometry(), mol->isRestricted(), mol->at(Utils::ElementType::E).basisSet);
  sad.computeInitialDensityMatrices();
  sad.runScf();

  C = sad.getMOs();
  auto energies = sad.getEnergies();

  if (mol->isRestricted()) {
    auto tmp = HueckelGuessCoefficients(C.restrictedMatrix(), data->S[type], data->X[type], energies);
    C = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Matrix>(std::move(tmp));
  }
  else {
    auto tmp = HueckelGuessCoefficients(C.alphaMatrix(), data->S[type], data->X[type], energies);
    auto tmpcp = tmp;
    C = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Matrix>(std::move(tmp), std::move(tmpcp));
  }
  D = HartreeFockUtils::makeDensity(type, mol, C);
}

template<>
auto InitialGuess::perform<InitialGuessSCF::BO>(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data,
                                                const Utils::ElementType type, Utils::DensityMatrix& D,
                                                Utils::MolecularOrbitals& C, bool verbose) -> void {
  if (type != Utils::ElementType::E) {
    throw std::runtime_error("BO guess was called for nuclei!");
  }

  if (verbose) {
    std::cout << "-----------------\n";
    std::cout << "  BO-guess\n";
    std::cout << "-----------------\n\n";
  }

  auto boMol = std::make_shared<Molecule>(mol->getBOMolecule());
  auto boData = std::make_shared<Kiwi::Data>(boMol);
  boData->integralDirect = data->integralDirect;
  boData->oneBodyIntegrals();
  boData->twoBodyIntegrals();
  boData->makeExchange();
  boData->X[Utils::ElementType::E] = data->X[Utils::ElementType::E];
  Kiwi::HartreeFockSettings settings;
  auto trahSettings = std::make_shared<Kiwi::TRAHSettings>();
  Utils::DensityMatrix D_bo;
  Utils::MolecularOrbitals C_bo;

  perform<InitialGuessSCF::SAD>(boMol, boData, type, D_bo, C_bo);

  finalize(boMol, boData, type, D_bo, C_bo);

  Kiwi::SerialScf<Kiwi::SymmetryType::None> scf(boMol, boData, settings, ProjectionParameters(), true);
  scf.singleIteration();

  Kiwi::SecondOrder<Kiwi::SymmetryType::None> arh(boMol, boData, settings, ProjectionParameters(), verbose);
  arh.setTrahSettings(trahSettings);
  arh.run();

  if (verbose) {
    std::cout << "\n\n";
  }

  C = boData->C.at(Utils::ElementType::E);
  D = boData->D.at(Utils::ElementType::E);
}

auto InitialGuess::getDensity() const -> Utils::DensityMatrix {
  return D_;
}

auto InitialGuess::getMOs() const -> Utils::MolecularOrbitals {
  return C_;
}

auto InitialGuess::setVerbose(bool verbose) -> void {
  _verbose = verbose;
}

auto InitialGuess::CoreGuessCoefficients(const Eigen::MatrixXd& H, const Eigen::MatrixXd& X) -> Eigen::MatrixXd {
  auto F_guess = X.transpose() * H * X;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(F_guess);
  // C              = X             * eigvec
  // fullDim x fin_dim  = fullDim x fin_dim * fin_dim x fin_dim
  return X * es.eigenvectors();
}

auto InitialGuess::HueckelGuessCoefficients(const Eigen::MatrixXd& coefficients, const Eigen::MatrixXd& S,
                                            const Eigen::MatrixXd& X, const Eigen::VectorXd& ionizationPotentials)
    -> Eigen::MatrixXd {
  Eigen::MatrixXd Hueckel;
  Hueckel.resizeLike(S);

  Eigen::MatrixXd SC = S * coefficients;
  Eigen::MatrixXd overlapHueckel = coefficients.transpose() * SC;

  static constexpr double K = 1.75;

  for (auto i = 0; i < Hueckel.rows(); ++i) {
    for (auto j = 0; j < Hueckel.cols(); ++j) {
      if (i == j) {
        Hueckel(i, i) = ionizationPotentials(i);
      }
      else {
        // Hueckel(i, j) = 0.5 * S(i, j) * K * (H(i, i) + H(j, j));
        Hueckel(i, j) = 0.5 * overlapHueckel(i, j) * K * (ionizationPotentials(i) + ionizationPotentials(j));
      }
    }
  }

  // convert Hueckel matrix to AO basis
  Hueckel = SC * Hueckel * SC.transpose();

  auto F_guess = X.transpose() * Hueckel * X;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(F_guess);
  // C              = X             * eigvec
  // fullDim x fin_dim  = fullDim x fin_dim * fin_dim x fin_dim
  return X * es.eigenvectors();
}

/**
 * @brief This is a Hueckel guess, similar to the one implemented in Serenity. (Thanks Jan!)
 * @note It can only be used for electrons. The other Hueckel guess can be used for nuclei, too, although I am not
 *       sure if it makes sense.
 * @param mol
 * @param dataa
 * @return
 */
auto InitialGuess::HueckelSTO3gGuessCoefficients(const std::shared_ptr<Molecule>& mol,
                                                 const std::shared_ptr<Data>& dataa, bool verbose) -> Eigen::MatrixXd {
  Integrals::LibintIntegrals eval;

  if (verbose) {
    auto& clock = Clock::getInstance();
    clock.time("Hueckel");
    std::cout << "\tEvaluating integrals in minimal basis." << std::endl;
  }

  auto basisA = eval.initializeBasisSet("sto-3g", mol->getGeometry());

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Kinetic;
  auto resultMapKin = Integrals::LibintIntegrals::evaluate(specifier, basisA, basisA);
  const auto& Kin = resultMapKin[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  specifier.op = Utils::Integrals::Operator::PointCharges;
  specifier.atoms = mol->getGeometry();
  auto resultMapPC = Integrals::LibintIntegrals::evaluate(specifier, basisA, basisA);
  const auto& PC = resultMapPC[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  specifier.op = Utils::Integrals::Operator::Overlap;
  auto resultMapS = Integrals::LibintIntegrals::evaluate(specifier, basisA, basisA);
  const auto& S = resultMapS[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  Eigen::MatrixXd Hcore = Kin + PC;
  Eigen::MatrixXd Hueckel = Eigen::MatrixXd::Zero(S.rows(), S.cols());

  auto Parameters = EHTParameters();

  auto atom2shell = basisA.atomToShell(basisA.getAtoms());
  auto shell2bf = basisA.shell2bf();
  if (verbose) {
    std::cout << "\tDiagonalizing EHT matrix." << std::endl;
  }

  /*
   * Fill diagonal elements with tabulated values.
   */
  auto const& atoms = mol->getGeometry();
  for (int i = 0; i < atoms.size(); ++i) {
    if (Parameters.find(atoms.getElement(i)) == Parameters.end()) {
      throw std::runtime_error("Extended Hueckel paremeters not available for " +
                               Utils::ElementInfo::symbol(atoms.getElement(i)) + ".");
    }
    // This will throw an error if no parameters are available for that atom
    const auto& paramsForAtom = Parameters.at(atoms.getElement(i));
    // Loop over basis function shells of this atom; j is an index for the basis function vector
    auto Max = atom2shell[i][atom2shell[i].size() - 1] + 1;
    for (unsigned int j = atom2shell[i][0]; j < Max; ++j) {
      /*
       * In this loop cycle we fill the matrix starting at the first index belonging to this shell
       * of basis functions
       */
      const auto& firstIndexShell = shell2bf[j];
      /*
       * The atom parameter vector of course only has entries for basis functions of that atom, thus
       * we need to subtract the first index we use for accessing the basis function vector.
       */
      const auto& thisParam = paramsForAtom[j - atom2shell[i][0]];
      auto N = basisA[j].size();
      // Loop over the basis functions within a shell
      for (unsigned int k = 0; k < N; ++k) {
        Hueckel(firstIndexShell + k, firstIndexShell + k) = thisParam;
      }
    }
  }

  static constexpr double K = 1.75;

  for (auto i = 1; i < Hueckel.rows(); ++i) {
    for (auto j = 0; j < i; ++j) {
      Hueckel(i, j) = 0.5 * K * S(i, j) * (Hcore(i, i) + Hcore(j, j));
      Hueckel(j, i) = Hueckel(i, j);
    }
  }

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(Hueckel, S);

  if (verbose) {
    std::cout << "\tPerforming basis set projection." << std::endl;
  }
  /*
   * Basis set projection:
   *
   *  phi = \sum_i c_i |i>
   *
   *  I = \sum_IJ <I|J>^-1 |I><J|
   *
   *  phi = \sum_iIJ c_i <I|J>^-1 |I> <J|i> = \sum_J C_J |J>
   *
   *  C_o = C * Q^-1/2,  Q = C^T * C
   *
   */
  const Eigen::MatrixXd& coefficientsA = es.eigenvectors();
  // The number of AOs in Basis B
  const unsigned int nOrbitalsA = basisA.nbf();
  // The number of AOs in Basis B
  const auto& basisB = mol->at(Utils::ElementType::E).basisSet;
  const unsigned int nOrbitalsB = basisB.nbf();
  Utils::Integrals::IntegralSpecifier specifierProj;
  specifierProj.op = Utils::Integrals::Operator::Overlap;
  auto resultMapSAB = Integrals::LibintIntegrals::evaluate(specifierProj, basisA, basisB);
  auto SAB = resultMapSAB[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];
  SAB.transposeInPlace();
  // Calculate inverse of AO overlap integrals in basis B. Use SVD,
  // i.e. S_B = U * D * V^T, since overlapB could be ill-conditioned
  Eigen::VectorXd singularValues;
  Eigen::MatrixXd V;
  Eigen::MatrixXd U;

  const auto& SB = dataa->S[Utils::ElementType::E];
  {
    // Execute the svd only on one thread, because it's not unique.
    // This makes the results reproducible.
    auto n = omp_get_num_threads();
    omp_set_num_threads(1);
    {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(SB, Eigen::ComputeThinU | Eigen::ComputeThinV);
      singularValues = svd.singularValues();
      V = svd.matrixV();
      U = svd.matrixU();
    }
    omp_set_num_threads(n);
  }
  for (unsigned int i = 0; i < singularValues.rows(); ++i) {
    // if singular value is too small, set it to zero
    if (singularValues(i) < 1.0e-10) {
      singularValues(i) = 0.0;
    }
    else {
      singularValues(i) = 1 / singularValues(i);
    }
  }
  // build inverse, i.e. V * D * U^T
  // Eigen::MatrixXd overlapBinverse = svd.matrixV() * singularValues.asDiagonal() * svd.matrixU().transpose();
  Eigen::MatrixXd overlapBinverse = V * singularValues.asDiagonal() * U.transpose();
  // Calculate projection operator P_{BA}, i.e. the projector for the
  // transformation from basis A to basis B
  Eigen::MatrixXd projectionOperator = overlapBinverse * SAB;
  // Calculate new MO coefficients for basis B
  Eigen::MatrixXd newCoefficientsB = projectionOperator * coefficientsA;
  // Orthonormalzation:
  for (unsigned int i = 0; i < nOrbitalsA; i++) {
    double factor = 0;
    for (unsigned int j = 0; j < nOrbitalsB; j++) {
      for (unsigned int k = 0; k < nOrbitalsB; k++) {
        factor += newCoefficientsB(j, i) * newCoefficientsB(k, i) * SB(j, k);
      }
    }
    factor = sqrt(factor);
    for (unsigned int j = 0; j < nOrbitalsB; j++) {
      newCoefficientsB(j, i) /= factor;
    }
  }

  if (verbose) {
    std::cout << "EHT guess took ";
    auto& clock = Clock::getInstance();
    clock.time("Hueckel");
  }

  return newCoefficientsB;
}

auto InitialGuess::EHTParameters() -> std::map<Utils::ElementType, std::vector<double>> {
  std::map<Utils::ElementType, std::vector<double>> parameters;
  std::vector<double> paramsForElement;

  // EHT parameters
  // using AO energies from http://www.few.vu.nl/~visscher/FiniteNuclei/Table3.html
  //   (Nov. 8, 2017) author Lucas Visscher

  // First period
  paramsForElement.resize(1);
  // Hydrogen
  paramsForElement[0] = -5.0000666E-01;
  parameters[Utils::ElementType::H] = paramsForElement;
  // Helium
  paramsForElement[0] = -9.1799069E-01;
  parameters[Utils::ElementType::He] = paramsForElement;

  // Second period
  paramsForElement.resize(3);
  // Lithium
  paramsForElement[0] = -2.4779791E+00;
  paramsForElement[1] = -1.9633887E-01;
  paramsForElement[2] = 0.0;
  parameters[Utils::ElementType::Li] = paramsForElement;
  // Berylium
  paramsForElement[0] = -4.7334980E+00;
  paramsForElement[1] = -3.0932208E-01;
  paramsForElement[2] = 0.0;
  parameters[Utils::ElementType::Be] = paramsForElement;
  // Bor
  paramsForElement[0] = -7.6976454E+00;
  paramsForElement[1] = -4.9491914E-01;
  paramsForElement[2] = -3.0972112E-01;
  parameters[Utils::ElementType::B] = paramsForElement;
  // Carbon
  paramsForElement[0] = -1.1343601E+01;
  paramsForElement[1] = -7.1260267E-01;
  paramsForElement[2] = -4.0664220E-01;
  parameters[Utils::ElementType::C] = paramsForElement;
  // Nitrogen
  paramsForElement[0] = -1.5676414E+01;
  paramsForElement[1] = -9.6478811E-01;
  paramsForElement[2] = -5.0818404E-01;
  parameters[Utils::ElementType::N] = paramsForElement;
  // Oxygen
  paramsForElement[0] = -2.0698561E+01;
  paramsForElement[1] = -1.2524124E+00;
  paramsForElement[2] = -6.1537304E-01;
  parameters[Utils::ElementType::O] = paramsForElement;
  // Fluorine
  paramsForElement[0] = -2.6411757E+01;
  paramsForElement[1] = -1.5759831E+00;
  paramsForElement[2] = -7.2866290E-01;
  parameters[Utils::ElementType::F] = paramsForElement;
  // Neon
  paramsForElement[0] = -3.2817472E+01;
  paramsForElement[1] = -1.9358461E+00;
  paramsForElement[2] = -8.4826678E-01;
  parameters[Utils::ElementType::Ne] = paramsForElement;
  return parameters;
}

InitialGuess::InitialGuess(const InitialGuessSCF guess, const std::shared_ptr<Molecule>& mol,
                           const std::shared_ptr<Data>& data, const Utils::ElementType type, bool verbose)
  : _type(type) {
  switch (guess) {
    case InitialGuessSCF::Core:
      perform<InitialGuessSCF::Core>(mol, data, type, D_, C_, verbose);
      break;
    case InitialGuessSCF::Hueckel:
      perform<InitialGuessSCF::Hueckel>(mol, data, type, D_, C_, verbose);
      break;
    case InitialGuessSCF::BO:
      perform<InitialGuessSCF::BO>(mol, data, type, D_, C_, verbose);
      break;
    case InitialGuessSCF::SAD:
      perform<InitialGuessSCF::SAD>(mol, data, type, D_, C_, verbose);
      break;
    case InitialGuessSCF::SadHueckel:
      perform<InitialGuessSCF::SadHueckel>(mol, data, type, D_, C_, verbose);
      break;
    case InitialGuessSCF::SADNO:
      perform<InitialGuessSCF::SADNO>(mol, data, type, D_, C_, verbose);
      break;
    case InitialGuessSCF::SND:
      perform<InitialGuessSCF::SND>(mol, data, type, D_, C_, verbose);
      break;
    case InitialGuessSCF::Read:
      throw std::runtime_error("The SCF guess 'Read' must be handled elsewhere.!");
      break;
  }

  finalize(mol, data, type, D_, C_);
}

auto InitialGuess::finalize(const std::shared_ptr<Molecule>& mol, const std::shared_ptr<Data>& data,
                            const Utils::ElementType type, const Utils::DensityMatrix& D,
                            const Utils::MolecularOrbitals& C) -> void {
  Eigen::MatrixXd X_inv = data->X[type].completeOrthogonalDecomposition().pseudoInverse();
  data->D[type] = D;
  data->C[type] = C;
  if (mol->at(type).isRestricted) {
    Eigen::MatrixXd Crest_OAO = X_inv * data->C[type].restrictedMatrix();
    auto C_OAO = Utils::MolecularOrbitals::createFromRestrictedCoefficients<Eigen::MatrixXd>(std::move(Crest_OAO));
    data->C_OAO[type] = C_OAO;
    data->D_OAO[type] = HartreeFockUtils::makeDensity(type, mol, C_OAO);
  }
  else {
    Eigen::MatrixXd Cup_OAO = X_inv * data->C[type].alphaMatrix();
    Eigen::MatrixXd Cdown_OAO = X_inv * data->C[type].betaMatrix();
    auto C_OAO = Utils::MolecularOrbitals::createFromUnrestrictedCoefficients<Eigen::MatrixXd>(std::move(Cup_OAO),
                                                                                               std::move(Cdown_OAO));
    data->C_OAO[type] = C_OAO;
    data->D_OAO[type] = HartreeFockUtils::makeDensity(type, mol, C_OAO);
  }
}

} // namespace Kiwi
} // namespace Scine
