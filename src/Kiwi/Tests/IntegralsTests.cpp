/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <LibintIntegrals/LibintIntegrals.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Settings.h>
#include <gmock/gmock.h>

namespace Scine {

using namespace testing;

class IntegralsTest : public Test {};

TEST_F(IntegralsTest, canCalculate) {
  Eigen::MatrixXd pyScfOverlapH2Def2SVP;
  pyScfOverlapH2Def2SVP.resize(10, 10);
  pyScfOverlapH2Def2SVP << 1.00000000e+00, 6.84799825e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      2.53303682e-01, 4.14114249e-01, -2.86782459e-01, 0.00000000e+00, 0.00000000e+00, 6.84799825e-01, 1.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 4.14114249e-01, 7.30845790e-01, -1.73676647e-01, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.86782459e-01,
      1.73676647e-01, -3.98093618e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.27845428e-01, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 1.27845428e-01, 2.53303682e-01, 4.14114249e-01, 2.86782459e-01, 0.00000000e+00,
      0.00000000e+00, 1.00000000e+00, 6.84799825e-01, -8.49505555e-18, 0.00000000e+00, 0.00000000e+00, 4.14114249e-01,
      7.30845790e-01, 1.73676647e-01, 0.00000000e+00, 0.00000000e+00, 6.84799825e-01, 1.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, -2.86782459e-01, -1.73676647e-01, -3.98093618e-01, 0.00000000e+00, 0.00000000e+00,
      -8.49505555e-18, 0.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
      0.00000000e+00, 1.27845428e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.27845428e-01, 0.00000000e+00,
      0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00;

  std::stringstream h2("2\n\n"
                       "H 0 0 0\n"
                       "H 1.2 0 0");
  auto scineAtoms = Utils::XyzStreamHandler::read(h2);

  Integrals::LibintIntegrals eval;

  Utils::Settings& settings = eval.settings();

  settings.modifyBool("use_pure_spherical", false);

  auto basis = eval.initializeBasisSet("def2-svp", scineAtoms);

  Utils::Integrals::IntegralSpecifier specifier;
  specifier.op = Utils::Integrals::Operator::Overlap;

  auto result_map = Integrals::LibintIntegrals::evaluate(specifier, basis, basis);

  const auto& result = result_map[{Utils::Integrals::Component::none, Utils::Integrals::DerivKey::value, 0}];

  for (int row = 0; row < 10; ++row) {
    for (int col = 0; col < 10; ++col) {
      EXPECT_THAT(result(row, col), DoubleNear(pyScfOverlapH2Def2SVP(row, col), 1e-8));
    }
  }
}

} // namespace Scine
