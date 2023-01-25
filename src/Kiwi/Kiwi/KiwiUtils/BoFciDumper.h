/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_BOFCIDUMPER_H
#define KIWI_BOFCIDUMPER_H

#include <LibintIntegrals/TwoBodyIntegrals/SymmetryHelper.h>
#include <Eigen/Core>
#include <string>

namespace Scine {
namespace Kiwi {

/**
 * @brief Struct containing the information to be parsed from/written into a FciDump file.
 */
struct BoFciDumpData {
  //! @brief NORB
  int nOrbitals;
  //! @brief NELEC
  int nElectrons;
  //! @brief \f$E_{core}\f$
  double coreEnergy;
  //! @brief MS2, spin multiplicity
  int spinPolarization = 0;
  //! @brief UHF
  bool unrestrictedReference = false;
  //! @brief ISYM
  int wfSymIfLowestOccupied = 1;
  //! @brief ORBSYM
  std::vector<int> orbitalSymmetry;
};

/**
 * @class FciDumper @file FciDumper.h
 * @brief Class for input and output of FCIDump files for one- and two-bodies integrals.
 */
class BoFciDumper {
 public:
  /**
   * @brief Enum class for possible input and output formats.
   */
  enum class format { fci };
  /**
   * @brief Static variable to control the threshold value for integrals to be printed.
   */
  static double integralThreshold;
  /**
   * @brief Write one- and two-bodies integrals to a file.
   *
   * @param fileName The file path.
   * @param hCore The one-body integrals matrix.
   * @param eris The two-bodies integrals matrix.
   * @param fciDumpData Struct containing auxiliary infos as number of orbitals, number of electrons, symmetries...
   * @param f The format to be used/expected.
   * @throws std::runtime_error If the file could not be created.
   */
  static void write(const std::string& filename, const Eigen::MatrixXd& hCore, const Eigen::MatrixXd& eris,
                    const BoFciDumpData& fciDumpData, format f = format::fci);
  /**
   * @brief Write one- and two-bodies integrals to a stream.
   *
   * @param out  The output stream.
   * @param hCore The one-body integrals matrix.
   * @param eris The two-bodies integrals matrix.
   * @param fciDumpData Struct containing auxiliary infos as number of orbitals, number of electrons, symmetries...
   * @param f The format to be used/expected.
   */
  static void write(std::ostream& out, const Eigen::MatrixXd& hCore, const Eigen::MatrixXd& eris,
                    const BoFciDumpData& fciDumpData, format f = format::fci);
  /**
   * @brief Read one- and two-bodies integrals from a file.
   *
   * @param fileName The file path.
   * @param f The format to be used/expected.
   * @throws std::runtime_error If the file could not be opened.
   * @return Returns a std::tuple<Eigen::MatrixXd, EriType> with the one-, two-bodies integrals and auxiliary infos.
   */
  static std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, BoFciDumpData> read(const std::string& filename, format f = format::fci);
  /**
   * @brief Read one- and two-bodies integrals from a stream.
   *
   * @param in The input stream.
   * @param f The format to be used/expected.
   * @throws std::runtime_error If the format isn't supported.
   * @return Returns a std::tuple<Eigen::MatrixXd, EriType, BoFciDumpData> with one-, two-bodies integrals and auxiliary
   * infos.
   */
  static std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, BoFciDumpData> read(std::istream& in, format f = format::fci);

 private:
  static void writeFci(std::ostream& out, const Eigen::MatrixXd& hCore, const Eigen::MatrixXd& eris,
                       const BoFciDumpData& fciDumpData);
  static void writeHeaderFci(std::ostream& out, const BoFciDumpData& fciDumpData);
  static void writeOneElectronIntegralsFci(std::ostream& out, const Eigen::MatrixXd& hCore);
  static void writeTwoElectronIntegralsFci(std::ostream& out, const Eigen::MatrixXd& eris);
  static std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, BoFciDumpData> readFci(std::istream& in);
  static BoFciDumpData parseFciHeader(std::istream& in);
  static std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> parseFciIntegrals(std::istream& in, BoFciDumpData& fciDumpData);
};

} // namespace Kiwi
} // namespace Scine

#endif // KIWI_BOFCIDUMPER_H
