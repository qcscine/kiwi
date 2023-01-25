/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef KIWI_KEYWORDS_H
#define KIWI_KEYWORDS_H

#include <string>
#include <vector>

namespace Scine {
namespace Kiwi {
namespace Keywords {

static const std::vector<std::string> referenceCalculationSettings = {"type",
                                                                      "print keywords",
                                                                      "scf",
                                                                      "trah settings",
                                                                      "max iter",
                                                                      "guess",
                                                                      "nuclear guess",
                                                                      "integral direct",
                                                                      "reset incremental",
                                                                      "disable incremental",
                                                                      "max iter nested",
                                                                      "thresh",
                                                                      "loewdin thresh",
                                                                      "thresh fock",
                                                                      "mix angle",
                                                                      "mixes",
                                                                      "perturbations",
                                                                      "accelerator",
                                                                      "diis thresh",
                                                                      "nuclear diis thresh",
                                                                      "ediis thresh",
                                                                      "nuclear ediis thresh",
                                                                      "verbose",
                                                                      "stability analysis",
                                                                      "scf type"};

static const std::vector<std::string> trahSettings = {"optimizer",
                                                      "print keywords",
                                                      "initial trust radius",
                                                      "max trust radius",
                                                      "max davidson iterations",
                                                      "max davidson subspace dim",
                                                      "grad scaling",
                                                      "min thresh",
                                                      "local thresh",
                                                      "max arh dim",
                                                      "2nd start vector",
                                                      "2nd start vector noise",
                                                      "3rd start vector",
                                                      "3rd start vector noise"};

static const std::vector<std::string> ao2moSettings = {"type", "print keywords", "write", "thresh", "format"};

static const std::vector<std::string> particleDensitySettings = {
    "format",   "print keywords", "mode",           "particle type", "ms index", "mo index", "xmin", "xmax",
    "ymin",     "ymax",           "zmin",           "zmax",          "NX",       "NY",       "NZ",   "N",
    "rdm file", "type",           "rdm file alpha", "rdm file beta"};

static const std::vector<std::string> naturalOrbitalSettings = {"type", "print keywords", "rdm files"};

} // namespace Keywords
} // namespace Kiwi
} // namespace Scine

#endif // KIWI_KEYWORDS_H
