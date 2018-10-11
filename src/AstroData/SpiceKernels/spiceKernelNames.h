/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/

#ifndef SMART_ASTRO_SPICELSK_H
#define SMART_ASTRO_SPICELSK_H

#include <string>

namespace smartastro
{
namespace spiceKernels
{

    // Leap seconds
    const std::string       leap    = "SpiceKernels/lsk/naif0012.tls";

    // Planetary kernels
    const std::string       planets = "SpiceKernels/spk/de432s.bsp";

} // namespace spiceKernels
} // namespace smartastro

#endif //SMART_ASTRO_SPICELSK_H
