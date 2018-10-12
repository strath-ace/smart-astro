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
    const std::string       leap                    = "SpiceKernels/lsk/naif0012.tls";

    // Planetary kernels
    const std::string       planetsEph              = "SpiceKernels/spk/de432s.bsp";

    // Planetary Orientation
    const std::string       planetsData             = "SpiceKernels/pck/pck00010.tpc";

    // Earth high accuracy orientation data
    const std::string       earthHighAccOrientation = "SpiceKernels/pck/earth_720101_070426.bpc";

    // Ground station ephemerides
    const std::string       groundStationEph        = "SpiceKernels/spk/earthstns_itrf93_050714.bsp";

    // Ground station topocentric frame
    const std::string       groundStationTopo       = "SpiceKernels/fk/earth_topo_050714.tf";


} // namespace spiceKernels
} // namespace smartastro

#endif //SMART_ASTRO_SPICELSK_H
