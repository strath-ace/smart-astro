/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#include "Astro-Core/spice_general_functions.h"

using namespace smartastro;
using namespace smartastro::astrocore;

void spice_general_functions::subslr(const std::string &method, const std::string &target, double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr, std::vector<double> &spoint, double &trgepc, std::vector<double> &srfvec)
{
  return subslr_c(&method.at(0), &target.at(0), et, &fixref.at(0), &abcorr.at(0), &obsrvr.at(0), &spoint[0], &trgepc, &srfvec[0]);
}

void spice_general_functions::subpnt(const std::string &method, const std::string &target, double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> srfvec)
{
  return subpnt_c(&method.at(0), &target.at(0), et, &fixref.at(0), &abcorr.at(0), &obsrvr.at(0), &spoint[0], &trgepc, &srfvec[0]);
}

double spice_general_functions::lspcn(const std::string &body, double &et, const std::string &abcorr)
{
  return lspcn_c(&body.at(0), et, &abcorr.at(0));
}

void spice_general_functions::furnsh(const std::string &file)
{
  return furnsh_c(&file.at(0));
}

void spice_general_functions::unload(const std::string &file)
{
  return unload_c(&file.at(0));
}
