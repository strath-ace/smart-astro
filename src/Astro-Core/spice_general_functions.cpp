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

void spice_general_functions::bodc2n(int &code, int &lenout, std::string &name, int &found)
{
  name.push_back('\0');
  
  return bodc2n_c(code, lenout, &name.at(0), &found);

  name.pop_back();
}

void spice_general_functions::bodn2c(std::string &name, int &code, int found)
{
  name.push_back('\0');

  return bodn2c_c(&name.at(0), &code, &found);
  
  name.pop_back();
}


void spice_general_functions::subslr(std::string &method, std::string &target, double &et, std::string &fixref, std::string &abcorr, std::string &obsrvr, std::vector<double> &spoint, double &trgepc, std::vector<double> &srfvec)
{
  method.push_back('\0');
  target.push_back('\0');
  fixref.push_back('\0');
  abcorr.push_back('\0');
  obsrvr.push_back('\0');
  
  return subslr_c(&method.at(0), &target.at(0), et, &fixref.at(0), &abcorr.at(0), &obsrvr.at(0), &spoint[0], &trgepc, &srfvec[0]);

  method.pop_back();
  target.pop_back();
  fixref.pop_back();
  abcorr.pop_back();
  obsrvr.pop_back();
}

void spice_general_functions::subpnt(std::string &method, std::string &target, double &et, std::string &fixref, std::string &abcorr, std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> srfvec)
{
  method.push_back('\0');
  target.push_back('\0');
  fixref.push_back('\0');
  abcorr.push_back('\0');
  obsrvr.push_back('\0');
  
  return subpnt_c(&method.at(0), &target.at(0), et, &fixref.at(0), &abcorr.at(0), &obsrvr.at(0), &spoint[0], &trgepc, &srfvec[0]);

  method.pop_back();
  target.pop_back();
  fixref.pop_back();
  abcorr.pop_back();
  obsrvr.pop_back();
}

double spice_general_functions::lspcn(std::string &body, double &et, std::string &abcorr)
{
  body.push_back('\0');
  abcorr.push_back('\0');
  
  return lspcn_c(&body.at(0), et, &abcorr.at(0));

  body.pop_back();
  abcorr.pop_back();
}

void spice_general_functions::furnsh(std::string &file)
{
  file.push_back('\0');
  
  return furnsh_c(&file.at(0));

  file.pop_back();
}

void spice_general_functions::unload(std::string &file)
{
  file.push_back('\0');
  
  return unload_c(&file.at(0));

  file.pop_back();
}
