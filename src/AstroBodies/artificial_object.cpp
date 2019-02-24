/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#include "AstroBodies/artificial_object.h"

using namespace smartastro;
using namespace smartastro::astrobodies;

Artificial_Object::Artificial_Object(std::string givenName) : Astro_Body(givenName)
{
}

Artificial_Object::Artificial_Object(int givenId) : Astro_Body(givenId)
{
}


void Artificial_Object::gfrfov(const std::string &inst, const std::vector<double> &raydir, const std::string &rframe, const std::string &abcorr, const double &step, SpiceCell  &cnfine,  SpiceCell &result)
{
  return gfrfov_c(inst.c_str(), &raydir[0], rframe.c_str(), abcorr.c_str(), name.c_str(), step, &cnfine, &result);
}

void Artificial_Object::gftfov(const std::string &inst, const std::string &target, const std::string &tshape, const std::string &tframe, const std::string &abcorr, const double &step,  SpiceCell &cnfine,  SpiceCell &result)
{
  return gftfov_c(inst.c_str(), target.c_str(), tshape.c_str(), tframe.c_str(), abcorr.c_str(), name.c_str(), step, &cnfine, &result);
}

void Artificial_Object::fovray(const std::string &inst, const std::vector<double> &raydir, const std::string &rframe, const std::string &abcorr, double &et, int &visible)
{
  return fovray_c(inst.c_str(), &raydir[0], rframe.c_str(), abcorr.c_str(), name.c_str(), &et, &visible);
}

void Artificial_Object::fovtrg(const std::string &inst, const std::string &target, const std::string &tshape, const std::string &tframe, const std::string &abcorr, double &et, int &visible)
{ 
  return fovtrg_c(inst.c_str(), target.c_str(), tshape.c_str(), tframe.c_str(), abcorr.c_str(), name.c_str(), &et, &visible);
}

void Artificial_Object::getfov(const int &instid, const int &room, const int &shapelen, const int &framelen, std::string &shape, std::string &frame,  std::vector<double> &bsight, int &n, std::vector<double> &bounds)
{
  if(bounds.size() % 3 == 0){
    int noOfRows = bounds.size() / 3;
    double boundsA[noOfRows][3];
    smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(bounds, boundsA);
    shape.push_back('\0');
    frame.push_back('\0');

    getfov_c(instid, room, shapelen, framelen, &shape.at(0), &frame.at(0), &bsight[0], &n, boundsA);

    smartastro::astrocore::spice_general_functions::double2dArray3ToVector(boundsA, bounds);
    shape.pop_back();
    frame.pop_back();
  } else {
    smartastro_throw("bounds vector is not the correct size for the array.");
  }
}

void Artificial_Object::gffove(const std::string &inst, const std::string &tshape, const std::vector<double> &raydir, const std::string &target, const std::string &tframe, const std::string &abcorr, const double &tol, const int &rpt, const int &bail, std::vector<double> &cnfine, std::vector<double> &result)
{
  SPICEDOUBLE_CELL (cnfineCell, 100);
  SPICEDOUBLE_CELL (resultCell, 100);

  for(int i = 0; i < cnfine.size(); ++i){
    SPICE_CELL_SET_D (cnfine[i], i, &cnfineCell);
  }

  for(int i = 0; i < result.size(); ++i){
    SPICE_CELL_SET_D (result[i], i, &resultCell);
  }
  
  gffove_c(inst.c_str(), tshape.c_str(), &raydir[0], target.c_str(), tframe.c_str(), abcorr.c_str(), name.c_str(), tol, gfstep_c, gfrefn_c, rpt, gfrepi_c, gfrepu_c, gfrepf_c, bail, gfbail_c, &cnfineCell, &resultCell);

  for(int i = 0; i < cnfine.size(); ++i){
    SPICE_CELL_GET_D (&cnfineCell, i, &cnfine[i]);
  }

  for(int i = 0; i < result.size(); ++i){
    SPICE_CELL_GET_D (&resultCell, i, &result[i]);
  }
}
