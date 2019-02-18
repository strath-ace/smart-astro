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

Artificial_Object::Artificial_Object(std::string givenName, int givenId) : Astro_Body(givenName, givenId)
{

}

Artificial_Object::~Artificial_Object()
{
}

void Artificial_Object::gfrfov(const std::string &inst,  std::vector<double> &raydir, const std::string &rframe, const std::string &abcorr, double &step, SpiceCell  &cnfine,  SpiceCell &result)
{  
  return gfrfov_c(&inst.at(0), &raydir[0], &rframe.at(0), &abcorr.at(0), &name.at(0), step, &cnfine, &result);
}

void Artificial_Object::gftfov(const std::string &inst, const std::string &target, const std::string &tshape, const std::string &tframe, const std::string &abcorr, double &step,  SpiceCell &cnfine,  SpiceCell &result)
{
  return gftfov_c(&inst.at(0), &target.at(0), &tshape.at(0), &tframe.at(0), &abcorr.at(0), &name.at(0), step, &cnfine, &result);
}

void Artificial_Object::fovray(const std::string &inst,  std::vector<double> &raydir, const std::string &rframe, const std::string &abcorr, double &et,  int &visibl)
{
  return fovray_c(&inst.at(0), &raydir[0], &rframe.at(0), &abcorr.at(0), &name.at(0), &et, &visibl);
}

void Artificial_Object::fovtrg( const std::string &inst, const std::string &target, const std::string &tshape, const std::string &tframe, const std::string &abcorr, double &et, int &visibl)
{
  return fovtrg_c(&inst.at(0), &target.at(0), &tshape.at(0), &tframe.at(0), &abcorr.at(0), &name.at(0), &et, &visibl);
}

/*void Artificial_Object::getfov(int &instid, int &room, int &shapelen, int &framelen, const std::string &shape, const std::string &frame,  std::vector<double> &bsight,  std::vector<int> &n,  std::vector<vector<double>> &bounds)
{
  return getfov_c(instid, room, shapelen, framelen, &shape.at(0), &frame.at(0), &bsight[0], &n, &bounds[0][0]);
  }*/

/*int Artificial_Object::gffove_( std::string &inst,  std::string &tshape,  std::vector<double> &raydir,  std::string &target,  std::string &tframe,  std::string &abcorr,  std::vector<double> &tol, int udstep, int udrefn,  std::vector<int> &rpt, int udrepi, int udrepu, int udrepf,  std::vector<int> &bail, L_fp udbail,  std::vector<double> &cnfine,  std::vector<double> &result, int inst_len, int tshape_len, int target_len, int tframe_len, int abcorr_len, int obsrvr_len)
{
  char *instA = smartastro::astrocore::spice_general_functions::stringToCharArray(inst);
  char *tshapeA = smartastro::astrocore::spice_general_functions::stringToCharArray(tshape);
  double *raydirA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(raydir);
  char *targetA = smartastro::astrocore::spice_general_functions::stringToCharArray(target);
  char *tframeA = smartastro::astrocore::spice_general_functions::stringToCharArray(tframe);
  char *abcorrA = smartastro::astrocore::spice_general_functions::stringToCharArray(abcorr);
  double *tolA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(tol);
  int *rptA = smartastro::astrocore::spice_general_functions::intVectorToArray(rpt);
  int *bailA = smartastro::astrocore::spice_general_functions::intVectorToArray(bail);
  double *cnfineA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(cnfine);
  double *resultA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(result);
  
  return gffove_(instA, tshapeA, raydirA, targetA, tframeA, abcorrA, nameA, tolA, udstep, udrefn, rptA, udrepi, udrepu, udrepf, bailA, udbail, cnfineA, resultA, inst_len, tshape_len, target_len, tframe_len, abcorr_len, obsrvr_len);
  }*/
