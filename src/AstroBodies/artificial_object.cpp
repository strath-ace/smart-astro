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

Artificial_Object::Artificial_Object( std::string givenName, int givenId)
{
  Astro_Body(givenName, givenId);
}

Artificial_Object::~Artificial_Object()
{
}

int Artificial_Object::gfrfov_( std::string &inst,  std::vector<double> &raydir,  std::string &rframe,  std::string &abcorr,  std::vector<double> &step,  std::vector<double> &cnfine,  std::vector<double> &result, int inst_len, int rframe_len, int abcorr_len, int obsrvr_len)
{
  char *instA = smartastro::astrcore::spice_general_functions::stringToCharArray(inst);
  double *raydirA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(raydir);
  char *rframeA = smartastro::astrcore::spice_general_functions::stringToCharArray(rframe);
  char *abcorrA = smartastro::astrcore::spice_general_functions::stringToCharArray(abcorr);
  double *stepA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(step);
  double *cnfineA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(cnfine);
  double *resultA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(result);
  
  return gfrfov_(instA, raydirA, rframeA, abcorrA, name, stepA, cnfineA, resultA, inst_len, rframe_len, abcorr_len, obsrvr_len);
}

int Artificial_Object::gftfov_( std::string &inst,  std::string &target,  std::string &tshape,  std::string &tframe,  std::string &abcorr,  std::vector<double> &step,  std::vector<double> &cnfine,  std::vector<double> &result, int inst_len, int target_len, int tshape_len, int tframe_len, int abcorr_len, int obsrvr_len)
{
  char *instA = smartastro::astrcore::spice_general_functions::stringToCharArray(inst);
  char *targetA = smartastro::astrcore::spice_general_functions::stringToCharArray(target);
  char *tshapeA = smartastro::astrcore::spice_general_functions::stringToCharArray(tshape);
  char *tframeA = smartastro::astrcore::spice_general_functions::stringToCharArray(tframe);
  char *abcorrA = smartastro::astrcore::spice_general_functions::stringToCharArray(abcorr);
  double *stepA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(step);
  double *cnfineA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(cnfine);
  double *resultA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(result);
  
  return gftfov_(instA, targetA, tshapeA, tframeA, abcorrA, nameA, stepA, cnfineA, resultA, inst_len, target_len, tshape_len, tframe_len, abcorr_len, obsrvr_len);
}

int Artificial_Object::fovray_( std::string &inst,  std::vector<double> &raydir,  std::string &rframe,  std::string &abcorr,  std::vector<double> &et,  std::vector<int> &visibl, int inst_len, int rframe_len, int abcorr_len, int obsrvr_len)
{
  char *instA = smartastro::astrcore::spice_general_functions::stringToCharArray(inst);
  double *raydirA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(raydir);
  char *rframeA = smartastro::astrcore::spice_general_functions::stringToCharArray(rframe);
  char *abcorrA = smartastro::astrcore::spice_general_functions::stringToCharArray(abcorr);
  double *etA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(et);
  int *visiblA = smartastro::astrcore::spice_general_functions::intVectorToArray(visibl);
  
  return fovray_(instA, raydirA, rframeA, abcorrA, name, etA, visiblA, inst_len, rframe_len, abcorr_len, obsrvr_len);
}

int Artificial_Object::fovtrg_( std::string &inst,  std::string &target,  std::string &tshape,  std::string &tframe,  std::string &abcorr,  std::vector<double> &et,  std::vector<int> &visibl, int inst_len, int target_len, int tshape_len, int tframe_len, int abcorr_len, int obsrvr_len)
{
  char *instA = smartastro::astrcore::spice_general_functions::stringToCharArray(inst);
  char *targetA = smartastro::astrcore::spice_general_functions::stringToCharArray(target);
  char *tshapeA = smartastro::astrcore::spice_general_functions::stringToCharArray(tshape);
  char *tframeA = smartastro::astrcore::spice_general_functions::stringToCharArray(tframe);
  char *abcorrA = smartastro::astrcore::spice_general_functions::stringToCharArray(abcorr);
  double *etA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(et);
  int *visiblA = smartastro::astrcore::spice_general_functions::intVectorToArray(visibl);
  
  return fovtrg_(instA, targetA, tshapeA, tframeA, abcorrA, nameA, etA, visiblA, inst_len, target_len, tshape_len, tframe_len, abcorr_len, obsrvr_len);
}

int Artificial_Object::getfov_( std::vector<int> &instid,  std::vector<int> &room,  std::string &shape,  std::string &frame,  std::vector<double> &bsight,  std::vector<int> &n,  std::vector<double> &bounds, int shape_len, int frame_len)
{
  int *instidA = smartastro::astrcore::spice_general_functions::intVectorToArray(instid);
  int *roomA = smartastro::astrcore::spice_general_functions::intVectorToArray(room);
  char *shapeA = smartastro::astrcore::spice_general_functions::stringToCharArray(shape);
  char *frameA = smartastro::astrcore::spice_general_functions::stringToCharArray(frame);
  double *bsightA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(bsight);
  int *nA = smartastro::astrcore::spice_general_functions::intVectorToArray(n);
  double *etA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(bounds);
  
  return getfov_(instidA, roomA, shapeA, frameA, bsightA, nA, boundsA, shape_len, frame_len);
}

int Artificial_Object::gffove_( std::string &inst,  std::string &tshape,  std::vector<double> &raydir,  std::string &target,  std::string &tframe,  std::string &abcorr,  std::vector<double> &tol, int udstep, int udrefn,  std::vector<int> &rpt, int udrepi, int udrepu, int udrepf,  std::vector<int> &bail, L_fp udbail,  std::vector<double> &cnfine,  std::vector<double> &result, int inst_len, int tshape_len, int target_len, int tframe_len, int abcorr_len, int obsrvr_len)
{
  char *instA = smartastro::astrcore::spice_general_functions::stringToCharArray(inst);
  char *tshapeA = smartastro::astrcore::spice_general_functions::stringToCharArray(tshape);
  double *raydirA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(raydir);
  char *targetA = smartastro::astrcore::spice_general_functions::stringToCharArray(target);
  char *tframeA = smartastro::astrcore::spice_general_functions::stringToCharArray(tframe);
  char *abcorrA = smartastro::astrcore::spice_general_functions::stringToCharArray(abcorr);
  double *tolA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(tol);
  int *rptA = smartastro::astrcore::spice_general_functions::intVectorToArray(rpt);
  int *bailA = smartastro::astrcore::spice_general_functions::intVectorToArray(bail);
  double *cnfineA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(cnfine);
  double *resultA = smartastro::astrcore::spice_general_functions::doubleVectorToArray(result);
  
  return gffove_(instA, tshapeA, raydirA, targetA, tframeA, abcorrA, name, tolA, udstep, udrefn, rptA, udrepi, udrepu, udrepf, bailA, udbail, cnfineA, resultA, inst_len, tshape_len, target_len, tframe_len, abcorr_len, obsrvr_len);
}
