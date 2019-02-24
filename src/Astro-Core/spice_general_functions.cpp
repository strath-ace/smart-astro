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

void spice_general_functions::vectorToSpiceDSKDescr(const std::vector<int> &vector, SpiceDSKDescr &descr)
{
  if(vector.size() == 24){
    
    descr.surfce = vector[0];
    descr.center = vector[1];
    descr.dclass = vector[2];
    descr.dtype = vector[3];
    descr.frmcde = vector[4];
    descr.corsys = vector[5];
    descr.co1min = vector[16];
    descr.co1max = vector[17];
    descr.co2min = vector[18];
    descr.co2max = vector[19];
    descr.co3min = vector[20];
    descr.co3max = vector[21];
    descr.start = vector[22];
    descr.stop = vector[23];

    for(int i = 6; i < 16; ++i){
      descr.corpar[i - 6] = vector[i];
    }
  } else {
    smartastro_throw("Vector is not the correct size for a SpiceDSKDescr.");
  }
}

void spice_general_functions::spiceDSKDescrToVector(const SpiceDSKDescr &descr, std::vector<int> &vector)
{
  if(vector.size() == 24){
    vector[0] = descr.surfce;
    vector[1] = descr.center;
    vector[2] = descr.dclass;
    vector[3] = descr.dtype;
    vector[4] = descr.frmcde;
    vector[5] = descr.corsys;
    vector[16] = descr.co1min;
    vector[17] = descr.co1max;
    vector[18] = descr.co2min;
    vector[19] = descr.co2max;
    vector[20] = descr.co3min;
    vector[21] = descr.co3max;
    vector[22] = descr.start;
    vector[23] = descr.stop;

    for(int i = 6; i < 16; ++i){
      vector[i] = descr.corpar[i - 6];
    }
  } else {
    smartastro_throw("Vector is not the correct size for a SpiceDSKDescr.");
  }
}

void spice_general_functions::vectorToSpiceDLADescr(const std::vector<int> &vector, SpiceDLADescr &descr)
{
  if(vector.size() == 8){
    descr.bwdptr = vector[0];
    descr.fwdptr = vector[1];
    descr.ibase = vector[2];
    descr.isize = vector[3];
    descr.dbase = vector[4];
    descr.dsize = vector[5];
    descr.cbase = vector[6];
    descr.csize = vector[7];
  } else {
    smartastro_throw("Vector is not the correct size for a SpiceDLADescr.");
  }
}

void spice_general_functions::spiceDLADescrToVector(const SpiceDLADescr &descr, std::vector<int> &vector)
{
  if(vector.size() == 8){
    vector[0] = descr.bwdptr;
    vector[1] = descr.fwdptr;
    vector[2] = descr.ibase;
    vector[3] = descr.isize;
    vector[4] = descr.dbase;
    vector[5] = descr.dsize;
    vector[6] = descr.cbase;
    vector[7] = descr.csize;
  } else {
    smartastro_throw("Vector is not the correct size for a SpiceDLADescr.");
  }
}

void spice_general_functions::vectorToSpicePlaneNVC(const std::vector<double> &vector, SpicePlane &plane)
{
  if(vector.size() == 4){
    double normal[3] = {vector[0], vector[1], vector[2]};
    double constant = vector[4];

    nvc2pl_c(&normal[0], constant, &plane);
  } else {
    smartastro_throw("Vector is not the correct size for a plane.");
  }
}

void spice_general_functions::vectorToSpicePlaneNVP(const std::vector<double> &vector, SpicePlane &plane)
{
  if(vector.size() == 6){
    double normal[3] = {vector[0], vector[1], vector[2]};
    double point[3] = {vector[0], vector[1], vector[2]};

    nvp2pl_c(&normal[0], &point[0], &plane);
  } else {
    smartastro_throw("Vector is not the correct size for a plane.");
  }
}

void spice_general_functions::vectorToSpicePlanePSV(const std::vector<double> &vector, SpicePlane &plane)
{
  if(vector.size() == 9){
    double point[3] = {vector[0], vector[1], vector[2]};
    double span1[3] = {vector[3], vector[4], vector[5]};
    double span2[3] = {vector[6], vector[7], vector[8]};

    psv2pl_c(&point[0], &span1[0], &span2[0], &plane);
  } else {
    smartastro_throw("Vector is not the correct size for a plane.");
  }
}

void spice_general_functions::vectorToSpiceEllipse(const std::vector<double> &vector, SpiceEllipse &ellipse)
{
  if(vector.size() == 9){
    double center[3] = {vector[0], vector[1], vector[2]};
    double v1[3] = {vector[3], vector[4], vector[5]};
    double v2[3] = {vector[6], vector[7], vector[8]};

    cgv2el_c (center, v1, v2, &ellipse);
  } else {
    smartastro_throw("Vector is not the correct size for an ellipse.");
  }
}

void spice_general_functions::spiceEllipseToVector(const SpiceEllipse &ellipse, std::vector<double> &vector)
{
  if(vector.size() == 9){
    double center[3];
    double v1[3];
    double v2[3];
    
    el2cgv_c (&ellipse, center, v1, v2);

    vector[0] = center[0];
    vector[1] = center[1];
    vector[2] = center[2];
    vector[3] = v1[0];
    vector[4] = v1[1];
    vector[5] = v1[2];
    vector[6] = v2[0];
    vector[7] = v2[1];
    vector[8] = v2[2];
    
  } else {
    smartastro_throw("Vector is not the correct size for an ellipse.");
  }
}

void spice_general_functions::doubleVectorTo2dArray2(const std::vector<double> &vector, double (*array)[2])
{
  int rows =  sizeof array / sizeof array[0];
  int size = vector.size() - 1;
  int vectorIndexCnt = 0;

  if(vector.size() == rows * 2){
  for(int i = 0; i < rows; ++i){
      for(int j = 0; j < 2; ++j){
        array[i][j] = vector[vectorIndexCnt];
	++vectorIndexCnt;
      }
    }
  } else {
    smartastro_throw("Vector is not the correct size for the array.");
  }
}

void spice_general_functions::doubleVectorTo2dArray3(const std::vector<double> &vector, double (*array)[3])
{
  int rows =  sizeof array / sizeof array[0];
  int vectorIndexCnt = 0;

  if(vector.size() == rows * 3){
  for(int i = 0; i < rows; ++i){
      for(int j = 0; j < 3; ++j){
        array[i][j] = vector[vectorIndexCnt];
	++vectorIndexCnt;
      }
    }
  } else {
    smartastro_throw("Vector is not the correct size for the array.");
  }
}
  
void spice_general_functions::double2dArray3ToVector(const double (*array)[3], std::vector<double> &vector)
{
  int rows =  sizeof array / sizeof array[0];
  int vectorIndexCnt = 0;

  if(vector.size() == rows * 3){
  for(int i = 0; i < rows; ++i){
      for(int j = 0; j < 3; ++j){
	vector[vectorIndexCnt] = array[i][j];
	++vectorIndexCnt;
      }
    }
  } else {
    smartastro_throw("Vector is not the correct size for the array.");
  }
}

void spice_general_functions::double2dArray2ToVector(const double (*array)[2], std::vector<double> &vector)
{
  int rows =  sizeof array / sizeof array[0];
  int vectorIndexCnt = 0;

  if(vector.size() == rows * 2){
  for(int i = 0; i < rows; ++i){
      for(int j = 0; j < 3; ++j){
	vector[vectorIndexCnt] = array[i][j];
        ++vectorIndexCnt;
      }
    }
  } else {
    smartastro_throw("Vector is not the correct size for the array.");
  }
}

void spice_general_functions::bodc2n(const int &code, const int &lenout, std::string &name, int &found)
{
  name.push_back('\0');
  
  bodc2n_c(code, lenout, &name.at(0), &found);

  name.pop_back();
}

void spice_general_functions::bodn2c(const std::string &name, int &code, int found)
{
  return bodn2c_c(name.c_str(), &code, &found);
}

void spice_general_functions::subslr(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr, std::vector<double> &spoint, double &trgepc, std::vector<double> &srfvec)
{
  return subslr_c(method.c_str(), target.c_str(), et, fixref.c_str(), abcorr.c_str(), obsrvr.c_str(), &spoint[0], &trgepc, &srfvec[0]);
}

void spice_general_functions::subpnt(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> srfvec)
{
  return subpnt_c(method.c_str(), target.c_str(), et, fixref.c_str(), abcorr.c_str(), obsrvr.c_str(), &spoint[0], &trgepc, &srfvec[0]);
}

double spice_general_functions::lspcn(const std::string &body, const double &et, const std::string &abcorr)
{
  return lspcn_c(body.c_str(), et, abcorr.c_str());
}

void spice_general_functions::furnsh(const std::string &file)
{
  return furnsh_c(file.c_str());
}

void spice_general_functions::unload(const std::string &file)
{
  return unload_c(file.c_str());
}
