/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#include "AstroBodies/celestial_object.h"

using namespace smartastro;
using namespace smartastro::astrobodies;

Celestial_Object::Celestial_Object(std::string givenName, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &givenPositn, double givenMu) : Astro_Body(givenName, spiceEphemeridesParams)
{
  positn = givenPositn;
  mu = givenMu;
}

Celestial_Object::Celestial_Object(std::string givenName, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>, std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &givenPositn, double givenMu) : Astro_Body(givenName, referenceFrame, integratedEphemeridesFunction, integratedEphemeridesParams)
{
  positn = givenPositn;
  mu = givenMu;
}

Celestial_Object::Celestial_Object(int givenId, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &givenPositn, double givenMu) : Astro_Body(givenId, spiceEphemeridesParams)
{
  positn = givenPositn;
  mu = givenMu;
}

Celestial_Object::Celestial_Object(int givenId, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>, std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &givenPositn, double givenMu) : Astro_Body(givenId, referenceFrame, integratedEphemeridesFunction, integratedEphemeridesParams)
{
  positn = givenPositn;
  mu = givenMu;
}

void Celestial_Object::latsrf(const std::string &method, const std::string &target, double &et, const std::string &fixref, const int &npts, const std::vector<double> &lonlat, std::vector<double> &srfpts)
{
  if(lonlat.size() % 2 == 0 && srfpts.size() % 3 == 0){
     int noOfRowsl = lonlat.size() / 2;
     double lonlatA[noOfRowsl][2];
     int noOfRowss = srfpts.size() / 3;
     double srfptsA[noOfRowss][3];
     
     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray2(lonlat, lonlatA);
     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(srfpts, srfptsA);

     latsrf_c(method.c_str(), target.c_str(), et, fixref.c_str(), npts, lonlatA, srfptsA);

     smartastro::astrocore::spice_general_functions::double2dArray3ToVector(srfptsA, srfpts);
  } else {
     smartastro_throw("Vectors are not the correct size for the arrays.");
  }
}

void Celestial_Object::srfnrm(const std::string &method, const std::string &target, double &et, const std::string &fixref, const int &npts, const std::vector<double> &srfpts, std::vector<double> &normls)
{
  if(srfpts.size() % 3 == 0 && normls.size() % 3 == 0){
     int noOfRowss = srfpts.size() / 3;
     double srfptsA[noOfRowss][3];
     int noOfRowsn = normls.size() / 3;
     double normlsA[noOfRowsn][3];

     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(srfpts, srfptsA);
     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(normls, normlsA);
     
     srfnrm_c(method.c_str(), target.c_str(), et, fixref.c_str(), npts, srfptsA, normlsA);

     smartastro::astrocore::spice_general_functions::double2dArray3ToVector(normlsA, normls);
  } else {
     smartastro_throw("Vectors are not the correct size for the arrays.");
  }
}

void Celestial_Object::nearpt(const double &a, const double &b, const double &c, std::vector<double> &npoint, double &alt)
{
  return nearpt_c(&positn[0], a, b, c, &npoint[0], &alt);
}

void Celestial_Object::surfpt(const std::vector<double> &u, const double &a, const double &b, const double &c,  std::vector<double> &point, int &found)
{
  return surfpt_c(&positn.at(0), &u[0], a, b, c, &point[0], &found);
}

void Celestial_Object::surfnm(const double &a, const double &b, const double &c, const std::vector<double> &point,  std::vector<double> &normal)
{
  return surfnm_c(a, b, c, &point[0], &normal[0]);
}

void Celestial_Object::edlimb(const double &a, const double &b, const double &c, const std::vector<double> &viewpt,  std::vector<double> &limb)
{
  SpiceEllipse limbE;
  smartastro::astrocore::spice_general_functions::vectorToSpiceEllipse(limb, limbE);
  
  edlimb_c(a, b, c, &viewpt[0], &limbE);

  smartastro::astrocore::spice_general_functions::spiceEllipseToVector(limbE, limb);
}

void Celestial_Object::inelpl(const std::vector<double> &ellips, const std::vector<double> &plane, const int spicePlaneType, int &nxpts, std::vector<double> &xpt1, std::vector<double> &xpt2)
{
  SpiceEllipse ellipsE;
  SpicePlane planeP;

  smartastro::astrocore::spice_general_functions::vectorToSpiceEllipse(ellips, ellipsE);

  switch(spicePlaneType){
  case 1: smartastro::astrocore::spice_general_functions::vectorToSpicePlaneNVC(plane, planeP);
    break;
  case 2: smartastro::astrocore::spice_general_functions::vectorToSpicePlaneNVP(plane, planeP);
    break;
  case 3: smartastro::astrocore::spice_general_functions::vectorToSpicePlanePSV(plane, planeP);
    break;
  }
  
  return inelpl_c(&ellipsE, &planeP, &nxpts, &xpt1[0], &xpt2[0]);
}

void Celestial_Object::npedln(const double &a, const double &b, const double &c, const std::vector<double> &linept, const std::vector<double> &linedr, std::vector<double> &pnear, double &dist)
{
  npedln_c(a, b, c, &linept[0], &linedr[0], &pnear[0], &dist);
}

void Celestial_Object::inrypl(const std::vector<double> &vertex, const std::vector<double> &dir, const std::vector<double> &plane, const int spicePlaneType, int &nxpts, std::vector<double> &xpt)
{
  SpicePlane planeP;
  
  switch(spicePlaneType){
  case 1: smartastro::astrocore::spice_general_functions::vectorToSpicePlaneNVC(plane, planeP);
    break;
  case 2: smartastro::astrocore::spice_general_functions::vectorToSpicePlaneNVP(plane, planeP);
    break;
  case 3: smartastro::astrocore::spice_general_functions::vectorToSpicePlanePSV(plane, planeP);
    break;
  }
  
  return inrypl_c(&vertex[0], &dir[0], &planeP, &nxpts, &xpt[0]);
}

void Celestial_Object::pjelpl(const std::vector<double> &elin, const std::vector<double> &plane, const int spicePlaneType, std::vector<double> &elout)
{
  SpiceEllipse elinE;
  SpiceEllipse eloutE;
  SpicePlane planeP;

  smartastro::astrocore::spice_general_functions::vectorToSpiceEllipse(elin, elinE);
  smartastro::astrocore::spice_general_functions::vectorToSpiceEllipse(elout, eloutE);
  
  switch(spicePlaneType){
  case 1: smartastro::astrocore::spice_general_functions::vectorToSpicePlaneNVC(plane, planeP);
    break;
  case 2: smartastro::astrocore::spice_general_functions::vectorToSpicePlaneNVP(plane, planeP);
    break;
  case 3: smartastro::astrocore::spice_general_functions::vectorToSpicePlanePSV(plane, planeP);
    break;
  }
  
  pjelpl_c(&elinE, &planeP, &eloutE);

  smartastro::astrocore::spice_general_functions::spiceEllipseToVector(eloutE, elout);
}

void Celestial_Object::saelgv(const std::vector<double> &vec1, const std::vector<double> &vec2,  std::vector<double> &smajor,  std::vector<double> &sminor)
{
  return saelgv_c(&vec1[0], &vec2[0], &smajor[0], &sminor[0]);
}

void Celestial_Object::nplnpt(const std::vector<double> &linpt, const std::vector<double> &lindir, const std::vector<double> &point,  std::vector<double> &pnear, double &dist)
{
  return nplnpt_c(&linpt[0], &lindir[0], &point[0], &pnear[0], &dist);
}

void Celestial_Object::sincpt(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &dref, const std::vector<double> &dvec,  std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, int &found)
{
  return sincpt_c(method.c_str(), target.c_str(), et, fixref.c_str(), abcorr.c_str(), name.c_str(), dref.c_str(), &dvec[0], &spoint[0], &trgepc, &srfvec[0], &found);
}

void Celestial_Object::dskxv(const int &pri, const std::string &target, const int &nsurf, const std::vector<int> &srflst, const double &et, const std::string &fixref, const int &nrays, std::vector<double> &vtxarr, const std::vector<double> &dirarr, std::vector<double> &xptarr, std::vector<int> &fndarr)
{
   if(vtxarr.size() % 3 == 0 && dirarr.size() % 3 == 0 && xptarr.size() % 3 == 0){
     int noOfRowsv = vtxarr.size() / 3;
     double vtxarrA[noOfRowsv][3];
     int noOfRowsd = dirarr.size() / 3;
     double dirarrA[noOfRowsd][3];
     int noOfRowsx = xptarr.size() / 3;
     double xptarrA[noOfRowsx][3];

     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(vtxarr, vtxarrA);
     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(dirarr, dirarrA);
     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(xptarr, xptarrA);

     dskxv_c(pri, target.c_str(), nsurf, &srflst[0], et, fixref.c_str(), nrays, vtxarrA, dirarrA, xptarrA, &fndarr[0]);

     smartastro::astrocore::spice_general_functions::double2dArray3ToVector(vtxarrA, vtxarr);
   } else {
     smartastro_throw("Vectors are not the correct size for the arrays.");
   }
}

void Celestial_Object::dskxsi(const int &pri, const std::string &target, const int &nsurf, const std::vector<int> &srflst, const double &et, const std::string &fixref, const std::vector<double> &vertex, const std::vector<double> &raydir, const int &maxd, const int &maxi, std::vector<double> &xpt, int &handle, std::vector<int> &dladsc, std::vector<double> &dskdsc, std::vector<double> &dc, std::vector<int> &ic, int &found)
{
  SpiceDLADescr dladscD;
  SpiceDSKDescr dskdscD;

  smartastro::astrocore::spice_general_functions::vectorToSpiceDLADescr(dladsc, dladscD);
  smartastro::astrocore::spice_general_functions::vectorToSpiceDSKDescr(dskdsc, dskdscD);
  
  dskxsi_c(pri, target.c_str(), nsurf, &srflst[0], et, fixref.c_str(), &vertex[0], &raydir[0], maxd, maxi, &xpt[0], &handle, &dladscD, &dskdscD, &dc[0], &ic[0], &found);

  smartastro::astrocore::spice_general_functions::spiceDLADescrToVector(dladscD, dladsc);
  smartastro::astrocore::spice_general_functions::spiceDSKDescrToVector(dskdscD, dskdsc);
}

void Celestial_Object::ilumin(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn)
{
  return ilumin_c(method.c_str(), target.c_str(), et, fixref.c_str(), abcorr.c_str(), name.c_str(), &spoint[0], &trgepc, &srfvec[0], &phase, &incdnc, &emissn);
}

void Celestial_Object::illumg(const std::string &method, const std::string &target, const std::string &ilusrc, const double &et, const std::string &fixref, const std::string &abcorr, const std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn)
{
  return illumg_c(method.c_str(), target.c_str(), ilusrc.c_str(), et, fixref.c_str(), abcorr.c_str(), name.c_str(), &spoint[0], &trgepc, &srfvec[0], &phase, &incdnc, &emissn);
}

void Celestial_Object::illumf(const std::string &method, const std::string &target, const std::string &ilusrc, const double &et, const std::string &fixref, const std::string &abcorr, const std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn, int &visibl, int &lit)
{
  return illumf_c(method.c_str(), target.c_str(), ilusrc.c_str(), et, fixref.c_str(), abcorr.c_str(), name.c_str(), &spoint[0], &trgepc, &srfvec[0], &phase, &incdnc, &emissn, &visibl, &lit);
}

void Celestial_Object::limbpt(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &corloc, const std::vector<double> &refvec, const double &rolstp, const int &ncuts, const double &schstp, const double &soltol, const int &maxn,  std::vector<int> &npts,  std::vector<double> &points, std::vector<double> &epochs, std::vector<double> &tangts)
{
  if(points.size() % 3 == 0 && tangts.size() % 3 == 0){
     int noOfRowsp = points.size() / 3;
     double pointsA[noOfRowsp][3];
     int noOfRowst = tangts.size() / 3;
     double tangtsA[noOfRowst][3];

     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(points, pointsA);
     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(tangts, tangtsA);
     
     limbpt_c(method.c_str(), target.c_str(), et, fixref.c_str(), abcorr.c_str(), corloc.c_str(), name.c_str(), &refvec[0], rolstp, ncuts, schstp, soltol, maxn, &npts[0], pointsA, &epochs[0], tangtsA);
  
     smartastro::astrocore::spice_general_functions::double2dArray3ToVector(pointsA, points);
     smartastro::astrocore::spice_general_functions::double2dArray3ToVector(tangtsA, tangts);
  } else {
     smartastro_throw("Vectors are not the correct size for the arrays.");
  }
}

void Celestial_Object::termpt(const std::string &method, const std::string &ilusrc, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &corloc, const std::vector<double> &refvec, const double &rolstp, const int &ncuts, const double &schstp, const double &soltol, const int &maxn,  std::vector<int> &npts, std::vector<double> &points, std::vector<double> &epochs, std::vector<double> &trmvcs)
{
  if(points.size() % 3 == 0 && trmvcs.size() % 3 == 0){
     int noOfRowsp = points.size() / 3;
     double pointsA[noOfRowsp][3];
     int noOfRowst = trmvcs.size() / 3;
     double trmvcsA[noOfRowst][3];

     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(points, pointsA);
     smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(trmvcs, trmvcsA);
     
     termpt_c(method.c_str(), ilusrc.c_str(), target.c_str(), et, fixref.c_str(), abcorr.c_str(), corloc.c_str(), name.c_str(), &refvec[0], rolstp, ncuts, schstp, soltol, maxn, &npts[0], pointsA, &epochs[0], trmvcsA);

     smartastro::astrocore::spice_general_functions::double2dArray3ToVector(pointsA, points);
     smartastro::astrocore::spice_general_functions::double2dArray3ToVector(trmvcsA, trmvcs);
  } else {
     smartastro_throw("Vectors are not the correct size for the arrays.");
  }
}

void Celestial_Object::conics(const std::vector<double> &elts, const double &et, std::vector<double> &state)
{
  return conics_c(&elts[0], et, &state[0]);
}

void Celestial_Object::oscelt(const std::vector<double> &state, const double &et, std::vector<double> &elts)
{
  return oscelt_c(&state[0], et, mu, &elts[0]);
}

void Celestial_Object::occult(const std::string &targ1, const std::string &shape1, const std::string &frame1, const std::string &targ2, const std::string &shape2, const std::string &frame2, const std::string &abcorr, double &et, int &ocltid)
{
  return occult_c(targ1.c_str(), shape1.c_str(), frame1.c_str(), targ2.c_str(), shape2.c_str(), frame2.c_str(), abcorr.c_str(), name.c_str(), et, &ocltid);
}

void Celestial_Object::gfoclt(const std::string &occtyp, const std::string &front, const std::string &fshape, const std::string &fframe, const std::string &back, const std::string &bshape, const std::string &bframe, const std::string &abcorr, const double &step, std::vector<double> &cnfine, std::vector<double> &result)
{
  SPICEDOUBLE_CELL (cnfineCell, 100);
  SPICEDOUBLE_CELL (resultCell, 100);

  for(int i = 0; i < cnfine.size(); ++i){
    SPICE_CELL_SET_D (cnfine[i], i, &cnfineCell);
  }

  for(int i = 0; i < result.size(); ++i){
    SPICE_CELL_SET_D (result[i], i, &resultCell);
  }
  
  gfoclt_c(occtyp.c_str(), front.c_str(), fshape.c_str(), fframe.c_str(), back.c_str(), bshape.c_str(), bframe.c_str(), abcorr.c_str(), name.c_str(), step, &cnfineCell, &resultCell);

  for(int i = 0; i < cnfine.size(); ++i){
    SPICE_CELL_GET_D (&cnfineCell, i, &cnfine[i]);
  }

  for(int i = 0; i < result.size(); ++i){
    SPICE_CELL_GET_D (&resultCell, i, &result[i]);
  }
}

void Celestial_Object::srfcss(const int &code, const int srflen, std::string &srfstr, int &isname)
{
  return srfcss_c(code, name.c_str(), srflen, &srfstr.at(0), &isname);
}

void Celestial_Object::srfs2c(const std::string &srfstr, int &code, int &found)
{
  return srfs2c_c(srfstr.c_str(), name.c_str(), &code, &found);
}

void Celestial_Object::srfc2s(const int &code, const int srflen, std::string &srfstr, int &isname)
{
  srfstr.push_back('\0');
  
  srfc2s_c(code, id, srflen, &srfstr.at(0), &isname);

  srfstr.pop_back();
}

void Celestial_Object::srfscc(const std::string &srfstr, int &code, int &found)
{
  return srfscc_c(srfstr.c_str(), id, &code, &found);
}

void Celestial_Object::srfrec(const double &longitude, const double &latitude,  std::vector<double> &rectan)
{
  return srfrec_c(id, longitude, latitude, &rectan[0]);
}
