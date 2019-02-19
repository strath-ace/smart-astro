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

Celestial_Object::Celestial_Object(std::string givenName, std::vector<double> &givenPositn, double givenMu) : Astro_Body(givenName)
{
  positn = givenPositn;
  mu = givenMu;
}

Celestial_Object::Celestial_Object(int givenId,  std::vector<double> &givenPositn, double givenMu) : Astro_Body(givenId)
{
  positn = givenPositn;
  mu = givenMu;
}



/*int Celestial_Object::latsrf_(const std::string &method, const std::string &target, double &et, const std::string &fixref, int &npts, const std::vector<std::vector::<double>> &lonlat, std::vector<std::vector::<double>> &srfpts)
{
  return latsrf_(&method.at(0), &target.at(0), et, &fixref.at(0), npts, setupHMM(lonlat, lonlat.size(), lonlat[0].size()), setupHMM(srfpts, srfpts.size(), srfpts[0].size()));
}

int Celestial_Object::srfnrm_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &npts,  std::vector<double> &srfpts,  std::vector<double> &normls, int method_len, int target_len, int fixref_len)
{
  return srfnrm_(methodA, targetA, etA, fixrefA, nptsA, srfptsA, normlsA,  method_len, target_len, fixref_len);
  }*/

void Celestial_Object::nearpt(double &a, double &b, double &c, std::vector<double> &npoint, double &alt)
{
  return nearpt_c(&positn[0], a, b, c, &npoint[0], &alt);
}

void Celestial_Object::surfpt(const std::vector<double> &u, double &a, double &b, double &c,  std::vector<double> &point, int &found)
{
  return surfpt_c(&positn.at(0), &u[0], a, b, c, &point[0], &found);
}

void Celestial_Object::surfnm(double &a, double &b, double &c, const std::vector<double> &point,  std::vector<double> &normal)
{
  return surfnm_c(a, b, c, &point[0], &normal[0]);
}

/*int Celestial_Object::edlimb_(double &a, double &b, double &c, const std::vector<double> &viewpt,  std::vector<double> &limb)
{
  double *aA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(a);
  double *bA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(b);
  double *c__A = smartastro::astrocore::spice_general_functions::doubleVectorToArray(c__);
  double *viewptA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(viewpt);
  double *limbA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(limb);
  
  return edlimb_(a, b, c__, viewpt, limb);
}

int Celestial_Object::inelpl_( std::vector<double> &ellips,  std::vector<double> &plane,  std::vector<int> &nxpts,  std::vector<double> &xpt1,  std::vector<double> &xpt2)
{
  double *ellipsA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(ellips);
  double *planeA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(plane);
  double *c__A = smartastro::astrocore::spice_general_functions::doubleVectorToArray(nxpts);
  double *viewptA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(xpt1);
  double *limbA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(xpt2);
  
  return inelpl_(ellipsA, planeA, nxptsA, xpt1A, xpt2A);
  }*/

void Celestial_Object::npedln(double &a, double &b, double &c, const std::vector<double> &linept, const std::vector<double> &linedr, std::vector<double> &pnear, double &dist)
{
  return npedln_c(a, b, c, &linept[0], &linedr[0], &pnear[0], &dist);
}

/*int Celestial_Object::inrypl_( std::vector<double> &vertex,  std::vector<double> &dir,  std::vector<double> &plane,  std::vector<int> &nxpts,  std::vector<double> &xpt)
{
  double *vertexA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(vertex);
  double *dirA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(dir);
  double *planeA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(plane);
  double *nxptsA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(nxpts);
  double *xptA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(xpt);
  
  return inrypl_(vertexA, dirA, planeA, nxptsA, xptA);
}

int Celestial_Object::pjelpl_( std::vector<double> &elin,  std::vector<double> &plane,  std::vector<double> &elout)
{
  double *elinA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(elin);
  double *planeA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(plane);
  double *eloutA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(elout);
  
  return pjelpl_(elinA, planeA, eloutA);
  }*/

void Celestial_Object::saelgv(const std::vector<double> &vec1, const std::vector<double> &vec2,  std::vector<double> &smajor,  std::vector<double> &sminor)
{
  return saelgv_c(&vec1[0], &vec2[0], &smajor[0], &sminor[0]);
}

void Celestial_Object::nplnpt(const std::vector<double> &linpt, const std::vector<double> &lindir, const std::vector<double> &point,  std::vector<double> &pnear, double &dist)
{
  return nplnpt_c(&linpt[0], &lindir[0], &point[0], &pnear[0], &dist);
}

void Celestial_Object::sincpt(std::string &method, std::string &target, double &et, std::string &fixref, std::string &abcorr, std::string &dref, const std::vector<double> &dvec,  std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, int &found)
{
  return sincpt_c(&method.at(0), &target.at(0), et, &fixref.at(0), &abcorr.at(0), &name.at(0), &dref.at(0), &dvec[0], &spoint[0], &trgepc, &srfvec[0], &found);
}

/*int Celestial_Object::dskxv_( std::vector<int> &pri,  std::string &target,  std::vector<int> &nsurf,  std::vector<int> &srflst,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &nrays,  std::vector<double> &vtxarr,  std::vector<double> &dirarr,  std::vector<double> &xptarr,  std::vector<int> &fndarr, int target_len, int fixref_len)
{
  int *priA = smartastro::astrocore::spice_general_functions::intVectorToArray(pri);
  char *targetA = smartastro::astrocore::spice_general_functions::stringToCharArray(target);
  int *nsurfA = smartastro::astrocore::spice_general_functions::intVectorToArray(nsurf);
  int *srflstA = smartastro::astrocore::spice_general_functions::intVectorToArray(srflst);
  double *etA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(et);
  char *fixrefA = smartastro::astrocore::spice_general_functions::stringToCharArray(fixref);
  int *nraysA = smartastro::astrocore::spice_general_functions::intVectorToArray(nrays);
  double *vtxarrA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(vtxarr);
  double *dirarrA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(dirarr);
  double *xptarrA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(xptarr);
  int *fndarrA = smartastro::astrocore::spice_general_functions::intVectorToArray(fndarr);

  return dskxv_(priA, targetA, nsurfA, srflstA, etA, fixrefA, nraysA, vtxarrA, dirarrA, xptarrA, fndarrA, target_len, fixref_len);
}

int Celestial_Object::dskxsi_( std::vector<int> &pri,  std::string &target,  std::vector<int> &nsurf,  std::vector<int> &srflst,  std::vector<double> &et,  std::string &fixref,  std::vector<double> &vertex,  std::vector<double> &raydir,  std::vector<int> &maxd,  std::vector<int> &maxi,  std::vector<double> &xpt,  std::vector<int> &handle,  std::vector<int> &dladsc,  std::vector<double> &dskdsc,  std::vector<double> &dc,  std::vector<int> &ic,  std::vector<int> &found, int target_len, int fixref_len)
{
  int *priA = smartastro::astrocore::spice_general_functions::intVectorToArray(pri);
  char *targetA = smartastro::astrocore::spice_general_functions::stringToCharArray(target);
  int *nsurfA = smartastro::astrocore::spice_general_functions::intVectorToArray(nsurf);
  int *srflstA = smartastro::astrocore::spice_general_functions::intVectorToArray(srflst);
  double *etA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(et);
  char *fixrefA = smartastro::astrocore::spice_general_functions::stringToCharArray(fixref);
  double *vertexA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(vertex);
  double *raydirA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(raydir);
  int *maxdA = smartastro::astrocore::spice_general_functions::intVectorToArray(maxd);
  int *maxiA = smartastro::astrocore::spice_general_functions::intVectorToArray(maxi);
  double *xptA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(xpt);
  int *handleA = smartastro::astrocore::spice_general_functions::intVectorToArray(handle);
  int *dladscA = smartastro::astrocore::spice_general_functions::intVectorToArray(dladsc);
  double *dskdscA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(dskdsc);
  double *dcA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(dc);
  int *icA = smartastro::astrocore::spice_general_functions::intVectorToArray(ic);
  int *foundA = smartastro::astrocore::spice_general_functions::intVectorToArray(found);

  return dskxsi_(priA, targetA, nsurfA, srflstA, etA, fixrefA, vertexA, raydirA, maxdA, maxiA, xptA, handleA, dladscA, dskdscA, dcA, icA, foundA, target_len, fixref_len);
  }*/

void Celestial_Object::ilumin(std::string &method, std::string &target, double &et, std::string &fixref, std::string &abcorr, const std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn)
{
  return ilumin_c(&method.at(0), &target.at(0), et, &fixref.at(0), &abcorr.at(0), &name.at(0), &spoint[0], &trgepc, &srfvec[0], &phase, &incdnc, &emissn);
}

void Celestial_Object::illumg(std::string &method, std::string &target, std::string &ilusrc, double &et, std::string &fixref, std::string &abcorr, std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn)
{
  return illumg_c(&method.at(0), &target.at(0), &ilusrc.at(0), et, &fixref.at(0), &abcorr.at(0), &name.at(0), &spoint[0], &trgepc, &srfvec[0], &phase, &incdnc, &emissn);
}

void Celestial_Object::illumf(std::string &method, std::string &target, std::string &ilusrc, double &et, std::string &fixref, std::string &abcorr, std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn, int &visibl, int &lit)
{
  return illumf_c(&method.at(0), &target.at(0), &ilusrc.at(0), et, &fixref.at(0), &abcorr.at(0), &name.at(0), &spoint[0], &trgepc, &srfvec[0], &phase, &incdnc, &emissn, &visibl, &lit);
}

/*int Celestial_Object::limbpt_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &corloc,  std::vector<double> &refvec,  std::vector<double> &rolstp,  std::vector<int> &ncuts,  std::vector<double> &schstp,  std::vector<double> &soltol,  std::vector<int> &maxn,  std::vector<int> &npts,  std::vector<double> &points,  std::vector<double> &epochs,  std::vector<double> &tangts, int method_len, int target_len, int fixref_len, int abcorr_len, int corloc_len, int obsrvr_len)
{
  char *methodA = smartastro::astrocore::spice_general_functions::stringToCharArray(method);
  char *targetA = smartastro::astrocore::spice_general_functions::stringToCharArray(target);
  double *etA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(et);
  char *fixrefA = smartastro::astrocore::spice_general_functions::stringToCharArray(fixref);
  char *abcorrA = smartastro::astrocore::spice_general_functions::stringToCharArray(abcorr);
  char *corlocA = smartastro::astrocore::spice_general_functions::stringToCharArray(corloc);
  double *refvecA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(refvec);
  double *rolstpA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(rolstp);
  int *ncutsA = smartastro::astrocore::spice_general_functions::intVectorToArray(ncuts);
  double *schstpA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(schstp);
  double *soltolA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(soltol);
  int *maxnA = smartastro::astrocore::spice_general_functions::intVectorToArray(maxn);
  int *nptsA = smartastro::astrocore::spice_general_functions::intVectorToArray(npts);
  double *pointsA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(points);
  double *epochsA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(epochs);
  double *tangtsA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(tangts);
  
  return limbpt_(methodA, targetA, etA, fixrefA, abcorrA, corlocA, name, refvecA, rolstpA, ncutsA, schstpA, soltolA, maxnA, nptsA, pointsA, epochsA, tangtsA, method_len, target_len, fixref_len, abcorr_len, corloc_len, obsrvr_len);
}

int Celestial_Object::termpt_( std::string &method,  std::string &ilusrc,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &corloc,  std::vector<double> &refvec,  std::vector<double> &rolstp,  std::vector<int> &ncuts,  std::vector<double> &schstp,  std::vector<double> &soltol,  std::vector<int> &maxn,  std::vector<int> &npts,  std::vector<double> &points,  std::vector<double> &epochs,  std::vector<double> &trmvcs, int method_len, int ilusrc_len, int target_len, int fixref_len, int abcorr_len, int corloc_len, int obsrvr_len)
{
  char *methodA = smartastro::astrocore::spice_general_functions::stringToCharArray(method);
  char *ilusrcA = smartastro::astrocore::spice_general_functions::stringToCharArray(ilusrc);
  char *targetA = smartastro::astrocore::spice_general_functions::stringToCharArray(target);
  double *etA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(et);
  char *fixrefA = smartastro::astrocore::spice_general_functions::stringToCharArray(fixref);
  char *abcorrA = smartastro::astrocore::spice_general_functions::stringToCharArray(abcorr);
  char *corlocA = smartastro::astrocore::spice_general_functions::stringToCharArray(corloc);
  double *refvecA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(refvec);
  double *rolstpA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(rolstp);
  int *ncutsA = smartastro::astrocore::spice_general_functions::intVectorToArray(ncuts);
  double *schstpA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(schstp);
  double *soltolA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(soltol);
  int *maxnA = smartastro::astrocore::spice_general_functions::intVectorToArray(maxn);
  int *nptsA = smartastro::astrocore::spice_general_functions::intVectorToArray(npts);
  double *pointsA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(points);
  double *epochsA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(epochs);
  double *trmvcsA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(trmvcs);
  
  return termpt_(methodA, ilusrcA, targetA, etA, fixrefA, abcorrA, corlocA, name, refvecA, rolstpA, ncutsA, schstpA, soltolA, maxnA, nptsA, pointsA, epochsA, trmvcsA, method_len, ilusrc_len, target_len, fixref_len, abcorr_len, corloc_len, obsrvr_len);
  }*/

void Celestial_Object::conics(const std::vector<double> &elts, double &et,  std::vector<double> &state)
{
  return conics_c(&elts[0], et, &state[0]);
}

void Celestial_Object::oscelt(const std::vector<double> &state, double &et, std::vector<double> &elts)
{
  return oscelt_c(&state[0], et, mu, &elts[0]);
}

void Celestial_Object::occult(std::string &targ1, std::string &shape1, std::string &frame1, std::string &targ2, std::string &shape2, std::string &frame2, std::string &abcorr, double &et, int &ocltid)
{
  return occult_c(&targ1.at(0), &shape1.at(0), &frame1.at(0), &targ2.at(0), &shape2.at(0), &frame2.at(0), &abcorr.at(0), &name.at(0), et, &ocltid);
}

/*int Celestial_Object::gfoclt_( std::string &occtyp,  std::string &front,  std::string &fshape,  std::string &fframe,  std::string &back,  std::string &bshape,  std::string &bframe,  std::string &abcorr,  std::vector<double> &step,  std::vector<double> &cnfine,  std::vector<double> &result, int occtyp_len, int front_len, int fshape_len, int fframe_len, int back_len, int bshape_len, int bframe_len, int abcorr_len, int obsrvr_len)
{
  char *occtypA = smartastro::astrocore::spice_general_functions::stringToCharArray(occtyp);
  char *frontA = smartastro::astrocore::spice_general_functions::stringToCharArray(front);
  char *fshapeA = smartastro::astrocore::spice_general_functions::stringToCharArray(fshape);
  char *fframeA = smartastro::astrocore::spice_general_functions::stringToCharArray(fframe);
  char *backA = smartastro::astrocore::spice_general_functions::stringToCharArray(back);
  char *bshapeA = smartastro::astrocore::spice_general_functions::stringToCharArray(bshape);
  char *bframeA = smartastro::astrocore::spice_general_functions::stringToCharArray(bframe);
  char *abcorrA = smartastro::astrocore::spice_general_functions::stringToCharArray(abcorr);
  double *stepA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(step);
  double *cnfineA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(cnfine);
  double *resultA = smartastro::astrocore::spice_general_functions::doubleVectorToArray(result);
  
  return gfoclt_(occtypA, frontA, fshapeA, fframeA, backA, bshapeA, bframeA, abcorrA, name, stepA, cnfineA, resultA, occtyp_len, front_len, fshape_len, fframe_len, back_len, bshape_len, bframe_len, abcorr_len, obsrvr_len);
  }*/

void Celestial_Object::srfcss(int &code, int srflen, std::string &srfstr, int &isname)
{
  return srfcss_c(code, &name.at(0), srflen, &srfstr.at(0), &isname);
}

void Celestial_Object::srfs2c(std::string &srfstr, int &code, int &found)
{
  return srfs2c_c(&srfstr.at(0), &name.at(0), &code, &found);
}

void Celestial_Object::srfc2s(int &code, int srflen, std::string &srfstr, int &isname)
{
  return srfc2s_c(code, id, srflen, &srfstr.at(0), &isname);
}

void Celestial_Object::srfscc(std::string &srfstr, int &code, int &found)
{
  return srfscc_c(&srfstr.at(0), id, &code, &found);
}

void Celestial_Object::srfrec(double &longitude, double &latitude,  std::vector<double> &rectan)
{
  return srfrec_c(id, longitude, latitude, &rectan[0]);
}
