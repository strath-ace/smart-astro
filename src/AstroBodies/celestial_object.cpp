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

Celestial_Object::Celestial_Object( std::string givenName, int givenId,  std::vector<double> &givenPositn,  std::vector<double> &givenMu)
{
  Astro_Body(givenName, givenId);
  positn = givenPositn;
  mu = givenMu;
}

Celestial_Object::~Celestial_Object()
{
}

int Celestial_Object::latsrf_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &npts,  std::vector<double> &lonlat,  std::vector<double> &srfpts, int method_len, int target_len, int fixref_len)
{
  return latsrf_(method, target, et, fixref, npts, lonlat, srfpts, method_len, target_len, fixref_len);
}

int Celestial_Object::srfnrm_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &npts,  std::vector<double> &srfpts,  std::vector<double> &normls, int method_len, int target_len, int fixref_len)
{
  return srfnrm_(method, target, et, fixref, npts, srfpts, normls,  method_len, target_len, fixref_len);
}

int Celestial_Object::nearpt_( std::vector<double> &positn,  std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &npoint,  std::vector<double> &alt)
{
  return nearpt_(positn, a, b, c__, npoint, alt);
}

int Celestial_Object::dnearp_( std::vector<double> &state,  std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &dnear,  std::vector<double> &dalt,  std::vector<int> &found)
{
  return dnearp_(state, a, b, c__, dnear, dalt, found);
}

int Celestial_Object::surfpt_( std::vector<double> &u,  std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &point,  std::vector<int> &found)
{
  return surfpt_(positn, u, a, b, c__, point, found);
}

int Celestial_Object::surfnm_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &point,  std::vector<double> &normal)
{
  return surfnm_(a, b, c__, point, normal);
}

int Celestial_Object::ednmpt_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &normal,  std::vector<double> &point)
{
  return ednmpt_(a, b, c__, normal, point);
}

int Celestial_Object::edpnt_( std::vector<double> &p,  std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &ep)
{
  return edpnt_(p, a, b, c__, ep);
}

int Celestial_Object::edlimb_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &viewpt,  std::vector<double> &limb)
{
  return edlimb_(a, b, c__, viewpt, limb);
}

int Celestial_Object::inelpl_( std::vector<double> &ellips,  std::vector<double> &plane,  std::vector<int> &nxpts,  std::vector<double> &xpt1,  std::vector<double> &xpt2)
{
  return inelpl_(ellips, plane, nxpts, xpt1, xpt2);
}

int Celestial_Object::npedln_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &linept,  std::vector<double> &linedr,  std::vector<double> &pnear,  std::vector<double> &dist)
{
  return npedln_(a, b, c__, linept, linedr, pnear, dist);
}

int Celestial_Object::inrypl_( std::vector<double> &vertex,  std::vector<double> &dir,  std::vector<double> &plane,  std::vector<int> &nxpts,  std::vector<double> &xpt)
{
  return inrypl_(vertex, dir, plane, nxpts, xpt);
}

int Celestial_Object::pjelpl_( std::vector<double> &elin,  std::vector<double> &plane,  std::vector<double> &elout)
{
  return pjelpl_(elin, plane, elout);
}

int Celestial_Object::saelgv_( std::vector<double> &vec1,  std::vector<double> &vec2,  std::vector<double> &smajor,  std::vector<double> &sminor)
{
  return saelgv_(vec1, vec2, smajor, sminor);
}

int Celestial_Object::nplnpt_( std::vector<double> &linpt,  std::vector<double> &lindir,  std::vector<double> &point,  std::vector<double> &pnear,  std::vector<double> &dist)
{
  return nplnpt_(linpt, lindir, point, pnear, dist);
}

int Celestial_Object::npsgpt_( std::vector<double> &ep1,  std::vector<double> &ep2,  std::vector<double> &point,  std::vector<double> &pnear,  std::vector<double> &dist)
{
  return npsgpt_(ep1, ep2, point, pnear, dist);
}

int Celestial_Object::sincpt_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &dref,  std::vector<double> &dvec,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec,  std::vector<int> &found, int method_len, int target_len, int fixref_len, int abcorr_len, int obsrvr_len, int dref_len)
{
  return sincpt_(method, target, et, fixref, abcorr, name, dref, dvec, spoint, trgepc, srfvec, found, method_len, target_len, fixref_len, abcorr_len, obsrvr_len, dref_len);
}

int Celestial_Object::dskxv_( std::vector<int> &pri,  std::string &target,  std::vector<int> &nsurf,  std::vector<int> &srflst,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &nrays,  std::vector<double> &vtxarr,  std::vector<double> &dirarr,  std::vector<double> &xptarr,  std::vector<int> &fndarr, int target_len, int fixref_len)
{
  return dskxv_(pri, target, nsurf, srflst, et, fixref, nrays, vtxarr, dirarr, xptarr, fndarr, target_len, fixref_len);
}

int Celestial_Object::dskxsi_( std::vector<int> &pri,  std::string &target,  std::vector<int> &nsurf,  std::vector<int> &srflst,  std::vector<double> &et,  std::string &fixref,  std::vector<double> &vertex,  std::vector<double> &raydir,  std::vector<int> &maxd,  std::vector<int> &maxi,  std::vector<double> &xpt,  std::vector<int> &handle,  std::vector<int> &dladsc,  std::vector<double> &dskdsc,  std::vector<double> &dc,  std::vector<int> &ic,  std::vector<int> &found, int target_len, int fixref_len)
{
  return dskxsi_(pri, target, nsurf, srflst, et, fixref, vertex, raydir, maxd, maxi, xpt, handle, dladsc, dskdsc, dc, ic, found, target_len, fixref_len);
}

int Celestial_Object::ilumin_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec,  std::vector<double> &phase,  std::vector<double> &solar,  std::vector<double> &emissn, int method_len, int target_len, int fixref_len, int abcorr_len, int obsrvr_len)
{
  return ilumin_(method, target, et, fixref, abcorr, name, spoint, trgepc, srfvec, phase, solar, emissn, method_len, target_len, fixref_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::illumg_( std::string &method,  std::string &target,  std::string &ilusrc,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec,  std::vector<double> &phase,  std::vector<double> &incdnc,  std::vector<double> &emissn, int method_len, int target_len, int ilusrc_len, int fixref_len, int abcorr_len, int obsrvr_len)
{
  return illumg_(method, target, ilusrc, et, fixref, abcorr, name, spoint, trgepc, srfvec, phase, incdnc, emissn, method_len, target_len, ilusrc_len, fixref_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::illumf_( std::string &method,  std::string &target,  std::string &ilusrc,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec,  std::vector<double> &phase,  std::vector<double> &incdnc,  std::vector<double> &emissn,  std::vector<int> &visibl,  std::vector<int> &lit, int method_len, int target_len, int ilusrc_len, int fixref_len, int abcorr_len, int obsrvr_len)
{
  return illumf_(method, target, ilusrc, et, fixref, abcorr, name, spoint, trgepc, srfvec, phase, incdnc, emissn, visibl, lit, method_len, target_len, ilusrc_len, fixref_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::limbpt_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &corloc,  std::vector<double> &refvec,  std::vector<double> &rolstp,  std::vector<int> &ncuts,  std::vector<double> &schstp,  std::vector<double> &soltol,  std::vector<int> &maxn,  std::vector<int> &npts,  std::vector<double> &points,  std::vector<double> &epochs,  std::vector<double> &tangts, int method_len, int target_len, int fixref_len, int abcorr_len, int corloc_len, int obsrvr_len)
{
  return limbpt_(method, target, et, fixref, abcorr, corloc, name, refvec, rolstp, ncuts, schstp, soltol, maxn, npts, points, epochs, tangts, method_len, target_len, fixref_len, abcorr_len, corloc_len, obsrvr_len);
}

int Celestial_Object::termpt_( std::string &method,  std::string &ilusrc,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &corloc,  std::vector<double> &refvec,  std::vector<double> &rolstp,  std::vector<int> &ncuts,  std::vector<double> &schstp,  std::vector<double> &soltol,  std::vector<int> &maxn,  std::vector<int> &npts,  std::vector<double> &points,  std::vector<double> &epochs,  std::vector<double> &trmvcs, int method_len, int ilusrc_len, int target_len, int fixref_len, int abcorr_len, int corloc_len, int obsrvr_len)
{
  return termpt_(method, ilusrc, target, et, fixref, abcorr, corloc, name, refvec, rolstp, ncuts, schstp, soltol, maxn, npts, points, epochs, trmvcs, method_len, ilusrc_len, target_len, fixref_len, abcorr_len, corloc_len, obsrvr_len);
}

int Celestial_Object::conics_( std::vector<double> &elts,  std::vector<double> &et,  std::vector<double> &state)
{
  return conics_(elts, et, state);
}

int Celestial_Object::oscelt_( std::vector<double> &state,  std::vector<double> &et,  std::vector<double> &elts)
{
  return oscelt_(state, et, mu, elts);
}

int Celestial_Object::occult_( std::string &targ1,  std::string &shape1,  std::string &frame1,  std::string &targ2,  std::string &shape2,  std::string &frame2,  std::string &abcorr,  std::vector<double> &et,  std::vector<int> &ocltid, int targ1_len, int shape1_len, int frame1_len, int targ2_len, int shape2_len, int frame2_len, int abcorr_len, int obsrvr_len)
{
  return occult_(targ1, shape1, frame1, targ2, shape2, frame2, abcorr, name, et, ocltid, targ1_len, shape1_len, frame1_len, targ2_len, shape2_len, frame2_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::gfoclt_( std::string &occtyp,  std::string &front,  std::string &fshape,  std::string &fframe,  std::string &back,  std::string &bshape,  std::string &bframe,  std::string &abcorr,  std::vector<double> &step,  std::vector<double> &cnfine,  std::vector<double> &result, int occtyp_len, int front_len, int fshape_len, int fframe_len, int back_len, int bshape_len, int bframe_len, int abcorr_len, int obsrvr_len)
{
  return gfoclt_(occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, name, step, cnfine, result, occtyp_len, front_len, fshape_len, fframe_len, back_len, bshape_len, bframe_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::srfcss_( std::vector<int> &code,  std::string &srfstr,  std::vector<int> &isname, int bodstr_len, int srfstr_len)
{
  return srfcss_(code, name, srfstr, isname, bodstr_len, srfstr_len);
}

int Celestial_Object::srfs2c_( std::string &srfstr,  std::vector<int> &code,  std::vector<int> &found, int srfstr_len, int bodstr_len)
{
  return srfs2c_(srfstr, name, code, found, srfstr_len, bodstr_len);
}

int Celestial_Object::srfc2s_( std::vector<int> &code,  std::string &srfstr,  std::vector<int> &isname, int srfstr_len)
{
  return srfc2s_(code, id, srfstr, isname, srfstr_len);
}

int Celestial_Object::srfscc_( std::string &srfstr,  std::vector<int> &code,  std::vector<int> &found, int srfstr_len)
{
  return srfscc_(srfstr, id, code, found, srfstr_len);
}

int Celestial_Object::bodc2n_( std::string &name__,  std::vector<int> &found, int name_len)
{
  return bodc2n_(id, name__, found, name_len);
}

int Celestial_Object::srfrec_( std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &rectan)
{
  return srfrec_(id, long__, lat, rectan);
}
