#include "AstroBodies/celestial_object.h"

using namespace smartastro;
using namespace smartastro::astrobodies;

Celestial_Object::Celestial_Object()
{
}

Celestial_Object::~Celestial_Object()
{
}

int Celestial_Object::latsrf_(char *method, char *target, doublereal *et, char *fixref, integer *npts, doublereal *lonlat, doublereal *srfpts, ftnlen method_len, ftnlen target_len, ftnlen fixref_len)
{
  return latsrf_(method, target, et, fixref, npts, lonlat, srfpts, method_len, target_len, fixref_len);
}

int Celestial_Object::srfnrm_(char *method, char *target, doublereal *et, char *fixref, integer *npts, doublereal *srfpts, doublereal *normls, ftnlen method_len, ftnlen target_len, ftnlen fixref_len)
{
  return srfnrm_(method, target, et, fixref, npts, srfpts, normls,  method_len, target_len, fixref_len);
}

int Celestial_Object::nearpt_(doublereal *positn, doublereal *a, doublereal *b, doublereal *c__, doublereal *npoint, doublereal *alt)
{
  return nearpt_(positn, a, b, c__, npoint, alt);
}

int Celestial_Object::dnearp_(doublereal *state, doublereal *a, doublereal *b, doublereal *c__, doublereal *dnear, doublereal *dalt, logical *found)
{
  return dnearp_(state, a, b, c__, dnear, dalt, found);
}

int Celestial_Object::surfpt_(doublereal *positn, doublereal *u, doublereal *a, doublereal *b, doublereal *c__, doublereal *point, logical *found)
{
  return surfpt_(positn, u, a, b, c__, point, found);
}

int Celestial_Object::surfnm_(doublereal *a, doublereal *b, doublereal *c__, doublereal *point, doublereal *normal)
{
  return surfnm_(a, b, c__, point, normal);
}

int Celestial_Object::ednmpt_(doublereal *a, doublereal *b, doublereal *c__, doublereal *normal, doublereal *point)
{
  return ednmpt_(a, b, c__, normal, point);
}

int Celestial_Object::edpnt_(doublereal *p, doublereal *a, doublereal *b, doublereal *c__, doublereal *ep)
{
  return edpnt_(p, a, b, c__, ep);
}

int Celestial_Object::edlimb_(doublereal *a, doublereal *b, doublereal *c__, doublereal *viewpt, doublereal *limb)
{
  return edlimb_(a, b, c__, viewpt, limb);
}

int Celestial_Object::inelpl_(doublereal *ellips, doublereal *plane, integer *nxpts, doublereal *xpt1, doublereal *xpt2)
{
  return inelpl_(ellips, plane, nxpts, xpt1, xpt2);
}

int Celestial_Object::npedln_(doublereal *a, doublereal *b, doublereal *c__, doublereal *linept, doublereal *linedr, doublereal *pnear, doublereal *dist)
{
  return npedln_(a, b, c__, linept, linedr, pnear, dist);
}

int Celestial_Object::inrypl_(doublereal *vertex, doublereal *dir, doublereal *plane, integer *nxpts, doublereal *xpt)
{
  return inrypl_(vertex, dir, plane, nxpts, xpt);
}

int Celestial_Object::pjelpl_(doublereal *elin, doublereal *plane, doublereal *elout)
{
  return pjelpl_(elin, plane, elout);
}

int Celestial_Object::saelgv_(doublereal *vec1, doublereal *vec2, doublereal *smajor, doublereal *sminor)
{
  return saelgv_(vec1, vec2, smajor, sminor);
}

int Celestial_Object::nplnpt_(doublereal *linpt, doublereal *lindir, doublereal *point, doublereal *pnear, doublereal *dist)
{
  return nplnpt_(linpt, lindir, point, pnear, dist);
}

int Celestial_Object::npsgpt_(doublereal *ep1, doublereal *ep2, doublereal *point, doublereal *pnear, doublereal *dist)
{
  return npsgpt_(ep1, ep2, point, pnear, dist);
}

int Celestial_Object::sincpt_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *obsrvr, char *dref, doublereal *dvec, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, logical *found, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len, ftnlen dref_len)
{
  return sincpt_(method, target, et, fixref, abcorr, obsrvr, dref, dvec, spoint, trgepc, srfvec, found, method_len, target_len, fixref_len, abcorr_len, obsrvr_len, dref_len);
}

int Celestial_Object::dskxv_(logical *pri, char *target, integer *nsurf, integer *srflst, doublereal *et, char *fixref, integer *nrays, doublereal *vtxarr, doublereal *dirarr, doublereal *xptarr, logical *fndarr, ftnlen target_len, ftnlen fixref_len)
{
  return dskxv_(pri, target, nsurf, srflst, et, fixref, nrays, vtxarr, dirarr, xptarr, fndarr, target_len, fixref_len);
}

int Celestial_Object::dskxsi_(logical *pri, char *target, integer *nsurf, integer *srflst, doublereal *et, char *fixref, doublereal *vertex, doublereal *raydir, integer *maxd, integer *maxi, doublereal *xpt, integer *handle, integer *dladsc, doublereal *dskdsc, doublereal *dc, integer *ic, logical *found, ftnlen target_len, ftnlen fixref_len)
{
  return dskxsi_(pri, target, nsurf, srflst, et, fixref, vertex, raydir, maxd, maxi, xpt, handle, dladsc, dskdsc, dc, ic, found, target_len, fixref_len);
}

int Celestial_Object::ilumin_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, doublereal *phase, doublereal *solar, doublereal *emissn, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return ilumin_(method, target, et, fixref, abcorr, obsrvr, spoint, trgepc, srfvec, phase, solar, emissn, method_len, target_len, fixref_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::illumg_(char *method, char *target, char *ilusrc, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, doublereal *phase, doublereal *incdnc, doublereal *emissn, ftnlen method_len, ftnlen target_len, ftnlen ilusrc_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return illumg_(method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint, trgepc, srfvec, phase, incdnc, emissn, method_len, target_len, ilusrc_len, fixref_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::illumf_(char *method, char *target, char *ilusrc, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, doublereal *phase, doublereal *incdnc, doublereal *emissn, logical *visibl, logical *lit, ftnlen method_len, ftnlen target_len, ftnlen ilusrc_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return illumf_(method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint, trgepc, srfvec, phase, incdnc, emissn, visibl, lit, method_len, target_len, ilusrc_len, fixref_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::limbpt_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *corloc, char *obsrvr, doublereal *refvec, doublereal *rolstp, integer *ncuts, doublereal *schstp, doublereal *soltol, integer *maxn, integer *npts, doublereal *points, doublereal *epochs, doublereal *tangts, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen corloc_len, ftnlen obsrvr_len)
{
  return limbpt_(method, target, et, fixref, abcorr, corloc, obsrvr, refvec, rolstp, ncuts, schstp, soltol, maxn, npts, points, epochs, tangts, method_len, target_len, fixref_len, abcorr_len, corloc_len, obsrvr_len);
}

int Celestial_Object::termpt_(char *method, char *ilusrc, char *target, doublereal *et, char *fixref, char *abcorr, char *corloc, char *obsrvr, doublereal *refvec, doublereal *rolstp, integer *ncuts, doublereal *schstp, doublereal *soltol, integer *maxn, integer *npts, doublereal *points, doublereal *epochs, doublereal *trmvcs, ftnlen method_len, ftnlen ilusrc_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen corloc_len, ftnlen obsrvr_len)
{
  return termpt_(method, ilusrc, target, et, fixref, abcorr, corloc, obsrvr, refvec, rolstp, ncuts, schstp, soltol, maxn, npts, points, epochs, trmvcs, method_len, ilusrc_len, target_len, fixref_len, abcorr_len, corloc_len, obsrvr_len);
}

int Celestial_Object::conics_(doublereal *elts, doublereal *et, doublereal *state)
{
  return conics_(elts, et, state);
}

int Celestial_Object::oscelt_(doublereal *state, doublereal *et, doublereal *mu, doublereal *elts)
{
  return oscelt_(state, et, mu, elts);
}

int Celestial_Object::occult_(char *targ1, char *shape1, char *frame1, char *targ2, char *shape2, char *frame2, char *abcorr, char *obsrvr, doublereal *et, integer *ocltid, ftnlen targ1_len, ftnlen shape1_len, ftnlen frame1_len, ftnlen targ2_len, ftnlen shape2_len, ftnlen frame2_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return occult_(targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, et, ocltid, targ1_len, shape1_len, frame1_len, targ2_len, shape2_len, frame2_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::gfoclt_(char *occtyp, char *front, char *fshape, char *fframe, char *back, char *bshape, char *bframe, char *abcorr, char *obsrvr, doublereal *step, doublereal *cnfine, doublereal *result, ftnlen occtyp_len, ftnlen front_len, ftnlen fshape_len, ftnlen fframe_len, ftnlen back_len, ftnlen bshape_len, ftnlen bframe_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return gfoclt_(occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr, step, cnfine, result, occtyp_len, front_len, fshape_len, fframe_len, back_len, bshape_len, bframe_len, abcorr_len, obsrvr_len);
}

int Celestial_Object::srfcss_(integer *code, char *bodstr, char *srfstr, logical *isname, ftnlen bodstr_len, ftnlen srfstr_len)
{
  return srfcss_(code, bodstr, srfstr, isname, bodstr_len, srfstr_len);
}

int Celestial_Object::srfs2c_(char *srfstr, char *bodstr, integer *code, logical *found, ftnlen srfstr_len, ftnlen bodstr_len)
{
  return srfs2c_(srfstr, bodstr, code, found, srfstr_len, bodstr_len);
}

int Celestial_Object::srfc2s_(integer *code, integer *bodyid, char *srfstr, logical *isname, ftnlen srfstr_len)
{
  return srfc2s_(code, bodyid, srfstr, isname, srfstr_len);
}

int Celestial_Object::srfscc_(char *srfstr, integer *bodyid, integer *code, logical *found, ftnlen srfstr_len)
{
  return srfscc_(srfstr, bodyid, code, found, srfstr_len);
}

int Celestial_Object::bodc2n_(integer *code, char *name__, logical *found, ftnlen name_len)
{
  return bodc2n_(code, name__, found, name_len);
}

int Celestial_Object::srfrec_(integer *body, doublereal *long__, doublereal *lat, doublereal *rectan)
{
  return srfrec_(body, long__, lat, rectan);
}
