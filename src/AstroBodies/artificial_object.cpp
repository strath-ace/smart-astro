#include "AstroBodies/artificial_object.h"

using namespace smartastro;
using namespace smartastro::astrobodies;

Artificial_Object::Artificial_Object()
{
}

Artificial_Object::~Artificial_Object()
{
}

int Artificial_Object::gfrfov_(char *inst, doublereal *raydir, char *rframe, char *abcorr, char *obsrvr, doublereal *step, doublereal *cnfine, doublereal *result, ftnlen inst_len, ftnlen rframe_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return gfrfov_(inst, raydir, rframe, abcorr, obsrvr, step, cnfine, result, inst_len, rframe_len, abcorr_len, obsrvr_len);
}

int Artificial_Object::gftfov_(char *inst, char *target, char *tshape, char *tframe, char *abcorr, char *obsrvr, doublereal *step, doublereal *cnfine, doublereal *result, ftnlen inst_len, ftnlen target_len, ftnlen tshape_len, ftnlen tframe_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return gftfov_(inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine, result, inst_len, target_len, tshape_len, tframe_len, abcorr_len, obsrvr_len);
}

int Artificial_Object::fovray_(char *inst, doublereal *raydir, char *rframe, char *abcorr, char *obsrvr, doublereal *et, logical *visibl, ftnlen inst_len, ftnlen rframe_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return fovray_(inst, raydir, rframe, abcorr, obsrvr, et, visibl, inst_len, rframe_len, abcorr_len, obsrvr_len);
}

int Artificial_Object::fovtrg_(char *inst, char *target, char *tshape, char *tframe, char *abcorr, char *obsrvr, doublereal *et, logical *visibl, ftnlen inst_len, ftnlen target_len, ftnlen tshape_len, ftnlen tframe_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return fovtrg_(inst, target, tshape, tframe, abcorr, obsrvr, et, visibl, inst_len, target_len, tshape_len, tframe_len, abcorr_len, obsrvr_len);
}

int Artificial_Object::getfov_(integer *instid, integer *room, char *shape, char *frame, doublereal *bsight, integer *n, doublereal *bounds, ftnlen shape_len, ftnlen frame_len)
{
  return getfov_(instid, room, shape, frame, bsight, n, bounds, shape_len, frame_len);
}

int Artificial_Object::gffove_(char *inst, char *tshape, doublereal *raydir, char *target, char *tframe, char *abcorr, char *obsrvr, doublereal *tol, U_fp udstep, U_fp udrefn, logical *rpt, S_fp udrepi, U_fp udrepu, S_fp udrepf, logical *bail, L_fp udbail, doublereal *cnfine, doublereal *result, ftnlen inst_len, ftnlen tshape_len, ftnlen target_len, ftnlen tframe_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
  return gffove_(inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol, udstep, udrefn, rpt, udrepi, udrepu, udrepf, bail, udbail, cnfine, result, inst_len, tshape_len, target_len, tframe_len, abcorr_len, obsrvr_len);
}
