#ifndef SMARTASTO_CELESTIAL_OBJECT_H
#define SMARTASTO_CELESTIAL_OBJECT_H

extern "C" {
#include "../../../CSpice/cspice/include/SpiceUsr.h"
#include "../../../CSpice/cspice/include/SpiceZfc.h"
}

#include "astro_body.h"

namespace smartastro
{
	namespace astrobodies 
	{
		class Celestial_Object : public Astro_Body
		{
			public:
				Celestial_Object();

				~Celestial_Object();

				int latsrf_(char *method, char *target, doublereal *et, char *fixref, integer *npts, doublereal *lonlat, doublereal *srfpts, ftnlen method_len, ftnlen target_len, ftnlen fixref_len);

				int srfnrm_(char *method, char *target, doublereal *et, char *fixref, integer *npts, doublereal *srfpts, doublereal *normls, ftnlen method_len, ftnlen target_len, ftnlen fixref_len);

				int nearpt_(doublereal *positn, doublereal *a, doublereal *b, doublereal *c__, doublereal *npoint, doublereal *alt);

				int dnearp_(doublereal *state, doublereal *a, doublereal *b, doublereal *c__, doublereal *dnear, doublereal *dalt, logical *found);

				int surfpt_(doublereal *positn, doublereal *u, doublereal *a, doublereal *b, doublereal *c__, doublereal *point, logical *found);

				int surfnm_(doublereal *a, doublereal *b, doublereal *c__, doublereal *point, doublereal *normal);

				int ednmpt_(doublereal *a, doublereal *b, doublereal *c__, doublereal *normal, doublereal *point);

				int edpnt_(doublereal *p, doublereal *a, doublereal *b, doublereal *c__, doublereal *ep);

				int edlimb_(doublereal *a, doublereal *b, doublereal *c__, doublereal *viewpt, doublereal *limb);

				int inelpl_(doublereal *ellips, doublereal *plane, integer *nxpts, doublereal *xpt1, doublereal *xpt2);

				int npedln_(doublereal *a, doublereal *b, doublereal *c__, doublereal *linept, doublereal *linedr, doublereal *pnear, doublereal *dist);

				int inrypl_(doublereal *vertex, doublereal *dir, doublereal *plane, integer *nxpts, doublereal *xpt);

				int pjelpl_(doublereal *elin, doublereal *plane, doublereal *elout);

				int saelgv_(doublereal *vec1, doublereal *vec2, doublereal *smajor, doublereal *sminor);

				int nplnpt_(doublereal *linpt, doublereal *lindir, doublereal *point, doublereal *pnear, doublereal *dist);

				int npsgpt_(doublereal *ep1, doublereal *ep2, doublereal *point, doublereal *pnear, doublereal *dist);

				int sincpt_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *obsrvr, char *dref, doublereal *dvec, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, logical *found, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len, ftnlen dref_len);

				int dskxv_(logical *pri, char *target, integer *nsurf, integer *srflst, doublereal *et, char *fixref, integer *nrays, doublereal *vtxarr, doublereal *dirarr, doublereal *xptarr, logical *fndarr, ftnlen target_len, ftnlen fixref_len);

				int dskxsi_(logical *pri, char *target, integer *nsurf, integer *srflst, doublereal *et, char *fixref, doublereal *vertex, doublereal *raydir, integer *maxd, integer *maxi, doublereal *xpt, integer *handle, integer *dladsc, doublereal *dskdsc, doublereal *dc, integer *ic, logical *found, ftnlen target_len, ftnlen fixref_len);

				int ilumin_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, doublereal *phase, doublereal *solar, doublereal *emissn, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int illumg_(char *method, char *target, char *ilusrc, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, doublereal *phase, doublereal *incdnc, doublereal *emissn, ftnlen method_len, ftnlen target_len, ftnlen ilusrc_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int illumf_(char *method, char *target, char *ilusrc, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, doublereal *phase, doublereal *incdnc, doublereal *emissn, logical *visibl, logical *lit, ftnlen method_len, ftnlen target_len, ftnlen ilusrc_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int limbpt_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *corloc, char *obsrvr, doublereal *refvec, doublereal *rolstp, integer *ncuts, doublereal *schstp, doublereal *soltol, integer *maxn, integer *npts, doublereal *points, doublereal *epochs, doublereal *tangts, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen corloc_len, ftnlen obsrvr_len);

				int termpt_(char *method, char *ilusrc, char *target, doublereal *et, char *fixref, char *abcorr, char *corloc, char *obsrvr, doublereal *refvec, doublereal *rolstp, integer *ncuts, doublereal *schstp, doublereal *soltol, integer *maxn, integer *npts, doublereal *points, doublereal *epochs, doublereal *trmvcs, ftnlen method_len, ftnlen ilusrc_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen corloc_len, ftnlen obsrvr_len);

				int conics_(doublereal *elts, doublereal *et, doublereal *state);

				int oscelt_(doublereal *state, doublereal *et, doublereal *mu, doublereal *elts);

				int occult_(char *targ1, char *shape1, char *frame1, char *targ2, char *shape2, char *frame2, char *abcorr, char *obsrvr, doublereal *et, integer *ocltid, ftnlen targ1_len, ftnlen shape1_len, ftnlen frame1_len, ftnlen targ2_len, ftnlen shape2_len, ftnlen frame2_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int gfoclt_(char *occtyp, char *front, char *fshape, char *fframe, char *back, char *bshape, char *bframe, char *abcorr, char *obsrvr, doublereal *step, doublereal *cnfine, doublereal *result, ftnlen occtyp_len, ftnlen front_len, ftnlen fshape_len, ftnlen fframe_len, ftnlen back_len, ftnlen bshape_len, ftnlen bframe_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int srfcss_(integer *code, char *bodstr, char *srfstr, logical *isname, ftnlen bodstr_len, ftnlen srfstr_len);

				int srfs2c_(char *srfstr, char *bodstr, integer *code, logical *found, ftnlen srfstr_len, ftnlen bodstr_len);

				int srfc2s_(integer *code, integer *bodyid, char *srfstr, logical *isname, ftnlen srfstr_len);

				int srfscc_(char *srfstr, integer *bodyid, integer *code, logical *found, ftnlen srfstr_len);

				int bodc2n_(integer *code, char *name__, logical *found, ftnlen name_len);

				int srfrec_(integer *body, doublereal *long__, doublereal *lat, doublereal *rectan);
		};
	}
}

#endif // SMARTASTRO_CELESTIAL_OBJECT_H
