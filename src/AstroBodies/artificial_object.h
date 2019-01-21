#ifndef SMARTASTO_MAN_MADE_OBJECT_H
#define SMARTASTO_MAN_MADE_OBJECT_H

extern "C" {
#include "../../../CSpice/cspice/include/SpiceUsr.h"
#include "../../../CSpice/cspice/include/SpiceZfc.h"
}

#include "astro_body.h"

namespace smartastro
{
	namespace astrobodies
	{
		class Artificial_Object : public Astro_Body
		{
			public:
				Artificial_Object();

				~Artificial_Object();

				int gfrfov_(char *inst, doublereal *raydir, char *rframe, char *abcorr, char *obsrvr, doublereal *step, doublereal *cnfine, doublereal *result, ftnlen inst_len, ftnlen rframe_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int gftfov_(char *inst, char *target, char *tshape, char *tframe, char *abcorr, char *obsrvr, doublereal *step, doublereal *cnfine, doublereal *result, ftnlen inst_len, ftnlen target_len, ftnlen tshape_len, ftnlen tframe_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int fovray_(char *inst, doublereal *raydir, char *rframe, char *abcorr, char *obsrvr, doublereal *et, logical *visibl, ftnlen inst_len, ftnlen rframe_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int fovtrg_(char *inst, char *target, char *tshape, char *tframe, char *abcorr, char *obsrvr, doublereal *et, logical *visibl, ftnlen inst_len, ftnlen target_len, ftnlen tshape_len, ftnlen tframe_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				int getfov_(integer *instid, integer *room, char *shape, char *frame, doublereal *bsight, integer *n, doublereal *bounds, ftnlen shape_len, ftnlen frame_len);

				int gffove_(char *inst, char *tshape, doublereal *raydir, char *target, char *tframe, char *abcorr, char *obsrvr, doublereal *tol, U_fp udstep, U_fp udrefn, logical *rpt, S_fp udrepi, U_fp udrepu, S_fp udrepf, logical *bail, L_fp udbail, doublereal *cnfine, doublereal *result, ftnlen inst_len, ftnlen tshape_len, ftnlen target_len, ftnlen tframe_len, ftnlen abcorr_len, ftnlen obsrvr_len);
		};
	}
}

#endif // SMARTASTRO_ARTIFICIAL_OBJECT_H
