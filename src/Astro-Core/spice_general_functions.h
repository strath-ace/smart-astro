#ifndef SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H
#define SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H

extern "C" {
#include "../../../CSpice/cspice/include/SpiceUsr.h"
#include "../../../CSpice/cspice/include/SpiceZfc.h"
}

namespace smartastro
{
	namespace astrocore
	{

		class spice_general_functions
		{
			public:

				static int subslr_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				static int subpnt_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len);

				static doublereal lspcn_(char *body, doublereal *et, char *abcorr, ftnlen body_len, ftnlen abcorr_len);

				static int furnsh_(char *file, ftnlen file_len);

				static int unload_(char *file, ftnlen file_len);

		};
	}
}

#endif /* SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H */
