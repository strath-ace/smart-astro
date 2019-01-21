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

				static int radrec_(doublereal *range, doublereal *ra, doublereal *dec, doublereal *rectan);

				static int recrad_(doublereal *rectan, doublereal *range, doublereal *ra, doublereal *dec);

				static int cylrec_(doublereal *r__, doublereal *long__, doublereal *z__, doublereal *rectan);

				static int reccyl_(doublereal *rectan, doublereal *r__, doublereal *long__, doublereal *z__);

				static int sphrec_(doublereal *r__, doublereal *colat, doublereal *long__, doublereal *rectan);

				static int recsph_(doublereal *rectan, doublereal *r__, doublereal *colat, doublereal *long__);

				static int sphcyl_(doublereal *radius, doublereal *colat, doublereal *slong, doublereal *r__, doublereal *long__, doublereal *z__);

				static int cylsph_(doublereal *r__, doublereal *longc, doublereal *z__, doublereal *radius, doublereal *colat, doublereal *long__);

				static int sphlat_(doublereal *r__, doublereal *colat, doublereal *longs, doublereal *radius, doublereal *long__, doublereal *lat);

				static int latsph_(doublereal *radius, doublereal *long__, doublereal *lat, doublereal *rho, doublereal *colat, doublereal *longs);

				static int latcyl_(doublereal *radius, doublereal *long__, doublereal *lat, doublereal *r__, doublereal *longc, doublereal *z__);

				static int cyllat_(doublereal *r__, doublereal *longc, doublereal *z__, doublereal *radius, doublereal *long__, doublereal *lat);

				static int furnsh_(char *file, ftnlen file_len);

				static int unload_(char *file, ftnlen file_len);

				static int reclat_(doublereal *rectan, doublereal *radius, doublereal *long__, doublereal *lat);

				static int latrec_(doublereal *radius, doublereal *long__, doublereal *lat, doublereal *rectan);

				static int recgeo_(doublereal *rectan, doublereal *re, doublereal *f, doublereal *long__, doublereal *lat, doublereal *alt);

				static int georec_(doublereal *long__, doublereal *lat, doublereal *alt, doublereal *re, doublereal *f, doublereal *rectan);

				static int recpgr_(char *body, doublereal *rectan, doublereal *re, doublereal *f, doublereal *lon, doublereal *lat, doublereal *alt, ftnlen body_len);

				static int pgrrec_(char *body, doublereal *lon, doublereal *lat, doublereal *alt, doublereal *re, doublereal *f, doublereal *rectan, ftnlen body_len);

				static int xfmsta_(doublereal *istate, char *icosys, char *ocosys, char *body, doublereal *ostate, ftnlen icosys_len, ftnlen ocosys_len, ftnlen body_len);

				static int drdlat_(doublereal *r__, doublereal *long__, doublereal *lat, doublereal *jacobi);

				static int dlatdr_(doublereal *x, doublereal *y, doublereal *z__, doublereal *jacobi);

				static int drdpgr_(char *body, doublereal *lon, doublereal *lat, doublereal *alt, doublereal *re, doublereal *f, doublereal *jacobi, ftnlen body_len);

				static int dpgrdr_(char *body, doublereal *x, doublereal *y, doublereal *z__, doublereal *re, doublereal *f, doublereal *jacobi, ftnlen body_len);

				static int drdgeo_(doublereal *long__, doublereal *lat, doublereal *alt, doublereal *re, doublereal *f, doublereal *jacobi);

				static int dgeodr_(doublereal *x, doublereal *y, doublereal *z__, doublereal *re, doublereal *f, doublereal *jacobi);

				static int drdcyl_(doublereal *r__, doublereal *long__, doublereal *z__, doublereal *jacobi);

				static int dcyldr_(doublereal *x, doublereal *y, doublereal *z__, doublereal *jacobi);

				static int drdsph_(doublereal *r__, doublereal *colat, doublereal *long__, doublereal *jacobi);

				static int dsphdr_(doublereal *x, doublereal *y, doublereal *z__, doublereal *jacobi);

		};
	}
}

#endif /* SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H */
