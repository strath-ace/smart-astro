#include "Astro-Core/spice_general_functions.h"

using namespace smartastro;
using namespace smartastro::astrocore;

int spice_general_functions::subslr_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
	return subslr_(method, target, et, fixref, abcorr, obsrvr, spoint, trgepc, srfvec, method_len, target_len, fixref_len, abcorr_len, obsrvr_len);
}

int spice_general_functions::subpnt_(char *method, char *target, doublereal *et, char *fixref, char *abcorr, char *obsrvr, doublereal *spoint, doublereal *trgepc, doublereal *srfvec, ftnlen method_len, ftnlen target_len, ftnlen fixref_len, ftnlen abcorr_len, ftnlen obsrvr_len)
{
	return subpnt_(method, target, et, fixref, abcorr, obsrvr, spoint, trgepc, srfvec, method_len, target_len, fixref_len, abcorr_len, obsrvr_len);
}

doublereal spice_general_functions::lspcn_(char *body, doublereal *et, char *abcorr, ftnlen body_len, ftnlen abcorr_len)
{
	return lspcn_(body, et, abcorr, body_len, abcorr_len);
}

int spice_general_functions::radrec_(doublereal *range, doublereal *ra, doublereal *dec, doublereal *rectan)
{
		return radrec_(range, ra, dec, rectan);
}

int spice_general_functions::recrad_(doublereal *rectan, doublereal *range, doublereal *ra, doublereal *dec)
{
	return recrad_(rectan, range, ra, dec);
}

int spice_general_functions::cylrec_(doublereal *r__, doublereal *long__, doublereal *z__, doublereal *rectan)
{
	return cylrec_(r__, long__, z__, rectan);
}

int spice_general_functions::reccyl_(doublereal *rectan, doublereal *r__, doublereal *long__, doublereal *z__)
{
	return reccyl_(rectan, r__, long__, z__);
}

int spice_general_functions::sphrec_(doublereal *r__, doublereal *colat, doublereal *long__, doublereal *rectan)
{
	return sphrec_(r__, colat, long__, rectan);
}

int spice_general_functions::recsph_(doublereal *rectan, doublereal *r__, doublereal *colat, doublereal *long__)
{
	return recsph_(rectan, r__, colat, long__);
}

int spice_general_functions::sphcyl_(doublereal *radius, doublereal *colat, doublereal *slong, doublereal *r__, doublereal *long__, doublereal *z__)
{
		return sphcyl_(radius, colat, slong, r__, long__, z__);
}

int spice_general_functions::cylsph_(doublereal *r__, doublereal *longc, doublereal *z__, doublereal *radius, doublereal *colat, doublereal *long__)
{
	return cylsph_(r__, longc, z__, radius, colat, long__);
}

int spice_general_functions::sphlat_(doublereal *r__, doublereal *colat, doublereal *longs, doublereal *radius, doublereal *long__, doublereal *lat)
{
	return sphlat_(r__, colat, longs, radius, long__, lat);
}

int spice_general_functions::latsph_(doublereal *radius, doublereal *long__, doublereal *lat, doublereal *rho, doublereal *colat, doublereal *longs)
{
	return latsph_(radius, long__, lat, rho, colat, longs);
}

int spice_general_functions::latcyl_(doublereal *radius, doublereal *long__, doublereal *lat, doublereal *r__, doublereal *longc, doublereal *z__)
{
	return latcyl_(radius, long__, lat, r__, longc, z__);
}

int spice_general_functions::cyllat_(doublereal *r__, doublereal *longc, doublereal *z__, doublereal *radius, doublereal *long__, doublereal *lat)
{
	return cyllat_(r__, longc, z__, radius, long__, lat);
}

int spice_general_functions::furnsh_(char *file, ftnlen file_len)
{
	return furnsh_(file, file_len);
}

int spice_general_functions::unload_(char *file, ftnlen file_len)
{
	return unload_(file, file_len);
}

int spice_general_functions::reclat_(doublereal *rectan, doublereal *radius, doublereal *long__, doublereal *lat)
{
	return reclat_(rectan, radius, long__, lat);
}

int spice_general_functions::latrec_(doublereal *radius, doublereal *long__, doublereal *lat, doublereal *rectan)
{
	return latrec_(radius, long__, lat, rectan);
}

int spice_general_functions::recgeo_(doublereal *rectan, doublereal *re, doublereal *f, doublereal *long__, doublereal *lat, doublereal *alt)
{
  return recgeo_(rectan, re, f, long__, lat, alt);
}

int spice_general_functions::georec_(doublereal *long__, doublereal *lat, doublereal *alt, doublereal *re, doublereal *f, doublereal *rectan)
{
  return georec_(long__, lat, alt, re, f, rectan);
}

int spice_general_functions::recpgr_(char *body, doublereal *rectan, doublereal *re, doublereal *f, doublereal *lon, doublereal *lat, doublereal *alt, ftnlen body_len)
{
  return recpgr_(body, rectan, re, f,lon, lat, alt, body_len);
}

int spice_general_functions::pgrrec_(char *body, doublereal *lon, doublereal *lat, doublereal *alt, doublereal *re, doublereal *f, doublereal *rectan, ftnlen body_len)
{
  return pgrrec_(body, lon, lat, alt, re, f, rectan, body_len);
}

int spice_general_functions::xfmsta_(doublereal *istate, char *icosys, char *ocosys, char *body, doublereal *ostate, ftnlen icosys_len, ftnlen ocosys_len, ftnlen body_len)
{
  return xfmsta_(istate, icosys, ocosys, body, ostate, icosys_len, ocosys_len, body_len);
}

int spice_general_functions::drdlat_(doublereal *r__, doublereal *long__, doublereal *lat, doublereal *jacobi)
{
  return drdlat_(r__, long__, lat, jacobi);
}

int spice_general_functions::dlatdr_(doublereal *x, doublereal *y, doublereal *z__, doublereal *jacobi)
{
  return dlatdr_(x, y, z__, jacobi);
}

int spice_general_functions::drdpgr_(char *body, doublereal *lon, doublereal *lat, doublereal *alt, doublereal *re, doublereal *f, doublereal *jacobi, ftnlen body_len)
{
  return drdpgr_(body, lon, lat, alt, re, f, jacobi, body_len);
}

int spice_general_functions::dpgrdr_(char *body, doublereal *x, doublereal *y, doublereal *z__, doublereal *re, doublereal *f, doublereal *jacobi, ftnlen body_len)
{
  return dpgrdr_(body, x, y, z__, re, f, jacobi, body_len);
}

int spice_general_functions::drdgeo_(doublereal *long__, doublereal *lat, doublereal *alt, doublereal *re, doublereal *f, doublereal *jacobi)
{
  return drdgeo_(long__, lat, alt, re, f, jacobi);
}

int spice_general_functions::dgeodr_(doublereal *x, doublereal *y, doublereal *z__, doublereal *re, doublereal *f, doublereal *jacobi)
{
  return dgeodr_(x, y, z__, re, f, jacobi);
}

int spice_general_functions::drdcyl_(doublereal *r__, doublereal *long__, doublereal *z__, doublereal *jacobi)
{
  return drdcyl_(r__, long__, z__, jacobi);
}

int spice_general_functions::dcyldr_(doublereal *x, doublereal *y, doublereal *z__, doublereal *jacobi)
{
  return dcyldr_(x, y, z__, jacobi);
}

int spice_general_functions::drdsph_(doublereal *r__, doublereal *colat, doublereal *long__, doublereal *jacobi)
{
  return drdsph_(r__, colat, long__, jacobi);
}

int spice_general_functions::dsphdr_(doublereal *x, doublereal *y, doublereal *z__, doublereal *jacobi)
{
  return dsphdr_(x, y, z__, jacobi);
}
