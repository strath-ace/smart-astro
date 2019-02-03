/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

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

int spice_general_functions::furnsh_(char *file, ftnlen file_len)
{
	return furnsh_(file, file_len);
}

int spice_general_functions::unload_(char *file, ftnlen file_len)
{
	return unload_(file, file_len);
}
