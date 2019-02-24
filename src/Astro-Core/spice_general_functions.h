/* Thi Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#ifndef SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H
#define SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H

#include <vector>
#include <string>
#include <cspice/SpiceUsr.h>
#include "../exception.h"

namespace smartastro
{
	namespace astrocore
	{

		class spice_general_functions
		{
			public:	

				/**
				* @brief Convert a one-dimensional vector to a SpiceDSKDescr
				*
				* @param[in] vector 1D vector to be converted from. Format: [surfce, center, dclass, dtype, 
				* frmcde, corsys, corpar[10], co1min, co1max, co2min, co2max, co3min, co3max, start, stop].
				* @param[out] descr SpiceDSKDescr to be converted to. 
				*
				*/

				static void vectorToSpiceDSKDescr(const std::vector<int> &vector, SpiceDSKDescr &descr);

				
				/**
				* @brief Convert a SpiceDSKDescr to a one-dimensional vector 
				*
				* @param[in] descr SpiceDSKDescr to be converted from.
				* @param[out] vector 1D vector to be converted to. Format: [surfce, center, dclass, dtype, 
				* frmcde, corsys, corpar[10], co1min, co1max, co2min, co2max, co3min, co3max, start, stop].
				*
				*/

				static void spiceDSKDescrToVector(const SpiceDSKDescr &descr, std::vector<int> &vector);

				/**
				* @brief Convert a one-dimensional vector to a SpiceDLADescr
				*
				* @param[in] vector 1D vector to be converted from. Format: [bwdptr, fwdptr, ibase, isize, 
				* dbase, dsize, cbase, csize].
				*
				*/

				static void vectorToSpiceDLADescr(const std::vector<int> &vector, SpiceDLADescr &descr);

				/**
				* @brief Convert a SpiceDLADescr to a one-dimensional vector 
				*
				* @param[in] descr SpiceDLADescr to be converted from.
				* @param[out] vector 1D vector to be converted to. Format: [bwdptr, fwdptr, ibase, isize, 
				* dbase, dsize, cbase, csize].
				*
				*/
				
				static void spiceDLADescrToVector(const SpiceDLADescr &descr, std::vector<int> &vector);

				/**
				* @brief Convert a one-dimensional vector with a normal vector and constant to a SpicePlane 
				*
				* @param[in] vector 1D vector to be converted from. Format: [bwdptr, fwdptr, ibase, isize, 
				* dbase, dsize, cbase, csize].
				* @param[out] plane SpicePlane to be converted to.
				*
				*/

				static void vectorToSpicePlaneNVC(const std::vector<double> &vector, SpicePlane &plane);

				/**
				* @brief Convert a one-dimensional vector with a normal vector and point to a SpicePlane 
				*
				* @param[in] vector 1D vector to be converted from. Format: [normal[3], constant].
				* @param[out] plane SpicePlane to be converted to.
				*
				*/

				static void vectorToSpicePlaneNVP(const std::vector<double> &vector, SpicePlane &plane);

				/**
				* @brief Convert a one-dimensional vector with a point and spanning vectors to a SpicePlane 
				*
				* @param[in] vector 1D vector to be converted from. Format: [normal[3], point[3]].
				* @param[out] plane SpicePlane to be converted to.
				*
				*/

				static void vectorToSpicePlanePSV(const std::vector<double> &vector, SpicePlane &plane);

				/**
				* @brief Convert a one-dimensional vector to a SpiceEllipse 
				*
				* @param[in] vector 1D vector to be converted from. Format: [point[3], span1[3], span2[3]].
				* @param[out] ellispe SpiceEllipse to be converted to.
				*
				*/

				static void vectorToSpiceEllipse(const std::vector<double> &vector, SpiceEllipse &ellipse);

				/**
				* @brief Convert a SpiceEllipse to a one-dimensional vector  
				*
				* @param[in] ellispe SPiceEllipse to be converted from
				* @param[out] vector 1D vector to be converted to. Format: [center[3], v1[3], v2[3]].
				*
				*/

				static void spiceEllipseToVector(const SpiceEllipse &ellipse, std::vector<double> &vector);

				/**
				* @brief Convert a one-dimensional vector to a two-dimensional array with a breadth of 2  
				*
				* @param[in] vector 1D vector to be converted from. Format: [center[3], v1[3], v2[3]].
				* @param[out] array 2D array to be converted to.
				*
				*/

				static void doubleVectorTo2dArray2(const std::vector<double> &vector, double (*array)[2]);

				/**
				* @brief Convert a one-dimensional vector to a two-dimensional array with a breadth of 3  
				*
				* @param[in] vector 1D vector to be converted from.
				* @param[out] array 2D array to be converted to.
				*
				*/

				static void doubleVectorTo2dArray3(const std::vector<double> &vector, double (*array)[3]);

				/**
				* @brief Convert a two-dimensional array with a breadth of 3 to a one-dimensional vector 
				*
				* @param[in] array 2D array to be converted from.
				* @param[out] vector 1D vector to be converted to.
				*
				*/

				static void double2dArray3ToVector(const double (*array)[3], std::vector<double> &vector);

				/**
				* @brief Convert a two-dimensional array with a breadth of 2 to a one-dimensional vector 
				*
				* @param[in] array 2D array to be converted from.
				* @param[out] vector 1D vector to be converted to.
				*
				*/

				static void double2dArray2ToVector(const double (*array)[2], std::vector<double> &vector);

				/**
				* @brief Translate the SPICE integer code of a body into a common name
   				* for that body. 
				*
				* @param[in] code Integer ID code to be translated into a name.
				* @param[in] lenout Maximum length of output name.
				* @param[out] name A common name for the body identified by code.
				* @param[out] found True if translated, otherwise false.
				*
				*/

				static void bodc2n(const int &code, const int &lenout, std::string &name, int &found);

				/**
				* @brief Translate the name of a body or object to the corresponding SPICE
   				* integer ID code.
				*
				* @param[in] name Body name to be translated into a SPICE ID code.
				* @param[out] code SPICE integer ID code for the named body.
				* @param[out] found SPICETRUE if translated, otherwise SPICEFALSE.
				*
				*/

				static void bodn2c(const std::string &name, int &code, int found);

				/**
				* @brief Compute the rectangular coordinates of the sub-solar point on 
   				* a target body at a specified epoch, optionally corrected for 
   				* light time and stellar aberration. The surface of the target body may 
				* be represented by a triaxial ellipsoid or by topographic data provided 
				* by DSK files. This routine supersedes subsol_c. 
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] et Epoch in ephemeris seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] abcorr Aberration correction flag.
				* @param[out] spoint Sub-solar point on the target body. 
				* @param[out] trgepc Sub-solar point epoch.
				* @param[out] srfvec Vector from observer to sub-solar point.
				* @return void
				*
				*/
				
				static void subslr(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr, std::vector<double> &spoint, double &trgepc, std::vector<double> &srfvec);

				/**
				* @brief Compute the rectangular coordinates of the sub-observer point on 
   				* a target body at a specified epoch, optionally corrected for 
   				* light time and stellar aberration. The surface of the target body may be 
				* represented by a triaxial ellipsoid or by topographic data provided by DSK files. 
				* This routine supersedes subpt_c.
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] et Epoch in ephemeris seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] abcorr Aberration correction flag.
				* @param[out] spoint Sub-solar point on the target body. 
				* @param[out] trgepc Sub-solar point epoch.
				* @param[out] srfvec Vector from observer to sub-solar point.
				* @return void
				*
				*/

				static void subpnt(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> srfvec);

				/**
				* @brief Compute L_s, the planetocentric longitude of the sun, as seen 
   				* from a specified body.
				*
				* @param[in] body Name of central body. 
				* @param[in] et Epoch in seconds past J2000 TDB. 
				* @param[in] abcorr Aberration correction. 
				* @return double
				*
				*/

				static double lspcn(const std::string &body, const double &et, const std::string &abcorr);

				/**
				* @brief Load one or more SPICE kernels into a program.
				*
				* @param[in] file Name of SPICE kernel file (text or binary). 
				* @return void
				*
				*/

				static void furnsh(const std::string &file);

				/**
				* @brief Unload a SPICE kernel.
				*
				* @param[in] file The name of a kernel to unload.
				* @return void
				*
				*/

				static void unload(const std::string &file);
		};
	}
}

#endif /* SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H */
