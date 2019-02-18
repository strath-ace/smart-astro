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

namespace smartastro
{
	namespace astrocore
	{

		class spice_general_functions
		{
			public:

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

				static void bodc2n(int &code, int &lenout, std::string &name, int &found);

				/**
				* @brief Translate the name of a body or object to the corresponding SPICE
   				* integer ID code.
				*
				* @param[in] name Body name to be translated into a SPICE ID code.
				* @param[out] code SPICE integer ID code for the named body.
				* @param[out] found SPICETRUE if translated, otherwise SPICEFALSE.
				*
				*/

				static void bodn2c(std::string &name, int &code, int found);

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
				
				static void subslr(std::string &method, std::string &target, double &et, std::string &fixref, std::string &abcorr, std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec);

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

				static void subpnt(std::string &method, std::string &target, double &et, std::string &fixref, std::string &abcorr, std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> srfvec);

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

				static double lspcn(std::string &body, double &et, std::string &abcorr);

				/**
				* @brief Load one or more SPICE kernels into a program.
				*
				* @param[in] file Name of SPICE kernel file (text or binary). 
				* @return void
				*
				*/

				static void furnsh(std::string &file);

				/**
				* @brief Unload a SPICE kernel.
				*
				* @param[in] file The name of a kernel to unload.
				* @return void
				*
				*/

				static void unload(std::string &file);
		};
	}
}

#endif /* SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H */
