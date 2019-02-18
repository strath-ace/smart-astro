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
				
				static void subslr(const std::string &method, const std::string &target, double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec);

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

				static void subpnt(const std::string &method, const std::string &target, double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> srfvec);

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

				static double lspcn(const std::string &body, double &et, const std::string &abcorr);

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
