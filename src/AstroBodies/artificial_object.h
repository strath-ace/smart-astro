/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#ifndef SMARTASTO_MAN_MADE_OBJECT_H
#define SMARTASTO_MAN_MADE_OBJECT_H

#include "astro_body.h"

namespace smartastro
{
	namespace astrobodies
	{
		class Artificial_Object : public Astro_Body
		{
			public:
				Artificial_Object( std::string givenName, int givenId);

				~Artificial_Object();

				/**
				* @brief Determine time intervals when a specified ray intersects the 
   				* space bounded by the field-of-view (FOV) of a specified 
   				* instrument.
				*
				* @param[in] inst Name of the instrument.
				* @param[in] raydir Ray's direction vector.
				* @param[in] rframe Reference frame of ray's direction vector.
				* @param[in] abcorr Aberration correction flag.
				* @param[in] step Step size in seconds for finding FOV events.
				* @param[in/out] cnfine SPICE window to which the search is restricted.
				* @param[out] SPICE window containing results.
				* @return void
				*
				*/

				void gfrfov(const std::string &inst,  std::vector<double> &raydir, const std::string &rframe, const std::string &abcorr, double &step, SpiceCell &cnfine, SpiceCell &result);

				void gftfov(const std::string &inst, const std::string &target, const std::string &tshape, const std::string &tframe, const std::string &abcorr, double &step, SpiceCell &cnfine, SpiceCell &result);

				void fovray(const std::string &inst,  std::vector<double> &raydir, const std::string &rframe, const std::string &abcorr, double &et,  int &visibl);

				void fovtrg( const std::string &inst, const std::string &target, const std::string &tshape, const std::string &tframe, const std::string &abcorr, double &et, int &visibl);

				/*void getfov(int &instid, int &room, int &shapelen, int &framelen, const std::string &shape, const std::string &frame,  std::vector<double> &bsight,  std::vector<int> &n,  std::vector<std::vector<double>> &bounds);

				int gffove_(std::string &inst, std::string &tshape,  std::vector<double> &raydir,  std::string &target, std::string &tframe, std::string &abcorr,   std::vector<double> &tol, int udstep, int udrefn, std::vector<int> &rpt, int udrepi, int udrepu, int udrepf, std::vector<int> &bail, L_fp udbail, std::vector<double> &cnfine,  std::vector<double> &result, int inst_len, int tshape_len, int target_len, int tframe_len, int abcorr_len, int obsrvr_len);*/
		};
	}
}

#endif // SMARTASTRO_ARTIFICIAL_OBJECT_H
