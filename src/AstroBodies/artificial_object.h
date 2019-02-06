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
				Artificial_Object(std::string givenName, int givenId);

				~Artificial_Object();

				int gfrfov_(std::string &inst, std::vector<double> &raydir,  std::string &rframe, std::string &abcorr, std::vector<double> &step, std::vector<double> &cnfine, std::vector<double> &result, int inst_len, int rframe_len, int abcorr_len, int obsrvr_len);

				int gftfov_(std::string &inst, std::string &target,  std::string &tshape, std::string &tframe,  std::string &abcorr, std::string &obsrvr,  std::vector<double> &step, std::vector<double> &cnfine, std::vector<double> &result, int inst_len, int target_len, int tshape_len, int tframe_len, int abcorr_len, int obsrvr_len);

				int fovray_( std::string &inst, std::vector<double> &raydir,  std::string &rframe,  std::string &abcorr,  std::string &obsrvr, std::vector<double> &et,  std::vector<int> &visibl, int inst_len, int rframe_len, int abcorr_len, int obsrvr_len);

				int fovtrg_(std::string &inst, std::string &target,  std::string &tshape, std::string &tframe, std::string &abcorr, std::string &obsrvr,  std::vector<double> &et, std::vector<int> &visibl, int inst_len, int target_len, int tshape_len, int tframe_len, int abcorr_len, int obsrvr_len);

				int getfov_(std::vector<int> &instid, std::vector<int> &room,  std::string &shape, std::string &frame,  std::vector<double> &bsight, std::vector<int> &n,  std::vector<double> &bounds, int shape_len, int frame_len);

				int gffove_(std::string &inst, std::string &tshape,  std::vector<double> &raydir,  std::string &target, std::string &tframe, std::string &abcorr,  std::string &obsrvr, std::vector<double> &tol, int udstep, int udrefn, std::vector<int> &rpt, int udrepi, int udrepu, int udrepf, std::vector<int> &bail, L_fp udbail, std::vector<double> &cnfine,  std::vector<double> &result, int inst_len, int tshape_len, int target_len, int tframe_len, int abcorr_len, int obsrvr_len);
		};
	}
}

#endif // SMARTASTRO_ARTIFICIAL_OBJECT_H
