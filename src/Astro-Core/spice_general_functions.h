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
				
				static void subslr(const std::string &method, const std::string &target, double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec);

				static void subpnt(const std::string &method, const std::string &target, double &et, const std::string &fixref, const std::string &abcorr, const std::string &obsrvr,  std::vector<double> &spoint, double &trgepc,  std::vector<double> srfvec);

				static double lspcn(const std::string &body, double &et, const std::string &abcorr);

				static void furnsh(const std::string &file);

				static void unload(const std::string &file);
		};
	}
}

#endif /* SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H */
