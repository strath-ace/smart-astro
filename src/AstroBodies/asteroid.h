/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#ifndef ASTEROID_H
#define ASTEROID_H

#include "celestial_object.h"

namespace smartastro
{
	namespace astrobodies
	{
		class Asteroid : public Celestial_Object
		{
			public:
				Asteroid(std::string givenName, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);

				Asteroid(std::string givenName, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);

				Asteroid(int givenId, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);

				Asteroid(int givenId, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);
				
				virtual ~Asteroid() = default;
		};
	}
}

#endif // SMARTASTRO_ASTEROID_H
