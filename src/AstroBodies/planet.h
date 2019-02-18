/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#ifndef PLANET_H
#define PLANET_H

#include "celestial_object.h"

namespace smartastro
{
	namespace astrobodies
	{

		class Planet : public Celestial_Object
		{
			public:
				Planet(std::string givenName, int givenId, std::vector<double> &Givenpositn, double givenMu);

				~Planet();
		};
	}
}

#endif // SMARTASTRO_PLANET_H
