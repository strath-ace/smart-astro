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
				Asteroid( std::string givenName, int givenId,  std::vector<double> &Givenpositn,  std::vector<double> &givenMu);
				
				~Asteroid();
		};
	}
}

#endif // SMARTASTRO_ASTEROID_H
