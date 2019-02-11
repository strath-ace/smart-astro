/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#ifndef SMARTASTO_ASTRO_BODY_H
#define SMARTASTO_ASTRO_BODY_H

#include <vector>
#include <cspice/SpiceUsr.h>
#include "../Astro-Core/spice_general_functions.h"

namespace smartastro
{
	namespace astrobodies
	{
		class Astro_Body
		{
			public:
				Astro_Body(std::string givenName, int givenId);

				~Astro_Body();

			protected:
				std::string name;
				int id;
		};
	}
}

#endif // SMARTASTRO_ASTRO_BODY_H
