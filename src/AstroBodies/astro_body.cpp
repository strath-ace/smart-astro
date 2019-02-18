/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#include "AstroBodies/astro_body.h"

using namespace smartastro;
using namespace smartastro::astrobodies;

Astro_Body::Astro_Body(std::string givenName)
{
  int found;
  
  name = givenName;
  smartastro::astrocore::spice_general_functions::bodn2c(name, id, found);
}

Astro_Body::Astro_Body(int givenId)
{
  int found;
  int lenout = 100;
  
  id = givenId;
  smartastro::astrocore::spice_general_functions::bodc2n(id, lenout, name, found);
}
				 
Astro_Body::~Astro_Body()
{
  delete &name;
  delete &id;
}
