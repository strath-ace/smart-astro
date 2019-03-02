/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#include "AstroBodies/planet.h"

using namespace smartastro;
using namespace smartastro::astrobodies;

Planet::Planet(std::string givenName, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &givenPosition, double givenMu) : Celestial_Object(givenName, spiceEphemeridesParams, givenPosition, givenMu)
{
}

Planet::Planet(std::string givenName, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>, std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &givenPosition, double givenMu) : Celestial_Object(givenName, referenceFrame, integratedEphemeridesFunction, integratedEphemeridesParams, givenPosition, givenMu)
{
}

Planet::Planet(int givenId, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &givenPosition, double givenMu) : Celestial_Object(givenId, spiceEphemeridesParams, givenPosition, givenMu)
{
}

Planet::Planet(int givenId, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &givenPosition, double givenMu) : Celestial_Object(givenId, referenceFrame, integratedEphemeridesFunction, integratedEphemeridesParams, givenPosition, givenMu)
{
}
