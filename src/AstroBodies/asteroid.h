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

				/**
         			* Constructor for providing a name and a spiceEphemeris
         			*
         			* @param givenName Name of astro body
				* @param spiceEphemeridesParams Parameters to create spiceEphemeris. Format: referenceFrame, kernelsToLoad, target, abberrationCorrection
				* @param Givenpositn Position of the observer in body-fixed frame.
				* @param givenMu Gravitational parameter (GM)
         			*/

				Asteroid(std::string givenName, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);

				/**
         			* Constructor for providing a name and an integratedEphemeris
         			*
         			* @param givenName Name of astro body
				* @param referenceFrame Reference frame for integratedEphemeris
				* @param integratedEphemeridesFunction Integration function for integratedEphemeris
				* @param integratedEphemeridesParams  Parameters to create integratedEphemeris. Format: ti, xi, hstep
				* @param Givenpositn Position of the observer in body-fixed frame.
				* @param givenMu Gravitational parameter (GM)
         			*/

				Asteroid(std::string givenName, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);

				/**
         			* Constructor for providing an ID and a spiceEphemeris
         			*
         			* @param givenId ID of astro body
				* @param spiceEphemeridesParams Parameters to create spiceEphemeris. Format: referenceFrame, kernelsToLoad, target, abberrationCorrection
				* @param Givenpositn Position of the observer in body-fixed frame.
				* @param givenMu Gravitational parameter (GM)
         			*/

				Asteroid(int givenId, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);

				/**
         			* Constructor for providing an ID and an integratedEphemeris
         			*
         			* @param givenId ID of astro body
				* @param referenceFrame Reference frame for integratedEphemeris
				* @param integratedEphemeridesFunction Integration function for integratedEphemeris
				* @param integratedEphemeridesParams  Parameters to create integratedEphemeris. Format: ti, xi, hstep
				* @param Givenpositn Position of the observer in body-fixed frame.
				* @param givenMu Gravitational parameter (GM)
         			*/

				Asteroid(int givenId, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);
				
				virtual ~Asteroid() = default;
		};
	}
}

#endif // SMARTASTRO_ASTEROID_H
