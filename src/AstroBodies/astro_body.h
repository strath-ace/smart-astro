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
#include "../exception.h"
#include "../Ephemerides/smartastro_ephemerides.h"

namespace smartastro
{
	namespace astrobodies
	{
		class Astro_Body
		{
			protected:
				std::string name;
				int id;
				bool ephemerisType;
				smartastro::ephemerides::integratedEphemeris *integratedephemeris;
				smartastro::ephemerides::spiceEphemeris *spiceephemeris;
				smartastro::ephemerides::spiceEphemeris::spiceEphemerisParams spiceEphemerisParamsStruct;
				smartastro::ephemerides::integratedEphemeris::integratedEphemerisParams integratedEphemerisParamsStruct;

			public:

				enum StateType {CARTESIAN, KEPLERIAN, EQUINOCTIAL};

				/**
         			* Constructor for providing a name and a spiceEphemeris
         			*
         			* @param givenName Name of astro body
				* @param spiceEphemeridesParams Parameters to create spiceEphemeris. Format: referenceFrame, kernelsToLoad, target, abberrationCorrection
         			*/

				Astro_Body(std::string givenName, std::vector<std::string> spiceEphemeridesParams);

				/**
         			* Constructor for providing a name and an integratedEphemeris
         			*
         			* @param givenName Name of astro body
				* @param referenceFrame Reference frame for integratedEphemeris
				* @param integratedEphemeridesFunction Integration function for integratedEphemeris
				* @param integratedEphemeridesParams  Parameters to create integratedEphemeris. Format: ti, xi, hstep
         			*/

				Astro_Body(std::string givenName, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>, std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams);

				/**
         			* Constructor for providing an ID and a spiceEphemeris
         			*
         			* @param givenId ID of astro body
				* @param spiceEphemeridesParams Parameters to create spiceEphemeris. Format: referenceFrame, kernelsToLoad, target, abberrationCorrection
         			*/

				Astro_Body(int givenId, std::vector<std::string> spiceEphemeridesParams);

				/**
         			* Constructor for providing an ID and an integratedEphemeris
         			*
         			* @param givenId ID of astro body
				* @param referenceFrame Reference frame for integratedEphemeris
				* @param integratedEphemeridesFunction Integration function for integratedEphemeris
				* @param integratedEphemeridesParams  Parameters to create integratedEphemeris. Format: ti, xi, hstep
         			*/

				Astro_Body(int givenId, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams);

				virtual ~Astro_Body() = default;

				/**
				* @brief Function that returns object position-velocity for time time
				*
				* @param stateType Type of state to be returned, CARTESIAN, KEPLERIAN or EQUINOCTIAL.
				* @param time Time at which the state is desired
				* @return Cartesian state at time t
				*
				*/

				std::vector<double> getState(const StateType& stateType, const double& time);		
		};
	}
}

#endif // SMARTASTRO_ASTRO_BODY_H
