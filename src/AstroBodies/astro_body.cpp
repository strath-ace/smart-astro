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

Astro_Body::Astro_Body(std::string givenName, std::vector<std::string> spiceEphemeridesParams)
{
  int found;
  
  name = givenName;
  smartastro::astrocore::spice_general_functions::bodn2c(name, id, found);

  if (found==0){
    smartastro_throw("Astro body name not found in SPICE.");
  } else {
    
    std::vector<std::string> kernel;
    
    for(int i = 1; i < spiceEphemeridesParams.size() - 3; ++i){
      kernel.push_back(spiceEphemeridesParams[i]);
    }
    
    smartastro::ephemerides::spiceEphemeris::spiceEphemerisParams spiceEphemerisParamsStruct = {
      spiceEphemeridesParams[0],
      kernel,
      spiceEphemeridesParams[spiceEphemeridesParams.size() - 2],
      name,
      spiceEphemeridesParams[spiceEphemeridesParams.size() - 1],
    };

    *spiceEphemeris = smartastro::ephemerides::spiceEphemeris(&spiceEphemerisParamsStruct);
  }
}

Astro_Body::Astro_Body(std::string givenName, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams)
{
  int found;

  name = givenName;
  smartastro::astrocore::spice_general_functions::bodn2c(name, id, found);

  if (found==0){
    smartastro_throw("Astro body name not found in SPICE.");
  } else {

    std::vector<double> xiVector;

    for(int i = 1; i < integratedEphemeridesParams.size() - 2; ++i){
      xiVector.push_back(integratedEphemeridesParams[i]);
    }

    smartastro::ephemerides::integratedEphemeris::integratedEphemerisParams integratedEphemerisParamsStruct = {
      referenceFrame,
      integratedEphemeridesFunction,
      integratedEphemeridesParams[0],
      xiVector,
      integratedEphemeridesParams[integratedEphemeridesParams.size() - 1],
    };

    *integratedEphemeris = smartastro::ephemerides::integratedEphemeris(&integratedEphemerisParamsStruct);
  }
}

Astro_Body::Astro_Body(int givenId, std::vector<std::string> spiceEphemeridesParams)
{
  int found;
  int lenout = 100;
  
  id = givenId;
  smartastro::astrocore::spice_general_functions::bodc2n(id, lenout, name, found);

  if (found==0){
    smartastro_throw("Astro body code not found in SPICE.");
  } else {

    std::vector<std::string> kernal;

    for(int i = 1; i < spiceEphemeridesParams.size() - 3; ++i){
      kernal.push_back(spiceEphemeridesParams[i]);
    }

    smartastro::ephemerides::spiceEphemeris::spiceEphemerisParams spiceEphemerisParamsStruct = {
      spiceEphemeridesParams[0],
      kernal,
      spiceEphemeridesParams[spiceEphemeridesParams.size() - 2],
      name,
      spiceEphemeridesParams[spiceEphemeridesParams.size() - 1],
    };

    *spiceEphemeris = smartastro::ephemerides::spiceEphemeris(&spiceEphemerisParamsStruct);
  }
}

Astro_Body::Astro_Body(int givenId, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams)
{
  int found;
  int lenout = 100;

  id = givenId;
  smartastro::astrocore::spice_general_functions::bodc2n(id, lenout, name, found);

  if (found==0){
    smartastro_throw("Astro body code not found in SPICE.");
  } else {

    std::vector<double> xiVector;

    for(int i = 1; i < integratedEphemeridesParams.size() - 2; ++i){
      xiVector.push_back(integratedEphemeridesParams[i]);
    }

    smartastro::ephemerides::integratedEphemeris::integratedEphemerisParams integratedEphemerisParamsStruct = {
      referenceFrame,
      integratedEphemeridesFunction,
      integratedEphemeridesParams[0],
      xiVector,
      integratedEphemeridesParams[integratedEphemeridesParams.size() - 1],
    };

    *integratedEphemeris = smartastro::ephemerides::integratedEphemeris(&integratedEphemerisParamsStruct);
  }
}

smartastro::ephemerides::base_ephemeris * Astro_Body::getEphemerides(const int ephemeridesType, const smartastro::ephemerides::base_ephemeris::ephemerisParams &pParams)
{
  switch(ephemeridesType) {
    case 1:
      if(integratedEphemeris == NULL){
	*integratedEphemeris = smartastro::ephemerides::integratedEphemeris(&pParams);
      }
      return integratedEphemeris;
      break;
    case 2:
      if(spiceEphemeris == NULL){
       *spiceEphemeris = smartastro::ephemerides::spiceEphemeris(&pParams);
      }
      return spiceEphemeris;
      break;
  }
}

void Astro_Body::getState(const StateType& stateType, const double& time, std::vector<double>& state)
{
  switch(stateType){
  case CARTESIAN:
    if(integratedEphemeris != NULL){
      state = integratedEphemeris->getCartesianState(time); 
    } else  {
      state = spiceEphemeris->getCartesianState(time);
    }
    break;
  case KEPLERIAN:
    break;
  case EQUINOCTIAL:
    break;
  }
}
