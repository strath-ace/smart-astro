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
  ephemerisType = false;
  int found;
  std::vector<std::string> kernel;
  
  name = givenName;
  smartastro::astrocore::spice_general_functions::bodn2c(name, id, found);

  if (found==0){
    smartastro_throw("Astro body name not found in SPICE.");
  } else {

    for(int i = 1; i < spiceEphemeridesParams.size() - 2; ++i){
	kernel.push_back(spiceEphemeridesParams[i]);
    }

    spiceEphemerisParamsStruct.referenceFrame = spiceEphemeridesParams[0];
    spiceEphemerisParamsStruct.kernelToLoad = kernel;
    spiceEphemerisParamsStruct.target = spiceEphemeridesParams[spiceEphemeridesParams.size() - 2];
    spiceEphemerisParamsStruct.observer = name;
    spiceEphemerisParamsStruct.abberrationCorrection = spiceEphemeridesParams[spiceEphemeridesParams.size() - 1];

    spiceephemeris = new smartastro::ephemerides::spiceEphemeris(&spiceEphemerisParamsStruct);
  }
}

Astro_Body::Astro_Body(std::string givenName, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams)
{
  int found;
  ephemerisType = true;
  name = givenName;
  smartastro::astrocore::spice_general_functions::bodn2c(name, id, found);

  if (found==0){
    smartastro_throw("Astro body name not found in SPICE.");
  } else {

    std::vector<double> xiVector;

    for(int i = 1; i < integratedEphemeridesParams.size() - 1; ++i){
      xiVector.push_back(integratedEphemeridesParams[i]);
    }

     integratedEphemerisParamsStruct.referenceFrame = referenceFrame;
     integratedEphemerisParamsStruct.integrate = integratedEphemeridesFunction;
     integratedEphemerisParamsStruct.ti = integratedEphemeridesParams[0];
     integratedEphemerisParamsStruct.xi = xiVector;
     integratedEphemerisParamsStruct.hstep = integratedEphemeridesParams[integratedEphemeridesParams.size() - 1];

    integratedephemeris = new smartastro::ephemerides::integratedEphemeris(&integratedEphemerisParamsStruct);
  }
}

Astro_Body::Astro_Body(int givenId, std::vector<std::string> spiceEphemeridesParams)
{
  ephemerisType = false;
  int found =1;
  int lenout = 100;
  std::vector<std::string> kernel;
  
  id = givenId;

  smartastro::astrocore::spice_general_functions::bodc2n(id, lenout, name, found);

  if (found==0){
    smartastro_throw("Astro body name not found in SPICE.");
  } else {

    for(int i = 1; i < spiceEphemeridesParams.size() - 2; ++i){
	kernel.push_back(spiceEphemeridesParams[i]);
    }

    spiceEphemerisParamsStruct.referenceFrame = spiceEphemeridesParams[0];
    spiceEphemerisParamsStruct.kernelToLoad = kernel;
    spiceEphemerisParamsStruct.target = spiceEphemeridesParams[spiceEphemeridesParams.size() - 2];
    spiceEphemerisParamsStruct.observer = name;
    spiceEphemerisParamsStruct.abberrationCorrection = spiceEphemeridesParams[spiceEphemeridesParams.size() - 1];

    spiceephemeris = new smartastro::ephemerides::spiceEphemeris(&spiceEphemerisParamsStruct);
  }
}

Astro_Body::Astro_Body(int givenId, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams)
{
  ephemerisType = true;
  int found;
  int lenout = 100;

  id = givenId;
  smartastro::astrocore::spice_general_functions::bodc2n(id, lenout, name, found);

  if (found==0){
    smartastro_throw("Astro body code not found in SPICE.");
  } else {

    std::vector<double> xiVector;

    for(int i = 1; i < integratedEphemeridesParams.size() - 1; ++i){
      xiVector.push_back(integratedEphemeridesParams[i]);
    }

   integratedEphemerisParamsStruct.referenceFrame = referenceFrame;
   integratedEphemerisParamsStruct.integrate = integratedEphemeridesFunction;
   integratedEphemerisParamsStruct.ti = integratedEphemeridesParams[0];
   integratedEphemerisParamsStruct.xi = xiVector;
   integratedEphemerisParamsStruct.hstep = integratedEphemeridesParams[integratedEphemeridesParams.size() - 1];

  integratedephemeris = new  smartastro::ephemerides::integratedEphemeris(&integratedEphemerisParamsStruct);
  }
}

std::vector<double> Astro_Body::getState(const StateType& stateType, const double& time)
{
  switch(stateType){
  case CARTESIAN:
    if(ephemerisType){
      return integratedephemeris->getCartesianState(time); 
    } else if(!ephemerisType) {
		return spiceephemeris->getCartesianState(time);
	} else {
		smartastro_throw("This body has no epheneris data.");   		
	}
    break;
  case KEPLERIAN:
    break;
  case EQUINOCTIAL:
    break;
  }
}
