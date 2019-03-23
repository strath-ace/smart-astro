#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include "smartastro.h"
#include "AstroData/SpiceKernels/spiceKernelNames.h"
#include "Ephemerides/spiceEphemeris.h"
#include "Observations/smartastro_observations.h"
#include "../src/AstroBodies/astro_body.h"

#include <cspice/SpiceUsr.h>

#include <stdio.h>
#include <string>
#include "catch.hpp"

using namespace std;
using namespace smartastro;
using namespace ephemerides;
using namespace spiceKernels;
using namespace observations;

int main() {

    	string SPK (smartastro::spiceKernels::planetsEph);
	
	std::vector<std::string> params{"J2000", SPK, "moon", "NONE"};

	smartastro::astrobodies::Astro_Body astroBody("Earth", params);

	double mjd2000 = 0.0;
	vector<double> stateV;
	smartastro::astrobodies::Astro_Body::StateType stateType = smartastro::astrobodies::Astro_Body::CARTESIAN;

	astroBody.getState(stateType, mjd2000, stateV);

	//cout << "state = ";
    	//	for (unsigned int i = 0 ; i < 6 ; i++)
        //		cout << stateV[i] << " ";
    	//cout << endl;

	return 0;
}
