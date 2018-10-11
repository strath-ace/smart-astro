/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/* */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/


#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include "smartastro.h"
#include "AstroData/SpiceKernels/spiceKernelNames.h"
#include "Ephemerides/spiceEphemeris.h"


#include <cspice/SpiceUsr.h>

#include <stdio.h>
#include <string>

using namespace std;

int main()
{

    string  SPK (smartastro::spiceKernels::planets);


    /*
    Local variables
    */
    SpiceDouble    et = 0.0;
    SpiceDouble    lt;
    SpiceDouble    state [6];


    /*
     * Load the spk file.
     */
    furnsh_c ( SPK.c_str() );
    furnsh_c ( smartastro::spiceKernels::leap.c_str() );

    // Convert time
    str2et_c("2000 Jan 01 12:00:00",&et);

    // find position
    spkezr_c ( "moon",    et,     "J2000",  "NONE",
               "earth",  state,  &lt             );

    cout << "state = ";
    for (unsigned int i = 0 ; i < 6 ; i++)
        cout << state[i] << " ";
    cout << endl;

    /**
     * Now try with ephemerides
     */
    smartastro::ephemerides::spiceEphemeris::spiceEphemerisParams epParams;
    epParams.observer = "earth";
    epParams.target   = "moon";
    epParams.abberrationCorrection = "NONE";
    epParams.referenceFrame = "J2000";
    epParams.kernelToLoad = vector<string>(1,SPK);

    smartastro::ephemerides::spiceEphemeris spiceEp (&epParams);

    double mjd2000 = 0.0;
    vector<double> epState = spiceEp.getCartesianState(mjd2000);

    cout << "state = ";
    for (unsigned int i = 0 ; i < 6 ; i++)
        cout << epState[i] << " ";
    cout << endl;


    return 0;
}