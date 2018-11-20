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
#include "Observations/smartastro_observations.h"


#include <cspice/SpiceUsr.h>

#include <stdio.h>
#include <string>

using namespace std;
using namespace smartastro;
using namespace ephemerides;
using namespace spiceKernels;
using namespace observations;

int moonFromEarth()
{

    string  SPK (smartastro::spiceKernels::planetsEph);

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

int sunAzimuthElevationFromStation()
{

    std::string  META ("spkcpo.tm");
    std::string  TIMFMT ("YYYY MON DD HR:MN:SC.###### UTC");
    std::string  TIMFM2 ("YYYY MON DD HR:MN:SC.###### TDB ::TDB");
    int          TIMLEN =  41;

    /*
     *   Local variables
     */
    ConstSpiceChar             * abcorr;
    SpiceChar               emitim  [ TIMLEN ];
    SpiceChar               epcstr  [ TIMLEN ];
    ConstSpiceChar             * refloc;
    ConstSpiceChar             * obsctr;
    ConstSpiceChar             * obsref;
    ConstSpiceChar             * obstim;
    ConstSpiceChar             * outref;
    ConstSpiceChar             * target;

    SpiceDouble             az;
    SpiceDouble             el;
    SpiceDouble             et;
    SpiceDouble             lat;
    SpiceDouble             lon;
    SpiceDouble             lt0;
    SpiceDouble             obsepc;
    SpiceDouble             obspos [ 3 ];
    SpiceDouble             r;
    SpiceDouble             state0 [ 6 ];


    /*
    Load SPICE kernels.
                       */
    std::vector<std::string> kernelToLoad ;
    kernelToLoad.push_back(leap);
    kernelToLoad.push_back(planetsEph);
    kernelToLoad.push_back(planetsData);
    kernelToLoad.push_back(earthHighAccOrientation);
    kernelToLoad.push_back(groundStationEph);
    kernelToLoad.push_back(groundStationTopo);

    for (unsigned int i = 0 ; i < kernelToLoad.size(); i++)
        furnsh_c ( kernelToLoad[i].c_str() );

    /*
    Convert the observation time to seconds past J2000 TDB.
                                                               */
    obstim = "2003 OCT 13 06:00:00.000000 UTC";
    str2et_c ( obstim, &et );

    /*
    Set the target, observer center, and observer frame.
                                                          */
    target = "MOON";
    obsctr = "EARTH";
    obsref = "ITRF93";

    /*
    Set the position of DSS-14 relative to the earth's
    center at the J2000 epoch, expressed in the
    ITRF93 reference frame. Values come from the
    earth station SPK specified in the meta-kernel.

            The actual station velocity is non-zero due
    to tectonic plate motion; we ignore the motion
    in this example. See the routine SPKCVO for an
    example in which the plate motion is accounted for.
    */
    obsepc    =  0.0;

    obspos[0] =  -2353.6213656676991;
    obspos[1] =  -4641.3414911499403;
    obspos[2] =   3677.0523293197439;

    /*
        compute geometric state of observer relative to Earth center
    */
    SpiceDouble   state[6];
    SpiceDouble   lt;
    spkezr_c ( "DSS-14", et, "ITRF93", "NONE", "EARTH", state,
               &lt );

    cout << "DSS-14 state = " ;
    for (unsigned int i = 0 ; i < 6 ; i++)
        cout << state[i] << " " ;
    cout << endl;

    /*
    Find the apparent state of the sun relative
    to the station in the DSS-14_TOPO reference frame.
            Evaluate the output frame's orientation, that is the
    orientation of the DSS-14_TOPO frame relative to the
    J2000 frame, at the observation epoch. This
    correction is obtained by setting `refloc' to
    "OBSERVER".
            */

    outref = "DSS-14_TOPO";
    abcorr = "CN+S";

    refloc = "OBSERVER";

    /*
    Compute the observer-target state.
                                        */
    spkcpo_c ( target, et,     outref, refloc,
               abcorr, obspos, obsctr,
               obsref, state0, &lt0            );

    /*
    Compute planetocentric coordinates of the
    observer-target position in the local
    topocentric reference frame DSS-14_TOPO.
                                            */
    reclat_c ( state0, &r, &lon, &lat );

    /*
    Compute solar azimuth. The latitude we've
    already computed is the elevation. Express
    both angles in degrees.
                           */
    el =   lat * dpr_c();
    az = - lon * dpr_c();

    if ( az < 0.0 )
    {
        az +=  360.0;
    }

    /*
    Display the computed state, light time. and angles.
                                                        */
    timout_c ( et-lt0, TIMFMT.c_str(), TIMLEN, emitim );
    timout_c ( obsepc, TIMFM2.c_str(), TIMLEN, epcstr );

    printf ( "\n"
             " Frame evaluation locus:     %s\n"
             "\n"
             " Target:                     %s\n"
             " Observation time:           %s\n"
             " Observer center:            %s\n"
             " Observer-center state time: %s\n"
             " Observer frame:             %s\n"
             " Emission time:              %s\n"
             " Output reference frame:     %s\n"
             " Aberration correction:      %s\n"
             "\n"
             " Observer-target position (km):\n"
             "   %20.8f %20.8f %20.8f\n"
             " Observer-target velocity (km/s):\n"
             "   %20.8f %20.8f %20.8f\n"
             " Light time (s):        %20.8f\n",

             refloc,    target,    obstim,    obsctr,
             epcstr,    obsref,    emitim,    outref,
             abcorr,    state0[0], state0[1], state0[2],
             state0[3], state0[4], state0[5], lt0   );

    printf ( "\n"
             " Solar azimuth (deg):   %20.8f\n"
             " Solar elevation (deg): %20.8f\n",
             az, el                             );


    /*
     * UnLoad SPICE kernels.
     */
    for (unsigned int i = 0 ; i < kernelToLoad.size(); i++)
        unload_c ( kernelToLoad[i].c_str() );

    return 0;
}



int sunAzimuthElevationFromObservation()
{

    std::string  TIMFMT ("YYYY MON DD HR:MN:SC.###### UTC");
    std::string  TIMFMT3 ("JD .############");
    std::string  TIMFM2 ("YYYY MON DD HR:MN:SC.###### TDB ::TDB");
    int          TIMLEN =  41;

    /*
     *   Local variables
     */
    ConstSpiceChar             * abcorr;
    SpiceChar               emitim  [ TIMLEN ];
    SpiceChar               epcstr  [ TIMLEN ];
    ConstSpiceChar             * refloc;
    ConstSpiceChar             * obsctr;
    ConstSpiceChar             * obsref;
    ConstSpiceChar             * obstim;
    ConstSpiceChar             * outref;
    ConstSpiceChar             * target;

    SpiceDouble             az;
    SpiceDouble             el;
    SpiceDouble             et;
    SpiceDouble             lat;
    SpiceDouble             lon;
    SpiceDouble             lt0;
    SpiceDouble             obsepc;
    SpiceDouble             obspos [ 3 ];
    SpiceDouble             r;
    SpiceDouble             state0 [ 6 ];


    /*
    Load SPICE kernels.
                       */
    std::vector<std::string> kernelToLoad ;
    kernelToLoad.push_back(leap);
    kernelToLoad.push_back(planetsEph);
    kernelToLoad.push_back(planetsData);
    kernelToLoad.push_back(earthHighAccOrientation);
    kernelToLoad.push_back(groundStationEph);
    kernelToLoad.push_back(groundStationTopo);

    for (unsigned int i = 0 ; i < kernelToLoad.size(); i++)
        furnsh_c ( kernelToLoad[i].c_str() );

    /*
    Convert the observation time to seconds past J2000 TDB.
                                                               */
    obstim = "2451545.000000 JD";
    str2et_c ( obstim, &et );

    double jd = unitim_c(et,"ET","JED");



    /*
    Set the target, observer center, and observer frame.
                                                          */
    target = "MOON";
    obsctr = "EARTH";
    obsref = "ITRF93";

    /*
    Set the position of DSS-14 relative to the earth's
    center at the J2000 epoch, expressed in the
    ITRF93 reference frame. Values come from the
    earth station SPK specified in the meta-kernel.

            The actual station velocity is non-zero due
    to tectonic plate motion; we ignore the motion
    in this example. See the routine SPKCVO for an
    example in which the plate motion is accounted for.
    */

    obsepc    =  et;

    obspos[0] =  -2353.6213656676991;
    obspos[1] =  -4641.3414911499403;
    obspos[2] =   3677.0523293197439;

    /*
    compute geometric state of observer relative to Earth center
        */


    /*
    Find the apparent state of the sun relative
    to the station in the DSS-14_TOPO reference frame.
            Evaluate the output frame's orientation, that is the
    orientation of the DSS-14_TOPO frame relative to the
    J2000 frame, at the observation epoch. This
    correction is obtained by setting `refloc' to
    "OBSERVER".
            */

    outref = "DSS-14_TOPO";
    abcorr = "CN+S";

    refloc = "OBSERVER";

    /*
    Compute the observer-target state.
                                        */
    spkcpo_c ( target, et,     outref, refloc,
               abcorr, obspos, obsctr,
               obsref, state0, &lt0            );

    /*
    Compute planetocentric coordinates of the
    observer-target position in the local
    topocentric reference frame DSS-14_TOPO.
                                            */
    reclat_c ( state0, &r, &lon, &lat );

    /*
    Compute solar azimuth. The latitude we've
    already computed is the elevation. Express
    both angles in degrees.
                           */
    el =   lat * dpr_c();
    az = - lon * dpr_c();

    if ( az < 0.0 )
    {
        az +=  360.0;
    }

    /*
    Display the computed state, light time. and angles.
                                                        */
    timout_c ( et-lt0, TIMFMT.c_str(), TIMLEN, emitim );
    timout_c ( obsepc, TIMFM2.c_str(), TIMLEN, epcstr );

    printf ( "\n"
             " Frame evaluation locus:     %s\n"
             "\n"
             " Target:                     %s\n"
             " Observation time:           %s\n"
             " Observer center:            %s\n"
             " Observer-center state time: %s\n"
             " Observer frame:             %s\n"
             " Emission time:              %s\n"
             " Output reference frame:     %s\n"
             " Aberration correction:      %s\n"
             "\n"
             " Observer-target position (km):\n"
             "   %20.8f %20.8f %20.8f\n"
             " Observer-target velocity (km/s):\n"
             "   %20.8f %20.8f %20.8f\n"
             " Light time (s):        %20.8f\n",

             refloc,    target,    obstim,    obsctr,
             epcstr,    obsref,    emitim,    outref,
             abcorr,    state0[0], state0[1], state0[2],
             state0[3], state0[4], state0[5], lt0   );

    printf ( "\n"
             " Solar azimuth (deg):   %20.8f\n"
             " Solar elevation (deg): %20.8f\n",
             az, el                             );







    /**
     * Use observations and ephemerides now
     *
     */



    /**
     * Ephemerides for Ground station
     */
    spiceEphemeris::spiceEphemerisParams groundStationParams;
    groundStationParams.referenceFrame        = "DSS-14_TOPO";
    groundStationParams.abberrationCorrection = "CN+S";
    groundStationParams.observer              = "EARTH";
    groundStationParams.target                = "DSS-14";

    spiceEphemeris groundStationEphemeris (&groundStationParams);

    // Link function to compute station position
    function<vector<double>(double)> getStationState = bind(&spiceEphemeris::getCartesianState,
                                                            groundStationEphemeris,
                                                            placeholders::_1);

    /**
     * Ephemerides for Sun
     */
    spiceEphemeris::spiceEphemerisParams sunParams;
    sunParams.referenceFrame        = "DSS-14_TOPO";
    sunParams.abberrationCorrection = "CN+S";
    sunParams.observer              = "EARTH";
    sunParams.target                = "MOON";

    spiceEphemeris sunEphemeris (&sunParams);

    // Link function to compute station position
    function<vector<double>(double)> getSunState = bind(&spiceEphemeris::getCartesianState,
                                                        sunEphemeris,
                                                        placeholders::_1);




    // Check
    vector<double> sunState = getSunState(jd);
    vector<double> staState = getStationState(jd);
    vector<double> relState (6);
    for (unsigned int i = 0 ;i < 6 ; i++)
        relState[i] = sunState[i]-staState[i];

    cout << "Sun State = ";
    for (unsigned int i = 0 ;i < 6 ; i++)
        cout << setprecision(15) << sunState[i] << " ";
    cout << endl;
    cout << "Sta State = ";
    for (unsigned int i = 0 ;i < 6 ; i++)
        cout << setprecision(15) << staState[i] << " ";
    cout << endl;
    cout << "Rel State = ";
    for (unsigned int i = 0 ;i < 6 ; i++)
        cout << setprecision(15) << relState[i] << " ";
    cout << endl << endl;



    /**
     * Create observation
     */
    rangeObservation::observationParams sParams;
    std::vector<double> sensorState = getStationState(jd);
    std::vector<double> targetState = getSunState(jd);
    

    //range
    rangeObservation range(&sParams);

    cout << "JD = " << jd << "; Range Ephemeris = " << setprecision(15) << range.getPerfectObservation(sensorState,targetState) << endl;
    cout << "JD = " << jd << "; Range Cspice    = " << setprecision(15) << sqrt(pow(state0[0],2.0)+pow(state0[1],2.0)+pow(state0[2],2.0)) << endl;
    cout << endl;


    // Elevation-azimuth
    azimuthElevationObservation azEl(&sParams);

    vector<double> azElMeas = azEl.getPerfectObservation(sensorState,targetState);

    cout << "JD = " << jd << "; Az = " << setprecision(15) << azElMeas[0]*180.0/M_PI << "; El = " << azElMeas[1]*180.0/M_PI << endl;
    cout << "JD = " << jd << "; Az = " << setprecision(15) << az  << "; El = " << el  << endl;




















    /*
     * UnLoad SPICE kernels.
     */
    for (unsigned int i = 0 ; i < kernelToLoad.size(); i++)
        unload_c ( kernelToLoad[i].c_str() );

    return 0;
}



int main()
{

    // Compute Moon position from Earth using ephemerides
//    moonFromEarth();


    // Compute Sun azimuth and elevation from spice routine
//    sunAzimuthElevationFromStation();

    sunAzimuthElevationFromObservation();

    return 0;
}
