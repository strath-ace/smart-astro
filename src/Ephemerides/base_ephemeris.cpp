/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2018 University of Strathclyde and Authors ------
-------------------- e-mail: c.greco@strath.ac.uk --------------------
----------------------- Author: Cristian Greco -----------------------
*/

#include "base_ephemeris.h"


using namespace smartastro::ephemerides;
using namespace std;


/**
 * Default constructor
 *
 * @param referenceFrame: name of reference frame to which ephemerides are referred
 */
base_ephemeris::base_ephemeris(const ephemerisParams* pParams ) :
        m_pParams(pParams)
{
    // Nothing to initialise
}


/**
 * Default destructor
 *
 */
base_ephemeris::~base_ephemeris()
{
    // Nothing to delete
}


//
///**
// * getCartesianState: Function that returns object position-velocity for time t
// *
// * @param jd: time at which the Cartesian State is desired in jd
// * @param jd: Reference frame name of output state
// * @param jd: Reference frame center of output state
// * @return Cartesian state at time jd
// *
// */
//vector<double> base_ephemeris::getCartesianState( const double& jd,
//                                                  const string& outputReferenceFrame,
//                                                  const vector<double>& outputReferenceFrameCenter )
//{
//    // Cstate
//    SpiceDouble       x_or[6], x_tr[6];
//
//    // Cartesian state
//    vector<double> car = getCartesianState(jd);
//    for (unsigned int i = 0 ; i < 6 ; i++)
//        x_or[i] = car[i];
//
//    // Change origin
//    if (!outputReferenceFrame.empty())
//    {
//        if (outputReferenceFrame.size()!=3)
//            smartastro_throw("Center shall be 3 dimensional");
//
//        for (unsigned int i = 0 ; i < 6 ; i++)
//            x_or[i] = car[i];
//    }
//
//
//
//    // Compute rotation matrix
//    //
//    // Rotation matrix
//    SpiceDouble       xform[6][6];
//    sxform_c(m_pParams->referenceFrame.c_str(),outputReferenceFrame.c_str(),jd,xform);
//
//    // Rotate vector
//    mxvg_c ( xform,x_or,6,6,x_tr);
//
//
//
//}