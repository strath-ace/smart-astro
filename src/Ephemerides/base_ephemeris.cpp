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

