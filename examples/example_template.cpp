/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
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

#include <cspice/SpiceUsr.h>

#include <stdio.h>
#include <string>

int main()
{


    double   et = 0.0;
    SpiceChar     utc[32];

    /*
       load kernels: LSK, PCK, CASSINI SCLK, FK and CK
    */
    furnsh_c ( "naif0012.tls" );

    /*
       convert UTC to ET
    */
    ConstSpiceChar * string = "Oct 1, 1996 09:12:32";

    /*
        add 1 day to ET and convert it back to UTC
    */
    timout_c ( et+1.0, "YYYY-DOYTHR:MN:SC.### ::RND", 32, utc );


    return 0;
}
