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

using namespace std;

int main()
{

    double   et = 0.0;
    SpiceChar     utc[32];

    /*
       load LSK file
    */
    furnsh_c  ( "naif0012.tls" );

    /*
       convert UTC to ET
    */
    str2et_c  ( "2005 DEC 31 12:00", &et );

    /*
       add 1 day to ET and convert it back to UTC
    */
    timout_c ( et+spd_c(), "YYYY-DOYTHR:MN:SC.### ::RND", 32, utc );

    cout << endl;
    cout << "et = " << et << endl;
    cout << "utc = " ;
    for (unsigned int i = 0 ; i < 32; i++)
        cout << utc[i] ;
    cout << endl;

    cout << endl << "Welcome to template example" << endl;

    return 0;
}
