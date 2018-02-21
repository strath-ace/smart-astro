/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/
 
#include <catch.hpp>

#include "Astro-Core/astrocore.hpp"


namespace smart-sim
{
namespace tests
{

TEST_CASE( "The kep2car has been sussesfully tested", "[Conversion]" )
{

     double epsilon = 10e-6;
	
    // Test 1: Mercury Ephemerides 

    // Keplerian Reference
    Vector<double> in_kep(5.7909,  // semimajor axis in km
                          0.2056,  // eccentricity
                          0.1222,  // inclination in rad
	                  0.8461,  // right ascension of the ascending node in rad
                          0.5091,  // argument of perigee in rad
                         -2.6073); // true anomaly in rad

    // Planetary constant of the central body.	
    double in_mu = 1.3272448769E11;

    //
    Vector<double> com_car;
			    
    // Run and check that kep2car works without any exception
    REQUIRE( kep2car(in_kep, in_mu, com_car) );
	 	
    // Check the Position Vector
    REQUIRE( com_car[0] ==  34153299.6856891);   
    REQUIRE( com_car[1] == -54793665.566915);    
    REQUIRE( com_car[2] ==  -7604827.16037849);  
	 	
    // Check the Velocity Vector 
    REQUIRE( com_car[3] ==  31.6072051729721);
    REQUIRE( com_car[4] ==  28.1372492902702); 
    REQUIRE( com_car[5] ==  -0.616300097438105);
	 			 	
    // Test 2: Venus Ephemerides 

    // Orbit Type = KEPLERIAN
    in_kep[0] =  1.08208; // semimajor axis in km
    in_kep[1] =  0.00677; // eccentricity
    in_kep[2] =  0.05924; // inclination in rad
    in_kep[3] =  1.33674; // right ascension of the ascending node in rad
    in_kep[4] =  0.95714; // argument of perigee in rad
    in_kep[5] = -0.69183; // true anomaly in rad

    // Run and check that kep2car works without any exception
    REQUIRE( kep2car(in_kep, in_mu, com_car) );

    // Check the Position Vector
    REQUIRE( com_car[0] ==  -4241125.92693356); 
    REQUIRE( com_car[1] == 107540694.29033);
    REQUIRE( com_car[2] ==   1724073.54453162);
	 	
    // Check the Velocity Vector
    REQUIRE( com_car[3] == -35.1154203806911);
    REQUIRE( com_car[4] ==  -1.56698938124944);
    REQUIRE( com_car[5] ==   2.00453335564705);

}

} // namespace tests
} // namespace smart-sim
