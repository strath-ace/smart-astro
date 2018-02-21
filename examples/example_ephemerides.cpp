/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example shows how to use the analytical (thus low-precision) ephemerides for some bodies of the solar system */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../include/smartastro.h"

using namespace std;


int main(){

  std::ofstream myfile;
  double jd0, r, lambda;

  /** scaling parameters */
  double R_ref=smartastro::constants::R_earth; // distance in meters
  double T_ref=smartastro::constants::T_Earth; // time in seconds

  /** initial gregorian date */
  std::vector<double> date0(6);

  date0[0] = 2005.0;
  date0[1] = 2.0;
  date0[2] = 16.0;
  date0[3] = 8.0;
  date0[4] = 59.0;
  date0[5] = 24.0;

  std::vector<double> date = date0;
  smartastro::astrocore::conversion_time::date2jd(date0, jd0);

  // double eps = constants::obliquity_ecliptic_2000;

  smartastro::astrodynamics::dearth_3dof<double> *dyn = new smartastro::astrodynamics::dearth_3dof<double>(R_ref, T_ref);

  dyn->Sun(jd0, r, lambda);
  r *= R_ref * 1.0e-3;

  std::cout << std::setprecision(16) << - r * cos(lambda) << " " << - r * sin(lambda) << " " << 0.0 << std::endl;

  delete dyn;

  double mjd2000;
  smartastro::astrocore::conversion_time::date2mjd2000(date, mjd2000);
  std::vector<double> kep(6), car(6);
  smartastro::ephemerides::analytical_planets::get_orbel(mjd2000, 3, kep);
  smartastro::astrocore::conversion_coordinates::kep2car(kep, smartastro::constants::mu_earth, car);

  std::cout << std::setprecision(16) << car[0] << " " << car[1] << " " << car[2] << std::endl;

  return 0;

}
