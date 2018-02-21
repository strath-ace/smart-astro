/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example is for the semi-analytical case of a submitted paper to AIAA JGCD on the effects of the Lorentz force */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../include/smartastro.h"

using namespace std;

int main(){

  std::ofstream semianalytical, semianalytical2;

  double pi=smartastro::constants::pi;

  int n_max = 3;

  /** scaling parameters */
  double R_ref=smartastro::constants::R_earth; // distance in meters
  double T_ref=smartastro::constants::T_Earth; // time in seconds

  std::vector<bool> flag(4, false);
  flag[0] = true; // drag
  flag[1] = true; // lunisolar
  flag[2] = true; // SRP

  std::vector<bool> flag2 = flag;
  flag2[3] = true; // magnetic

  std::vector<double> params(5);
  params[0] = 0.1; // mass in kg
  params[1] = 1.0; // cross-section in square meters
  params[2] = 2.0; // drag coefficient
  params[3] = 187.0; // mean solar flux
  params[4] = 1.0; // geomagnetic index

  std::vector<double> other(2);
  other[0] = params[1];
  other[1] = 1.3;

  std::vector<double> x0(6), x00(6), x0bis(6), xfbis(6), x0prime(6), xfprime(6), xf(6), yf(6), yfbis(6), kep(6), car(6), modeq(6);// z0(5), zf(5);
  x0[0] = 8000.0e3 / R_ref;
  x0[1] = 1.0e-3;
  x0[2] = 89.0 * pi / 180.0;
  x0[3] = 60.0 * pi / 180.0;
  x0[4] = 30.0 * pi / 180.0;
  x0[5] = 100.0 * pi / 180.0;

  smartastro::astrocore::conversion_coordinates::oscul2mean(x0, 1.0, 1.0, kep);
  smartastro::astrocore::conversion_coordinates::kep2modeq(kep, x0bis);
  x0prime = x0bis; // initial mean modified equinoctial elements

  smartastro::astrodynamics::dearth_semianalytical *dyn = new smartastro::astrodynamics::dearth_semianalytical(flag, n_max, params, other);
  smartmath::integrator::rk4<double> prop_semianalytical(dyn);

  smartastro::astrodynamics::dearth_semianalytical *dyn2 = new smartastro::astrodynamics::dearth_semianalytical(flag2, n_max, params, other);
  smartmath::integrator::rk4<double> prop_semianalytical2(dyn2);  

  double mjd = 51624.5;
  double t = (mjd + smartastro::constants::mjd_zero_in_jd) * 3600.0 * 24.0 / T_ref;
  double dt = 24.0 * 3600.0 / T_ref;
  int steps = 4; // + floor(dt * T_ref / (12.0 * 3600.0)); 

  if(flag[2])
  {
    semianalytical.open("semianalytical_magnetic_polar_SRP.txt");
    semianalytical2.open("semianalytical_magnetic_polar2_SRP.txt");
  }
  else
  {
    semianalytical.open("semianalytical_magnetic_polar.txt");
    semianalytical2.open("semianalytical_magnetic_polar2.txt");
  }

  semianalytical << setprecision(16) << t * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi << " " << kep[5] * 180.0 /pi << " ";
  semianalytical << "\n";
  semianalytical2 << setprecision(16) << t * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi << " " << kep[5] * 180.0 /pi << " ";
  semianalytical2 << "\n";  

  int propagations = 10 * 365.25 * 1;
  for(int k = 0; k < propagations; k++)
  {
    prop_semianalytical.integrate(t, t + dt, steps, x0bis, xfbis);
    x0bis = xfbis;
    smartastro::astrocore::conversion_coordinates::modeq2kep(x0bis, kep);
    semianalytical << setprecision(16) << (t + dt) * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi << " " << kep[5] * 180.0 /pi << " ";
    semianalytical << "\n";

    prop_semianalytical2.integrate(t, t + dt, steps, x0prime, xfprime);
    x0prime = xfprime;
    smartastro::astrocore::conversion_coordinates::modeq2kep(x0prime, kep);
    semianalytical2 << setprecision(16) << (t + dt) * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi << " " << kep[5] * 180.0 /pi << " ";
    semianalytical2 << "\n";

    t += dt;
  }

  semianalytical.close();
  semianalytical2.close();

  delete dyn;
  delete dyn2;

  return 0;

}