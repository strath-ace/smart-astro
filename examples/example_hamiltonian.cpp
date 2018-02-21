/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example shows how to use symplectic integrators for Hamiltonian Earth orbits (standard and with mixed variables) */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../include/smartastro.h"
#include "smartmath.h"

using namespace std;


 int main(){

    double R_ref;
    double pi=smartastro::constants::pi;

    int steps1, steps2;
    int n_max = 2;
    //bool drift_kick_drift = true;

    /** scaling parameters */
    R_ref = smartastro::constants::R_earth; // distance in meters

    std::vector<double> y0(6), x0(6), xf(6), yf(6), xf2(6), yf2(6), xf3(6), yf3(6);
     y0[0]= 8000.0e3 / R_ref;
     y0[1]= 1.0e-2;
     y0[2]= 40.0 * pi / 180.0;
     y0[3]= 100.0 * pi / 180.0;
     y0[4]= 60.0 * pi / 180.0;
     y0[5]= 200.0 * pi / 180.0; 
    smartastro::astrocore::conversion_coordinates::kep2car(y0, 1.0, x0);
    std::vector<double> x02 = x0;

    smartastro::astrodynamics::dearth_hamiltonian<double> *dyn = new smartastro::astrodynamics::dearth_hamiltonian<double>(n_max);
    smartmath::integrator::ABM<double> prop(dyn,6);

    smartastro::astrodynamics::dearth_mixedvar<double> *dyn2 = new smartastro::astrodynamics::dearth_mixedvar<double>(n_max);
    smartmath::integrator::forest_mixedvar<double> prop2(dyn2);

    double t = 0.0;
    double dt = 1.0e2;
    steps1 = floor(dt * smartastro::constants::T_Earth / (60.0)) + 1;
    steps2 = steps1 / 10;
    for(int k = 0; k < 10; k++){
      prop.integrate(t, t + dt, steps1, x0, xf);
      prop2.integrate(t, t + dt, steps2, x02, xf2);
      t += dt;
      x0 = xf;
      x02 = xf2; 
      smartastro::astrocore::conversion_coordinates::car2kep(xf, 1.0, yf);
      smartastro::astrocore::conversion_coordinates::car2kep(xf2, 1.0, yf2);
      std::cout << "Step " << k + 1 << std::endl;
      std::cout << std::setprecision(10) << "Propagation V1: " << R_ref * yf[0] / 1000.0 << " " << yf[1] << " " << yf[2] << " "  << yf[3] << " " << yf[4] << " " << yf[5] << std::endl;
      std::cout << std::setprecision(10) << "Propagation V2: " << R_ref * yf2[0] / 1000.0 << " " << yf2[1] << " " << yf2[2] << " "  << yf2[3] << " " << yf2[4] << " " << yf2[5] << std::endl;
    }

    delete dyn;
    delete dyn2;

    return 0;

}

