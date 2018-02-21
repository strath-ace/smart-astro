/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example compares first-order analytical propagation of Earth orbits with symplectic integration */


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

  std::ofstream analytical, numerical;
  clock_t begin, end;

  double pi=smartastro::constants::pi;

  int steps;
  int n_max = 0;

  /** scaling parameters */
  double R_ref=smartastro::constants::R_earth; // distance in meters
  double T_ref=smartastro::constants::T_Earth; // time in seconds

  std::vector<double> x0(6), x00(6), x0bis(6), xfbis(6), xf(6), y0(6), yf(6), yfbis(6), y0bis(6), kep(6), car(6);
  x0[0] = 7000.0e3 / R_ref;
  x0[1] = 5.0e-4;
  x0[2] = 80.0 * pi / 180.0;
  x0[3] = 10.0 * pi / 180.0;
  x0[4] = 60.0 * pi / 180.0;
  x0[5] = 200.0 * pi / 180.0;
  // double T_period = (2.0 * smartastro::constants::pi) * pow(x0[0], 1.5);
  smartastro::astrocore::conversion_coordinates::kep2car(x0, 1.0, y0);
  kep = x0;
  y0bis = y0;
  smartastro::astrocore::conversion_coordinates::kep2modeq(x0, xf);
  x0 = xf;
  x00 = x0;
  x0bis = x0;
  double L0 = x0[5], L = L0;
  std::vector<double> modeq(x0);

  std::vector<bool> flag(1, true); // lunisolar ON or OFF

  smartastro::astrodynamics::dearth_hamiltonian<double> *dyn = new smartastro::astrodynamics::dearth_hamiltonian<double>(n_max, flag);
  smartmath::integrator::yoshida6<double> prop_numerical(dyn);    

  smartastro::propagator::perturbation_propagator prop_perturbed(flag, n_max); 

  double jd = smartastro::constants::mjd_zero_in_jd + 53000.0;
  prop_perturbed.set_timing(L0, jd);

  analytical.open ("analytical.txt");
  numerical.open ("numerical.txt");

  double t = jd * 3600.0 * 24.0 / T_ref;

  // numerical << setprecision(16) << t * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi  << " " << kep[5] * 180.0 /pi << " ";
  // analytical << setprecision(16) << t * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi  << " " << kep[5] * 180.0 /pi << " ";
  numerical << setprecision(16) << t * T_ref / (3600.0 * 24.0) << " " << modeq[0] * R_ref << " " << modeq[1] << " " << modeq[2] << " " << modeq[3]  << " " << modeq[4]  << " " << modeq[5] * 180.0 /pi << " ";
  analytical << setprecision(16) << t * T_ref / (3600.0 * 24.0) << " " << modeq[0] * R_ref << " " << modeq[1] << " " << modeq[2] << " " << modeq[3]  << " " << modeq[4]  << " " << modeq[5] * 180.0 /pi << " ";  

  numerical << "\n";
  analytical << "\n";

  double time_elapsed1 = 0.0, time_elapsed2 = 0.0;

  int counts_L = 0, counts_L0 = 0;

  double dt = 1.0e-1;
  steps = 1 + int(floor(fabs(dt) * T_ref / 60.0));
  int propagations = int(floor(1.0e2 / dt)); // 1e2;
  for(int k = 0; k < propagations; k++)
  {
    begin = clock();
    prop_numerical.integrate(t, t + dt, steps, y0, yf);
    end = clock();
    time_elapsed1 += (double (end-begin)) / CLOCKS_PER_SEC;

    // prop_numerical.integrate(t, t + dt, steps, y0bis, yfbis);

    smartastro::astrocore::conversion_coordinates::car2kep(y0, 1.0, kep);
    smartastro::astrocore::conversion_coordinates::kep2modeq(kep, x0);
    smartastro::astrocore::conversion_coordinates::car2kep(yf, 1.0, kep);
    smartastro::astrocore::conversion_coordinates::kep2modeq(kep, xf); 

    /* business to keep track of actual true longitudes (without the modulo) */
    if(dt > 0.0)
    {
      if(xf[5] + 2.0 * pi * double(counts_L) < L)
        counts_L++;
      L = xf[5] + 2.0 * pi * double(counts_L);

      if(x0[5] + 2.0 * pi * double(counts_L0) < L0)
        counts_L0++;
      L0 = x0[5] + 2.0 * pi * double(counts_L0);  
    }
    else
    {
      if(xf[5] - 2.0 * pi * double(counts_L) > L)
        counts_L++;
      L = xf[5] - 2.0 * pi * double(counts_L);

      if(x0[5] - 2.0 * pi * double(counts_L0) > L0)
        counts_L0++;
      L0 = x0[5] - 2.0 * pi * double(counts_L0);       
    }    
    xf[5] = L;
    x0[5] = L0;

    // numerical << setprecision(16) << (t + dt) * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi  << " " << kep[5] * 180.0 /pi << " ";

    smartastro::astrocore::conversion_coordinates::kep2modeq(kep, modeq);

    // std::cout << "Equinoctial elements propagation number " << k + 1 << " compared between mixed-variables, standard numerical and analytical" << std::endl;
    // std::cout << std::setprecision(10) << R_ref * xf[0] / 1000.0 << " " << xf[1] << " " << xf[2] << " "  << xf[3] << " " << xf[4] << " " << xf[5] << std::endl;

    // smartastro::astrocore::conversion_coordinates::car2kep(yfbis, 1.0, kep);
    // smartastro::astrocore::conversion_coordinates::kep2modeq(kep, xf); 

    // std::cout << std::setprecision(10) << R_ref * xf[0] / 1000.0 << " " << xf[1] << " " << xf[2] << " "  << xf[3] << " " << xf[4] << " " << xf[5] << std::endl;

    jd += dt * T_ref / (3600.0 * 24.0);
    prop_perturbed.set_timing(L0, jd);
    x0bis[5] = L0;
    begin = clock();
    prop_perturbed.propagate(L0, L, x0bis, xfbis); // prop_perturbed.propagate(L00, L, x00, xfbis)
    end = clock();
    time_elapsed2 += (double (end-begin)) / CLOCKS_PER_SEC;
    x0bis = xfbis;

    // std::cout << std::setprecision(10) << R_ref * xfbis[0] / 1000.0 << " " << xfbis[1] << " " << xfbis[2] << " "  << xfbis[3] << " " << xfbis[4] << " " << xfbis[5] << std::endl;

    smartastro::astrocore::conversion_coordinates::modeq2kep(xfbis, kep); 
    smartastro::astrocore::conversion_coordinates::kep2car(kep, 1.0, car);
    // analytical << setprecision(16) << (t + dt) * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi  << " " << kep[5] * 180.0 /pi << " ";
    
    numerical << setprecision(16) << t * T_ref / (3600.0 * 24.0) << " " << modeq[0] * R_ref << " " << modeq[1] << " " << modeq[2] << " " << modeq[3]  << " " << modeq[4]  << " " << modeq[5] * 180.0 / pi << " ";
    analytical << setprecision(16) << t * T_ref / (3600.0 * 24.0) << " " << xfbis[0] * R_ref << " " << xfbis[1] << " " << xfbis[2] << " " << xfbis[3]  << " " << xfbis[4]  << " " << (xfbis[5] - 2.0 * pi * floor(xfbis[5] / (2.0 * pi))) * 180.0 / pi << " ";  

    numerical << "\n";
    analytical << "\n";

    t += dt;
    y0 = yf;
    y0bis = yfbis;

    // std::cout << std::endl;
  }

  cout << "Time to obtain final state with numerical : " << time_elapsed1 << " VS analytical propagator: " << time_elapsed2 << endl;

  numerical.close();
  analytical.close();

  delete dyn;

  return 0;

}

