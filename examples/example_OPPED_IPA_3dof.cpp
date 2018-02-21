/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example shows how to use the 3 degrees-of-freedom version of the Orbit Propagation with PrecisE Dynamics - and Intrusive Polynomial Approximation */

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

  /* declarations */
  std::ofstream output_file;
  int saving_steps;
  double T_ref, R_ref, jd0, tof, dt, t0, t, m;
  std::vector<double> x0(7), xf(7), p(1), kep0(6), kep(6), car(6), other(2, 0.0);
  double pi = smartastro::constants::pi;

  /** scaling parameters */
  R_ref = smartastro::constants::R_earth; // distance in meters
  T_ref = smartastro::constants::T_Earth; // time in seconds

  tof = 10.0 * 24.0 * 3600.0 / T_ref; // scaled time of flight
  dt = 120.0 / T_ref; // data saving frequency
  saving_steps = floor(tof / dt); // number of saved intermediate states

  int n_max = smartastro::constants::n_Earth_gravity; // degree of EGM96 harmonics used 
  std::vector<int> flags(4, 0); // flags for other orbital perturbations 
  flags[0] = 1; // drag
  flags[1] = 1; // lunisolar
  flags[2] = 1; // SRP with conical shadow

  /** initial gregorian date */
  std::vector<double> date0(6);
  date0[0] = 2013.0; // year
  date0[1] = 10.0; // month
  date0[2] = 22.0; // day
  date0[3] = 3.0; // hour
  date0[4] = 0.0; // min
  date0[5] = 0.0; // s

  /** initial Keplerian coordinates */
  kep0[0]=8000.0e3 / R_ref; // scaled semi-major axis
  kep0[1]=1.0e-3; // eccentricity
  kep0[2]=89.0*pi/180.0; // inclination
  kep0[3]=60.0*pi/180.0; // right ascension of the ascending node
  kep0[4]=30.0*pi/180.0; // argument of periapsis
  kep0[5]=100.0*pi/180.0; // true anomaly

  m = 1100.0; // initial mass in kilograms

  /** drag-related parameters */
  p[0] = 15.0; // drag cross-section in square meters 
  p.push_back(2.0); // drag coefficient
  p.push_back(106.4); // mean solar flux in solar flux units
  p.push_back(3.85); // logarithmic geomagnetic index

  other[0] = 1.625; // SRP cross-section in square meters 
  other[1] = 1.3; // SRP coefficient

  /** conversions **/
  smartastro::astrocore::conversion_time::date2jd(date0, jd0);
  t0 = jd0 * 3600.0 * 24.0 / T_ref;
  t = t0+dt;
  smartastro::astrocore::conversion_coordinates::kep2car(kep0,1.0,car);   
  for(unsigned int i = 0; i < 6; i++)
    x0[i] = car[i];
  x0[6] = m;

  /** initialization of integrator */
  smartastro::astrodynamics::dearth_3dof<double> *dyn = new smartastro::astrodynamics::dearth_3dof<double>(R_ref, T_ref, n_max, flags, p, other);
  int order_integration = 8;
  smartmath::integrator::ABM<double> prop(dyn, order_integration, true); // the boolean is true a Bulirsh-Stoer of same order is used to initialize the predictor-corrector, false if RK4 is used instead
  std::vector<std::vector<double> > f; // vector to save the multiple steps
  int n_mini_step = floor(dt * T_ref / 30.0); // integration step
  prop.initialize(order_integration, t0, dt / double(n_mini_step), x0, f); // initializing f

  /** storing initial conditions */
  output_file.open("OPPED_IPA_3dof.txt");
  output_file << setprecision(16) << 0.0 << " " << kep0[0]*R_ref <<  " " << kep0[1] << " " << kep0[2]*180.0/pi << " " << kep0[3]*180.0/pi << " " << kep0[4]*180.0/pi << " " << kep0[5]*180.0/pi << "\n" ; 

  for(int j = 0; j < saving_steps; j++)
  {
    prop.integrate(t,t+dt,n_mini_step,x0,xf,f); // integration re-using previously computed multi-steps (gain of time w.r.t. prop.integrate(t,t+dt,n_mini_step,x0,xf))

    x0=xf;    
    t+=dt;

    for(unsigned int i=0; i<6; i++)
        car[i]=xf[i];
    smartastro::astrocore::conversion_coordinates::car2kep(car,1.0,kep);     
    output_file << setprecision(16) << (t - t0) * T_ref << " " << kep[0]*R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
  }

  delete dyn;
  output_file.close();

  return 0;
}