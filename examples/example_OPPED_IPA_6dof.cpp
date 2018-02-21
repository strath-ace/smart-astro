/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example shows how to use the 6 degrees-of-freedom version of OPPED-IPA */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../include/smartastro.h"
#include "smartmath.h"

using namespace std;

int main()
{

  /* declarations */
  std::ofstream output_file;
  int saving_steps;
  double T_ref, R_ref, jd0, tof, dt, t0, t, m;
  std::vector<double> x0(14), xf(14), p(2), kep0(6), kep(6), car(6), other(2, 0.0);
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
  flags[2] = 2; // SRP with cylindrical shadow

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

  /** initial attitude */
  double beta=50.0*smartastro::constants::pi/180.0;
  x0[7]=0.0; // quaternion 1
  x0[8]=0.5*sqrt(3.0)*sin(0.5*beta); // quaternion 2
  x0[9]=0.5*sin(0.5*beta); // quaternion 3
  x0[10]=cos(0.5*beta); // quaternion 4
  x0[11]=1.0e-2*T_ref; // angular rate 1
  x0[12]=-2.0e-2*T_ref; // angular rate 2
  x0[13]=3.0e-2*T_ref; // angular rate 3

  /* mass parameters */
  m = 1100.0; // initial mass in kilograms
  int geometry = 0; // 0 for flag for square flat plate (1 for GOCE-like cylinder and 2 for sphere)
  double a = 1.0; // length in meters
  double w = a / 100.0; // width in meters
  std::vector<double> I(3);
  I[0]=(m/6.0)*a*a; // moment of inertia 1
  I[1]=(m/12.0)*(a*a+w*w); // moment of inertia 2
  I[2]=(m/12.0)*(a*a+w*w); // moment of inertia 3
  std::vector<double> CoS(3,0.0); // displacement vector for center of symmetry w.r.t. barycenter

  p[0] = 4.0; // cross-section in square meters 
  p[1] = 2.0; // reference length

  /** atmospheric parameters */
  p.push_back(106.4); // mean solar flux in solar flux units
  p.push_back(3.85); // logarithmic geomagnetic index

  other[0] = 0.3; // specular reflectivity
  other[1] = 0.0; // diffusive reflexion

  /** conversions **/
  smartastro::astrocore::conversion_time::date2jd(date0, jd0);
  t0 = jd0 * 3600.0 * 24.0 / T_ref;
  t = t0+dt;
  smartastro::astrocore::conversion_coordinates::kep2car(kep0,1.0,car);   
  for(unsigned int i = 0; i < 6; i++)
    x0[i] = car[i];
  x0[6] = m;

  /** initialization of integrator */
  smartastro::astrodynamics::dearth_6dof<double> *dyn = new smartastro::astrodynamics::dearth_6dof<double>(I, R_ref, T_ref, n_max, flags, p, geometry, CoS, other);
  int order_integration = 8;
  smartmath::integrator::ABM<double> prop(dyn, order_integration, true); // the boolean is true a Bulirsh-Stoer of same order is used to initialize the predictor-corrector, false if RK4 is used instead
  std::vector<std::vector<double> > f; // vector to save the multiple steps
  int n_mini_step = floor(dt * T_ref / 5.0); // integration step
  prop.initialize(order_integration, t0, dt / double(n_mini_step), x0, f); // initializing f

  /** storing initial conditions */
  output_file.open("OPPED_IPA_6dof.txt");
  output_file << setprecision(16) << 0.0 << " " << x0[0]*R_ref <<  " " << x0[1]*R_ref << " " << x0[2]*R_ref << " " << x0[3]*R_ref/T_ref << " " << x0[4]*R_ref/T_ref << " " << x0[5]*R_ref/T_ref << " "; 
  output_file << setprecision(16) << x0[7] <<  " " << x0[8] << " " << x0[9] << " " << x0[10] << " " << x0[11]/T_ref << " " << x0[12]/T_ref << " " << x0[13]/T_ref << "\n" ;

  for(int j = 0; j < saving_steps; j++)
  {
    prop.integrate(t,t+dt,n_mini_step,x0,xf,f); // integration re-using previously computed multi-steps (gain of time w.r.t. prop.integrate(t,t+dt,n_mini_step,x0,xf))

    x0=xf;    
    t+=dt;
    
    output_file << setprecision(16) << (t - t0) * T_ref << xf[0]*R_ref <<  " " << xf[1]*R_ref << " " << xf[2]*R_ref << " " << xf[3]*R_ref/T_ref << " " << xf[4]*R_ref/T_ref << " " << xf[5]*R_ref/T_ref << " "; 
    output_file << setprecision(16) << xf[7] <<  " " << xf[8] << " " << xf[9] << " " << xf[10] << " " << xf[11]/T_ref << " " << xf[12]/T_ref << " " << xf[13]/T_ref << "\n" ;
  }

  delete dyn;
  output_file.close();

  return 0;
}