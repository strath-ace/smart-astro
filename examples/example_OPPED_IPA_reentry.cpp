/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example shows how to use OPPED-IPA until re-entry occurs (within a maximum time of flight) */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../include/smartastro.h"
#include "smartmath.h"

using namespace std;

/** scaling parameters */
const double R_ref=smartastro::constants::R_earth; // distance in meters
const double T_ref=smartastro::constants::T_Earth; // time in seconds
const double alt_min = 120.0; // in km

void alti(const double &jd, const std::vector<double> &x, double &h){
    std::vector<double> M(9), pos(3);
    double rd,rk,phi,phi0,C,C0,e2,Re;

    smartastro::astrocore::conversion_frames::inertial_to_tod(jd,M); // rotation matrix
    /** computing position in true of date */
    pos[0]=M[0]*x[0]+M[3]*x[1]+M[6]*x[2];
    pos[1]=M[1]*x[0]+M[4]*x[1]+M[7]*x[2];
    pos[2]=M[2]*x[0]+M[5]*x[1]+M[8]*x[2];

    /** computing altitude */
    e2=1.0-(1.0-smartastro::constants::f_earth)*(1.0-smartastro::constants::f_earth);
    Re=smartastro::constants::R_earth/R_ref;
    rd=sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    rk=pos[2];
    phi=atan(pos[2]/rd);
    C=Re/sqrt(1.0-e2*sin(phi)*sin(phi));
    for(int i=0; i<2; i++)
    {
        C0=C;
        phi0=phi;
        C=Re/sqrt(1.0-e2*sin(phi0)*sin(phi0));
        phi=atan((rk+C0*e2*sin(phi0))/rd);
    }
    h=rd/cos(phi)-C; // altitude in kilometers
    h*=R_ref/1000.0; // scaling
}

std::vector<int> reentry(std::vector<double> x, double jd){
    std::vector<int> v(1, 0);
    
    double h;
    alti(jd, x, h); // computing altitude
  
    if(h < alt_min) 
      v[0] = 1; // output is 1 if critical altitude is reached, 0 otherwise

    return v;
}

int main(){

  /* declarations */
  std::ofstream output_file;
  int saving_steps;
  double T_ref, R_ref, jd0, tof, dt, t0, t, tf, m;
  std::vector<double> x0(7), xf(7), p(1), kep0(6), kep(6), car(6), other(2);
  double pi = smartastro::constants::pi;

  tof = 100.0 * 24.0 * 3600.0 / T_ref; // scaled time of flight
  dt = 2.0 * 3600.0 / T_ref; // data saving frequency
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
  kep0[0]=6800.0e3 / R_ref; // scaled semi-major axis
  kep0[1]=1.0e-3; // eccentricity
  kep0[2]=89.0*pi/180.0; // inclination
  kep0[3]=60.0*pi/180.0; // right ascension of the ascending node
  kep0[4]=30.0*pi/180.0; // argument of periapsis
  kep0[5]=100.0*pi/180.0; // true anomaly

  m = 100.0; // initial mass in kilograms

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
  t = t0;
  tf = t+dt;
  smartastro::astrocore::conversion_coordinates::kep2car(kep0,1.0,car);   
  for(unsigned int i = 0; i < 6; i++)
    x0[i] = car[i];
  x0[6] = m;

  /** initialization of integrator */
  smartastro::astrodynamics::dearth_3dof<double> *dyn = new smartastro::astrodynamics::dearth_3dof<double>(R_ref, T_ref, n_max, flags, p, other);
  smartmath::integrator::rkf45<double> prop(dyn, 1.0e-9, 5.0, 1.0e-2 / T_ref, 60.0 / T_ref); /* second input is integration tolerance, third is max. step growth, fourth is minimum stepsize and fifth is max. stepsize */
  int n_mini_step = floor(dt * T_ref / 60.0); // guess for integrations step

  /** storing initial conditions */
  output_file.open("OPPED_IPA_reentry.txt");
  output_file << setprecision(16) << 0.0 << " " << kep0[0]*R_ref <<  " " << kep0[1] << " " << kep0[2]*180.0/pi << " " << kep0[3]*180.0/pi << " " << kep0[4]*180.0/pi << " " << kep0[5]*180.0/pi << "\n" ; 

  int j = 0;
  while(j < saving_steps)
  {
    prop.integrate(t,tf,n_mini_step,x0,xf,reentry); // integration with even-handler for re-entry

    if(tf == t+dt)
    {
        x0=xf;    
        t=tf;
        tf+=dt;
        j++;
    }
    else
        j = saving_steps;

    for(unsigned int i=0; i<6; i++)
        car[i]=xf[i];
    smartastro::astrocore::conversion_coordinates::car2kep(car,1.0,kep);     
    output_file << setprecision(16) << (t - t0) * T_ref << " " << kep[0]*R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
  }

  delete dyn;
  output_file.close();

  return 0;
}