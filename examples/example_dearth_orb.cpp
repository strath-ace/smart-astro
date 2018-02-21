/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example shows how to perform orbit propagation using the modified equinoctial elements */

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

double alti(std::vector<double> x, double jd){
    std::vector<double> pos(3), M(9), car(6);
    double rd, rk, phi, phi0, C, C0, e, h;

    smartastro::astrocore::conversion_coordinates::kep2car(x, 1.0, car);
    smartastro::astrocore::conversion_frames::inertial_to_tod(jd, M); // rotation matrix
    /** computing position in true of date */
    pos[0] = M[0] * car[0] + M[3] * car[1] + M[6] * car[2];
    pos[1] = M[1] * car[0] + M[4] * car[1] + M[7] * car[2];
    pos[2] = M[2] * car[0] + M[5] * car[1] + M[8] * car[2];
    // pos[0] = car[0]; pos[1] = car[1]; pos[2] = car[2];

    // std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;

    /** computing altitude */
    e = sqrt(1.0 - (1.0 - smartastro::constants::f_earth) * (1.0 - smartastro::constants::f_earth));
    rd = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    rk = pos[2];
    phi = atan(pos[2] / rd);
    C = (smartastro::constants::R_earth / R_ref) / sqrt(1.0 - e * e * sin(phi) * sin(phi));
    for(int i = 0; i < 2; i++){
        C0 = C;
        phi0 = phi;
        C = (smartastro::constants::R_earth / R_ref) / sqrt(1.0 - e * e * sin(phi0) * sin(phi0));
        phi = atan((rk + C0 * e * e * sin(phi0)) / rd);
    }
    h = rd / cos(phi) - C;
    h *= R_ref / 1000.0;
  
    return h;
}

int main(){

    std::ofstream orbital;
    double dt,s,t,jd0,jd,tt;
    int yr,mo,d,hr,min;
    double pi=smartastro::constants::pi;
    int stop;

    /** initial gregorian date */
    yr=2010;mo=1;d=1;hr=0;min=0;s=0.0;

    dt=2.0*3600.0/T_ref; // saving frequency

    jd0=367*yr - floor(7*(yr+floor((mo+9)/12))/4)+floor(275*mo/9)+d+1721013.5+((s/60+min)/60+hr)/24;
    jd=jd0;

    std::vector<double> modeq(6),modeq0(6),modeqf(6),p(4),kep0(6),kep(6);

    kep0[0]=7000.0e3 / R_ref; //7378
    kep0[1]=1.0e-2;
    kep0[2]=80.0*pi/180.0;
    kep0[3]=60.0*pi/180.0;
    kep0[4]=30.0*pi/180.0;
    kep0[5]=0.0*pi/180.0;
    smartastro::astrocore::conversion_coordinates::kep2modeq(kep0,modeq0);           

    double m0=1000.0; // mass

    p[0]=1.0e-2 * m0; // drag cross-section in square meters 
    p[1]=2.2; // drag coefficient
    p[2]=150.0; // F10.7
    p[3]=3.0; // Kp

    double t0=(3600.0*24.0*jd0)/T_ref;

    std::vector<int> flags(4,0);
    flags[0]=1; // drag

    std::vector<double> other(2,0.0);

    int n_gravity=0;

    int order_integration = 8;
    bool initializer_bulirschstoer = true;

    orbital.open ("orbital.txt");
    smartastro::astrodynamics::dearth_orb<double> dyn3dof(R_ref,T_ref,n_gravity,m0,flags,p,other);
    smartmath::integrator::ABM<double> prop3dof(&dyn3dof, order_integration, initializer_bulirschstoer);

    orbital << setprecision(16) << 0.0 << " " << kep0[0] * R_ref << " " << kep0[1] << " " << kep0[2]*180.0/pi << " " << kep0[3]*180.0/pi << " "  << kep0[4]*180.0/pi << " " << kep0[5]*180.0/pi << "\n" ;

    double delta_jd=30.0;
    int n_mini_step=floor((T_ref*dt)/30.0) + 1;
    double h = dt / double(n_mini_step);
    std::vector<std::vector<double> > f;

    // std::cout << alti(kep0, jd0) << std::endl;
    //std::vector<double> dx;
    //dyn3dof.evaluate(jd/T_ref*(3600.0*24.0), modeq0, dx);
    //return 0;

    stop=0;
    t=t0;
    jd=jd0;
    kep=kep0;       
    modeq=modeq0;
    prop3dof.initialize(order_integration, t, h, modeq, f);
    while(stop==0)
    {
        tt=t+dt;
        prop3dof.integrate(t,tt,n_mini_step,modeq,modeqf,f);
      
        modeq=modeqf;
        if(tt<t+dt)
            stop=1;       
        t=tt;
        jd=t*T_ref/(3600.0*24.0);
        modeqf[5] -= 2.0 * pi * double(floor(modeqf[5] / (2.0 * pi)));  
        smartastro::astrocore::conversion_coordinates::modeq2kep(modeqf,kep);

        orbital << setprecision(16) << (tt - t0) * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
        if(jd-jd0>delta_jd)
            stop=1;
    }

    orbital.close();

    return 0;
}