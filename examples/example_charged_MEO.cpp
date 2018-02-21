/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example is for the MEO case of a submitted paper to AIAA JGCD on the effects of the Lorentz force */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../include/smartastro.h"
#include "smartmath.h"

using namespace smartmath;
using namespace std;

/** scaling parameters */
const double R_ref=smartastro::constants::R_earth; // distance in meters
const double T_ref=smartastro::constants::T_Earth; // time in seconds
const double alt_min = 120.0; // in km

void alti(const double &jd, const std::vector<double> &modeq, double &h){
    std::vector<double> M(9), pos(3), kep(6), x(6);
    double rd,rk,phi,phi0,C,C0,e2,Re;

    smartastro::astrocore::conversion_coordinates::modeq2kep(modeq,kep);
    smartastro::astrocore::conversion_coordinates::kep2car(kep,1.0,x);
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
    for(int i=0; i<2; i++){
        C0=C;
        phi0=phi;
        C=Re/sqrt(1.0-e2*sin(phi0)*sin(phi0));
        phi=atan((rk+C0*e2*sin(phi0))/rd);
    }
    h=rd/cos(phi)-C;
    h*=R_ref/1000.0;
}

std::vector<int> reentry(std::vector<double> x, double jd){
    std::vector<int> v(1, 0);
    
    double h;
    alti(jd, x, h);
  
    if(h < alt_min)
      v[0] = 1;

    return v;
}

int main(){

    std::ofstream altitudes, orbital;
    std::ofstream altitudes2, orbital2;
    std::ofstream electrostatic2;
    double dt,s,t,jd0,jd,h,tt;
    int yr,mo,d,hr,min;
    double pi=smartastro::constants::pi;
    int stop;

    /** initial gregorian date */
    yr=2000;mo=3;d=21;hr=12;min=0;s=0.0;

    dt=24.0*3600.0/T_ref; // saving frequency

    jd0=367*yr - floor(7*(yr+floor((mo+9)/12))/4)+floor(275*mo/9)+d+1721013.5+((s/60+min)/60+hr)/24;
    jd=jd0;

    std::vector<double> x(6),x0(6),xf(6),p(4),p2(4),p3(4),p4(4),kep0(6),kep(6);

    kep0[0]=20000.0e3 / R_ref; 
    kep0[1]=1.0e-3;
    kep0[2]=89.0*pi/180.0;
    kep0[3]=60.0*pi/180.0;
    kep0[4]=30.0*pi/180.0;
    kep0[5]=100.0*pi/180.0;
          
    p[0]=pi * 1.0e-12; // 1.0
    p[1]=2.0;
    p[2]=187.0;
    p[3]=1.0;
   
    p2=p;

    double m0= (4.0 / 3.0) * pi * pow(1.0e-6, 3) * 1000.0;// 0.1;

    double t0=(3600.0*24.0*jd0)/T_ref;

    std::vector<int> flags(4,0);
    flags[0]=1; // drag
    flags[1]=1; // lunisolar
    flags[2]=1; // SRP
    flags[3]=0;
    std::vector<int> flags2=flags;
    flags2[3]=1;

    std::vector<double> other(2,0.0);
    other[1]=p[0];
    other[1]=1.3;

    double C, V; // charging quantities
    int n_gravity=14;
    int n_abm = 8;

    altitudes.open ("altitudes_charged_polarMEO_dust.txt");
    orbital.open ("orbital_charged_polarMEO_dust.txt");
    smartastro::astrodynamics::dearth_orb<double> dyn3dof(R_ref,T_ref,n_gravity,m0,flags,p,other);
    smartmath::integrator::rkf45<double> prop3dof(&dyn3dof, 1.0e-9, 5.0, 1.0e-2 / T_ref, 90.0 / T_ref);

    altitudes2.open ("altitudes2_charged_polarMEO_dust.txt");
    orbital2.open ("orbital2_charged_polarMEO_dust.txt");
    electrostatic2.open("electrostatic2_charged_polarMEO_dust.txt");
    smartastro::astrodynamics::dearth_orb<double> dyn3dof2(R_ref,T_ref,n_gravity,m0,flags2,p2,other);
    //smartmath::integrator::ABM<double> prop3dof2(&dyn3dof2, n_abm, true);
    smartmath::integrator::rkf45<double> prop3dof2(&dyn3dof2, 1.0e-9, 5.0, 1.0e-2 / T_ref, 90.0 / T_ref);

    double delta_jd=100.0*365.25;
    int n_mini_step=floor((T_ref*dt)/90.0) + 1;
    std::vector<std::vector<double> > f;

    stop=0;
    t=t0;
    jd=jd0;
    kep=kep0;
    smartastro::astrocore::conversion_coordinates::kep2modeq(kep,x0);
    x=x0;
    //prop3dof.initialize(n_abm, t0, dt / double(n_mini_step), x0, f);
    while(stop==0)
    {
        tt=t+dt;
        prop3dof.integrate(t,tt,n_mini_step,x,xf,reentry);
        x=xf;
        if(tt<t+dt)
            stop=1;       
        t=tt;
        jd=t*T_ref/(3600.0*24.0);
        smartastro::astrocore::conversion_coordinates::modeq2kep(x,kep);
        alti(jd, x, h);

        orbital << setprecision(16) << (tt - t0) * T_ref << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
        altitudes << h << "\n" ;

        if(jd-jd0>delta_jd)
            stop=1;
    }

    altitudes.close();
    orbital.close();

    std::cout << "simulation 1 over 2 done" << std::endl;

    stop=0;
    t=t0;
    jd=jd0;
    kep=kep0;
    smartastro::astrocore::conversion_coordinates::kep2modeq(kep,x0);
    x=x0;
    //prop3dof2.initialize(n_abm, t0, dt / double(n_mini_step), x0, f);
    while(stop==0)
    {
        tt=t+dt;
        prop3dof2.integrate(t,tt,n_mini_step,x,xf,reentry);
        x=xf;
        if(tt<t+dt)
            stop=1;       
        t=tt;
        jd=t*T_ref/(3600.0*24.0);
        smartastro::astrocore::conversion_coordinates::modeq2kep(x,kep);
        alti(jd, x, h);
        dyn3dof2.charging(jd,xf,C,V);

        orbital2 << setprecision(16) << (tt - t0) * T_ref << " "  << kep[0] * R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
        altitudes2 << h << "\n" ;
        electrostatic2 << V << " " << C*V << "\n" ;

        if(jd-jd0>delta_jd)
            stop=1;
    }

    altitudes2.close();
    orbital2.close();
    electrostatic2.close();

    return 0;
}