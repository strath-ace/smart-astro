/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example is for the 3dof case of a submitted paper to AIAA JGCD on the effects of the Lorentz force */

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

int main(){

    std::ofstream orbital;
    std::ofstream orbital2;
    std::ofstream orbital3;
    std::ofstream orbital4;
    std::ofstream electrostatic2, electrostatic4;
    double dt,s,t,jd0,jd,tt;
    int yr,mo,d,hr,min;
    double pi=smartastro::constants::pi;
    int stop;

    /** initial gregorian date */
    yr=2000;mo=3;d=21;hr=12;min=0;s=0.0;
    //yr=2024;mo=3;d=21;hr=12;min=0;s=0.0;
    //yr=2018;mo=3;d=21;hr=12;min=0;s=0.0;

    dt=2.0*3600.0/T_ref; // saving frequency

    jd0=367*yr - floor(7*(yr+floor((mo+9)/12))/4)+floor(275*mo/9)+d+1721013.5+((s/60+min)/60+hr)/24;
    jd=jd0;

    std::vector<double> modeq(6),modeq0(6),modeqf(6),p(4),kep0(6),kep(6);

    kep0[0]=8000.0e3 / R_ref; //7378
    kep0[1]=1.0e-3;
    kep0[2]=89.0*pi/180.0;
    kep0[3]=60.0*pi/180.0;
    kep0[4]=30.0*pi/180.0;
    kep0[5]=100.0*pi/180.0;
    smartastro::astrocore::conversion_coordinates::kep2modeq(kep0,modeq0);           

    p[0]=1.0;//23.6e0;
    p[1]=2.0;
    p[2]=187.0;
    p[3]=1.0;//2.33;
   
    double m0=p[0] / 10.0;//23.6; // mass

    double t0=(3600.0*24.0*jd0)/T_ref;

    std::vector<int> flags(4,0);
    flags[0]=1; // drag
    flags[1]=1; // lunisolar
    std::vector<int> flags2 = flags;
    flags2[3]=1; // magnetic
    std::vector<int> flags3 = flags;
    flags3[2]=1; // SRP
    std::vector<int> flags4 = flags3;
    flags4[3]=1;

    std::vector<double> other(2,0.0);
    other[0]=p[0];
    other[1]=1.3;

    double C, V; // charging quantities
    int n_gravity=3;

    int order_integration = 8;
    bool initializer_bulirschstoer = true;

    orbital.open ("orbital_charged_orb_LEO.txt");
    smartastro::astrodynamics::dearth_orb<double> dyn3dof(R_ref,T_ref,n_gravity,m0,flags,p,other);
    smartmath::integrator::ABM<double> prop3dof(&dyn3dof, order_integration, initializer_bulirschstoer);

    orbital2.open ("orbital2_charged_orb_LEO.txt");
    electrostatic2.open("electrostatic2_charged_orb_LEO.txt");
    smartastro::astrodynamics::dearth_orb<double> dyn3dof2(R_ref,T_ref,n_gravity,m0,flags2,p,other);
    smartmath::integrator::ABM<double> prop3dof2(&dyn3dof2, order_integration, initializer_bulirschstoer);

    orbital3.open ("orbital3_charged_orb_LEO.txt");
    smartastro::astrodynamics::dearth_orb<double> dyn3dof3(R_ref,T_ref,n_gravity,m0,flags3,p,other);
    smartmath::integrator::ABM<double> prop3dof3(&dyn3dof3, order_integration, initializer_bulirschstoer);

    orbital4.open ("orbital4_charged_orb_LEO.txt");
    electrostatic4.open("electrostatic4_charged_orb_LEO.txt");
    smartastro::astrodynamics::dearth_orb<double> dyn3dof4(R_ref,T_ref,n_gravity,m0,flags4,p,other);
    smartmath::integrator::ABM<double> prop3dof4(&dyn3dof4, order_integration, initializer_bulirschstoer);   

    orbital << setprecision(16) << 0.0 << " " << kep0[0] * R_ref << " " << kep0[1] << " " << kep0[2]*180.0/pi << " " << kep0[3]*180.0/pi << " "  << kep0[4]*180.0/pi << " " << kep0[5]*180.0/pi << "\n" ;
    orbital2 << setprecision(16) << 0.0 << " " << kep0[0] * R_ref << " " << kep0[1] << " " << kep0[2]*180.0/pi << " " << kep0[3]*180.0/pi << " "  << kep0[4]*180.0/pi << " " << kep0[5]*180.0/pi << "\n" ;
    orbital3 << setprecision(16) << 0.0 << " " << kep0[0] * R_ref << " " << kep0[1] << " " << kep0[2]*180.0/pi << " " << kep0[3]*180.0/pi << " "  << kep0[4]*180.0/pi << " " << kep0[5]*180.0/pi << "\n" ;
    orbital4 << setprecision(16) << 0.0 << " " << kep0[0] * R_ref << " " << kep0[1] << " " << kep0[2]*180.0/pi << " " << kep0[3]*180.0/pi << " "  << kep0[4]*180.0/pi << " " << kep0[5]*180.0/pi << "\n" ;

    dyn3dof2.charging(jd0,modeq0,C,V);
    electrostatic2 << V << " " << C * V << "\n" ;
    dyn3dof4.charging(jd0,modeq0,C,V);
    electrostatic4 << V << " " << C * V << "\n" ;

    double delta_jd=10.0*365.25;
    int n_mini_step=floor((T_ref*dt)/30.0) + 1;
    double h = dt / double(n_mini_step);
    std::vector<std::vector<double> > f;

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

        orbital << setprecision(16) << (tt - t0) * T_ref << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;

        if(jd-jd0>delta_jd)
            stop=1;
    }

    orbital.close();
    std::cout << "simulation 1 over 4 done" << std::endl;

    stop=0;
    t=t0;
    jd=jd0;
    kep=kep0;
    modeq=modeq0;
    prop3dof2.initialize(order_integration, t, h, modeq, f);
    while(stop==0)
    {
        tt=t+dt;
        prop3dof2.integrate(t,tt,n_mini_step,modeq,modeqf);
        modeq=modeqf;
        if(tt<t+dt)
            stop=1;       
        t=tt;
        jd=t*T_ref/(3600.0*24.0);
        modeqf[5] -= 2.0 * pi * double(floor(modeqf[5] / (2.0 * pi)));  
        smartastro::astrocore::conversion_coordinates::modeq2kep(modeqf,kep);

        dyn3dof2.charging(jd,modeqf,C,V);

        orbital2 << setprecision(16) << (tt - t0) * T_ref << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
        electrostatic2 << V << " " << C * V << "\n" ;

        if(jd-jd0>delta_jd)
            stop=1;
    }

    orbital2.close();
    electrostatic2.close();
    std::cout << "simulation 2 over 4 done" << std::endl;

    stop=0;
    t=t0;
    jd=jd0;
    kep=kep0;
    modeq=modeq0;
    prop3dof3.initialize(order_integration, t, h, modeq, f);
    while(stop==0)
    {
        tt=t+dt;
        prop3dof3.integrate(t,tt,n_mini_step,modeq,modeqf);
        modeq=modeqf;
        if(tt<t+dt)
            stop=1;       
        t=tt;
        jd=t*T_ref/(3600.0*24.0);
        modeqf[5] -= 2.0 * pi * double(floor(modeqf[5] / (2.0 * pi)));  
        smartastro::astrocore::conversion_coordinates::modeq2kep(modeqf,kep);

        orbital3 << setprecision(16) << (tt - t0) * T_ref << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;

        if(jd-jd0>delta_jd)
            stop=1;
    }

    orbital3.close();
    std::cout << "simulation 3 over 4 done" << std::endl;

    stop=0;
    t=t0;
    jd=jd0;
    kep=kep0;
    modeq=modeq0;
    prop3dof4.initialize(order_integration, t, h, modeq, f);
    while(stop==0)
    {
        tt=t+dt;
        prop3dof4.integrate(t,tt,n_mini_step,modeq,modeqf);
        modeq=modeqf;
        if(tt<t+dt)
            stop=1;       
        t=tt;
        jd=t*T_ref/(3600.0*24.0);
        modeqf[5] -= 2.0 * pi * double(floor(modeqf[5] / (2.0 * pi)));  
        smartastro::astrocore::conversion_coordinates::modeq2kep(modeqf,kep);
        dyn3dof4.charging(jd,modeqf,C,V);

        orbital4 << setprecision(16) << (tt - t0) * T_ref << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
        electrostatic4 << V << " " << C * V << "\n" ;

        if(jd-jd0>delta_jd)
            stop=1;
    }

    orbital4.close();
    electrostatic4.close();  
    std::cout << "simulation 4 over 4 done" << std::endl;

    return 0;
}
