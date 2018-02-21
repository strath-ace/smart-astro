/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example is for the 6dof case of a submitted paper to AIAA JGCD on the effects of the Lorentz force */

#include "smartmath.h"

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
const double R_ref=1000.0; // distance in meters
const double T_ref=1.0; // time in seconds

 int main(){

    std::ofstream orbital6dof, positions6dof, attitude;
    std::ofstream orbital6dof2, positions6dof2, attitude2;
    std::ofstream miscellaneous, miscellaneous2;
    double dt,s,t,jd0,jd,h,tt;
    int yr,mo,d,hr,min;
    double pi=smartastro::constants::pi;
    int stop;

    /** initial gregorian date */
    yr=2000;mo=3;d=21;hr=12;min=0;s=0.0;
    //yr=2017;mo=6;d=21;hr=12;min=0;s=0.0;
    //yr=2024;mo=3;d=21;hr=12;min=0;s=0.0;

    dt = 2.0 * 3600.0 / T_ref; // saving frequency

    jd0=367*yr - floor(7*(yr+floor((mo+9)/12))/4)+floor(275*mo/9)+d+1721013.5+((s/60+min)/60+hr)/24;
    //cout << setprecision(16) << "initial modified Julian date is " << jd0-smartastro::constants::mjd_zero_in_jd << endl;
    jd=jd0;

    std::vector<double> M(9),kep(6),car(6),other(2);
    std::vector<double> x(7),x0(7),xf(7),p(4),p2(4);

    // kep[0]=6900.0e3;
    // kep[1]=1.0e-3;
    // kep[2]=30.0*pi/180.0;
    // kep[3]=0.0*pi/180.0;
    // kep[4]=0.0*pi/180.0;
    // kep[5]=0.0*pi/180.0;

    kep[0]=8000.0e3; //7378
    kep[1]=1.0e-3;
    kep[2]=89.0*pi/180.0;
    kep[3]=60.0*pi/180.0;
    kep[4]=30.0*pi/180.0;
    kep[5]=100.0*pi/180.0;
    smartastro::astrocore::conversion_coordinates::kep2car(kep,smartastro::constants::mu_earth,car);
    x0[0]=car[0]/R_ref;
    x0[1]=car[1]/R_ref;
    x0[2]=car[2]/R_ref;
    x0[3]=car[3]*(T_ref/R_ref);
    x0[4]=car[4]*(T_ref/R_ref);
    x0[5]=car[5]*(T_ref/R_ref);
    x0[6]=1.0e-1; 

    double a=sqrt(1.0e0); // length in meters
    double w=a/100.0; // width in meters
    double ratio_m=0.7;

    p[0]=a*a;
    p[1]=a;
    p[2]=187.0;
    p[3]=1.0;

    p2=p;

    std::vector<double> I(3);
    // I[0]=(x0[6]/6.0)*a*a;
    // I[1]=(x0[6]/12.0)*(a*a+w*w);
    // I[2]=(x0[6]/12.0)*(a*a+w*w); 
    I[0]=(x0[6]/4.0)*a*a*(5.0/12.0+ratio_m*(1.0-ratio_m));
    I[1]=(x0[6]/12.0)*(a*a+w*w);
    I[2]=(x0[6]/12.0)*(a*a/4.0+w*w+3.0*a*a*ratio_m*(1.0-ratio_m));

    // I[0]=x0[6]*a*a*65.0/432.0;
    // I[1]=(x0[6]/144.0)*(a*a*29.0/3.0+9.0*w*w);
    // I[2]=(x0[6]/12.0)*(a*a+0.75*w*w);

    std::vector<double> CoS(3,0.0);
    CoS[1]=0.25*(2.0*ratio_m-1.0)*a; // to uncomment for displacement
    //CoS[2]=a/6.0;

    std::vector<double> I2=I;
    //I2[1]=I[2];I2[2]=I[1]; 

    std::vector<double> CoS2=CoS;

    other[0]=0.3;
    other[1]=0.0;
    //other.push_back(1.0e-8);

    double t0=(3600.0*24.0*jd0)/T_ref;

    std::vector<int> flags(4,0);
    flags[0]=1; // lunisolar
    flags[1]=1; // drag
    flags[2]=0; // SRP
    flags[3]=1; // magnetic
    std::vector<int> flags2=flags;
    flags2[3]=0;

    std::vector<double> y0(14),y(14),yf(14);
    for(int i=0;i<7;i++)
        y0[i]=x0[i];
    double beta=50.0*smartastro::constants::pi/180.0;
    y0[7]=0.0;
    y0[8]=0.5*sqrt(3.0)*sin(0.5*beta);
    y0[9]=0.5*sin(0.5*beta);
    y0[10]=cos(0.5*beta);
    y0[11]=1.0e-2*T_ref;
    y0[12]=-2.0e-2*T_ref;
    y0[13]=3.0e-2*T_ref; 

    int n_gravity = 29;
    int order_integration = 8;
    smartastro::astrodynamics::dearth_6dof<double> dyn6dof(I,R_ref,T_ref,n_gravity,flags,p,0,CoS,other);
    smartmath::integrator::ABM<double> prop6dof(&dyn6dof, order_integration, true);

    smartastro::astrodynamics::dearth_6dof<double> dyn6dof2(I2,R_ref,T_ref,n_gravity,flags2,p2,0,CoS2,other);
    smartmath::integrator::ABM<double> prop6dof2(&dyn6dof2, order_integration, true);

    std::vector<double> dy(14),magn(3);
    std::vector<std::vector<double> > f;
    double C, V;

    int n_mini_step=floor(dt*T_ref/0.5)+1;
    h = dt / double(n_mini_step);
    double delta_jd=0.1*365.25;

    orbital6dof.open ("orbital_charged_6dof_LEO.txt");
    positions6dof.open ("positions_charged_6dof_LEO.txt");
    attitude.open ("attitude_charged_LEO.txt");
    miscellaneous.open("miscellaneous_6dof_LEO.txt");

    orbital6dof2.open ("orbital_charged_6dof2_LEO.txt");
    positions6dof2.open ("positions_charged_6dof2_LEO.txt");
    attitude2.open ("attitude2_charged_LEO.txt");
    miscellaneous2.open("miscellaneous2_6dof_LEO.txt");

    orbital6dof << setprecision(16) << 0.0 << " " << kep[0] << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
    positions6dof << y0[0]*R_ref << " " << y0[1]*R_ref << " " << y0[2]*R_ref << " " << y0[3]*R_ref/T_ref << " " << y0[4]*R_ref/T_ref << " " << y0[5]*R_ref/T_ref << "\n" ;
    attitude << setprecision(16) << y0[7] << " " << y0[8] << " " << y0[9] << " " << y0[10] << " " << y0[11]/T_ref << " " << y0[12]/T_ref << " " << y0[13]/T_ref << "\n" ;
    dyn6dof.charging(jd,y0,C,V);
    miscellaneous << V << " " << C*V << " " << dyn6dof.shadow(jd,y) << "\n" ;

    orbital6dof2 << setprecision(16) << 0.0 << " " << kep[0] << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
    positions6dof2 << y[0]*R_ref << " " << y0[1]*R_ref << " " << y0[2]*R_ref << " " << y0[3]*R_ref/T_ref << " " << y0[4]*R_ref/T_ref << " " << y0[5]*R_ref/T_ref << "\n" ;
    attitude2 << setprecision(16) << y0[7] << " " << y0[8] << " " << y0[9] << " " << y0[10] << " " << y0[11]/T_ref << " " << y0[12]/T_ref << " " << y0[13]/T_ref << "\n" ;
    dyn6dof2.charging(jd,y0,C,V);
    miscellaneous2 << V << " " << C*V << " " << dyn6dof2.shadow(jd,y) << "\n" ;

    stop=0;
    y=y0;
    t=t0;
    jd=jd0;
    prop6dof.initialize(order_integration, t, h, y0, f);
    while(stop==0)
    {
        tt=t+dt;
        prop6dof.integrate(t,tt,n_mini_step,y,yf,f);
        y=yf;
        t+=dt;   
        if(tt<t)
           stop=1;
        jd+=dt*T_ref/(3600.0*24.0);  
        for(int i=0; i<3; i++)
            car[i]=yf[i]*R_ref;
        for(int i=3; i<6; i++)
            car[i]=yf[i]*R_ref/T_ref;
        smartastro::astrocore::conversion_coordinates::car2kep(car,smartastro::constants::mu_earth,kep);     
        dyn6dof.charging(jd,yf,C,V);

        orbital6dof << setprecision(16) << (tt - t0) * T_ref << " " << kep[0] << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
        positions6dof << y[0]*R_ref << " " << y[1]*R_ref << " " << y[2]*R_ref << " " << y[3]*R_ref/T_ref << " " << y[4]*R_ref/T_ref << " " << y[5]*R_ref/T_ref << "\n" ;
        attitude << setprecision(16) << y[7] << " " << y[8] << " " << y[9] << " " << y[10] << " " << y[11]/T_ref << " " << y[12]/T_ref << " " << y[13]/T_ref << "\n" ;
        miscellaneous << V << " " << C*V << " " << dyn6dof.shadow(jd,yf) << "\n" ;

        if(jd-jd0>delta_jd)
            stop=1;
        //std::cout << jd-jd0 << std::endl;
    }

    orbital6dof.close();
    positions6dof.close();
    attitude.close();   
    miscellaneous.close();

    std::cout << "simulation 1 over 2 done" << std::endl;

    stop=0;
    y=y0;
    t=t0;
    jd=jd0;
    prop6dof2.initialize(order_integration, t, h, y0, f);
    while(stop==0)
    {
        tt=t+dt;
        prop6dof2.integrate(t,tt,n_mini_step,y,yf,f);
        y=yf;
        t+=dt;
        if(tt<t)
           stop=1;
        jd+=dt*T_ref/(3600.0*24.0);
        for(int i=0; i<3; i++)
           car[i]=yf[i]*R_ref;
        for(int i=3; i<6; i++)
            car[i]=yf[i]*R_ref/T_ref;
        smartastro::astrocore::conversion_coordinates::car2kep(car,smartastro::constants::mu_earth,kep);
        dyn6dof2.charging(jd,yf,C,V);

        orbital6dof2 << setprecision(16) << (tt - t0) * T_ref << " " << kep[0] << " " << kep[1] << " " << kep[2]*180.0/pi << " " << kep[3]*180.0/pi << " "  << kep[4]*180.0/pi << " " << kep[5]*180.0/pi << "\n" ;
        positions6dof2 << y[0]*R_ref << " " << y[1]*R_ref << " " << y[2]*R_ref << " " << y[3]*R_ref/T_ref << " " << y[4]*R_ref/T_ref << " " << y[5]*R_ref/T_ref << "\n" ;
        attitude2 << setprecision(16) << y[7] << " " << y[8] << " " << y[9] << " " << y[10] << " " << y[11]/T_ref << " " << y[12]/T_ref << " " << y[13]/T_ref << "\n" ;
        miscellaneous2 << V << " " << C*V << " " << dyn6dof2.shadow(jd,yf) << "\n" ;

        if(jd-jd0>delta_jd)
            stop=1;
    }
  
    orbital6dof2.close();
    positions6dof2.close(); 
    attitude2.close(); 
    miscellaneous2.close();

    return 0;
}
