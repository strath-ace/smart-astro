/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example shows how to use OPPED-IPA with the Taylor Differential Algebra */

#include <fstream>
#include <typeinfo>
#include <string>

#include "../include/smartastro.h"
#include "smartmath.h"

const double pi = smartastro::constants::pi;

int main() 
{

    bool print_results_to_file = true;
    bool print_time_to_screen = true;

    //algebra
    int nvar = 7;
    int nparam  = 4;
    int poly_degree = 2;

    //polynomial allocation
    std::vector<smartuq::polynomial::taylor_polynomial> x0, param, xf, other(2);    

    /** scaling parameters */
    double R_ref = smartastro::constants::R_earth; // distance in meters
    double T_ref = smartastro::constants::T_Earth; // time in seconds
    double omega=smartastro::constants::omega_earth*T_ref;

    double tof = 24.0 * 3600.0 / T_ref; // scaled time of flight

    int n_max = 9; // degree of EGM96 harmonics used 
    std::vector<int> flags(4, 0); // flags for other orbital perturbations 
    flags[0] = 1; // drag
    flags[1] = 1; // lunisolar

    /** initial gregorian date */
    std::vector<double> date0(6);
    date0[0] = 2013.0; // year
    date0[1] = 10.0; // month
    date0[2] = 22.0; // day
    date0[3] = 3.0; // hour
    date0[4] = 0.0; // min
    date0[5] = 0.0; // s
    /* time conversions */
    double jd0;
    smartastro::astrocore::conversion_time::date2jd(date0, jd0);
    double t0 = jd0 * 3600.0 * 24.0 / T_ref;
    double t = t0+tof;

    /** initial Cartesian coordinates */
    std::vector<double> y(6); // initial conditions in Earth-fixed Earth-centered frame
    y[0]=2010038.0371203022/R_ref;
    y[1]=5868801.544770562/R_ref;
    y[2]=2277523.604184902/R_ref;
    y[3]=2185.0237988871622/(R_ref/T_ref);
    y[4]=2072.4192499635105/(R_ref/T_ref);
    y[5]=-7232.216463152104/(R_ref/T_ref);
    std::vector<double> M(9);
    smartastro::astrocore::conversion_frames::inertial_to_bf(jd0,M); // rotation matrix
    /** computing pos-vel in inertial frame */
    std::vector<double> x(7);
    x[0]=M[0]*y[0]+M[1]*y[1]+M[2]*y[2];
    x[1]=M[3]*y[0]+M[4]*y[1]+M[5]*y[2];
    x[2]=M[6]*y[0]+M[7]*y[1]+M[8]*y[2];
    x[3]=M[0]*(y[3]-omega*y[1])+M[1]*(y[4]+omega*y[0])+M[2]*y[5];
    x[4]=M[3]*(y[3]-omega*y[1])+M[4]*(y[4]+omega*y[0])+M[5]*y[5];
    x[5]=M[6]*(y[3]-omega*y[1])+M[7]*(y[4]+omega*y[0])+M[8]*y[5];

    x[6] = 1100.0; // initial mass in kilograms

    /** drag-related parameters */
    std::vector<double> p(4);
    p[0] = 5.0; // drag cross-section in square meters 
    p[1] = 2.0; // drag coefficient
    p[2] = 106.4; // mean solar flux in solar flux units
    p[3] = 3.85; // logarithmic geomagnetic index 

    //initialisation: uncertainty in initial states
    std::vector<double> unc_x(7);
    unc_x[0] = 1000.0 / R_ref;
    unc_x[1] = 1000.0 / R_ref;
    unc_x[2] = 1000.0 / R_ref;
    unc_x[3] = 0.1 / (R_ref/T_ref);
    unc_x[4] = 0.1 / (R_ref/T_ref);
    unc_x[5] = 0.1 / (R_ref/T_ref);
    unc_x[6] = 1.0;

    //initialisation: uncertainty in drag-related parameters
    std::vector<double> unc_p(4);
    unc_p[0] = 0.5;
    unc_p[1] = 0.2;    
    unc_p[2] = 2.0;
    unc_p[3] = 0.33;  

    //initialise state variables as Chebycheff/Taylor base of order 1 in the variable i mapped to [x-unc_x, x+unc_x]
    for(int i=0;i<nvar;i++)
        x0.push_back(smartuq::polynomial::taylor_polynomial(nvar+nparam, poly_degree, i, x[i]-unc_x[i], x[i]+unc_x[i]));

    //initialise parameter variables as Chebycheff/Taylor base of order 1 in the variable 7+i mapped to [p-unc_p, p+unc_p]
    for(int i=0;i<nparam;i++)
        param.push_back(smartuq::polynomial::taylor_polynomial(nvar+nparam, poly_degree, nvar+i, p[i]-unc_p[i], p[i]+unc_p[i]));

    //timing
    clock_t begin, end;
    begin=clock();

    //dynamical system
    smartastro::astrodynamics::dearth_3dof < smartuq::polynomial::taylor_polynomial > dyn(R_ref, T_ref, n_max, flags, param, other);
    smartmath::integrator::ABM< smartuq::polynomial::taylor_polynomial > integrator(&dyn,8,true); // false would mean that Adam-Bashforth-Moulton is initialized with RK4 rather than with a method of equivalent order
    int n_steps = floor(tof * T_ref / 60.0);

    //propagation 
    std::vector<std::vector<double> > coeffs_all;
    integrator.integrate(t0,t,n_steps,x0,xf);
    for(int j=0; j<6; j++)
    {
        std::vector<double> coeffs = xf[j].get_coeffs();
        coeffs_all.push_back(coeffs);
    }

    //deallocation of heavy structures
    x0[0].delete_M();

    //timing
    end=clock();
    double time_elapsed = (double (end-begin))/CLOCKS_PER_SEC;
    if (print_time_to_screen) cout << "Taylor " << poly_degree << " time elapsed : " << time_elapsed << endl;

    //printing
    if(print_results_to_file)
    {
        std::ofstream file;
        file.open ("OPPED_IPA_taylor.txt");
        for(unsigned int k=0; k<coeffs_all.size(); k++)
        {
            for(unsigned int kk=0; kk<coeffs_all[k].size(); kk++)
                file << setprecision(16) << coeffs_all[k][kk] << " ";
            file << "\n";
        }
        file << "\n\n\n\n";
        file.close();
    }
}