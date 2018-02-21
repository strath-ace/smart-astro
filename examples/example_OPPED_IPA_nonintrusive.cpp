/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This examples shows how to use OPPED-IPA with non-intrusive Chebyshev interpolation */

#include <fstream>
#include <typeinfo>
#include <string>


#include "../include/smartastro.h"
#include "smartmath.h"

const double pi = smartastro::constants::pi;

int main() 
{

    bool monomial_base = true; //true for output coefficients in monomial base, false for the chebyshev one
    bool print_results_to_file = true;
    bool print_time_to_screen = true;

    //algebra
    int nvar = 7;
    int nparam  = 6;
    int poly_degree = 2;

    //number of points in the sample
    int nsamples = smartmath::combination(nvar+nparam,poly_degree); // minimum number of samples required 
    nsamples += 100; // optional additional smamples

    /** scaling parameters */
    double R_ref = smartastro::constants::R_earth; // distance in meters
    double T_ref = smartastro::constants::T_Earth; // time in seconds

    double tof = 24.0 * 3600.0 / T_ref; // scaled time of flight
    int n_steps = floor(tof * T_ref / 60.0);

    int n_max = 9; // degree of EGM96 harmonics used 
    std::vector<int> flags(4, 0); // flags for other orbital perturbations 
    flags[0] = 1; // drag
    flags[1] = 1; // third-body
    flags[2] = 1; // SRP

    /** initial gregorian date */
    std::vector<double> date0(6);
    date0[0] = 2013.0; // year
    date0[1] = 10.0; // month
    date0[2] = 22.0; // day
    date0[3] = 3.0; // hour
    date0[4] = 0.0; // min
    date0[5] = 0.0; // s

    /** initial Keplerian coordinates */
    std::vector<double> kep0(6);
    kep0[0]=8000.0e3 / R_ref; // scaled semi-major axis
    kep0[1]=1.0e-3; // eccentricity
    kep0[2]=89.0*pi/180.0; // inclination
    kep0[3]=60.0*pi/180.0; // right ascension of the ascending node
    kep0[4]=30.0*pi/180.0; // argument of periapsis
    kep0[5]=100.0*pi/180.0; // true anomaly

    double m = 1100.0; // initial mass in kilograms

    /** drag-related parameters */
    std::vector<double> drag_param(4);
    drag_param[0] = 5.0; // drag cross-section in square meters 
    drag_param[1] = 2.0; // drag coefficient
    drag_param[2] = 106.4; // mean solar flux in solar flux units
    drag_param[3] = 3.85; // logarithmic geomagnetic index

    /** SRP parameters */
    std::vector<double> SRP_param(2);
    SRP_param[0] = 10.0; // SRP cross-section in square meters 
    SRP_param[1] = 1.3; // SRP coefficient   

    /** conversions **/
    double jd0; // initial Julian date
    smartastro::astrocore::conversion_time::date2jd(date0, jd0);
    double t0 = jd0 * 3600.0 * 24.0 / T_ref;
    double t = t0+tof;
    std::vector<double> car(6);
    smartastro::astrocore::conversion_coordinates::kep2car(kep0,1.0,car);   

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
    std::vector<double> unc_drag_param(4);
    unc_drag_param[0] = 0.5;
    unc_drag_param[1] = 0.2;    
    unc_drag_param[2] = 2.0;
    unc_drag_param[3] = 0.33;

    //initialisation: uncertainty in SRP parameters
    std::vector<double> unc_SRP_param(2);
    unc_SRP_param[0] = 0.5;
    unc_SRP_param[1] = 0.1;    

    //computing ranges for sampling
    std::vector<double> ranges_lb(nvar+nparam), ranges_ub(nvar+nparam);
    for(int i=0; i<6; i++)
    {
        ranges_lb[i] = car[i]-unc_x[i];
        ranges_ub[i] = car[i]+unc_x[i];
    }
    ranges_lb[6] = m-unc_x[6];
    ranges_ub[6] = m+unc_x[6];
    for(int i=0; i<4; i++)
    {
        ranges_lb[nvar+i] = drag_param[i]-unc_drag_param[i];
        ranges_ub[nvar+i] = drag_param[i]+unc_drag_param[i];
    }
    for(int i=0; i<2; i++)
    {
        ranges_lb[nvar+4+i] = SRP_param[i]-unc_SRP_param[i];
        ranges_ub[nvar+4+i] = SRP_param[i]+unc_SRP_param[i];
    }

    //timing
    clock_t begin, end;
    begin=clock();
    
    //construct LHS sampling
    smartuq::sampling::lhs lhs_gen(nvar+nparam,nsamples,ranges_lb, ranges_ub);
    std::vector<std::vector<double> > LHS, H;
    for(int i=0;i<nsamples;i++)
    {
        std::vector<double> sample=lhs_gen();
        LHS.push_back(sample);
    }

    //initialise polynomial for interpolation. the monomial flag will define the base
    smartuq::polynomial::chebyshev_polynomial poly(nvar+nparam,poly_degree, ranges_lb, ranges_ub, monomial_base);

    //propagation (MAIN LOOP)
    std::vector<std::vector<double> > coeffs_all;
    std::vector<std::vector<double> > z;
    std::vector<std::vector<double> > res_coeffs;

    // perform nsamples forward integrations
    for(int i=0;i<nsamples;i++)
    {
        std::vector<double> LHS_SRP, LHS_drag, LHS_x;
        // separate states from parameters in the LHS sampling
        for(int j=0;j<nvar+nparam; j++)
        {
            if(j<nvar)
                LHS_x.push_back(LHS[i][j]);
            else if(j<nvar+4)
                LHS_drag.push_back(LHS[i][j]);
            else
                LHS_SRP.push_back(LHS[i][j]);
        }

        //initialise dynamics and integrator with corresponding set of parameters
        smartastro::astrodynamics::dearth_3dof <double> dyn(R_ref, T_ref, n_max, flags, LHS_drag, LHS_SRP);
        smartmath::integrator::ABM<double> integrator(&dyn,8,true); // false would mean that Adam-Bashforth-Moulton is initialized with RK4 rather than with a method of equivalent order

        //perform integration and save results
        std::vector<double> z_tmp;
        integrator.integrate(t0,t,n_steps,LHS_x,z_tmp);
        z.push_back(z_tmp);

    }

    // perform interpolation. For efficiency reason the function that interpolate multiple outputs is used
    // poly will evaluate according to its base
    if(H.size()==0)
        poly.interpolation(LHS,z,H,res_coeffs);
    else
        poly.solve(H,z,res_coeffs);

    for(int i=0;i<nvar;i++)
        coeffs_all.push_back(res_coeffs[i]);

    //timing
    end=clock();
    double time = (double (end-begin))/CLOCKS_PER_SEC;
    if(print_time_to_screen) cout << "non-intrusive degree " << poly_degree << " time elapsed : " << time << endl << endl;

    //printing
    if(print_results_to_file)
    {
        std::ofstream file;
        file.open ("OPPED_IPA_nonintrusive_cheby.txt");
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