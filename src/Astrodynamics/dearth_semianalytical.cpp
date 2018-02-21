/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/dearth_semianalytical.h"
#include "Propagators/perturbation_propagator.h"
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

using namespace smartastro;
using namespace smartastro::astrodynamics;

dearth_semianalytical::dearth_semianalytical(const std::vector<bool> &flags, const int &zonal, const std::vector<double> &params_drag, const std::vector<double> &other, const int &subdivision) : 
    base_astrodynamics<double>("Semi-analytical Earth dynamics"), m_flags(flags), m_zonal(zonal), m_params(params_drag), m_subdivision(subdivision){

    /* Sanity checks */
    if(m_params.size() < 1)
        smartastro_throw("DEARTH_SEMIANALYTICAL: vector of parameters needs to have at least one component"); 
    if(m_params[0] <= 0.0)
        smartastro_throw("DEARTH_SEMIANALYTICAL: mass cannot be negative");  
    if(double(m_subdivision / 3) != double(m_subdivision) / 3.0)
        smartastro_throw("DEARTH_SEMIANALYTICAL: number of subdivisions for quadrature must be a multiple of 3");     
    if(m_flags.size() != 4)
        smartastro_throw("DEARTH_SEMIANALYTICAL: there must be four flags for perturbations");           
	if(m_flags[0])
	{
		if(params_drag.size() < 2)
			smartastro_throw("DEARTH_SEMIANALYTICAL: cross-section must be provided for atmospheric drag"); 
	}
	if(m_flags[2])
	{
		if(other.size() < 2)
			smartastro_throw("DEARTH_SEMIANALYTICAL: SRP area and coefficient must be provided"); 
	}

    std::vector<bool> flags_analytical(1, flags[1]); // flag for lunisolar
    m_prop = new smartastro::propagator::perturbation_propagator(flags_analytical, m_zonal); // temporary hack 

    std::vector<int> flags_cartesian(4, 0);
    std::vector<double> drag_param(1, 0.0);
    drag_param[0] = params_drag[1]; // cross-section
    if(m_flags[0]) // flag for drag
	{
		flags_cartesian[0] = 1;
		if(m_params.size() >= 3)
			drag_param.push_back(params_drag[2]); // drag coefficient
		if(m_params.size() >= 4)
			drag_param.push_back(params_drag[3]); // mean solar flux
		if(m_params.size() >= 5)
			drag_param.push_back(params_drag[4]); // geomagnetic index			
	}    	
	if(m_flags[2]) // flag for SRP
		flags_cartesian[2] = 2; // cylindrical model
	if(m_flags[3]) // flag for magnetic force
		flags_cartesian[3] = 1;			
	m_dyn = new smartastro::astrodynamics::dearth_3dof<double>(constants::R_earth, constants::T_Earth, m_zonal, flags_cartesian, drag_param, other);    

};

dearth_semianalytical::~dearth_semianalytical(){

}


int dearth_semianalytical::evaluate(const double &t, const std::vector<double> &x, std::vector<double> &dx) const{

	double pi = constants::pi;
	double e2 = x[1] * x[1] + x[2] * x[2]; // squared eccentricity
	double inter = sqrt(x[3] * x[3] + x[4] * x[4]); // tan(i/2)
	double si = 2.0 * inter / (1.0 + inter * inter); // sin(i)
	double n = pow(x[0] / (1.0 - e2), -1.5); // mean motion	
	if(m_zonal >= 2) // correction to n if J2 taken into account
		n *= 1.0 + 1.5 * m_J2 * (1.0 / x[0]) * (1.0 / x[0]) * sqrt(1.0 - e2) * (1.0 - 1.5 * si * si);
	double T = 2.0 * pi / n; // period

	std::vector<double> xp = x;
	m_prop->set_timing(x[5], t * constants::T_Earth / (3600.0 * 24.0));
	m_prop->propagate(x[5], x[5] + 2.0 * pi, x, xp); // computing state propagated over one revolution
	for(int j = 0; j < 6; j++)
		dx[j] = (xp[j] - x[j]) / T;

	std::vector<double> JDs;
	std::vector<std::vector<double> > KEs;
	if((m_flags[0])||(m_flags[2])||(m_flags[3])) // pre-computations for quadrature 
	{
		double jd = t * constants::T_Earth / (3600.0 * 24.0);
		double dt = T / double(m_subdivision);
		double M0, theta_0, theta_f, ecc;
		std::vector<double> oscul = x;
		if(m_zonal >= 2)
		{
			std::vector<double> kep_mean(6), kep_oscul(6);
			astrocore::conversion_coordinates::modeq2kep(x, kep_mean);
		    astrocore::conversion_coordinates::mean2oscul(kep_mean, 1.0, 1.0, kep_oscul);    
		    astrocore::conversion_coordinates::kep2modeq(kep_oscul, oscul);
		}
		else
			oscul = x;
		ecc = sqrt(oscul[1] * oscul[1] + oscul[2] * oscul[2]);
		theta_0 = oscul[5] - atan2(oscul[2], oscul[1]);
		theta_0 -= 2.0 * pi * floor(theta_0 / (2.0 * pi));
		astrocore::conversion_time::true2mean_anomaly(theta_0, ecc, M0);
		
		JDs.push_back(jd);
		KEs.push_back(oscul);
		for(int k = 1; k <= m_subdivision; k++) // computing quadrature points (date and osculating state)
		{
			JDs.push_back(jd + double(k) * dt * constants::T_Earth / (3600.0 * 24.0));
			KEs.push_back(oscul);
			astrocore::conversion_time::mean2true_anomaly(M0 + double(k) * n * dt, ecc, theta_f);
			if(theta_f > theta_0)
				KEs[k][5] += theta_f - theta_0;
			else
				KEs[k][5] += theta_f - theta_0 + 2.0 * pi;
		}
	}

	// std::vector<double> integral(5, 0.0);
	// std::vector<double> f0 = x, f1 = f0, f2 = f0, f3 = f0;
	// gravity_oscul(t * constants::T_Earth / (3600.0 * 24.0), KEs[0], f0);
	// for(int k = 0; k < m_subdivision / 3; k++)
	// {
	// 	gravity_oscul(JDs[3*k+1], KEs[3*k+1], f1);
	// 	gravity_oscul(JDs[3*k+2], KEs[3*k+2], f2);	
	// 	gravity_oscul(JDs[3*k+3], KEs[3*k+3], f3);	
	// 	for(int i = 0; i < 5; i++)
	// 		integral[i] += 3.0 * (f0[i] + 3.0 * f1[i] + 3.0 * f2[i] + f3[i]) / (8.0 * double(m_subdivision)); // composite Simpson's 3/8 rule
	// 	f0 = f3;
	// }
	// for(int i = 0; i < 5; i++) // no modification of the true longitude
	// 	dx[i] = integral[i];
	// dx[5] = 2.0 * pi / T;

	/* atmospheric drag */
	if(m_flags[0]) 
	{
		std::vector<double> integral(5, 0.0);
		std::vector<double> f0 = x, f1 = x, f2 = x, f3 = x;
		drag_oscul(JDs[0], KEs[0], f0);
		for(int k = 0; k < m_subdivision / 3; k++) 
		{
			drag_oscul(JDs[3*k+1], KEs[3*k+1], f1);
			drag_oscul(JDs[3*k+2], KEs[3*k+2], f2);	
			drag_oscul(JDs[3*k+3], KEs[3*k+3], f3);	
			for(int i = 0; i < 5; i++)
				integral[i] += 3.0 * (f0[i] + 3.0 * f1[i] + 3.0 * f2[i] + f3[i]) / (8.0 * double(m_subdivision)); // composite Simpson's 3/8 rule
	 		f0 = f3;
		}
		for(int i = 0; i < 5; i++) // no modification of the true longitude
			dx[i] += integral[i];
	}

	/* solar radiation pressure */
	if(m_flags[2]) 
	{
		std::vector<double> integral(5, 0.0);
		std::vector<double> f0 = x, f1 = x, f2 = x, f3 = x;
		SRP_oscul(JDs[0], KEs[0], f0);
		for(int k = 0; k < m_subdivision / 3; k++) 
		{
			SRP_oscul(JDs[3*k+1], KEs[3*k+1], f1);
			SRP_oscul(JDs[3*k+2], KEs[3*k+2], f2);	
			SRP_oscul(JDs[3*k+3], KEs[3*k+3], f3);	
			for(int i = 0; i < 5; i++)
				integral[i] += 3.0 * (f0[i] + 3.0 * f1[i] + 3.0 * f2[i] + f3[i]) / (8.0 * double(m_subdivision)); // composite Simpson's 3/8 rule
	 		f0 = f3;
		}
		for(int i = 0; i < 5; i++) // no modification of the true longitude
			dx[i] += integral[i];
	}	

	/* magnetic force */
	if(m_flags[3]) 
	{
		std::vector<double> integral(5, 0.0);
		std::vector<double> f0 = x, f1 = x, f2 = x, f3 = x;
		magnetic_oscul(JDs[0], KEs[0], f0);
		for(int k = 0; k < m_subdivision / 3; k++) 
		{
			magnetic_oscul(JDs[3*k+1], KEs[3*k+1], f1);
			magnetic_oscul(JDs[3*k+2], KEs[3*k+2], f2);	
			magnetic_oscul(JDs[3*k+3], KEs[3*k+3], f3);	
			for(int i = 0; i < 5; i++)
				integral[i] += 3.0 * (f0[i] + 3.0 * f1[i] + 3.0 * f2[i] + f3[i]) / (8.0 * double(m_subdivision)); // composite Simpson's 3/8 rule
	 		f0 = f3;
		}
		for(int i = 0; i < 5; i++) // no modification of the true longitude
			dx[i] += integral[i];
	}

  	return 0;
}

int dearth_semianalytical::drag_oscul(const double &jd, const std::vector<double> &x, std::vector<double> &dx) const{

	std::vector<double> car(6);
	// if(m_zonal >= 2)
	// {
	// 	std::vector<double> mean(6);
	// 	astrocore::conversion_coordinates::modeq2kep(x, mean);
	//     astrocore::conversion_coordinates::mean2oscul(mean, 1.0, 1.0, kep);    
	// }
	// else
	astrocore::conversion_coordinates::modeq2car(x, 1.0, car);   
    
    std::vector<double> posvelmass(car), acc(3), RTH(3);
    posvelmass.push_back(m_params[0]);
    m_dyn->aero(jd, posvelmass, acc);
    astrocore::conversion_coordinates::car2rth(acc, car, RTH);
	dearth_orb<double>::gauss_modeq(x, 1.0, RTH, dx);

	return 0;
}

int dearth_semianalytical::SRP_oscul(const double &jd, const std::vector<double> &x, std::vector<double> &dx) const{

	std::vector<double> car(6);
	astrocore::conversion_coordinates::modeq2car(x, 1.0, car);  
    
    std::vector<double> posvelmass(car), acc(3), RTH(3);
    posvelmass.push_back(m_params[0]);
    m_dyn->SRP(jd, posvelmass, acc);
    astrocore::conversion_coordinates::car2rth(acc, car, RTH);
	dearth_orb<double>::gauss_modeq(x, 1.0, RTH, dx);

	return 0;
}

int dearth_semianalytical::magnetic_oscul(const double &jd, const std::vector<double> &x, std::vector<double> &dx) const{

	// old implementation (analytical model from Hamilton)
	// double e2 = x[1] * x[1] + x[2] * x[2];
	// double inter = sqrt(x[3] * x[3] + x[4] * x[4]); // tan(i/2)
	// double si = 2.0 * inter / (1.0 + inter * inter);
	// double ci = (1.0 - inter * inter) / (1.0 + inter * inter);		
	// double n = pow(x[0] / (1.0 - e2), -1.5); 
	// if(m_zonal >= 2)  
	// 	n *= 1.0 + 1.5 * m_J2 * (1.0 / x[0]) * (1.0 / x[0]) * sqrt(1.0 - e2) * (1.0 - 1.5 * si * si);	
	// double rate = constants::T_Earth * constants::omega_earth;
	// double L = (m_q * constants::g10_2015) * constants::T_Earth * rate / m_params[0];
	// double a, e, Om, om;
	// double da, de, di, dOm, dom;

	// a = x[0] / (1.0 - e2);
	// e = sqrt(e2);
	// Om = atan2(x[4], x[3]);
	// om = atan2(x[2], x[1]) - Om;
	// if(Om < 0.0)
	// 	Om += 2.0 * constants::pi;
	// if(om < 0.0)
	// 	om += 2.0 * constants::pi;		

	// da = 0.0;
	// de = -0.25 * n * L * e * sqrt(1.0 - e2) * si * si * sin(2.0 * om);
	// di = 0.25 * n * L * e2 * si * ci * sin(2.0 * om) / sqrt(1.0 - e2); 
	// dOm = n * L * (ci - (n / rate) / (1.0 - e2)) / sqrt(1.0 - e2);
	// dom = n * L * ci * (-ci + 3.0 * (n / rate) / (1.0 - e2)) / sqrt(1.0 - e2);

	// dx[0] = da * (1.0 - e2) - 2.0 * a * e * de; 
	// dx[1] = de * cos(Om + om) - e * sin(Om + om) * (dOm + dom);
	// dx[2] = de * sin(Om + om) + e * cos(Om + om) * (dOm + dom);
	// dx[3] = -inter * sin(Om) * dOm + 0.5 * (1.0 + inter * inter) * di * cos(Om);
	// dx[4] = inter * cos(Om) * dOm + 0.5 * (1.0 + inter * inter) * di * sin(Om);

	std::vector<double> car(6);
	astrocore::conversion_coordinates::modeq2car(x, 1.0, car);
    
    std::vector<double> posvelmass(car), acc(3), RTH(3);
    posvelmass.push_back(m_params[0]);
    m_dyn->magnetism(jd, posvelmass, acc);
    astrocore::conversion_coordinates::car2rth(acc, car, RTH);
	dearth_orb<double>::gauss_modeq(x, 1.0, RTH, dx);

	return 0;
}

// int dearth_semianalytical::gravity_oscul(const double &jd, const std::vector<double> &x, std::vector<double> &dx) const{

// 	std::vector<double> car(6);
// 	astrocore::conversion_coordinates::modeq2car(x, 1.0, car);
    
//     std::vector<double> posvelmass(car), acc1(3), acc2(3, 0.0), acc(3), RTH(3);
//     posvelmass.push_back(m_params[0]);

// 	m_dyn->Earth_gravity(jd, posvelmass, acc1);
// 	if(m_flags[1])
//     	m_dyn->lunisolar(jd, posvelmass, acc2);
//  	for(int i = 0; i < 3; i++)
// 		acc[i] = acc1[i] + acc2[i];   

//     astrocore::conversion_coordinates::car2rth(acc, car, RTH);
// 	dearth_orb<double>::gauss_modeq(x, 1.0, RTH, dx);

// 	return 0;
// }