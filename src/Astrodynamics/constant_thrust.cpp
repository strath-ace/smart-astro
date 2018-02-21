/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/constant_thrust.h"
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

using namespace smartastro;
using namespace smartastro::astrodynamics;


template<class T>
constant_thrust<T>::constant_thrust(const double &L_scale, const double &t_scale, const double &mu, const std::vector<double> &thrust, const double &Isp, bool inertial):base_astrodynamics<T>("Keplerian motion with constant thrust in local/inertial orbital frame", L_scale, t_scale), m_mu(mu), m_thrust(thrust),m_Isp(Isp), m_inertial(inertial){

    /*Sanity check*/
    if(thrust.size()!=3)
        smartastro_throw("CONSTANT_THRUST: thrust vector must be 3-dimensional");
    if(mu<=0.0)
        smartastro_throw("CONSTANT_THRUST: gravitational constant cannot be negative");    
    if(Isp<=0.0)
        smartastro_throw("CONSTANT_THRUST: specific impulse cannot be negative"); 

    /* Computing constant rate of mass consumption */
    m_rate = -m_T_scale*sqrt(m_thrust[0]*m_thrust[0]+m_thrust[1]*m_thrust[1]+m_thrust[2]*m_thrust[2])/(constants::free_fall*m_Isp);

};


template < class T >
constant_thrust<T>::~constant_thrust()
{
      
}

template < class T >
void constant_thrust<T>::set_thrust(const std::vector<double> &thrust){
    /*Sanity check*/
    if(thrust.size()!=3)
        smartastro_throw("CONSTANT_THRUST: thrust vector must be 3-dimensional");
	
	m_thrust = thrust;
	/* Computing constant rate of mass consumption */
	m_rate = -m_T_scale*sqrt(m_thrust[0]*m_thrust[0]+m_thrust[1]*m_thrust[1]+m_thrust[2]*m_thrust[2])/(constants::free_fall*m_Isp);

}

template < class T >
int constant_thrust<T>::evaluate(const double &t, const std::vector<T> &x, std::vector<T> &dx) const{
      
    /*Sanity check*/
    if(x.size()!=7)
        smartastro_throw("CONSTANT_THRUST: state vector must be 7-dimensional");

	dx.clear();
	dx = x;

	/* Keplerian dynamics */
	dx[0]=x[3];
	dx[1]=x[4];
	dx[2]=x[5];
	T r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	T factor = m_mu/pow(r,3);
	dx[3]=-x[0]*factor;
	dx[4]=-x[1]*factor;
	dx[5]=-x[2]*factor;	
	dx[6]=m_rate;
	T factor2 = (1.0/x[6])*pow(m_T_scale,2)/m_L_scale;

	if (m_inertial){
		/* Adding thrust */
		dx[3] += (m_thrust[0])*factor2;
		dx[4] += (m_thrust[1])*factor2;
		dx[5] += (m_thrust[2])*factor2;	
	}
	else{
		/* Adding thrust */
		T inv_r = 1.0/r;
		T inv_v = 1.0/sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);
		/* computing radial direction */
		std::vector<T> dir_r; 
		for(int i=0;i<3;i++){
			dir_r.push_back(x[i]*inv_r);
		}
		/* computing velocity direction */
		std::vector<T> dir_v = dir_r, dir_h = dir_r, dir_t = dir_r;
		for(int i=0;i<3;i++){
			dir_v[i]=x[i+3]*inv_v;
		}
		/* computing out-of-plane direction */
		dir_h[0] = dir_r[1]*dir_v[2]-dir_r[2]*dir_v[1];;
		dir_h[1] = dir_r[2]*dir_v[0]-dir_r[0]*dir_v[2];
		dir_h[2] = dir_r[0]*dir_v[1]-dir_r[1]*dir_v[0];
		/* computing tranversal direction */
		dir_t[0] = dir_h[1]*dir_r[2]-dir_h[2]*dir_r[1];
		dir_t[1] = dir_h[2]*dir_r[0]-dir_h[0]*dir_r[2];
		dir_t[2] = dir_h[0]*dir_r[1]-dir_h[1]*dir_r[0];	

		T factor2 = (1.0/x[6])*pow(m_T_scale,2)/m_L_scale;
		dx[3] += (dir_r[0]*m_thrust[0] + dir_t[0]*m_thrust[1] + dir_h[0]*m_thrust[2])*factor2;
		dx[4] += (dir_r[1]*m_thrust[0] + dir_t[1]*m_thrust[1] + dir_h[1]*m_thrust[2])*factor2;
		dx[5] += (dir_r[2]*m_thrust[0] + dir_t[2]*m_thrust[1] + dir_h[2]*m_thrust[2])*factor2;		
	}
    return 0;
}


template class constant_thrust<double>;
template class constant_thrust<float>;
template class constant_thrust<long double>;
#ifdef ENABLE_SMARTUQ
template class constant_thrust<smartuq::polynomial::chebyshev_polynomial>;
template class constant_thrust<smartuq::polynomial::taylor_polynomial>;
#endif