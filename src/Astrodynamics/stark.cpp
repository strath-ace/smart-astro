/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/stark.h"
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

using namespace smartastro;
using namespace smartastro::astrodynamics;


template<class T>
stark<T>::stark(const double &L_scale, const double &t_scale, const double &mu, const std::vector<double> &acc):base_astrodynamics<T>("Accelerated Keplerian dynamics", L_scale, t_scale), m_mu(mu), m_acc(acc){

    /*Sanity check*/
    if(acc.size()!=3)
        smartastro_throw("STARK: thrust vector must be 3-dimensional");
};


template < class T >
stark<T>::~stark()
{
      
}


template < class T >
int stark<T>::evaluate(const double &t, const std::vector<T> &x, std::vector<T> &dx) const{
      
    /*Sanity check*/
    if(x.size()!=6)
        smartastro_throw("STARK: state vector (Cartesian coordinates and mass) must be 6-dimensional");

	dx.clear();
	dx = x;
	/* Computing time derivatives */
	/* Position */
	dx[0]=x[3];
	dx[1]=x[4];
	dx[2]=x[5];
	/* Velocity */
	T r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	T factor = m_mu/pow(r,3);
	dx[3]=-x[0]*factor+m_acc[0];
	dx[4]=-x[1]*factor+m_acc[1];
	dx[5]=-x[2]*factor+m_acc[2];	

    return 0;
}

template class stark<double>;
template class stark<float>;
template class stark<long double>;
#ifdef ENABLE_SMARTUQ
template class stark<smartuq::polynomial::chebyshev_polynomial>;
template class stark<smartuq::polynomial::taylor_polynomial>;
#endif