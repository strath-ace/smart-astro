/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/dearth_hamiltonian.h"
#include "Astrodynamics/base_dearth.h"
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

using namespace smartastro;
using namespace smartastro::astrodynamics;


template<class T>
dearth_hamiltonian<T>::dearth_hamiltonian(const int &n_max, const std::vector<bool> &flags): smartmath::dynamics::hamiltonian_momentum<T>("Hamiltonian dynamics for Earth orbits", 3, true), m_max_degree_Earth_gravity(n_max), m_flags(flags){

    if(m_flags.size() < 1)
        smartastro_throw("DEARTH_HAMILTONIAN: there must be at least one flag for perturbations");

    std::vector<double> C((m_max_degree_Earth_gravity + 1) * (m_max_degree_Earth_gravity + 1), 0.0), S((m_max_degree_Earth_gravity + 1) * (m_max_degree_Earth_gravity + 1), 0.0);
    m_C = C; m_S = S;
    /** retrieving zonal coefficients for zonal spherical harmonics **/
    for(int n = 0; n < m_max_degree_Earth_gravity + 1; n++)  
        m_C[n] = constants::C_Earth_norm[n] * sqrt((2.0 * double(n) + 1.0));
    
    std::vector<int> flags_3dof(4, 0);
    std::vector<T> params(2);
    m_dyn = new smartastro::astrodynamics::dearth_3dof<T>(constants::R_earth, constants::T_Earth, 0, flags_3dof, params, params); 

};


template < class T >
dearth_hamiltonian<T>::~dearth_hamiltonian()
{
      
}   

template < class T >
int dearth_hamiltonian<T>::DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{
    
    /* sanity checks */
    if(q.size() != m_dim)
        smartastro_throw("DHQ: the position must have dimension 3");

    std::vector<T> grad(3, 0.0 * q[0]);
    base_dearth<T>::harmonics(q, m_max_degree_Earth_gravity, 1.0, m_C, m_S, grad);
    dH.clear();
    for(unsigned int i = 0; i < m_dim; i++)
        dH.push_back(-grad[i]);

    if(m_flags[0])
    {
        std::vector<T> dH2(3, 0.0 * q[0]);
        m_dyn->lunisolar(t * constants::T_Earth / (3600.0 * 24.0), q, dH2);
        for(unsigned int i = 0; i < m_dim; i++)
            dH[i] -= dH2[i];
    }

    return 0;
}    

template class dearth_hamiltonian<double>;
#ifdef ENABLE_SMARTUQ
template class dearth_hamiltonian<smartuq::polynomial::chebyshev_polynomial>;
template class dearth_hamiltonian<smartuq::polynomial::taylor_polynomial>;
#endif