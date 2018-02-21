/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/

#include "Astrodynamics/d2bp_3dof.h"

using namespace smartastro;
using namespace smartastro::astrodynamics;

// Contructor
template < class T >
d2bp_3dof<T>::d2bp_3dof(const T &mass,
                        const double &L_scale,
                        const double &T_scale):
    base_astrodynamics<T>("Two-body dynamics with 3DoF", L_scale, T_scale),
    m_mu(mass*constants::gravity_constant*pow(L_scale,3)/pow(T_scale,2))
{
    // m_body = NULL;
    this->m_state_length = 6;
}


// Destructor
template < class T >
d2bp_3dof<T>::~d2bp_3dof()
{
    // NOTHING YET
}


// Evaluate
template < class T >
int d2bp_3dof<T>::evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const
{
    // sanity checks
    if(state.size() != (this->m_state_length*this->m_sc_number))
        smartastro_throw(this->m_name+": state must be a vector of length "+std::to_string(this->m_state_length*this->m_sc_number));

    T zero = state[0]*0.0;

    dstate = std::vector<T>(state.size(), zero);

    std::vector<T> pos(3, zero);
    T r;
    T factor;

    for(unsigned int index_sc = 0; index_sc < this->m_sc_number; ++index_sc)
    {
        /* xdot = v */
        for(unsigned int index = 0; index < 3; ++index)
        {
            dstate[index_sc*this->m_state_length+index] = state[index_sc*this->m_state_length+3+index];
        }

        pos[0] = state[index_sc*this->m_state_length];
        pos[1] = state[index_sc*this->m_state_length+1];
        pos[2] = state[index_sc*this->m_state_length+2];

        r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        factor = -m_mu/pow(r,3); //-mu/r^3

        /* vdot = a */
        for(unsigned int index = 0; index < 3; ++index)
        {
            dstate[index_sc*this->m_state_length+3+index] = factor*pos[index];
        }
    }


    return 0;
}

template class d2bp_3dof<double>;
template class d2bp_3dof<float>;
template class d2bp_3dof<long double>;
#ifdef ENABLE_SMARTUQ
template class d2bp_3dof<smartuq::polynomial::chebyshev_polynomial>;
template class d2bp_3dof<smartuq::polynomial::taylor_polynomial>;
#endif
