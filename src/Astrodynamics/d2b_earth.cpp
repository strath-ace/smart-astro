/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2016 University of Strathclyde and Authors ------
------------ Author: Annalisa Riccardi -------------------------------
------------ e-mail: annalisa.riccardi@strath.ac.uk ------------------
------------ Author: Carlos Ortega Absil -----------------------------
------------ e-mail: carlos.ortega@strath.ac.uk ----------------------
*/


#include "Astrodynamics/d2b_earth.h"


using namespace smartastro;
using namespace smartastro::astrodynamics;

template < class T >
d2b_earth<T>::d2b_earth(const std::vector<T> &param, const double &r_scale, const double &t_scale) :
    base_astrodynamics<T>("Two Body Problem around Earth",r_scale,t_scale),
    m_param(param), m_m_scale(1.0)
{
    if(m_param.size()!=10)
        smartastro_throw(m_name+": the parameters list need to be of size 10");

}

template < class T >
d2b_earth<T>::~d2b_earth()
{

}

template < class T >
int d2b_earth<T>::evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const
{
    //sanity checks
    if(t<0)
        smartastro_throw(m_name+": negative time supplied in evaluation of the dynamical system");
    if(state.size()!=7)
        smartastro_throw(m_name+": the state dimension needs to be 7");

    dstate.clear();

    //constant parameters
    double radius_earth = 6378.0*pow(10,3) / m_L_scale;
    double mu_earth = 398600.4415*pow(10,9) / (pow(m_L_scale,3)/pow(m_T_scale,2));
    double omega_earth = 7.2921150*pow(10,-5) * m_T_scale;
    double H0_atmosphere = 900000.0 / m_L_scale;

    //precomputations
    T r = sqrt(state[0]*state[0]+state[1]*state[1]+state[2]*state[2]);
    T tmp_3D =  mu_earth/pow(r,3);

    //atmospheric model
    T rho = m_param[4]*exp(-(r-radius_earth-H0_atmosphere)/m_param[5]);

    //relative velocity
    T rel_v_x = state[3]-omega_earth*state[1];
    T rel_v_y = state[4]+omega_earth*state[0];
    T mod_rel_v = sqrt(rel_v_x*rel_v_x+rel_v_y*rel_v_y+state[5]*state[5]);

    //drag computation
    T tmp_drag = 0.5*rho*m_param[6]*mod_rel_v/state[6];

    dstate.push_back(state[3]); //dx/dt
    dstate.push_back(state[4]); //dy/dt
    dstate.push_back(state[5]); //dz/dt
    dstate.push_back(-tmp_3D*state[0]+m_param[0]/state[6]+m_param[7]-tmp_drag*rel_v_x); //dv_x/dt
    dstate.push_back(-tmp_3D*state[1]+m_param[1]/state[6]+m_param[8]-tmp_drag*rel_v_y); //dv_y/dt
    dstate.push_back(-tmp_3D*state[2]+m_param[2]/state[6]+m_param[9]-tmp_drag*state[5]); //dv_z/dt
    dstate.push_back(-m_param[3]*sqrt(m_param[0]*m_param[0]+m_param[1]*m_param[1]+m_param[2]*m_param[2])); //dm/dt

    return 0;

}


template class d2b_earth<double>;
template class d2b_earth<float>;
template class d2b_earth<long double>;
#ifdef ENABLE_SMARTUQ
template class d2b_earth<smartuq::polynomial::chebyshev_polynomial>;
template class d2b_earth<smartuq::polynomial::taylor_polynomial>;
#endif
