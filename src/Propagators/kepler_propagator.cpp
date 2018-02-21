/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/

#include "../../include/Propagators/kepler_propagator.h"

using namespace smartastro;
using namespace propagator;

kepler_propagator::kepler_propagator(const double &mu, const double &tol, const int &maxiter, const double &R, const double &T):
    base_propagator("Kepler propagator")
{
    m_tol = tol;
    m_maxiter = maxiter;
    m_mu = mu;
    m_R = R; 
    m_T = T;

    if(mu<=0)
        smartastro_throw("KEPLER_PROPAGATOR: negative gravitational parameter");
    if(tol<=0)
        smartastro_throw("KEPLER_PROPAGATOR: negative tolerance");
    if(maxiter<0)
        smartastro_throw("KEPLER_PROPAGATOR: negative maximum number of iterations");
    if(R<=0)
        smartastro_throw("KEPLER_PROPAGATOR: negative scaling distance");
    if(T<=0)
        smartastro_throw("KEPLER_PROPAGATOR: negative scaling time");    

}

kepler_propagator::~kepler_propagator(){
}

int kepler_propagator::propagate(const double &t0, const double&tend, const std::vector<double> &initial_state, std::vector<double> &final_state) const{
    //Sanity checks
    if(initial_state.size()!=6)
        smartastro_throw("Kepler propagator: initial state needs to have 6 dimensions");
    if(final_state.size()!=6)
        smartastro_throw("Kepler propagator: final state needs to have 6 dimensions");

    // Declarations
    std::vector<double> r0(3), v0(3);
    double nr0 = 0, inv_mu = 0, alpha = 0, r0v0 = 0, chi0 = 0, h2 = 0, p = 0, s = 0, w = 0, a = 0, chin = 0, sqrt_mu = 0, numer = 0;
    double psi = 0, sqrt_psi = 0, c2 = 0, c3 = 0, r0v0_sqrt_mu = 0, r = 0, nr = 0, chin2 = 0, gd = 0, fd = 0, f = 0, g = 0;
    std::vector<double> h(3), r_vect(3), v_vect(3);
    int i;

    //Assignment
    double mu = m_mu*pow(m_T,2)/pow(m_R,3);
    double dt = tend-t0;

    r0[0] = initial_state[0]; r0[1] = initial_state[1]; r0[2] = initial_state[2];
    v0[0] = initial_state[3]; v0[1] = initial_state[4]; v0[2] = initial_state[5];

    nr0 = sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]);
    inv_mu = 1./mu;
    sqrt_mu = sqrt(mu);

    alpha = -(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]) * inv_mu + 2./nr0;
    r0v0 = r0[0]*v0[0] + r0[1]*v0[1] + r0[2]*v0[2];
    r0v0_sqrt_mu = r0v0 / sqrt_mu;

    double alpha_km = alpha / (m_R*1.0e-3);
    if(alpha_km > m_tol)
        chi0 = sqrt_mu*dt*alpha;
    else if (fabs(alpha_km) <= m_tol)
    {
	    h[0] = r0[1]*v0[2] - r0[2]*v0[1];
	    h[1] = r0[2]*v0[0] - r0[0]*v0[2];
	    h[2] = r0[0]*v0[1] - r0[1]*v0[0];

	    h2 = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
	    p = h2 * inv_mu;
	    s = 0.5 * atan(1. / (3.*sqrt(mu/pow(p,3))*dt));
	    w = atan(pow(tan(s), 1./3.));
	    chi0 = sqrt(p) * 2. / tan(2.*w);
    } 
    else // (alpha_km < -m_tol)
    {
        a = 1./alpha;
        if (dt >= 0.)
                chi0 =  sqrt(-a) * log(-2.0*mu*alpha*dt / (r0v0 + sqrt(-mu*a) * (1.0 - nr0*alpha)));
        else
            chi0 = -sqrt(-a) * log(-2.0*mu*alpha*dt / (r0v0 - sqrt(-mu*a) * (1.0 - nr0*alpha)));
    }

    i = 0;
    chin = chi0;
    numer = m_tol+1.;
    nr = 1.;
    while ((sqrt(m_R*1.0e-3)*fabs(numer)/nr > m_tol) && (i <= m_maxiter))
    {
        i += 1;
        psi = chin*chin*alpha;

        if (psi > m_tol)
        {
            sqrt_psi = sqrt(psi);
            c2 = (1.0 - cos(sqrt_psi)) / psi;
            c3 = (sqrt_psi-sin(sqrt_psi)) / (psi*sqrt_psi);
        } 
        else if (psi < -m_tol)
        {
            sqrt_psi = sqrt(-psi);
            c2 = (1.0 - cosh(sqrt_psi)) / psi;
            c3 = (sinh(sqrt_psi) - sqrt_psi) / (-psi*sqrt_psi);
        } 
        else
        {
            c2 = 0.5;
            c3 = 1./6.;
        }
        if(c2 > std::numeric_limits<double>::max() || c3 > std::numeric_limits<double>::max())
            return 3;
        
        r = chin*chin * c2 + r0v0_sqrt_mu * chin * (1.-psi*c3) + nr0*(1.-psi*c2);
        nr = fabs(r);
        numer = sqrt_mu*dt - chin*chin*chin*c3 - r0v0_sqrt_mu*chin*chin*c2 - nr0*chin*(1.-psi*c3);
        chin += numer/nr;
    }

    if (i > m_maxiter)
        return 2;

    chin2 = chin*chin;
    f = 1. - chin2*c2/nr0;
    g = dt - (chin2*chin)* c3 / sqrt_mu ;
    r_vect[0] = f*r0[0] + g*v0[0];
    r_vect[1] = f*r0[1] + g*v0[1];
    r_vect[2] = f*r0[2] + g*v0[2];
    nr = sqrt(r_vect[0]*r_vect[0]+r_vect[1]*r_vect[1]+r_vect[2]*r_vect[2]);
    gd = 1. - chin2*c2/nr;
    fd = sqrt_mu * chin * (psi*c3-1.) /(nr*nr0);
    v_vect[0] = fd*r0[0] + gd*v0[0];
    v_vect[1] = fd*r0[1] + gd*v0[1];
    v_vect[2] = fd*r0[2] + gd*v0[2];

    final_state[0]=r_vect[0]; final_state[1]=r_vect[1]; final_state[2]=r_vect[2];
    final_state[3]=v_vect[0]; final_state[4]=v_vect[1]; final_state[5]=v_vect[2];

    if (fabs(f*gd - fd*g - 1.) > 1.0e-2)
        return 1;

    return 0;
}

void kepler_propagator::set_scale_dist(const double &R){
    if(R<=0)
        smartastro_throw("Kepler propagator: negative scaling parameter for distance");
    m_R = R;
}

void kepler_propagator::set_scale_time(const double &T){
    if(T<=0)
        smartastro_throw("Kepler propagator: negative scaling parameter for time");
    m_T = T;
}

void kepler_propagator::set_tol(const double &tol){
    if(tol<=0)
        smartastro_throw("Kepler propagator: negative tolerance");
    m_tol = tol;
}

void kepler_propagator::set_maxiter(const int &maxiter){
    if(maxiter<0)
        smartastro_throw("Kepler propagator: negative maximum number of iterations");
    m_maxiter = maxiter;
}

double kepler_propagator::get_tol() const{
    return m_tol;
}

int kepler_propagator::get_maxiter() const{
    return m_maxiter;
}

double kepler_propagator::get_scale_dist() const{
    return m_R;
}

double kepler_propagator::get_scale_time() const{
    return m_T;
}
