/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/dearth_mixedvar.h"
#include "Astrodynamics/base_dearth.h"
#include "Astro-Core/conversion_coordinates.h"
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

using namespace smartastro;
using namespace smartastro::astrodynamics;


template<class T>
dearth_mixedvar<T>::dearth_mixedvar(const int &n_max, const std::vector<bool> &flags): smartmath::dynamics::hamiltonian_mixedvar<T>("Hamiltonian dynamics for Earth orbits with mixed variables", 3, true), m_max_degree_Earth_gravity(n_max), m_flags(flags){

    if(m_flags.size() < 1)
        smartastro_throw("DEARTH_MIXEDVAR: there must be at least one flag for perturbations");

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
dearth_mixedvar<T>::~dearth_mixedvar()
{
      
}   

template < class T >
int dearth_mixedvar<T>::evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const{

    /* sanity checks */
    if(state.size() != 2 * m_dim)
        smartastro_throw("EVALUATE: the Hamiltonian state must have a consistent dimension");          
   
    dstate = state;
    /* reconstituting the canonical variables q and p from the state vector */
    std::vector<T> q, p;
    for(unsigned int k = 0; k < m_dim; k++)
    {
        q.push_back(state[k]);
        p.push_back(state[k + m_dim]);
    }

    /* computing complete dynamics in Cartesian coordinates */
    std::vector<T> grad(3, 0.0 * q[0]);
    base_dearth<T>::harmonics(q, m_max_degree_Earth_gravity, 1.0, m_C, m_S, grad);

    /* computing the state derivative */
    for(unsigned int k = 0; k < m_dim; k++)
    {
        dstate[k] = p[k];
        dstate[k + m_dim] = grad[k];
    }

    if(m_flags[0])
    {
        std::vector<T> third_body(3, 0.0 * state[0]);
        m_dyn->lunisolar(t * constants::T_Earth / (3600.0 * 24.0), q, third_body);
        for(unsigned int i = 0; i < m_dim; i++)
            dstate[i + m_dim] += third_body[i];
    }

    return 0;
}

template < class T >
int dearth_mixedvar<T>::DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{
    
    /* sanity checks */
    if(q.size() != m_dim)
        smartastro_throw("DHQ: the position must have dimension 3");
    
    std::vector<T> grad(3, 0.0 * q[0]);
    base_dearth<T>::harmonics(q, m_max_degree_Earth_gravity, 1.0, m_C, m_S, grad);
    T r = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
    T factor = 1.0 / (r * r * r);
    dH.clear();
    for(unsigned int i = 0; i < m_dim; i++)
        dH.push_back(- grad[i] - q[i] * factor);

    if(m_flags[0])
    {
        std::vector<T> dH2(3, 0.0 * q[0]);
        m_dyn->lunisolar(t * constants::T_Earth / (3600.0 * 24.0), q, dH2);
        for(unsigned int i = 0; i < m_dim; i++)
            dH[i] -= dH2[i];
    }

    return 0;
}    

template < class T >
int dearth_mixedvar<T>::DHp(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{

    /* sanity checks */
    if(p.size() != m_dim)
        smartastro_throw("DHP: the momentum must have the correct dimension");
    
    dH.clear();
    for(unsigned int i = 0; i < m_dim; i++)
        dH.push_back(p[i]);

    return 0;

}

template < class T >
int dearth_mixedvar<T>::DHq2(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH2) const{
    
    /* sanity checks */
    if(q.size() != m_dim)
        smartastro_throw("DHQ2: the angle vector must have dimension 3");
    
    dH2.clear();
    for(unsigned int i = 0; i < m_dim; i++)
        dH2.push_back(0.0 * q[0]);

    return 0;
}    

template < class T >
int dearth_mixedvar<T>::DHp2(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const{

    /* sanity checks */
    if(p.size() != m_dim)
        smartastro_throw("DHP2: the action vector must have the correct dimension");
    
    dH.clear();
    dH.push_back(1.0 / (p[0] * p[0] * p[0]));
    dH.push_back(0.0 * p[0]);
    dH.push_back(0.0 * p[0]);

    return 0;

}

template < class T >
int dearth_mixedvar<T>::conversion(const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &q2, std::vector<T> &p2) const{
  
    q2 = q;
    p2 = p;
    std::vector<T> car;
    for(unsigned int i = 0; i < m_dim; i++)
        car.push_back(q[i]);
    for(unsigned int i = 0; i < m_dim; i++)
        car.push_back(p[i]);
    std::vector<T> elements(6, 0.0 * q[0]);

    /* easy implementation of Delaunay with doubles */
    // std::vector<T> kep = car;
    // astrocore::conversion_coordinates::car2kep(car, 1.0, kep);
    // astrocore::conversion_coordinates::kep2delaunay(kep, 1.0, elements);
    // std::cout << "From Cartesian to elements" << std::endl;
    // std::cout << elements[0] << " " << elements[1] << " " << elements[2] << " " << elements[3] << " "  << elements[4] << " "  << elements[5] << std::endl;

    /* polynomial-friendly version for Delaunay */
    // T r = sqrt(car[0] * car[0] + car[1] * car[1] + car[2] * car[2]);
    // T v2 = car[3] * car[3] + car[4] * car[4] + car[5] * car[5];
    // T energy = v2 / 2.0 - 1.0 / r; 
    // T a = - 1.0 / (2.0 * energy);
    // elements[3] = sqrt(a);
    // std::vector<T> h(3, 0.0 * q[0]);
    // h[0] = car[1] * car[5] - car[2] * car[4];
    // h[1] = car[2] * car[3] - car[0] * car[5];
    // h[2] = car[0] * car[4] - car[1] * car[3];
    // T h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
    // T e2 = 1.0 + 2.0 * energy * h2;
    // elements[4] = sqrt(h2);//elements[3] * sqrt(1.0 - e2);
    // T ci = h[2] / elements[4];
    // elements[5] = elements[4] * ci;
    // T n = sqrt(h[0] * h[0] + h[1] * h[1]);    
    // elements[2] = atan2(h[0], -h[1]); //acos(-h[1] / n);
    // if(-elements[2] > 0.0)
    //     elements[2] += 2.0 * smartastro::constants::pi;
    // T inter1 = v2 - 1.0 / r; 
    // T inter2 = car[0] * car[3] + car[1] * car[4] + car[2] * car[5];
    // std::vector<T> e(3, 0.0 * q[0]); // n0 = -h1; n1 = h0;
    // for(int i = 0; i < 3; i++)
    //     e[i] = inter1 * car[i] - inter2 * car[i + 3];
    // T eccentricity = sqrt(e2);       
    // elements[1] = acos((-e[0] * h[1] + e[1] * h[0]) / (eccentricity * n)); 
    // if(-e[2] > 0.0)
    //     elements[1] = 2.0 * smartastro::constants::pi - elements[1];
    // T theta = acos((e[0] * car[0] + e[1] * car[1] + e[2] * car[2]) / (eccentricity * r)); 
    // if(-inter2 > 0.0)
    //     theta = 2.0 * smartastro::constants::pi - theta;  
    // T E = 2.0 * atan(sqrt((1.0 - eccentricity) / (1.0 + eccentricity)) * tan(theta / 2.0)); 
    // if(-E > 0.0)
    //     E += 2.0 * smartastro::constants::pi;    
    // elements[0] = E - eccentricity * sin(E);

    /* polynomial-friendly version for Poincare */
    T r = sqrt(car[0] * car[0] + car[1] * car[1] + car[2] * car[2]);
    T v2 = car[3] * car[3] + car[4] * car[4] + car[5] * car[5];
    T energy = v2 / 2.0 - 1.0 / r; 
    T a = - 1.0 / (2.0 * energy);
    elements[3] = sqrt(a);
    std::vector<T> h(3, 0.0 * q[0]);
    h[0] = car[1] * car[5] - car[2] * car[4];
    h[1] = car[2] * car[3] - car[0] * car[5];
    h[2] = car[0] * car[4] - car[1] * car[3];
    T h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
    T e2 = 1.0 + 2.0 * energy * h2;
    T ci = h[2] / sqrt(h2);
    T n = sqrt(h[0] * h[0] + h[1] * h[1]);    
    T RAAN = atan2(h[0], -h[1]); //acos(-h[1] / n);
    if(-RAAN > 0.0)
        RAAN += 2.0 * smartastro::constants::pi;
    T inter1 = v2 - 1.0 / r; 
    T inter2 = car[0] * car[3] + car[1] * car[4] + car[2] * car[5];
    std::vector<T> e(3, 0.0 * q[0]); // n0 = -h1; n1 = h0;
    for(int i = 0; i < 3; i++)
        e[i] = inter1 * car[i] - inter2 * car[i + 3];
    T eccentricity = sqrt(e2);       
    T aop = acos((-e[0] * h[1] + e[1] * h[0]) / (eccentricity * n)); 
    if(-e[2] > 0.0)
        aop = 2.0 * smartastro::constants::pi - aop;
    T theta = acos((e[0] * car[0] + e[1] * car[1] + e[2] * car[2]) / (eccentricity * r)); 
    if(-inter2 > 0.0)
        theta = 2.0 * smartastro::constants::pi - theta;  
    T E = 2.0 * atan(sqrt((1.0 - eccentricity) / (1.0 + eccentricity)) * tan(theta / 2.0)); 
    if(-E > 0.0)
        E += 2.0 * smartastro::constants::pi;    
    T M = E - eccentricity * sin(E);
    T W = RAAN + aop;
    elements[0] = M + W;
    T sq = sqrt(1.0 - e2);
    T aux = sqrt(2.0 * elements[3] * (1.0 - sq));
    elements[1] = aux * cos(W);
    elements[4] = aux * sin(W); 
    T aux2 = sqrt(2.0 * elements[3] * sq * (1.0 - ci));
    elements[2] = aux2 * cos(RAAN);
    elements[5] = aux2 * sin(RAAN);     

    for(unsigned int i = 0; i < m_dim; i++)
    {
        q2[i] = elements[i];
        p2[i] = elements[i + m_dim];
    }
    // std::cout << "From Cartesian to elements" << std::endl;
    // std::cout << std::setprecision(10) << elements[0] << " " << elements[1] << " " << elements[2] << " " << elements[3] << " "  << elements[4] << " "  << elements[5] << std::endl;
    // std::cout << std::endl;

    return 0;
}; 


template < class T >
int dearth_mixedvar<T>::conversion2(const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &q, std::vector<T> &p) const{
    
    q = q2;
    p = p2;
    
    std::vector<T> elements;
    for(unsigned int i = 0; i < m_dim; i++)
        elements.push_back(q[i]);
    for(unsigned int i = 0; i < m_dim; i++)
        elements.push_back(p[i]);
    std::vector<T> car(6, 0.0 * q[0]);

    /* easy implementation of Delaunay with doubles */
    // double pi = smartastro::constants::pi;
    // std::vector<T> kep = elements;
    // astrocore::conversion_coordinates::delaunay2kep(elements, 1.0, kep);
    // kep[5] -= 2.0 * pi * double(floor(kep[5] / (2.0 * pi)));
    // astrocore::conversion_coordinates::kep2car(kep, 1.0, car);
    // std::cout << "From elements to Cartesian" << std::endl;
    // std::cout << car[0] << " " << car[1] << " " << car[2] << " " << car[3] << " "  << car[4] << " "  << car[5] << std::endl;

    /* polynomial-friendly version for Delaunay */
    // T M0 = elements[0];
    // T eccentricity = sqrt(1.0 - (elements[4] * elements[4]) / (elements[3] * elements[3]));
    // T E = M0, F_E = 0.0 * q[0], DF_E = 0.0 * q[0]; 
    // for(int k = 0; k < 10; k++){
    //     F_E = E - eccentricity * sin(E) - M0;
    //     DF_E = 1.0 - eccentricity * cos(E);
    //     E -= F_E / DF_E;      
    // }
    // T theta = 2.0 * atan(sqrt((1.0 + eccentricity) / (1.0 - eccentricity)) * tan(E / 2.0));
    // T a = elements[3] * elements[3];
    // T r = a * (1.0 - eccentricity * cos(E));
    // T v = sqrt(2.0 / r - 1.0 / a);
    // //T gamma = atan2(eccentricity * sin(theta), sqrt(1.0 - eccentricity * eccentricity));
    // T ci = elements[5] / elements[4];
    // T si = sqrt(1.0 - ci * ci);
    // T cO = cos(elements[2]);
    // T sO = sin(elements[2]);
    // T angle1 = theta + elements[1];
    // //T angle2 = angle1 - gamma;
    // T factor = 1.0 / sqrt(1.0 + eccentricity * (2.0 * cos(theta) + eccentricity));
    // T c1 = cos(angle1);
    // T s1 = sin(angle1);
    // T cGamma = (1.0 + eccentricity * cos(theta)) * factor;
    // T sGamma = eccentricity * sin(theta) * factor;
    // T c2 = c1 * cGamma + s1 * sGamma;
    // T s2 = s1 * cGamma - c1 * sGamma;   
    // car[0] = r * (c1 * cO - s1 * ci * sO);
    // car[1] = r * (c1 * sO + s1 * ci * cO);
    // car[2] = r * s1 * si;
    // car[3] = v * (-s2 * cO - c2 * ci * sO);
    // car[4] = v * (-s2 * sO + c2 * ci * cO);
    // car[5] = v * c2 * si;    

    /* polynomial-friendly version for Poincare */
    T W = atan2(elements[4], elements[1]);
    T cW = cos(W);
    T sW = sin(W);
    T RAAN = atan2(elements[5], elements[2]);
    T cO = cos(RAAN);
    T sO = sin(RAAN);    
    T a = elements[3] * elements[3];
    T aux = (elements[1] * elements[1] + elements[4] * elements[4]) / (2.0 * elements[3]);
    T sq = 1.0 - aux;
    T e2 = 1.0 - sq *sq;
    T eccentricity = sqrt(e2);
    T K = elements[0], F_K = 0.0 * q[0], DF_K = 0.0 * q[0]; 
    for(int k = 0; k < 10; k++)
    {
        F_K = K + eccentricity * (sW * cos(K) - cW * sin(K)) - elements[0];
        DF_K = 1.0 - eccentricity * (sW * sin(K) + cW * cos(K));
        K -= F_K / DF_K;      
    }
    T r = a * (1.0 - eccentricity * (sW * sin(K) + cW * cos(K)));
    T v = sqrt(2.0 / r - 1.0 / a);
    T ci = 1.0 - ((1.0 - sq) / sq) * (elements[2] * elements[2] + elements[5] * elements[5]) / (elements[1] * elements[1] + elements[4] * elements[4]);
    T si = sqrt(1.0 - ci * ci);
    T theta = 2.0 * atan(sqrt((1.0 + eccentricity) / (1.0 - eccentricity)) * tan((K - W) / 2.0));
    T angle1 = theta + W - RAAN;
    T factor = 1.0 / sqrt(1.0 + eccentricity * (2.0 * cos(theta) + eccentricity));
    T c1 = cos(angle1);
    T s1 = sin(angle1);
    T cGamma = (1.0 + eccentricity * cos(theta)) * factor;
    T sGamma = eccentricity * sin(theta) * factor;
    T c2 = c1 * cGamma + s1 * sGamma;
    T s2 = s1 * cGamma - c1 * sGamma;   
    car[0] = r * (c1 * cO - s1 * ci * sO);
    car[1] = r * (c1 * sO + s1 * ci * cO);
    car[2] = r * s1 * si;
    car[3] = v * (-s2 * cO - c2 * ci * sO);
    car[4] = v * (-s2 * sO + c2 * ci * cO);
    car[5] = v * c2 * si;   
    // std::cout << "From elements to Cartesian" << std::endl;
    // std::cout << car[0] << " " << car[1] << " " << car[2] << " " << car[3] << " "  << car[4] << " "  << car[5] << std::endl;

    for(unsigned int i = 0; i < m_dim; i++)
    {
        q[i] = car[i];
        p[i] = car[i + m_dim];
    }

    return 0;
}; 


template class dearth_mixedvar<double>;
#ifdef ENABLE_SMARTUQ
template class dearth_mixedvar<smartuq::polynomial::chebyshev_polynomial>;
template class dearth_mixedvar<smartuq::polynomial::taylor_polynomial>;
#endif