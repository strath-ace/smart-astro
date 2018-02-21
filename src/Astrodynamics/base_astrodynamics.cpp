/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2015 University of Strathclyde and Authors ------
-------------------- Author: Annalisa Riccardi -----------------------
-------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------
-------------------- Author: Romain Serra ----------------------------
-------------- e-mail: romain.serra@strath.ac.uk ---------------------
-------------------- Author: Francesco Torre -------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
*/

#include "Astrodynamics/base_astrodynamics.h"

using namespace smartastro;
using namespace smartastro::astrodynamics;

// Constructor
template < class T >
base_astrodynamics<T>::base_astrodynamics(const std::string &name,
                                const double &L_scale,
                                const double &T_scale) :
   				smartmath::dynamics::base_dynamics<T>(name)                          
{

    /*Sanity check*/
    if(L_scale<=0.0)
        smartastro_throw("BASE_ASTRODYNAMICS: scaling parameter for distances cannot be negative");
    if(T_scale<=0.0)
        smartastro_throw("BASE_ASTRODYNAMICS: scaling parameter for distances cannot be negative");    

    m_L_scale = L_scale;
    m_T_scale = T_scale;
    m_sc_number = 1;
    m_state_length = 0;

    m_comments = false;
}


// Destructor
template < class T >
base_astrodynamics<T>::~base_astrodynamics()
{
    // NOTHING YET    
}


// Methods for sc_number
template < class T >
unsigned int base_astrodynamics<T>::get_sc_number() const
{
    return m_sc_number;
}

template < class T >
void base_astrodynamics<T>::set_sc_number(const unsigned int &sc_number)
{
    m_sc_number = sc_number;
}


// Methods for state_length
template < class T >
unsigned int base_astrodynamics<T>::get_state_length() const
{
    return m_state_length;
}




// Methods for comments
template < class T >
bool base_astrodynamics<T>::get_comments()
{
    return m_comments;
}

template < class T >
void base_astrodynamics<T>::enable_comments()
{
    m_comments = true;
}

template < class T >
void base_astrodynamics<T>::disable_comments()
{
    m_comments = false;
}


template class base_astrodynamics<double>;
template class base_astrodynamics<float>;
template class base_astrodynamics<long double>;
#ifdef ENABLE_SMARTUQ
template class base_astrodynamics<smartuq::polynomial::chebyshev_polynomial>;
template class base_astrodynamics<smartuq::polynomial::taylor_polynomial>;
#endif
