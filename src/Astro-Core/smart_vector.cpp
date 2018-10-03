/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2016 University of Strathclyde and Authors ------
------------------- Author: Francesco Torre --------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
------------------- Author: Victor Rodriguez -------------------------
-------------- e-mail: victor.rodriguez@strath.ac.uk -----------------
*/

#include "Astro-Core/smart_vector.h"

using namespace smartastro;
using namespace astrocore;


// Default constructor
template < class T >
smart_vector<T>::smart_vector()
{
    m_type  = vector_types::MAX_VECTOR_TYPE;
    m_time  = -std::numeric_limits<double>::infinity();
    m_value = std::vector<T>();
}

template < class T >
smart_vector<T>::smart_vector(const vector_types::VECTOR_TYPE &type, const double &time, const int &length)
{
    m_type  = type;
    m_time  = time;
    m_value = std::vector<T>(length);
}

template< class T >
smart_vector<T>::smart_vector(const vector_types::VECTOR_TYPE &type, const double &time, const std::vector<T> &state)
{
    m_type   = type;
    m_time   = time;
    m_value  = state;
}

template < class T >
smart_vector<T>::~smart_vector()
{
    // NOTHING YET
}


// Methods for type
template < class T >
vector_types::VECTOR_TYPE smart_vector<T>::get_type() const
{
    return m_type;
}

template < class T >
std::string smart_vector<T>::get_type_str()
{
    return vector_types::vector_type_str[m_type];
}


// Methods for time
template < class T >
double smart_vector<T>::get_time() const
{
    return m_time;
}


// Methods for value
template < class T >
std::vector<T> smart_vector<T>::get_value() const
{
    return m_value;
}

template < class T >
const std::vector<T>* smart_vector<T>::get_pointer_to_value() const
{
    return &m_value;
}

template < class T >
T smart_vector<T>::get_value(const unsigned int &index) const
{
    return m_value[index];
}

template < class T >
void smart_vector<T>::set_new_value(const double &time, const std::vector<T> &state)
{
    m_time  = time;
    m_value = state;
}

template < class T >
unsigned int smart_vector<T>::get_length() const
{
    return m_value.size();
}


template class smart_vector<double>;
template class smart_vector<float>;
template class smart_vector<long double>;
#ifdef ENABLE_SMARTUQ
    template class smart_vector<smartuq::polynomial::chebyshev_polynomial>;
    template class smart_vector<smartuq::polynomial::taylor_polynomial>;
#endif
