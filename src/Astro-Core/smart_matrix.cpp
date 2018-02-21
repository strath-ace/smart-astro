/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2017 University of Strathclyde and Authors ------
------------------- Author: Francesco Torre --------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
*/

#include "Astro-Core/smart_matrix.h"

using namespace smartastro;
using namespace astrocore;


// Default constructor
template < class T >
smart_matrix<T>::smart_matrix()
{
    m_type  = matrix_types::MAX_MATRIX_TYPE;
    m_time  = -std::numeric_limits<double>::infinity();
    // m_value = std::vector<T>();
}

template < class T >
smart_matrix<T>::smart_matrix(const matrix_types::MATRIX_TYPE &type)
{
    m_type  = matrix_types::MAX_MATRIX_TYPE;
}

template < class T >
smart_matrix<T>::smart_matrix(const matrix_types::MATRIX_TYPE &type,
                              const double &time,
                              const unsigned int &rows,
                              const unsigned int &cols)
{
    m_type  = matrix_types::MAX_MATRIX_TYPE;
    m_time  = time;
    m_value = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows, cols);
}

template< class T >
smart_matrix<T>::smart_matrix(const matrix_types::MATRIX_TYPE &type,
                              const double &time,
                              const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix)
{
    m_type   = type;
    m_time   = time;
    m_value  = matrix;
}


// Destructor
template < class T >
smart_matrix<T>::~smart_matrix()
{
    // NOTHING YET
}


// Methods for type
template < class T >
matrix_types::MATRIX_TYPE smart_matrix<T>::get_type() const
{
    return m_type;
}

template < class T >
std::string smart_matrix<T>::get_type_str()
{
    return matrix_types::matrix_type_str[m_type];
}


// Methods for time
template < class T >
double smart_matrix<T>::get_time() const
{
    return m_time;
}


// Methods for value
template < class T >
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> smart_matrix<T>::get_value() const
{
    return m_value;
}

template < class T >
T smart_matrix<T>::get_value(const unsigned int &x, const unsigned int &y)
{
    return m_value(x,y);
}

template < class T >
void smart_matrix<T>::set_new_value(const double &time, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix)
{
    m_time  = time;
    m_value = matrix;
}

template < class T >
unsigned int smart_matrix<T>::rows()
{
    return m_value.rows();
}


template < class T >
unsigned int smart_matrix<T>::cols()
{
    return m_value.cols();
}


// // Operator()
// template < class T >
// T& smart_matrix<T>::operator()(const unsigned int &x, const unsigned int& y)
// {
//     if((x > m_value.rows()) || (y > m_value.cols()))
//     {
//         std::cout << "\nError in smart_matrix: invalid (x,y)" << std::endl;
//         std::cout << "matrix.size = (" << m_value.rows() << "," << m_value.cols() << ")" << " but (x,y) = (" << x << "," << y << ")" << std::endl;
//         smartastro_throw("Terminating");
//     }

//     return m_value(x,y);
// }


template class smart_matrix<double>;
template class smart_matrix<float>;
template class smart_matrix<long double>;