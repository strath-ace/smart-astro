/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2017 University of Strathclyde and Authors -------
------------------- Author: Víctor Rodríguez --------------------------
-------------- e-mail: victor.rodriguez@strath.ac.uk ------------------
------------------- Author: Francesco Torre ---------------------------
-------------- e-mail: francesco.torre@strath.ac.uk -------------------
*/

#include "Astro-Core/orbital_elements.h"

using namespace smartastro;
using namespace astrocore;

// Constructor: default and with given mu
template < class T >
orbital_elements<T>::orbital_elements(const double &time,
                                      const std::vector<T> &value,
                                      const T &mu):
    smart_vector<T>::smart_vector(vector_types::VECTOR_TYPE::KEPLERIAN_ELEMENTS, time, value)
{
    //Sanity checks
    if (value.size() != 6)
        smartastro_throw("The value of the keplerian elements vector must have 6 elements");
    
    
    m_mu = mu;

    m_period = smartastro::constants::pi2*sqrt(pow(get_semimajor_axis(), 3)/mu);

    m_mean_motion = sqrt(m_mu/pow(get_semimajor_axis(), 3));

    //Periapsis distance in m, based on the semimajor axis and the eccentricity
    m_periapsis_distance = (get_semimajor_axis(true))*(1.0 -get_eccentricity());

    double M0d;
    conversion_time::true2mean_anomaly(get_true_anomaly(), get_eccentricity(), M0d);
    m_M0 = (T)M0d;
}


// Constructor with given period
template < class T >
orbital_elements<T>::orbital_elements(const double &time,
                                      const double &period,
                                      const std::vector<T> &value):
    smart_vector<T>::smart_vector(vector_types::VECTOR_TYPE::KEPLERIAN_ELEMENTS, time, value)
{
    //Sanity checks
    if (value.size() != 6)
        smartastro_throw("The value of the keplerian elements vector must have 6 elements");
    
    m_period = period;
    
    m_mu = pow(smartastro::constants::pi2/period, 2)*pow(get_semimajor_axis(), 3);

    m_mean_motion = sqrt(m_mu/pow(get_semimajor_axis(),3));

    //Periapsis distance in m, based on the semimajor axis and the eccentricity
    m_periapsis_distance = (get_semimajor_axis(true))*(1.0 -get_eccentricity());

    double M0d;
    conversion_time::true2mean_anomaly(get_true_anomaly(), get_eccentricity(), M0d);
    m_M0 = (T)M0d;
}


// Destructor
template < class T >
orbital_elements<T>::~orbital_elements()
{
    // NOTHING YET
}


// Methods for semimajor axis
template < class T >
T orbital_elements<T>::get_semimajor_axis(const bool &in_meters) const
{
    return (in_meters   ?   this->m_value[0]*1000.0
                        :   this->m_value[0]);
}


// Methods for eccentricity
template < class T >
T orbital_elements<T>::get_eccentricity() const
{
    return this->m_value[1];
}


// Methods for inclunation
template < class T >
T orbital_elements<T>::get_inclination(const bool &in_degrees) const
{
    return (in_degrees  ?   this->m_value[2]*smartastro::constants::rad2deg
                        :   this->m_value[2]);
}


// Methods for right ascension
template < class T >
T orbital_elements<T>::get_right_ascension(const bool &in_degrees) const
{
    return (in_degrees  ?   this->m_value[3]*smartastro::constants::rad2deg
                        :   this->m_value[3]);
}


// Methods for argument of perigee
template < class T >
T orbital_elements<T>::get_argument_of_perigee(const bool &in_degrees) const
{
    return (in_degrees  ?   this->m_value[4]*smartastro::constants::rad2deg
                        :   this->m_value[4]);
}


// Methods for true_anomaly
template < class T >
T orbital_elements<T>::get_true_anomaly(const bool &in_degrees) const
{
    return (in_degrees  ?   this->m_value[5]*smartastro::constants::rad2deg
                        :   this->m_value[5]);
}


// Methods for period
template < class T >
double orbital_elements<T>::get_period() const
{
    return m_period;
}


// Methods for mu
template < class T >
T orbital_elements<T>::get_mu(bool in_meters) const
{
    return (in_meters   ?   m_mu*1.0e9
                        :   m_mu);
}


// Methods for mean_motion
template < class T >
T orbital_elements<T>::get_mean_motion() const
{
    return m_mean_motion;
}


// Methods for M0
template < class T >
T orbital_elements<T>::get_M0() const
{
    return m_M0;
}


// Methods for periapsis distance
template < class T >
T orbital_elements<T>::get_periapsis_distance() const
{
    return m_periapsis_distance;
}


// Get Keplerian elements
template < class T >
std::vector<T> orbital_elements<T>::get_keplerian_elements(const bool &in_meters, const bool &in_degrees) const
{
    std::vector<T> result(6);

    result[0] = get_semimajor_axis(in_meters);
    result[1] = get_eccentricity();
    result[2] = get_inclination(in_degrees);
    result[3] = get_right_ascension(in_degrees);
    result[4] = get_argument_of_perigee(in_degrees);
    result[5] = get_true_anomaly(in_degrees);

    return result;
}


// To Cartesian
template < class T >
std::vector<T> orbital_elements<T>::to_car(const bool &in_meters) const
{
    std::vector<double> result(6);
    conversion_coordinates::kep2car(std::vector<double>(this->m_value.begin(), this->m_value.end()), m_mu, result);

    if (in_meters)
    {
        for (size_t i=0; i<result.size(); i++)
        {
            result[i]*=1000.0;
        }
    }

    return std::vector<T>(result.begin(), result.end());
}


// Propagate
template < class T >
std::vector<T> orbital_elements<T>::propagate(const double &time)
{
    double delta_time = (time -this->m_time)/m_period;
    delta_time = (delta_time -(long)delta_time)*m_period;
    
    T Mf = m_mean_motion*delta_time +m_M0;
    if(Mf >= 2.0*constants::pi)
        Mf -= 2.0*constants::pi;
    
    double nif;
    astrocore::conversion_time::mean2true_anomaly(Mf, this->m_value[1], nif);

    if(nif < 0.0)
        nif += 2.0*constants::pi;

    std::vector<double> new_kep(this->m_value.begin(), this->m_value.end());
    // new_kep = this->m_value;
    new_kep[5] = nif;

    std::vector<double> result(6);
    astrocore::conversion_coordinates::kep2car(new_kep, m_mu, result);

    return std::vector<T>(result.begin(), result.end());
}


//JSON methods (namespace method)
#if defined(__GNUC__)
    #if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 9))
        //template < class T >
        void astrocore::to_json(json& j, const orbital_elements<double>& val)
        {
            j["time"] = val.get_time();
            j["semimajor_axis"] = val.get_semimajor_axis(false);
            j["eccentricity"] = val.get_eccentricity();
            j["inclination"] = val.get_inclination();
            j["right_ascension"] = val.get_right_ascension();
            j["argument_of_perigee"] = val.get_argument_of_perigee();
            j["true_anomaly"] = val.get_true_anomaly();
        }
    #endif
#endif


template class orbital_elements<double>;
template class orbital_elements<float>;
template class orbital_elements<long double>;
// #ifdef ENABLE_SMARTUQ
//     template class orbital_elements<smartuq::polynomial::chebyshev_polynomial>;
//     template class orbital_elements<smartuq::polynomial::taylor_polynomial>;
// #endif
