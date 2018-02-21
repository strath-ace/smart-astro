/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------- Copyright (C) 2017 University of Strathclyde and Authors ------
------------------- Author: Víctor Rodríguez --------------------------
-------------- e-mail: victor.rodriguez@strath.ac.uk ------------------
------------------- Author: Francesco Torre ---------------------------
-------------- e-mail: francesco.torre@strath.ac.uk -------------------
*/
#ifndef SMARTASTRO_ORBITAL_ELEMENTS_H
#define SMARTASTRO_ORBITAL_ELEMENTS_H

#include <cmath>

#include "exception.h"
#include "smart_vector.h"
#include "constants.h"
#include "conversion_time.h"
#include "conversion_coordinates.h"

namespace smartastro
{
    namespace astrocore
    {
        /** TODO Documentation: this constructor assumes Earth as central body
         * @brief The keplerian_elements class
         */
        template < class T >
        class orbital_elements : public smart_vector<T>
        {
        protected:
            
            double m_period;             // s
            T m_mu;                 // km^3kg^-1*s^-2
            T m_mean_motion;        // 1/s
            T m_periapsis_distance; // m
            T m_M0;                 // rad (initial mean anomaly - related to initial true anomaly)


        public:

            /** TODO Documentation
             * @brief keplerian_elements
             * @param time
             * @param value
             */
            orbital_elements(const double &time,
                             const std::vector<T> &value,
                             const T &mu = smartastro::constants::mu_earth*1e-9);
            

            /**
             *
             */
            orbital_elements(const double &time,
                             const double &period,
                             const std::vector<T> &value);
            

            /**
             *
             */
            ~orbital_elements();


            /**
             *
             */
            T get_semimajor_axis(const bool &in_meters = false) const; //km


            /**
             *
             */
            T get_eccentricity() const;


            /**
             *
             */
            T get_inclination(const bool &in_degrees = false) const;


            /**
             *
             */
            T get_right_ascension(const bool &in_degrees = false) const;


            /**
             *
             */
            T get_argument_of_perigee(const bool &in_degrees = false) const;


            /**
             *
             */
            T get_true_anomaly(const bool &in_degrees = false) const;


            /**
             *
             */
            double get_period() const;


            /**
             *
             */
            T get_mu(bool in_meters = false) const;


            /**
             *
             */
            T get_mean_motion() const;


            /**
             *
             */
            T get_M0() const;


            /**
             *
             */
            T get_periapsis_distance() const;


            /**
             * 
             */
            std::vector<T> get_keplerian_elements(const bool &in_meters=false,
                                                  const bool &in_degrees=false) const;


            /**
             * 
             */
            std::vector<T> to_car(const bool &in_meters = false) const;


            /**
             * 
             */
            std::vector<T> propagate(const double &time);
        };

        //JSON methods
        #if defined(__GNUC__)
            #if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 9))
                //template < class T >
                void to_json(json& j, const orbital_elements<double>& val);
            #endif
        #endif
    }
}

#endif // SMARTASTRO_ORBITAL_ELEMENTS_H
