/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_STARK_H
#define SMARTASTRO_STARK_H

#include "base_astrodynamics.h"

namespace smartastro
{
    namespace astrodynamics{

        /**
         * @brief stark is an implementation of the Stark problem for astrodynamics
         *
         * 
         * This dynamics simulates the so-called Stark problem a.k.a. augmented Keplerian motion i.e. Newtonian gravity from a central body plus a constant acceleration
         */
        template < class T >
        class stark: public base_astrodynamics<T>
        {

        protected:
            using base_astrodynamics<T>::m_name;
            using base_astrodynamics<T>::m_L_scale;
            using base_astrodynamics<T>::m_T_scale;
            /**
             * @brief m_mu scaled gravitiational constant
             */  
            double m_mu;
            /**
             * @brief m_acc constant acceleration defining the Stark problem
             */             
            std::vector<double> m_acc;

        public:
            /**
             * @brief stark constructor
             * @param L_scale scaling parameter for distances
             * @param T_scale scaling parameter for time
             * @param mu scaled gravitational constant of central body
             * @param acc vector of constant additional acceleration (scaled) in inertial frame 
             */
            stark(const double &L_scale, const double &t_scale, const double &mu, const std::vector<double> &acc);

            ~stark();
            
            /**
             * @brief stark evaluates dynamics for augmented Keplerian motion
             * @param[in] t time
             * @param[in] x state vector in scaled units
             * @param[in] dx state derivatives in scaled units
             */
            int evaluate(const double &t, const std::vector<T> &x, std::vector<T> &dx) const;

        };
    }
}


#endif // SMARTASTRO_STARK_H
