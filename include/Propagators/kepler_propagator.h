/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/


#ifndef SMARTASTRO_KEPLER_PROPAGATOR_H
#define SMARTASTRO_KEPLER_PROPAGATOR_H

#include "base_propagator.h"

namespace smartastro
{
    namespace propagator{
        /**
         * @brief Analytic orbit propagator for Keplerian orbits. Uses an algorithm from: D. A. Vallado, "Fundamentals of Astrodynamics and Applications, Second
         * Edition", Microcosm Press, pp. 101-102.
         *
         * @author Annalisa Riccardi, Luca Masi, Romain Serra
         */
        class kepler_propagator: public base_propagator<double>
        {
        public:
            /**
             * @brief kepler_propagator constructor
             * @param mu gravitational constant of central body in m^3/(s^2)
             * @param tol Tolerance on the time law. Uses the same number to check if the parameter alpha is zero. Default value 1e-06
             * @param maxiter Maximum number of iterations for convergence of the internal loop. Default value 100.      
             * @param R scaling parameter for distances in meters
             * @param T scaling parameter for times in seconds                             
             */
            kepler_propagator(const double &mu, const double &tol = 1.0e-6, const int &maxiter = 100, const double &R=1.0, const double &T=1.0);
            ~kepler_propagator();

            /**
             * @brief propagate kepler propagation routine
             * @param[in] t0 initial time in units consistent with scaling parameters
             * @param[in] tend final time in units consistent with scaling parameters
             * @param[in] initial_state initial state vector in units consistent with scaling parameters
             * @param[out] final_state final state vector in units consistent with scaling parameters
             * @return exit flag: 0 = success, 1 = check on (f*gd - fd*g == 1) failed, 2 = number of iterations exceeded
             */
            int propagate(const double &t0, const double&tend, const std::vector<double> &initial_state, std::vector<double> &final_state) const;

            //setter and getter
            void set_tol(const double &tol);
            void set_maxiter(const int &maxiter);
            void set_scale_dist(const double &R);
            void set_scale_time(const double &T);
            double get_tol() const;
            int get_maxiter() const;
            double get_scale_dist() const;
            double get_scale_time() const;

        private:
            //TODO Getters and setters for this reference?
            double m_mu;
            double m_tol;
            int m_maxiter;
            double m_R, m_T;

        };
    }
}



#endif // SMARTASTRO_KEPLER_PROPAGATOR_H
