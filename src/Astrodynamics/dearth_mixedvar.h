/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_DEARTH_MIXEDVAR_H
#define SMARTASTRO_DEARTH_MIXEDVAR_H

#include "Dynamics/hamiltonian_mixedvar.h"
#include "Astrodynamics/dearth_3dof.h"
#include "../exception.h"
#include "../config.h"

#ifdef ENABLE_SMARTUQ
#include "smartuq.h"
#endif

namespace smartastro
{
    namespace astrodynamics {
        /**
         * @brief dearth_mixedvar is a Hamiltonian formulation for the dynamics of Earth orbits using mixed-variables
         *
         * 
         * This Hamiltonian system with mixed-variables simulates Earth orbits and is suitable ONLY for modified symplectic integration
         * The mixed variables are: position-momentum vector and Poincare variables
         */
        template < class T >
        class dearth_mixedvar: public smartmath::dynamics::hamiltonian_mixedvar<T>
        {

        private:
            using smartmath::dynamics::hamiltonian_mixedvar<T>::m_dim;
            /**
            * @brief m_max_degree_Earth_gravity maximum degree used in geopotential
            */    
            int m_max_degree_Earth_gravity;
            /**
            * @brief m_flags list of booleans for orbital perturbations
            */                
            std::vector<bool> m_flags;
            /**
            * @brief m_C list of C Stokes coefficients
            */                
            std::vector<double> m_C;
            /**
            * @brief m_S list of S Stokes coefficients
            */ 
            std::vector<double> m_S;
            /**
            * @brief m_dyn pointer to 3dof dynamics with Cartesian coordinates to get lunisolar perturbations if necessary
            */                
            smartastro::astrodynamics::dearth_3dof<T> *m_dyn; 

        public:

            /**
             * @brief dearth_mixedvar constructor
             * @param n_max number of zonal harmonics considered in the geopotential 
             * @param flags list of booleans for orbital perturbations: lunisolar      
             */      
            dearth_mixedvar(const int &n_max = 0, const std::vector<bool> &flags = std::vector<bool>(1, false));

            ~dearth_mixedvar();

            /**
             * @brief DHq evaluate dummy implementation preventing use of classic integrator
             * @param[in] t time in TU
             * @param[in] state vector of primary canonical variables in DU, TU
             * @param[out] dstate vector of derivatives of primary canonical variables in DU, TU
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

            /**
             * @brief DHq evaluates partial derivative of Hamiltonian with respect to q
             * @param[in] t time in TU
             * @param[in] q vector of primary canonical position
             * @param[in] p vector of primary canonical momenta 
             * @param[out] dH partial derivative of H with respect to q in scaled units (DU, TU)
             * @return exit flag (0=success)
             */
			int DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const;

            /**
             * @brief DHp evaluates partial derivative of Hamiltonian with respect to p
             * @param[in] t time in TU
             * @param[in] q vector of primary canonical position
             * @param[in] p vector of primary canonical momenta 
             * @param[out] dH partial derivative of H with respect to p in scaled units (DU, TU)
             * @return exit flag (0=success)
             */
            int DHp(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const;

            /**
             * @brief DHq2 evaluates partial derivative of Hamiltonian with respect to q2
             * @param[in] t time in TU
             * @param[in] q2 vector of secondary set of canonical position
             * @param[in] p2 vector of secondary set of canonical momenta 
             * @param[out] dH partial derivative of H with respect to q2 in scaled units (DU, TU)
             * @return exit flag (0=success)
             */
            int DHq2(const double &t, const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &dH) const;

            /**
             * @brief DHp2 evaluates partial derivative of Hamiltonian with respect to p2
             * @param[in] t time in TU
             * @param[in] q2 vector of secondary set of canonical position
             * @param[in] p2 vector of secondary set of canonical momenta 
             * @param[out] dH partial derivative of H with respect to p2 in scaled units (DU, TU)
             * @return exit flag (0=success)
             */
            int DHp2(const double &t, const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &dH) const;

            /**
             * @brief conversion performs transformation from primary set of canonical variables (position-momentum) to secondary one (Poincare variables)
             * @param[in] q vector of primary canonical position
             * @param[in] p vector of primary canonical momenta 
             * @param[out] q2 vector of secondary set of canonical position
             * @param[out] p2 vector of secondary set of canonical momenta 
             * @return exit flag (0=success)
             */
            int conversion(const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &q2, std::vector<T> &p2) const; 

            /**
             * @brief conversion2 performs transformation from secondary set of canonical variables (Poincare variables) to primary one (position-momentum)
             * @param[in] q2 vector of secondary canonical position
             * @param[in] p2 vector of secondary canonical momenta 
             * @param[out] q vector of primary set of canonical position
             * @param[out] p vector of primary set of canonical momenta 
             * @return exit flag (0=success)
             */
            int conversion2(const std::vector<T> &q2, const std::vector<T> &p2, std::vector<T> &q, std::vector<T> &p) const; 
        };

    }
}

#endif // SMARTASTRO_DEARTH_MIXEDVAR_H
