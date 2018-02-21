/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/


#ifndef SMARTASTRO_DEARTH_HAMILTONIAN_H
#define SMARTASTRO_DEARTH_HAMILTONIAN_H

#include "Dynamics/hamiltonian_momentum.h"
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
         * @brief dearth_hamiltonian is a Hamiltonian formulation for the dynamics of Earth orbits
         *
         * 
         * This Hamiltonian dynamics simulates Earth orbits and is suitable for symplectic integration
         */
        template < class T >
        class dearth_hamiltonian: public smartmath::dynamics::hamiltonian_momentum<T>
        {

        private:
            using smartmath::dynamics::hamiltonian_momentum<T>::m_dim;
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
             * @brief dearth_hamiltonian constructor
             * @param n_max number of zonal harmonics considered in the geopotential     
             * @param flags list of booleans for orbital perturbations: lunisolar     
             */            
            dearth_hamiltonian(const int &n_max = 0, const std::vector<bool> &flags = std::vector<bool>(1, false));

            ~dearth_hamiltonian();

            /**
             * @brief DHq evaluates partial derivative of Hamiltonian with respect to q
             * @param[in] t time in TU
             * @param[in] q vector of canonical position
             * @param[in] p vector of canonical momenta 
             * @param[out] dH partial derivative of H with respect to q in scaled units (DU, TU)
             * @return exit flag (0=success)
             */
			int DHq(const double &t, const std::vector<T> &q, const std::vector<T> &p, std::vector<T> &dH) const;

        };

    }
}

#endif // SMARTASTRO_DEARTH_HAMILTONIAN_H
