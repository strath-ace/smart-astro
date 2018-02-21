/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_DEARTH_SEMIANALYTICAL_H
#define SMARTASTRO_DEARTH_SEMIANALYTICAL_H

#include "base_astrodynamics.h"
#include "Propagators/perturbation_propagator.h"
#include "Astrodynamics/dearth_3dof.h"
#include "Astrodynamics/dearth_orb.h"
#include "Astro-Core/conversion_coordinates.h"
#include "../exception.h"

namespace smartastro
{
    namespace astrodynamics{
        /**
         * @brief dearth_semianalytical is a dynamics for semi-analytical propagation of Earth orbits
         *
         * 
         * This semi-analytical dynamics is based on the analytical propagator of perturbations developped in CALISPO with some differences, including the first state variable 
         * The state variables are the modified equinoctial ones: p=a*(1-e^2), e*cos(\Omega+\omega), e*sin(\Omega+\omega), tan(i/2)*cos(\Omega), tan(i/2)*sin(\Omega), L=\Omega+\omega+\nu         
         */
        class dearth_semianalytical: public base_astrodynamics<double>
        {

        private:

            /**
             * @brief m_flags list of booleans for orbital perturbations: atmospheric drag, lunisolar effects, SRP (not yet implemented), magnetic force
             */  
            std::vector<bool> m_flags;
            /**
             * @brief m_zonal maximum degree of zonal harmonics taken into account
             */              
            int m_zonal;
            /**
             * @brief m_prop pointer to associated analytical propagator
             */              
            smartastro::propagator::perturbation_propagator *m_prop;
            /**
             * @brief m_params list of parameters
             */             
            std::vector<double> m_params;
            /**
             * @brief m_J2 value for J2 coefficient from EGM96
             */               
            double m_J2 = 0.001082626683553151; // 0.00108262617385222;
            /**
            * @brief m_dyn pointer to 3dof dynamics with Cartesian coordinates to get atmospheric drag if necessary
            */                
            smartastro::astrodynamics::dearth_3dof<double> *m_dyn; 
            /**
             * @brief m_quadrature number of quadrature points for orbital perturbations using Simpson's composite rule (must be even)
             */              
            int m_subdivision;            


        public:

            /**
             * @brief dearth_semianalytical constructor
             * @param flags list of booleans set to true if the corresponding perturbation is taken into account: drag, lunisolar, SRP, magnetic 
             * @param zonal number of zonal harmonics included if this type of perturbation is being considered      
             */
            dearth_semianalytical(const std::vector<bool> &flags = std::vector<bool>(2, false), const int &zonal = 2, const std::vector<double> &params = std::vector<double>(1, 1.0), const std::vector<double> &other = std::vector<double>(2, 0.0), const int &subdivision = 33);

            ~dearth_semianalytical();

            /**
             * @brief evaluates semi-analytical dynamics for Earth orbits
             * @param[in] time in TU
             * @param[in] state modified equinoctial elements in scaled units (DU, TU) 
             * @param[out] state derivative in scaled units (DU, TU)
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<double> &state, std::vector<double> &dstate) const;

            /**
             * @brief 
             * @param[in] jd Julian date
             * @param[in] x osculating orbital elements
             * @param[out] 
             * @return exit flag (0=success)
             */
            int drag_oscul(const double &jd, const std::vector<double> &x, std::vector<double> &dx) const;

            /**
             * @brief 
             * @param[in] jd Julian date
             * @param[in] x osculating orbital elements
             * @param[out] 
             * @return exit flag (0=success)
             */
            int SRP_oscul(const double &jd, const std::vector<double> &x, std::vector<double> &dx) const;

            /**
             * @brief 
             * @param[in] jd Julian date
             * @param[in] x osculating orbital elements
             * @param[out] 
             * @return exit flag (0=success)
             */
            int magnetic_oscul(const double &jd, const std::vector<double> &x, std::vector<double> &dx) const;

            // int gravity_oscul(const double &jd, const std::vector<double> &x, std::vector<double> &dx) const;

        };

    }
}


#endif // SMARTASTRO_DEARTH_SEMIANALYTICAL_H
