/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_PERTURBATION_PROPAGATOR_H
#define SMARTASTRO_PERTURBATION_PROPAGATOR_H

#include "base_propagator.h"
#include "Astro-Core/conversion_coordinates.h"
#include "Astro-Core/conversion_time.h"
#include "Astrodynamics/base_dearth.h"

namespace smartastro
{
    namespace propagator{
        /**
         * @brief perturbation_propagator is a class for first-order analytical propagation of Earth orbits under perturbations
         *
         * 
         * The code analytically propagates Earth orbits taking into account only the first-order effects of the orbital perturbations  
         * The state variables are the modified equinoctial ones: p=a*(1-e^2), e*cos(\Omega+\omega), e*sin(\Omega+\omega), tan(i/2)*cos(\Omega), tan(i/2)*sin(\Omega), L=\Omega+\omega+\nu          */
        class perturbation_propagator: public base_propagator<double>
        {
        public:
            /**
             * @brief perturbation_propagator constructor
             * @param flags list of booleans set to true if the corresponding perturbation is taken into account: lunisolar
             * @param zonal number of zonal harmonics included if this type of perturbation is being considered                          
             */
            perturbation_propagator(const std::vector<bool> &flags = std::vector<bool>(1, false), const int &zonal = 2);

            ~perturbation_propagator();

            /**
             * @brief propagate analytically approximates perturbed Keplerian motion
             * @param[in] L0 initial true longitude
             * @param[in] Lf final true longitude
             * @param[in] initial_state initial modified equinoctial elements in scaled units (DU, TU)
             * @param[out] final_state final modified equinoctial elements in scaled units (DU, TU)
             * @return exit flag (0=success)            
             */
            int propagate(const double &L0, const double&Lf, const std::vector<double> &initial_state, std::vector<double> &final_state) const;

            /**
             * @brief set_timing associates a time to a given initial true longitude for propagation
             * @param[in] L0 true longitude
             * @param[in] jd0 corresponding Julian date
             * @return exit flag (0=success)            
             */
            int set_timing(const double &L0, const double &jd0);

            /**
             * @brief get_timing retrieves initial time and true longitude of propagation
             * @param[out] L0 true longitude
             * @param[out] jd0 Julian date
             * @return exit flag (0=success)            
             */
            int get_timing(double &L0, double &jd0) const;

            /**
             * @brief zonal_harmonics computes contribution from few first zonal harmonics
             * @param[in] L0 initial true longitude
             * @param[in] Lf final true longitude
             * @param[in] initial_state initial modified equinoctial elements in scaled units (DU, TU)
             * @param[out] delta variation of orbital elements in scaled units (DU, TU), beware of: first one is semi-major axis instead of p and the inversion of some coordinates with respect to input
             * @return exit flag (0=success)
             */            
            int zonal_harmonics(const double &L0, const double &Lf, const std::vector<double> &initial_state, std::vector<double> &delta) const;

            /**
             * @brief lunisolar computes gravitational contribution from Moon and Sun
             * @param[in] L0 initial true longitude
             * @param[in] Lf final true longitude
             * @param[in] initial_state initial modified equinoctial elements in scaled units (DU, TU)
             * @param[out] delta variation of orbital elements in scaled units (DU, TU), beware of: first one is semi-major axis instead of p and the inversion of some coordinates with respect to input
             * @return exit flag (0=success)
             */            
            int lunisolar(const double &L0, const double &Lf, const std::vector<double> &initial_state, std::vector<double> &delta) const;
            
            /**
             * @brief integralsJ2 computes integrals needed for J2 contribution
             * @param[in] L0 initial true longitude
             * @param[in] Lf final true longitude
             * @return vector of integrals needed to compute the contribution of J2
             */              
            std::vector<double> integralsJ2(const double &L0, const double &Lf) const;

            /**
             * @brief integralsJ3 computes integrals needed for J2 contribution
             * @param[in] L0 initial true longitude
             * @param[in] Lf final true longitude
             * @param[in] P1 initial value for e*sin(Omega+omega)
             * @param[in] P2 initial value for e*cos(Omega+omega)
             * @param[in] Q1 initial value for tan(i/2)*sin(Omega)
             * @param[in] Q2 initial value for tan(i/2)*cos(Omega)             
             * @return vector of integrals needed to compute the contribution of J3
             */        
            std::vector<double> integralsJ3(const double &L0, const double &Lf, const double &P1, const double &P2, const double &Q1, const double &Q2) const;

            /**
             * @brief integralsJ4 computes integrals needed for J2 contribution
             * @param[in] L0 initial true longitude
             * @param[in] Lf final true longitude
             * @param[in] P1 initial value for e*sin(Omega+omega)
             * @param[in] P2 initial value for e*cos(Omega+omega)
             * @param[in] Q1 initial value for tan(i/2)*sin(Omega)
             * @param[in] Q2 initial value for tan(i/2)*cos(Omega) 
             * @return vector of integrals needed to compute the contribution of J4
             */        
            std::vector<double> integralsJ4(const double &L0, const double&Lf, const double &P1, const double &P2, const double &Q1, const double &Q2) const;

            /**
             * @brief integralsJ5 computes integrals needed for J2 contribution
             * @param[in] L0 initial true longitude
             * @param[in] Lf final true longitude
             * @param[in] P1 initial value for e*sin(Omega+omega)
             * @param[in] P2 initial value for e*cos(Omega+omega)
             * @param[in] Q1 initial value for tan(i/2)*sin(Omega)
             * @param[in] Q2 initial value for tan(i/2)*cos(Omega) 
             * @return vector of integrals needed to compute the contribution of J5
             */        
            std::vector<double> integralsJ5(const double &L0, const double&Lf, const double &P1, const double &P2, const double &Q1, const double &Q2) const;

            /**
             * @brief integrals_3rd_body computes integrals needed for 3rd body effects
             * @param[in] a initila value for semi-major axis
             * @param[in] P1 initial value for e*sin(Omega+omega)
             * @param[in] P2 initial value for e*cos(Omega+omega)
             * @param[in] Q1 initial value for tan(i/2)*sin(Omega)
             * @param[in] Q2 initial value for tan(i/2)*cos(Omega)    
             * @param[in] R distance to 3rd body
             * @param[in] dir vector of directions for 3rd body position          
             * @return vector of integrals needed to compute the contribution of 3rd body
             */        
            std::vector<double> integrals_3rd_body(const double &a, const double &P1, const double &P2, const double &Q1, const double &Q2, const double &R, const std::vector<double> &dir) const;

        private:

            /**
             * @brief m_flags list of booleans for orbital perturbations
             */  
            std::vector<bool> m_flags;
            /**
             * @brief m_zonal maximum degree of zonal harmonics taken into account
             */              
            int m_zonal;
            /**
             * @brief m_timing vector of dates and true longitudes for consistent analytical propagation
             */              
            std::vector<double> m_timing;
            /**
             * @brief m_J2 value for 2nd order zonal harmonics from EGM96
             */                       
            double m_J2 = 0.001082626683553151; // 0.00108262617385222; // constants::J2;
            /**
             * @brief m_J6 value for 3rd order zonal harmonics from EGM96
             */                 
            double m_J3 = -2.532656485353402e-06; // -2.5327e-6;
            /**
             * @brief m_J4 value for 4th order zonal harmonics from EGM96
             */                 
            double m_J4 = -1.6196215914e-06; // -1.6196e-6;
            /**
             * @brief m_J5 value for 5th order zonal harmonics from EGM96
             */                 
            double m_J5 = -2.272960829914134e-07; // -2.2730e-7;

        };
    }
}


#endif // SMARTASTRO_PERTURBATION_PROPAGATOR_H
