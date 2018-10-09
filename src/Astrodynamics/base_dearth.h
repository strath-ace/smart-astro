/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_BASE_DEARTH_H
#define SMARTASTRO_BASE_DEARTH_H

#include "base_astrodynamics.h"

namespace smartastro
{
    namespace astrodynamics{
        /**
         * @brief base_dearth is an abstract class for 3dof and 6dof propagators of Earth orbits
         *
         * 
         * The dynamics simulates Keplerian motion as well as a range of gravitational and non-graviational perturbations. 
         * The first six components of the state vector consist of position and velocity written with Cartesian coordinates (scaled) in an Earth-centered inertial frame. The seventh state is the mass in kg.
         */
        template < class T >
        class base_dearth: public base_astrodynamics<T>
        {

        protected:
            using base_astrodynamics<T>::m_L_scale;
            using base_astrodynamics<T>::m_T_scale;
            /**
             * @brief m_max_degree_Earth_gravity degree at which the geopotential is expanded
             */            
            int m_max_degree_Earth_gravity;
            /**
             * @brief m_max_degree_Earth_gravity degree at which the magnetic potential is expanded
             */            
            int m_max_degree_Earth_magnetic;
            /**
             * @brief m_flags list of booleans for perturbations switches
             */            
            std::vector<int> m_flags;
            /**
             * @brief m_drag list of parameters used in particular for atmospheric drag
             */            
            std::vector<T> m_drag;
            /**
             * @brief m_other list of additional parameters
             */                
            std::vector<T> m_other;  
            /**
             * @brief m_C Cmn Stockes coefficients for geopotential
             */                        
            std::vector<double> m_C;
            /**
             * @brief m_S Smn Stockes coefficients for geopotential
             */              
            std::vector<double> m_S;     
            /**
             * @brief m_F10dot7 set of parameters to simulate a periodic variation in solar flux: first vector is list of dates of, second is actual value for linear interpolation
             */  
            std::vector< std::vector<double> > m_F10dot7;

        public:
            /**
             * @brief base_dearth constructor
             * @param scaling factor in meters for distance
             * @param scaling factor in seconds for time
             * @param expansion order for geopotential
             * @param flags for perturbations: first one is drag, second is lunisolar, third is radiation pressure (1 for conical shadow, 2 for cylindrical one), fourth is magnetic force
             * @param drag parameters: cross-sectional area in square meters, drag coefficient, mean solar flux, geomagnetic index  
             * @param other parameters: see 3dof and 6dof instantiations
             */
            base_dearth(const std::string &name, const double &L_scale, const double &t_scale, const int &n, const std::vector<int> &flags, const std::vector<T> &p, const std::vector<T> &p2, const std::vector< std::vector<double> > &F10dot7);

            virtual ~base_dearth();

           /**
            * @brief computes gradient of potential in spherical harmonics
            * @param[in] position vector in scaled units
            * @param[in] n_max expansion order     
            * @param[in] Re scaled body radius          
            * @param[in] C un-normalized coefficients
            * @param[in] S un-normalized coefficients
            * @param[out] gradient of potential in spherical harmonics
            * @return exit flag (0=success)
            */
           static int harmonics(const std::vector<T> &pos, const int &n_max, const double &Re, const std::vector<double> &C, const std::vector<double> &S, std::vector<T> &grad);

            /**
             * @brief computes acceleration due to Earth's gravity
             * @param[in] Julian date
             * @return exit flag (0=success)
             */
            int tidal_corrections(const double &jd, std::vector<double> &C, std::vector<double> &S) const;

            /**
             * @brief computes acceleration due to Earth's gravity
             * @param[in] Julian date
             * @param[in] state vector in scaled units
             * @param[out] acceleration in scaled units
             * @return exit flag (0=success)
             */
            int Earth_gravity(const double &jd, const std::vector<T> &x, std::vector<T> &acc) const;

            /**
            * @brief computes altitude and geodetic latitude  
            * @param[in] position vector in Earth-fixed frame in scaled units
            * @param[out] altitude in km
            * @param[out] geodetic latitude in radians
            * @return exit flag (0=success)
            */
            int geodetic(const std::vector<T> &pos, T &Z, T &phi) const;

            /**
            * @brief computes atmospheric density according to Jacchia-Gill model 
            * @param[in] jd Julian date
            * @param[in] Z altitude in km
            * @param[in] phi geodetic latitude in rad
            * @param[in] alpha right ascension in Earth-fixed frame in rad
            * @param[out] rho atmospheric density in kg/m^3
            * @return exit flag (0=success)
            */
            int density(const double &jd, const T &Z, const T &phi, const T &alpha, T &rho) const;

            /**
            * @brief computes Sun's distance to Earth and right ascension in inertial frame  
            * @param[in] julian date             
            * @param[out] Sun's distance to Earth center in km 
            * @param[out] Sun's right ascension in inertial frame in radians
            * @return exit flag (0=success)
            */
            static int Sun(const double &jd, double &r, double &lambda);

            /**
            * @brief computes Moon's position from low-precision ephemerides in ecliptic plane
            * @param[in] julian date             
            * @param[out] Moon's position in km wrt ecliptic plane
            * @return exit flag (0=success)
            */
            static int moon_pos(const double &jd, std::vector<double> &pos_moon);

            /**
            * @brief computes Moon's position from low-precision ephemerides in reference inertial frame
            * @param[in] julian date             
            * @param[out] Moon's position in scaled units
            * @return exit flag (0=success)
            */
            int Moon(const double &jd, std::vector<double> &pos_moon) const;

            /**
             * @brief computes acceleration due to lunisolar perturbation
             * @param[in] Julian date
             * @param[in] state vector in scaled units
             * @param[out] acceleration in scaled units
             * @return exit flag (0=success)
             */
            int lunisolar(const double &jd, const std::vector<T> &x, std::vector<T> &acc) const;

            /**
            * @brief computes shadow function 
            * @param[in] Julian date             
            * @param[in] state vector in scaled units
            * @return shadow function (0=umbra, 1 otherwise)
            */
            T shadow(const double &jd, const std::vector<T> &x) const;

            /**
            * @brief computes Gaussian coefficients for Earth magnetic field
            * @param[in] Julian date             
            * @param[in] maximum degree of expansion
            * @param[out] g coefficients
            * @param[out] h coefficients
            * @return exit flag (0=success)
            */
            int Gaussian_coeff(const double &jd, const int &n_max, std::vector<double> &g, std::vector<double> &h) const;

            /**
            * @brief computes electron temperature according to partial implementation of IRI90
            * @param[in] Julian date             
            * @param[in] altitude in km
            * @param[in] geodetic latitude in rad
            * @param[in] local time in rad
            * @param[out] electron temperature in K
            * @return exit flag (0=success)
            */
            int IRI90(const double &jd, const T &Z, const T &phi, const T &local_time, T &Te) const;

            /**
            * @brief numerically solves for the surface potential in the current balance equation assuming only electron and proton 
            * @param[in] electron temperature in K           
            * @param[in] shape factor (0.0 for a flat plate, 0.5 for a cylinder and 1.0 for a sphere)     
            * @param[in] spacecraft velocity with respect to the ionosphere in m/s   
            * @param[in] cosinus of the angle relative to the SC velocity               
            * @return surface potential in volts
            */
            T surface_potential_LEO(const T &Te, const double &shape_factor, const T &Vsc, const T &cos_vel) const;

            /**
            * @brief computes electrostatic surface charge due to ambient plasma
            * @param[in] julian date     
            * @param[in] state (position-velocity-mass) vector in scaled units
            * @param[out] electrostatic capacitance in SI units
            * @param[out] electrostatic charge in SI units
            * @return exit flag (0=success)
            */
            virtual int charging(const double &jd, const std::vector<T> &x, T &C, T &V) const=0;

            /**
            * @brief computes solar flux as a function of time given periodic variation
            * @param[in] julian date     
            * @param[in] reference solar flux
            * @param[out] solar flux at given date 
            * @return exit flag (0=success)
            */
            int var_flux(const double &jd, const T &F0, T &F) const;

            /**
            * @brief logarithmic standard density according to Jacchia-Gill model
            * @param[in] altitude in km
            * @return logarithmic standard density
            */
            static double density_log(std::vector<float> X);
            static double density_log(std::vector<double> X);
            static double density_log(std::vector<long double> X);
            #ifdef ENABLE_SMARTUQ
            static smartuq::polynomial::taylor_polynomial density_log(const std::vector<smartuq::polynomial::taylor_polynomial> &X);
            static smartuq::polynomial::chebyshev_polynomial density_log(const std::vector<smartuq::polynomial::chebyshev_polynomial> &X);
            #endif

            static double density_log_sigmoid_Tinf(std::vector<double> X);

            /**
            * @brief logarithmic Helium density according to Jacchia-Gill model
            * @param[in] altitude in km
            * @return logarithmic standard density
            */
            static float He_number_density_log(std::vector<float> X);
            static double He_number_density_log(std::vector<double> X);
            static long double He_number_density_log(std::vector<long double> X);
            #ifdef ENABLE_SMARTUQ
            static smartuq::polynomial::taylor_polynomial He_number_density_log(const std::vector<smartuq::polynomial::taylor_polynomial> &X);
            static smartuq::polynomial::chebyshev_polynomial He_number_density_log(const std::vector<smartuq::polynomial::chebyshev_polynomial> &X);
            #endif

            /**
            * @brief computes the function x->sign(x)*x^2
            * @param[in] x
            * @return sign(x)*x^2
            */
            static float signed_square(float X);
            static double signed_square(double X);
            static long double signed_square(long double X);
            #ifdef ENABLE_SMARTUQ
            static smartuq::polynomial::taylor_polynomial signed_square(const smartuq::polynomial::taylor_polynomial &X);
            static smartuq::polynomial::chebyshev_polynomial signed_square(const smartuq::polynomial::chebyshev_polynomial &X);
            #endif


        };
    }
}


#endif // SMARTASTRO_BASE_DEARTH_H
