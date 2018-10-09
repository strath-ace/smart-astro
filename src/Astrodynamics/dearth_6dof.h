/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#ifndef SMARTASTRO_DEARTH_6DOF_H
#define SMARTASTRO_DEARTH_6DOF_H

#include "base_dearth.h"

namespace smartastro
{
    namespace astrodynamics{
        /**
         * @brief dearth_6dof is a 6 degrees of freedom propagator for Earth orbits
         *
         * 
         * This 6dof dynamics simulates Keplerian motion as well as a range of gravitational and non-graviational perturbations around the Earth. 
         * The state vector consists of Cartesian coordinates in the Earth-centered inertial frame, mass, quaternions and rotational rates in the principal frame of the rigid body
         */
        template < class T >
        class dearth_6dof: public base_dearth<T>
        {

        private:
            using base_dearth<T>::m_L_scale;
            using base_dearth<T>::m_T_scale;
            using base_dearth<T>::m_max_degree_Earth_gravity;
            using base_dearth<T>::m_max_degree_Earth_magnetic;
            using base_dearth<T>::m_flags;
            using base_dearth<T>::m_drag;
            using base_dearth<T>::m_other;
            using base_dearth<T>::m_C;
            using base_dearth<T>::m_S;
            /**
             * @brief m_I list of principal moments of inertia
             */              
            std::vector<T> m_I;
            /**
             * @brief m_object integer defining the geometry of the object 
             */              
            int m_object;
            /**
             * @brief m_CoS coordinates of the center of symmetry w.r.t. the barycenter
             */               
            std::vector<T> m_CoS;

        public:

            using base_dearth<T>::harmonics;
            using base_dearth<T>::Earth_gravity;
            using base_dearth<T>::geodetic;
            using base_dearth<T>::density;
            using base_dearth<T>::density_log;
            using base_dearth<T>::He_number_density_log;  
            using base_dearth<T>::density_log_sigmoid_Tinf;            
            using base_dearth<T>::Sun;
            using base_dearth<T>::lunisolar;
            using base_dearth<T>::shadow;
            using base_dearth<T>::Gaussian_coeff;
            using base_dearth<T>::IRI90;
            using base_dearth<T>::surface_potential_LEO;
            using base_dearth<T>::signed_square;

            /**
             * @brief dearth_3dof constructor
             * @param principal axis of inertia in kg*m2 around the center of gravity
             * @param scaling factor in meters for distance
             * @param scaling factor in seconds for time
             * @param expansion order for geopotential
             * @param flags for perturbations: first one is drag, second is lunisolar, third is radiation pressure, fourth is magnetic force
             * @param drag parameters: cross-sectional area in square meters, reference length in meters, mean solar flux and geomagnetic index 
             * @param flag for object's geometry: 0 for square flat plate (0.01x1x1m), 1 for cylinder (0.6m for radius, 4.93m for length) and 2 for sphere
             * @param Cartesian coordinates for center of symmetry in body-fixed frame in meters         
             * @param other parameters: specular reflectivity, diffuse reflection (optional), constant surface charge in Coulomb (optional) 
             */
            dearth_6dof(const std::vector<T> &I=std::vector<T>(3,1.0),const double &r_scale = 1.0, const double &t_scale = 1.0, const int &n=0, const std::vector<int> &flags=std::vector<int>(4,0), const std::vector<T> &p=std::vector<T>(1,0.0), const int &object=0, const std::vector<T> &coord=std::vector<T>(3,0.0), const std::vector<T> &m_other=std::vector<T>(1,0.3), const std::vector< std::vector<double> > &F10dot7 = std::vector< std::vector<double> >(2)) ;

            ~dearth_6dof();

            /**
             * @brief evaluate high-fidelity Earth-orbiting dynamics in celestial Earth-centered frame
             * @param[in] time in scaled units
             * @param[in] state (position-velocity-mass-quaternion-rate) vector in scaled units
             * @param[out] state derivative in scaled units
             * @return exit flag (0=success)
             */
            int evaluate(const double &t, const std::vector<T> &state, std::vector<T> &dstate) const;

            /**
             * @brief compute rotation matrix from body-fixed to inertial frame
             * @param[in] quaternions
             * @param[out] rotation matrix (returned in columns)
             * @return exit flag (0=success)
             */
            int orientation(const std::vector<T> &q,std::vector<T> &C) const;

            /**
             * @brief computes acceleration and torque due to atmospheric drag
             * @param[in] Julian date
             * @param[in] state (position-velocity-mass-quaternion-rate) vector in scaled units
             * @param[out] acceleration in scaled units
             * @param[out] torque in scaled units
             * @return exit flag (0=success)
             */
            int aero(const double &jd, const std::vector<T> &x,std::vector<T> &acc,std::vector<T> &tor) const;

            /**
            * @brief computes acceleration due to solar radiation pressure
            * @param[in] Julian date            
            * @param[in] state vector (position-velocity-mass-quaternion-rate) in scaled units
            * @param[out] acceleration in scaled units due to solar radiation pressure 
            * @param[out] torque in scaled units   
            * @return exit flag (0=success)
            */
            int SRP(const double &jd,const std::vector<T> &x, std::vector<T> &acc, std::vector<T> &tor) const; 

            /**
            * @brief computes acceleration due to Earth's magnetic field
            * @param[in] Julian date            
            * @param[in] state (position-velocity-mass) vector in scaled units
            * @param[out] acceleration in scaled units due to magnetic force
            * @param[out] torque in scaled units
            * @return exit flag (0=success)
            */
            int magnetism(const double &jd,const std::vector<T> &x, std::vector<T> &acc, std::vector<T> &tor) const; 
            
            /**
            * @brief computes electrostatic surface charge due to ambient plasma
            * @param[in] julian date     
            * @param[in] state (position-velocity-mass) vector in scaled units
            * @param[out] electrostatic capacitance in SI units
            * @param[out] electrostatic charge in SI units
            * @return exit flag (0=success)
            */
            int charging(const double &jd, const std::vector<T> &x, T &C, T &V) const;                     

        };
    }
}


#endif // SMARTASTRO_DEARTH_6DOF_H
