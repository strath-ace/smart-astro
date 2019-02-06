/******************************************************************************
 *                               CONVERSION_COORDINATES_H                     *
 *            List of conversion functions of the SMART-SIM toolbox           *
 ******************************************************************************/

#ifndef SMARTASTRO_CONVERSION_COORDINATES_H
#define SMARTASTRO_CONVERSION_COORDINATES_H

extern "C" {
#include <cspice/SpiceUsr.h>
#include <cspice/SpiceZfc.h>
}

#include <vector>
#include <cmath>
#include "LinearAlgebra/Eigen/Dense"
#include "../exception.h"
#include "../constants.h"
#include "conversion_time.h"

namespace smartastro{
    namespace astrocore{

        /**
         * @brief The conversion_space class implements a set of static methods for reference frame conversion.
         */
        class conversion_coordinates{
        public:

            /******************************************************************************
             *                    REFERENCE FRAME CONVERSION FUNCTIONS                    *
             ******************************************************************************/
            /**
             * @brief  Function to convert from Cartesian position and velocity to Keplerian elements.
             *
             * Function to convert from Cartesian position and velocity to Keplerian elements.
             * All the units have to be consistent, angles in radians.\n
             *  WARNINGS:
             *      - a circular orbit it considered for ecc < 1e-10,
             *      - an elliptic orbit for 1e-10 <= ecc <= 1 - 1e-10,
             *      - a parabola for 1 - 1e-10 < ecc < 1 + 1e-10
              *      - an hyperbola for ecc >= 1 + 1e-10.     *
             * @param[in] car vector of 6 doubles containing the position
             * and velocity vector in cartesian coordinates.
             * @param[in] mu planetary constant of the central body.
             * @param[out] kep vector of 6 doubles containing the mean
             * keplerian elements at date.\n
             * kep[0] = semimajor axis in km\n
             * kep[1] = eccentricity\n
             * kep[2] = inclination in rad (0<=kep[2]<=PI)\n
             * kep[3] = right ascension of the ascending
             *          node in rad (0<=kep[3]<PI2)\n
             * kep[4] = argument of perigee in rad
             *          (0<=kep[4]<PI2)\n
             * kep[5] = true anomaly in rad
             *          (0<=kep[5]<PI2)\n
             * @return Error code
             *
             * @see David A. Vallado, "Fundamentals of Astrodynamics and Applications",
             *      2007
             * @author Massimiliano Vasile 2002, Federico Zuiani 2008,
             */
            static bool car2kep(const std::vector<double> &car, const double &mu, std::vector<double> &kep);

            /**
             * @brief  Function to convert from Cartesian reference frame to B-plane reference frame.
             *
             * Function to convert from Cartesian reference frame to B-plane reference frame.\n
             * Cartesian reference frame: {x,y,z} inertial reference frame.
             * The b-plane is the plane perpendicular to the incoming relative velocity
             * of the small body at the planet arrival, containing the planet.
             * b-plane reference frame: {xi,eta,zeta} where
             * - eta-axis: as the incoming relative velocity of the small body on its
             *   arrival.
             * - zeta-axis: in the b-plane, direction opposite to the projection of
             *   the heliocentric velocity of the planet on the b-plane.
             * - xi-axis: in the b-plane, completes the reference frame.
             *
             * @param[in] x_car vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {x,y,z}.
             * @param[in] u_car vector of 3 doubles containing the
             *                  velocity of the small body relative to the
             *                  planet, expressed in {x,y,z}.
             * @param[in] vp_car vector of 3 doubles containing the
             *                   velocity of the planet, expressed in {x,y,z}.
             * @param[out] x_bpl vector of 3 doubles containing the
             *                   coordinates of the vector x_car expressed in
             *                   {xi,eta,zeta}.
             * @return Error code
             *
             * @author Massimiliano Vasile 2007, Luca Masi 2008,
             */
            static bool car2bpl(const std::vector<double> &x_car, const std::vector<double> &u_car, const std::vector<double> &vp_car, std::vector<double> &x_bpl);

            /**
             * @brief  Function to convert from rth reference frame to B-plane reference frame.
             *
             * Function to convert from rth reference frame to B-plane reference frame.\n
             * The b-plane is the plane perpendicular to the incoming relative velocity
             * of the small body at the planet arrival, containing the planet.
             * b-plane reference frame: {xi,eta,zeta} where
             * - eta-axis: as the incoming relative velocity of the small body on its
             *   arrival.
             * - zeta-axis: in the b-plane, direction opposite to the projection of
             *   the heliocentric velocity of the planet on the b-plane.
             * - xi-axis: in the b-plane, completes the reference frame.
             *
             * @param[in] x_rth vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {r,t,h}.
             * @param[in] s_car Position and Velocity [L/T] of the small body,
             *                  expressed in {x,y,z}.
             * @param[in] vp_car vector of 3 doubles containing the
             *                   velocity of the planet, expressed in {x,y,z}.
             * @param[out] x_bpl vector of 3 doubles containing the
             *                   coordinates of the vector x_car expressed in
             *                   {xi,eta,zeta}.
             * @return Error code
             *
             * @author Massimiliano Vasile 2014, Victor Rodriguez 2017
             */
            static bool rth2bpl(
                    const std::vector<double> &x_rth,
                    const std::vector<double> &s_car,
                    const std::vector<double> &vp_car,
                    std::vector<double> &x_bpl
                    );

            /**
             * @brief Function to transform  from Cartesian reference frame to Radial-Trasversal-h reference frame.
             *
             * Function to transform  from Cartesian reference frame to Radial-Trasversal-h reference frame.\n
             * Cartesian reference frame: {x,y,z} inertial reference frame.\n
             * RTh reference frame: {r,t,h}
             *      - r-axis: direction of the orbit radius
             *      - h-axis: direction of angular momentum
             *      - t-axis: in the orbit plane, completes the reference frame (inward)
             * @param[in] x_car vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {x,y,z}.
             * @param[in] s_car vector of 6 doubles containing the
             *                  position and velocity of the orbiting body,
             *                  expressed in {x,y,z}
             * @param[out] x_rth vector of 3 doubles containing the
             *                   coordinates of the vector x_car expressed in
             *                   {r,t,h}.
             * @return Error code:
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             */
            static bool car2rth(const std::vector<double> &x_car, const std::vector<double> &s_car, std::vector<double> &x_rth);

            /**
             * @brief Function to transform from Cartesian reference frame to Tangent-Normal-h reference frame.
             *
             * Function to transform from Cartesian reference frame to Tangent-Normal-h reference frame.\n
             * Cartesian reference frame: {x,y,z} inertial reference frame.\n
             * TNh reference frame: {t,n,h}
             *      - t-axis: tangent to the motion
             *      - h-axis: direction of angular momentum
             *      - n-axis: inward normal to t, in the orbit plane
             * @param[in] x_car vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {x,y,z}.
             * @param[in] s_car vector of 6 doubles containing the
             *                  position and velocity of the orbiting body,
             *                  expressed in {x,y,z}
             * @param[out] x_tnh  vector of 3 doubles containing the
             *                    coordinates of the vector x_car expressed in
             *                    {t,n,h}.
             * @return Error code
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             *
             */
            static bool car2tnh(const std::vector<double> &x_car, const std::vector<double> &s_car, std::vector<double> &x_tnh);

            /**
             * @brief Function to transform from Cartesian reference frame to Right Ascension
             * and Declination reference frame.
             *
             * Function to transform from Cartesian reference frame to Right Ascension
             * and Declination reference frame.
             * Cartesian reference frame: {x,y,z} inertial reference frame.\n
             * RADec reference frame: {r,alpha,delta} right ascension and declination
             * (spherical equatorial) reference frame.
             *      - r = modulus of the vector
             *      - alpha = in-plane right ascention angle, counted from the tangential
             *              direction to the projection of the vector on the orbital
             *              plane [rad]
             *      - delta = out-of-plane declination angle from the projection of the
             *              vector on the orbital plane up to the vector itself [rad]
             *
             * @param[in] x_car vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {x,y,z}.
             * @param[in] s_car vector of 6 doubles containing the
             *                  position and velocity of the orbiting body,
             *                  expressed in {x,y,z}
             * @param[out] x_radec vector of 3 doubles containing the
             *                     coordinates of the vector x_car expressed in
             *                     {r,alpha,delta}.
             * @return Error code
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             */
            static bool car2radec(const std::vector<double> &x_car, const std::vector<double> &s_car, std::vector<double> &x_radec);


            /**
             * @brief  Function to convert from Keplerian orbital elements to Cartesian coordinates.
             *
             * Function to convert from Keplerian orbital elements to Cartesian coordinates.
             * All units to be consistent each other. Angles in radians.\n
             * WARNINGS:
             *      - A circular orbit it considered for ecc < 1e-10,
             *        an elliptic orbit for 1e-10 <= ecc <= 1 - 1e-10,
             *        a parabola for 1 - 1e-10 < ecc < 1 + 1e-10
             *        an hyperbola for ecc >= 1 + 1e-10.
             *      - In case of a parabola, kep[0] must be the radius of pericentre,
             *        and not the semi-major axis (Inf for a parabola).
             *      - In the case of hyperbola, theta must be such that the point is on
             *        the physical leg of the hyperbola (the leg around the attracting
             *        body).
             *
             * @param[in] kep vector of 6 doubles containing the
             *                mean keplerian elements at date.\n
             *                kep[0] = semimajor axis in km\n
             *                kep[1] = eccentricity\n
             *                kep[2] = inclination in rad\n
             *                kep[3] = right ascension of the ascending node in rad\n
             *                kep[4] = argument of perigee in rad\n
             *                kep[5] = true anomaly in rad
             * @param[in] mu  planetary constant of the central body.
             * @param[out] car vector of 6 doubles containing the
             *                 position and velocity vector in cartesian
             *                 coordinates.
             * @return Error code
             *
             * @author Massimiliano Vasile 2002. Luca Masi 2008,
             */
            static bool kep2car(const std::vector<double> &kep, const double &mu, std::vector<double> &car);


            /**
             * @brief Function to transform from Radial-Trasversal-h reference frame to Cartesian reference frame.
             *
             * Function to transform from Radial-Trasversal-h reference frame to Cartesian reference frame.\n
             * Cartesian reference frame: {x,y,z} inertial reference frame.\n
             * RTh reference frame: {r,t,h}
             *      - r-axis: direction of the orbit radius
             *      - h-axis: direction of angular momentum
             *      - t-axis: in the orbit plane, completes the reference frame (inward)
             * @param[in] x_rth vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {r,t,h}.
             * @param[in] s_car vector of 6 doubles containing the
             *                  position and velocity of the orbiting body,
             *                  expressed in {x,y,z}
             * @param[out] x_car vector of 3 doubles containing the
             *                   coordinates of the vector x_rth expressed in
             *                   {x,y,z}.
             * @return Error code
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008,
             */
            static bool rth2car(const std::vector<double> &x_rth, const std::vector<double> &s_car, std::vector<double> &x_car);

            /**
             * @brief Function to transform from radial-transversal-h reference frame to
             * tangent-normal-h reference frame.
             *
             * Function to transform from radial-transversal-h reference frame to
             * tangent-normal-h reference frame.
             * This functions does not work for parabola.
             *
             * RTh reference frame: {r,t,h}
             *      - r-axis: direction of the orbit radius
             *      - h-axis: direction of angular momentum
             *      - t-axis: in the orbit plane, completes the reference frame (inward)
             *
             * TNh reference frame: {t,n,h}
             *      - t-axis: tangent to the motion
             *      - h-axis: direction of angular momentum
             *      - n-axis: inward normal to t, in the orbit plane
             *
             * @param[in] x_rth vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {r,t,h}.
             * @param[in] a double constaining the semi-major axis.
             * @param[in] e double constaining the eccentricity.
             * @param[in] f double constaining the true anomaly from the
             *              pericenter in rad.
             * @param[in] mu double constaining the gravitational constant of
             *               the central body.
             * @param[out] x_tnh vector of 3 doubles containing the
             *                   coordinates of the vector x_rth expressed in
             *                   {t,n,h}.
             * @return Error code (false if e is within 1 +/- 1e-10 (parabola),
             *                     true otherwise).
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             */
            static bool rth2tnh(const std::vector<double> &x_rth, const double &a, const double &e, const double &f, const double &mu, std::vector<double> &x_tnh);


            /**
             * @brief  Function to transform from Tangent-Normal-h reference frame to Cartesian reference frame.
             *
             * Function to transform from Tangent-Normal-h reference frame to Cartesian reference frame.\n
             * Cartesian reference frame: {x,y,z} inertial reference frame.\n
             * TNh reference frame: {t,n,h}
             *      - t-axis: tangent to the motion
             *      - h-axis: direction of angular momentum
             *      - n-axis: inward normal to t, in the orbit plane
             * @param[in] x_tnh vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {t,n,h}.
             * @param[in] s_car vector of 6 doubles containing the
             *                  position and velocity of the orbiting body,
             *                  expressed in {x,y,z}
             * @param[out] x_car vector of 3 doubles containing the
             *                   coordinates of the vector x_tnh expressed in
             *                   {x,y,z}.
             * @return Error code
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             */
            static bool tnh2car(const std::vector<double> &x_tnh, const std::vector<double> &s_car, std::vector<double> &x_car);

            /**
             * @brief Function to transform from Tangent-Normal-h reference frame to
             * Right Ascension and Declination reference frame.
             *
             * Function to transform from Tangent-Normal-h reference frame to
             * Right Ascension and Declination reference frame.\n
             *
             * RADec reference frame: {r,alpha,delta} right ascension and declination
             * (spherical equatorial) reference frame.
             *      - r = modulus of the vector
             *      - alpha = in-plane right ascention angle, counted from the tangential
             *              direction to the projection of the vector on the orbital
             *              plane [rad]
             *      - delta = out-of-plane declination angle from the projection of the
             *              vector on the orbital plane up to the vector itself [rad]
             *
             * TNh reference frame: {t,n,h}
             *      - t-axis: tangent to the motion
             *      - h-axis: direction of angular momentum
             *      - n-axis: inward normal to t, in the orbit plane
             * @param[in] x_tnh vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {t,n,h}.
             * @param[out] x_radec vector of 3 doubles containing the
             *                     coordinates of the vector x_tnh expressed in
             *                     {r,alpha,delta}.
             * @return Error code
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             */
            static bool tnh2radec(const std::vector<double> &x_tnh, std::vector<double> &x_radec);

            /**
             * @brief Function to transform from radial-transversal-h reference frame to
             * tangent-normal-h reference frame.
             *
             * Function to transform from radial-transversal-h reference frame to
             * tangent-normal-h reference frame.
             *
             * RTh reference frame: {r,t,h}
             *      - r-axis: direction of the orbit radius
             *      - h-axis: direction of angular momentum
             *      - t-axis: in the orbit plane, completes the reference frame (inward)
             *
             * TNh reference frame: {t,n,h}
             *      - t-axis: tangent to the motion
             *      - h-axis: direction of angular momentum
             *      - n-axis: inward normal to t, in the orbit plane
             *
             * @param[in] x_tnh vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {t,n,h}.
             * @param[in] a double constaining the semi-major axis.
             * @param[in] e double constaining the eccentricity.
             * @param[in] f double constaining the true anomaly from the
             *              pericenter in rad.
             * @param[in] mu double constaining the gravitational constant of
             *               the central body.
             * @param[out] x_rth vector of 3 doubles containing the
             *                   coordinates of the vector x_rth expressed in
             *                   {r,t,h}.
             * @return Error code (false if a is zero, true otherwise)
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             */
            static bool tnh2rth(const std::vector<double> &x_tnh, const double &a, const double &e, const double &f, const double &mu, std::vector<double> &x_rth);


            /**
             * @brief Function to transform from Right Ascension and Declination reference frame to
             * tangent-normal-h reference frame.
             *
             * Function to transform from Right Ascension and Declination reference frame to
             * tangent-normal-h reference frame.\n
             * RADec reference frame: {r,alpha,delta} right ascension and declination
             * (spherical equatorial) reference frame.
             *      - r = modulus of the vector
             *      - alpha = in-plane right ascention angle, counted from the tangential
             *              direction to the projection of the vector on the orbital
             *              plane [rad]
             *      - delta = out-of-plane declination angle from the projection of the
             *              vector on the orbital plane up to the vector itself [rad]
             *
             * TNh reference frame: {t,n,h}
             *     - t-axis: tangent to the motion
             *     - h-axis: direction of angular momentum
             *     - n-axis: inward normal to t, in the orbit plane
             *
             * @param[in] x_radec vector of 3 doubles containing the
             *                    coordinates of the vector to be transformed,
             *                    expressed in {r,alpha,delta}.
             *
             * @param[out] x_tnh vector of 3 doubles containing the
             *                   coordinates of the vector x_radec expressed in
             *                   {t,n,h}.
             * @return Error code
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             */
            static bool radec2tnh(const std::vector<double> &x_radec, std::vector<double> &x_tnh);

            /**
             * @brief Function to transform from Right Ascension and Declination reference frame to
             *        Cartesian reference frame.
             *
             * Function to transform from Right Ascension and Declination reference frame to
             * Cartesian reference frame.\n
             * Cartesian reference frame: {x,y,z} inertial reference frame.\n
             * RADec reference frame: {r,alpha,delta} right ascension and declination
             * (spherical equatorial) reference frame.
             *      - r = modulus of the vector
             *      - alpha = in-plane right ascention angle, counted from the tangential
             *              direction to the projection of the vector on the orbital
             *              plane [rad]
             *      - delta = out-of-plane declination angle from the projection of the
             *              vector on the orbital plane up to the vector itself [rad]
             * @param[in] x_radec vector of 3 doubles containing the
             *                    coordinates of the vector to be transformed,
             *                    expressed in {r,alpha,delta}.
             * @param[in] s_car vector of 6 doubles containing the
             *                  position and velocity of the orbiting body,
             *                  expressed in {x,y,z}
             * @param[out] x_car vector of 3 doubles containing the
             *                   coordinates of the vector x_radec expressed in
             *                   {x,y,z}.
             * @return Error code
             *
             * @author Luca Masi 2006, Massimiliano Vasile 2007, Federico Zuiani 2008
             */
            static bool radec2car(const std::vector<double> &x_radec, const std::vector<double> &s_car, std::vector<double> &x_car);

            /**
             * @brief Function converting Cartesian coordinates into spherical coordinates
             *
             * Function converting Cartesian coordinates into spherical coordinates
             *        
             * @param[in] vector of Cartesian coordinates (position-velocity)
             * @param[out] vector of spherical coordinates: radius, longitude (rad), latitude (rad), radial velocity, tranversal velocity and perpendicular velocity (?)
             *                   
             * @return Error code
             *
             * @author 
             */
            static bool car2spher(const std::vector<double> &car, std::vector<double> &spher);

            /**
             * @brief Function converting Keplerian coordinates into spherical coordinates
             *
             * Function converting Keplerian coordinates into spherical coordinates (scaling must be consistent)
             *        
             * @param[in] vector of Keplerian coordinates 
             * @param[in] gravitational constant of central body 
             * @param[out] vector of spherical coordinates: radius, longitude (rad), latitude (rad), radial velocity, tranversal velocity and perpendicular velocity (?)
             *                   
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool kep2spher(const std::vector<double> &car, const double &mu, std::vector<double> &spher);

            /**
             * @brief Function converting modified equinoctial elements into Keplerian coordinates
             *
             * Function converting modified equinoctial elements into Keplerian coordinates (scaling must be consistent)
             * Walker, M. J. H., Ireland, B., & Owens, J. (1985). A set modified equinoctial orbit elements. Celestial Mechanics and Dynamical Astronomy, 36(4), 409-419.             
             *        
             * @param[in] vector of modified equinoctial elements 
             * @param[out] vector of Keplerian coordinates
             *                   
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool modeq2kep(const std::vector<double> &modeq, std::vector<double> &kep);

            /**
             * @brief Function converting modified equinoctial elements into Cartesian coordinates
             *
             * Function converting modified equinoctial elements into Cartesian coordinates (scaling must be consistent)
             * Walker, M. J. H., Ireland, B., & Owens, J. (1985). A set modified equinoctial orbit elements. Celestial Mechanics and Dynamical Astronomy, 36(4), 409-419.             
             *     
             * @param[in] vector of modified equinoctial elements 
             * @param[in] gravitational constant of central body  
             * @param[out] vector of Cartesian coordinates
             *                   
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool modeq2car(const std::vector<double> &modeq, const double &mu, std::vector<double> &car);

            /**
             * @brief Function converting Keplerian coordinates into modified equinoctial elements 
             *
             * Function converting Keplerian coordinates into modified equinoctial elements (scaling must be consistent)
             * Walker, M. J. H., Ireland, B., & Owens, J. (1985). A set modified equinoctial orbit elements. Celestial Mechanics and Dynamical Astronomy, 36(4), 409-419.
             *        
             * @param[in] vector of Keplerian coordinates             
             * @param[out] vector of modified equinoctial elements 
             *                   
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool kep2modeq(const std::vector<double> &kep, std::vector<double> &modeq);

            /**
             * @brief Function converting Delaunay variables into Keplerian coordinates 
             *
             * Function converting Delaunay variables into Keplerian coordinates (scaling must be consistent)
             *        
             * @param[in] vector of Delaunay variables  
             * @param[in] gravitational constant of central body                          
             * @param[out] vector of Keplerian coordinates
             *                   
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool delaunay2kep(const std::vector<double> &delaunay, const double &mu, std::vector<double> &kep);

            /**
             * @brief Function converting Keplerian coordinates into Delaunay variables 
             *
             * Function converting Keplerian coordinates into Delaunay variables  (scaling must be consistent)
             *        
             * @param[in] vector of Keplerian coordinates 
             * @param[in] gravitational constant of central body                          
             * @param[out] vector of Delaunay variables  
             *                   
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool kep2delaunay(const std::vector<double> &kep, const double &mu, std::vector<double> &delaunay); 

            /**
             * @brief Function converting osculating Keplerian elements into mean-J2 Keplerian elements 
             *
             * Function converting osculating Keplerian elements into mean-J2 Keplerian elements (scaling must be consistent)
             *        
             * @param[in] vector of Keplerian coordinates 
             * @param[in] gravitational constant of central body  
             * @param[in] radius of central body                      
             * @param[out] vector of mean Keplerian elements (under J2 only) 
             *                   
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool oscul2mean(const std::vector<double> &oscul, const double &mu, const double & Re, std::vector<double> &mean);

            /**
             * @brief Function converting mean-J2 Keplerian elements into osculating Keplerian elements
             *
             * Function converting mean-J2 Keplerian elements into osculating Keplerian elements (scaling must be consistent)
             *        
             * @param[in] vector of mean Keplerian elements (under J2 only) 
             * @param[in] gravitational constant of central body  
             * @param[in] radius of central body                      
             * @param[out] vector of Keplerian coordinates  
             *                   
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool mean2oscul(const std::vector<double> &mean, const double &mu, const double & Re, std::vector<double> &oscul); 

                /*  EULER_AXIS_ANGLE - Rotates a vector about an axis of a given angle
                 *                     (Euler axis and angle rotation).
                 *
                 *  This function uses the right-hand rule.
                 *
                 *  PROTOTYPE
                 *      int euler_axis_angle(const double v[3],const double n[3],
                 *                           const double theta, double v1[3])
                 *
                 *  INPUTS
                 *      (double *) v            Pointer to a vector of 3 doubles containing
                 *                              the coordinates of the vector to be rotated.
                 *      (double *) n            Pointer to a vector of 3 doubles containing
                 *                              the coordinates of the axis of rotation.
                 *      (double) theta          Angle of rotation in radians.
                 *
                 *  OUTPUTS
                 *      (double *) vi           Pointer to a vector of 3 doubles containing
                 *                              the coordinates of the rotated vector.
                 *      (int) euler_axis_angle  Error code: 1 if the norm of n is null
                 *                                          0 otherwise.
                 *
                 *  NON-STANDARD LIBRARIES
                 *      math_utils
                 *
                 *  REFERENCES
                 *      (none)
                 *
                 *  ORIGINAL VERSION
                 *      Matteo Ceriotti, 11/01/2007, MATLAB, euler_axis_angle.m
                 *
                 *  AUTHOR
                 *      Nicolas Croisard, 19/09/2008
                 * 
                 *  CHANGELOG
                 *      14/05/2007, Revised by Camilla Colombo
                 *      19/09/2008, Nicolas Croisard: Conversion of the M-file into C
                 *
                 * ------------------------ - SpaceART Toolbox - ------------------------ */ 

            static bool euler_axis_angle(const std::vector<double> &v, const std::vector<double> &n, const double &theta, std::vector<double> &v1);           

/*
--------- SPICE functions implemented by: Scott Hurley (GitHub: Scott_James_Hurley)----
--------- e-mail: scott.james.hurley.97@gmail.com -------------------------------------
*/

/**
             * @brief  Function to convert from Cartesian reference frame to B-plane reference frame.
             *
             * Function to convert from Cartesian reference frame to B-plane reference frame.\n
             * Cartesian reference frame: {x,y,z} inertial reference frame.
             * The b-plane is the plane perpendicular to the incoming relative velocity
             * of the small body at the planet arrival, containing the planet.
             * b-plane reference frame: {xi,eta,zeta} where
             * - eta-axis: as the incoming relative velocity of the small body on its
             *   arrival.
             * - zeta-axis: in the b-plane, direction opposite to the projection of
             *   the heliocentric velocity of the planet on the b-plane.
             * - xi-axis: in the b-plane, completes the reference frame.
             *
             * @param[in] x_car vector of 3 doubles containing the
             *                  coordinates of the vector to be transformed,
             *                  expressed in {x,y,z}.
             * @param[in] u_car vector of 3 doubles containing the
             *                  velocity of the small body relative to the
             *                  planet, expressed in {x,y,z}.
             * @param[in] vp_car vector of 3 doubles containing the
             *                   velocity of the planet, expressed in {x,y,z}.
             * @param[out] x_bpl vector of 3 doubles containing the
             *                   coordinates of the vector x_car expressed in
             *                   {xi,eta,zeta}.
             * @return Error code
             *
             * @author Massimiliano Vasile 2007, Luca Masi 2008,
             */

	    static int radrec_( std::vector<double> &range,  std::vector<double> &ra,  std::vector<double> &dec,  std::vector<double> &rectan);

	    static int recrad_( std::vector<double> &rectan,  std::vector<double> &range,  std::vector<double> &ra,  std::vector<double> &dec);

	    static int cylrec_( std::vector<double> &r__,  std::vector<double> &long__,  std::vector<double> &z__,  std::vector<double> &rectan);

	    static int reccyl_( std::vector<double> &rectan,  std::vector<double> &r__,  std::vector<double> &long__,  std::vector<double> &z__);

	    static int sphrec_( std::vector<double> &r__,  std::vector<double> &colat,  std::vector<double> &long__,  std::vector<double> &rectan);

	    static int recsph_( std::vector<double> &rectan,  std::vector<double> &r__,  std::vector<double> &colat,  std::vector<double> &long__);

	    static int sphcyl_( std::vector<double> &radius,  std::vector<double> &colat,  std::vector<double> &slong,  std::vector<double> &r__,  std::vector<double> &long__,  std::vector<double> &z__);

	    static int cylsph_( std::vector<double> &r__,  std::vector<double> &longc,  std::vector<double> &z__,  std::vector<double> &radius,  std::vector<double> &colat,  std::vector<double> &long__);

	    static int sphlat_( std::vector<double> &r__,  std::vector<double> &colat,  std::vector<double> &longs,  std::vector<double> &radius,  std::vector<double> &long__,  std::vector<double> &lat);

	    static int latsph_( std::vector<double> &radius,  std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &rho,  std::vector<double> &colat,  std::vector<double> &longs);

	    static int latcyl_( std::vector<double> &radius,  std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &r__,  std::vector<double> &longc,  std::vector<double> &z__);

	    static int cyllat_( std::vector<double> &r__,  std::vector<double> &longc,  std::vector<double> &z__,  std::vector<double> &radius,  std::vector<double> &long__,  std::vector<double> &lat);
	
	    static int reclat_( std::vector<double> &rectan,  std::vector<double> &radius,  std::vector<double> &long__,  std::vector<double> &lat);

	    static int latrec_( std::vector<double> &radius,  std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &rectan);

	    static int recgeo_( std::vector<double> &rectan,  std::vector<double> &re,  std::vector<double> &f,  std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &alt);

	    static int georec_( std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &alt,  std::vector<double> &re,  std::vector<double> &f,  std::vector<double> &rectan);

	    static int recpgr_( std::string &body,  std::vector<double> &rectan,  std::vector<double> &re,  std::vector<double> &f,  std::vector<double> &lon,  std::vector<double> &lat,  std::vector<double> &alt, int body_len);

	    static int pgrrec_( std::string &body,  std::vector<double> &lon,  std::vector<double> &lat,  std::vector<double> &alt,  std::vector<double> &re,  std::vector<double> &f,  std::vector<double> &rectan, int body_len);

	    static int xfmsta_( std::vector<double> &istate,  std::string &icosys,  std::string &ocosys,  std::string &body,  std::vector<double> &ostate, int icosys_len, int ocosys_len, int body_len);
	
	    static int drdlat_( std::vector<double> &r__,  std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &jacobi);

	    static int dlatdr_( std::vector<double> &x,  std::vector<double> &y,  std::vector<double> &z__,  std::vector<double> &jacobi);

	    static int drdpgr_( std::string &body,  std::vector<double> &lon,  std::vector<double> &lat,  std::vector<double> &alt,  std::vector<double> &re,  std::vector<double> &f,  std::vector<double> &jacobi, int body_len);

	    static int dpgrdr_( std::string &body,  std::vector<double> &x,  std::vector<double> &y,  std::vector<double> &z__,  std::vector<double> &re,  std::vector<double> &f,  std::vector<double> &jacobi, int body_len);

	    static int drdgeo_( std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &alt,  std::vector<double> &re,  std::vector<double> &f,  std::vector<double> &jacobi);

	    static int dgeodr_( std::vector<double> &x,  std::vector<double> &y,  std::vector<double> &z__,  std::vector<double> &re,  std::vector<double> &f,  std::vector<double> &jacobi);

	    static int drdcyl_( std::vector<double> &r__,  std::vector<double> &long__,  std::vector<double> &z__,  std::vector<double> &jacobi);

	    static int dcyldr_( std::vector<double> &x,  std::vector<double> &y,  std::vector<double> &z__,  std::vector<double> &jacobi);

	    static int drdsph_( std::vector<double> &r__,  std::vector<double> &colat,  std::vector<double> &long__,  std::vector<double> &jacobi);

	    static int dsphdr_( std::vector<double> &x,  std::vector<double> &y,  std::vector<double> &z__,  std::vector<double> &jacobi);

		};

    }
}

#endif /* SMARTASTRO_CONVERSION_COORDINATES_H */

