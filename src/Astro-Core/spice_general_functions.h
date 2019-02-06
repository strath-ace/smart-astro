/* Thi Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#ifndef SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H
#define SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H

extern "C" {
#include <cspice/SpiceUsr.h>
#include <cspice/SpiceZfc.h>
}

#include <vector>

namespace smartastro
{
	namespace astrocore
	{

		class spice_general_functions
		{
			public:

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

				static std::vector<int> intArrayToVector(int array[]);

				static int* intVectorToArray(std::vector<int> vector);

				static std::vector<char> charArrayToVector(char array[]);

				static char* charVectorToArray(std::vector<char> vector);

				static std::vector<double> doubleArrayToVector(double array[]);

				static double* doubleVectorToArray(std::vector<double> vector);

				static char* stringToCharArray(std::string string);

				static int subslr_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec, int method_len, int target_len, int fixref_len, int abcorr_len, int obsrvr_len);

				static int subpnt_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec, int method_len, int target_len, int fixref_len, int abcorr_len, int obsrvr_len);

				static double lspcn_( std::string &body,  std::vector<double> &et,  std::string &abcorr, int body_len, int abcorr_len);

				static int furnsh_( std::string &file, int file_len);

				static int unload_( std::string &file, int file_len);

		};
	}
}

#endif /* SMARTASTRO_SPICE_GENERAL_FUNCTIONS_H */
