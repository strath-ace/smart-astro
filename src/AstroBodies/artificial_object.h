/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#ifndef SMARTASTO_MAN_MADE_OBJECT_H
#define SMARTASTO_MAN_MADE_OBJECT_H

#include "astro_body.h"

namespace smartastro
{
	namespace astrobodies
	{
		class Artificial_Object : public Astro_Body
		{
			public:
				Artificial_Object( std::string givenName, int givenId);

				~Artificial_Object();

				/**
				* @brief Determine time intervals when a specified ray intersects the 
   				* space bounded by the field-of-view (FOV) of a specified 
   				* instrument.
				*
				* @param[in] inst Name of the instrument.
				* @param[in] raydir Ray's direction vector.
				* @param[in] rframe Reference frame of ray's direction vector.
				* @param[in] abcorr Aberration correction flag.
				* @param[in] step Step size in seconds for finding FOV events.
				* @param[in/out] cnfine SPICE window to which the search is restricted.
				* @param[out] result SPICE window containing results.
				* @return void
				*
				*/

				void gfrfov(const std::string &inst, const std::vector<double> &raydir, const std::string &rframe, const std::string &abcorr, const double &step, SpiceCell &cnfine, SpiceCell &result);

				/**
				* @brief Determine time intervals when a specified ephemeris object 
   				* intersects the space bounded by the field-of-view (FOV) of a 
   				* specified instrument.
				*
				* @param[in] inst Name of the instrument.
				* @param[in] target Name of the target body.
				* @param[in] tshape Type of shape model used for target body.
				* @param[in] tframe Body-fixed, body-centered frame for target body.
				* @param[in] abcorr Aberration correction flag.
				* @param[in] step Step size in seconds for finding FOV events.
				* @param[in/out] cnfine SPICE window to which the search is restricted.
				* @param[out] result SPICE window containing results.
				* @return void
				*
				*/

				void gftfov(const std::string &inst, const std::string &target, const std::string &tshape, const std::string &tframe, const std::string &abcorr, const double &step, SpiceCell &cnfine, SpiceCell &result);

				/**
				* @brief Determine time intervals when a specified ephemeris object 
   				* intersects the space bounded by the field-of-view (FOV) of a 
   				* specified instrument.
				*
				* @param[in] inst Name of the instrument.
				* @param[in] raydir Ray's direction vector.
				* @param[in] rframe Reference frame of ray's direction vector.
				* @param[in] abcorr Aberration correction flag.
				* @param[in] et Time of the observation (seconds past J2000).
				* @param[out] visible Visibility flag (SPICETRUE/SPICEFALSE).
				* @return void
				*
				*/

				void fovray(const std::string &inst, const std::vector<double> &raydir, const std::string &rframe, const std::string &abcorr, const double &et, int &visibl);

				/**
				* @brief Determine if a specified ephemeris object is within the
   				* field-of-view (FOV) of a specified instrument at a given time.
				*
				* @param[in] inst Name of the instrument.
				* @param[in] target Name of the target body.
				* @param[in] tshape Type of shape model used for target body.
				* @param[in] tframe Body-fixed, body-centered frame for target body.
				* @param[in] abcorr Aberration correction flag.
				* @param[in] et Time of the observation (seconds past J2000).
				* @param[out] visible Visibility flag (SPICETRUE/SPICEFALSE).
				* @return void
				*
				*/

				void fovtrg(const std::string &inst, const std::string &target, const std::string &tshape, const std::string &tframe, const std::string &abcorr, const double &et, int &visibl);

				/**
				* @brief Return the field-of-view (FOV) parameters for a specified 
   				* instrument. The instrument is specified by its NAIF ID code.
				*
				* @param[in] instid NAIF ID of an instrument.
				* @param[in] room Maximum number of vectors that can be returned.
				* @param[in] shapelen Space available in the string `shape'.
				* @param[in] framelen Space available in the string `frame'.
				* @param[out] shape Instrument FOV shape. 
				* @param[out] frame Name of the frame in which FOV vectors are defined.
				* @param[out] bsight Boresight vector. 
				* @param[out] n Number of boundary vectors returned. 
				* @param[out] bounds FOV boundary vectors. 
				* @return void
				*
				*/

				/*void getfov(const int &instid, const int &room, const int &shapelen, const int &framelen, std::string &shape, std::string &frame,  std::vector<double> &bsight,  std::vector<int> &n,  std::vector<vector<double>> &bounds);

				/**
				* @brief Determine time intervals when a specified target body or ray
   				* intersects the space bounded by the field-of-view (FOV) of a
   				* specified instrument. Report progress and handle interrupts if so
   				* commanded.
				*
				* @param[in] inst Name of the instrument.
				* @param[in] tshape Type of shape model used for target body.
				* @param[in] raydir Ray's direction vector.
				* @param[in] target Name of the target body.
				* @param[in] tframe Body-fixed, body-centered frame for target body.
				* @param[in] abcorr Aberration correction flag.
				* @param[in] tol Convergence tolerance in seconds. 
				* @param[in] udstep Name of the routine returns a time step.
				* @param[in] udrefn Name of the routine that computes a refined time.
				* @param[in] rpt Progress report flag.
				* @param[in] udrepi Function that initializes progress reporting.
				* @param[in] udrepu Function that updates the progress report. 
				* @param[in] udrepf Function that finalizes progress reporting.
				* @param[in] bail Logical indicating program interrupt monitoring.
				* @param[in] udbail Name of a routine that signals a program interrupt. 
				* @param[in/out] cnfine SPICE window to which the search is restricted.
				* @param[out] result SPICE window containing results.  
				* @return int
				*
				*/

				/*int gffove_(const std::string &inst, const std::string &tshape, const std::vector<double> &raydir, const std::string &target, const std::string &tframe, const std::string &abcorr, const std::vector<double> &tol, const int udstep, const int udrefn, const std::vector<int> &rpt, const int udrepi, const int udrepu, const int udrepf, const std::vector<int> &bail, const L_fp udbail, const std::vector<double> &cnfine, std::vector<double> &result);*/
		};
	}
}

#endif // SMARTASTRO_ARTIFICIAL_OBJECT_H
