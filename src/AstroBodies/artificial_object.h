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
				Artificial_Object(std::string givenName);

				Artificial_Object(int givenId);

				virtual ~Artificial_Object() = default;

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

				void gfrfov(std::string &inst,  std::vector<double> &raydir, std::string &rframe, std::string &abcorr, double &step, SpiceCell &cnfine, SpiceCell &result);

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

				void gftfov(std::string &inst, std::string &target, std::string &tshape, std::string &tframe, std::string &abcorr, double &step, SpiceCell &cnfine, SpiceCell &result);

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

				void fovray(std::string &inst,  std::vector<double> &raydir, std::string &rframe, std::string &abcorr, double &et,  int &visibl);

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

				void fovtrg( std::string &inst, std::string &target, std::string &tshape, std::string &tframe, std::string &abcorr, double &et, int &visibl);

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

				/*void getfov(int &instid, int &room, int &shapelen, int &framelen, const std::string &shape, const std::string &frame,  std::vector<double> &bsight,  std::vector<int> &n,  std::vector<std::vector<double>> &bounds);

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

				/*int gffove_(std::string &inst, std::string &tshape,  std::vector<double> &raydir,  std::string &target, std::string &tframe, std::string &abcorr,   std::vector<double> &tol, int udstep, int udrefn, std::vector<int> &rpt, int udrepi, int udrepu, int udrepf, std::vector<int> &bail, L_fp udbail, std::vector<double> &cnfine,  std::vector<double> &result, int inst_len, int tshape_len, int target_len, int tframe_len, int abcorr_len, int obsrvr_len);*/
		};
	}
}

#endif // SMARTASTRO_ARTIFICIAL_OBJECT_H
