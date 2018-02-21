/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2016 University of Strathclyde and Authors ------
-------------------- Author: Francesco Torre -------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
*/

#ifndef SMARTASTRO_UTILS_H
#define SMARTASTRO_UTILS_H

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <cstdio>

#include "constants.h"
#include "LinearAlgebra/Eigen/Core"
#include "LinearAlgebra/Eigen/QR"

namespace smartastro
{
    class utilities
    {
        public:
        
        template < class T >
        static bool import_tetra_mesh(const std::string &file_name,
                                      std::vector< std::vector<T> > &mesh_vertices,
                                      std::vector< std::vector<unsigned int> > &mesh_faces);


        /**
         * @brief vectord2string This method "stringifies" a vector of doubles
         * @param vec vector of doubles
         * @return string containing each element of the vector separated by ' '
         */
        static std::string vectord2string(const std::vector<double> vec, int precision = 8);


        template < class T >
        /**
         * @brief scale_pos_vel_vector Scale a position-velocity vector,
         * in which the first components of the vector represent the
         * position and the [4-6] components represent the velocity. The rest
         * of elements in the vector are not altered.
         * @param v Vector to scale
         * @param pos_scale Scale factor for the position
         * @param vel_scale Scale factor for the velocity
         * @param v_scaled Outut argument. Scaled vector
         */
        static void scale_pos_vel_vector(const std::vector<T> &v,
                                         const double &pos_scale,
                                         const double &vel_scale,
                                         std::vector<T> &v_scaled)
        {
            v_scaled.clear();

            v_scaled.push_back(v[0] * pos_scale);
            v_scaled.push_back(v[1] * pos_scale);
            v_scaled.push_back(v[2] * pos_scale);
            v_scaled.push_back(v[3] * vel_scale);
            v_scaled.push_back(v[4] * vel_scale);
            v_scaled.push_back(v[5] * vel_scale);

            //Maintain the possible values of the state vector
            for (size_t i = 6; i < v.size(); i++)
                v_scaled.push_back(v[i]);
        }

        /**
         * @brief month_day Computes the day and month given a year day and a year
         * @param year positive integer (> 1583)
         * @param yearday day in year [1-365]
         * @param month Output argument containing the month [1-12]
         * @param monthday Output argument containig the day in the month (not higher than 31)
         */
        static void month_day(unsigned int year, int yearday, int &month, int &monthday)
        {
            char daytab[2][13] = {
             {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
             {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
            };

            int i, leap;
            ((leap = year%4 == 0) && (year%100 != 0))|| (year%400 == 0);

            for (i = 1; yearday > daytab[leap][i]; i++)
                yearday -= daytab[leap][i];

            monthday = yearday;
            month = i;
        }

        /**
         * @brief to_positive_angle Transform an angle, in degrees, to the range [0,360)
         * @param angle angle in degrees
         * @return transformed angle in degrees
         */
        static double to_positive_angle(double angle, const bool &in_radians = false)
        {
            double cycle = (in_radians ? smartastro::constants::pi2 : 360.0);

            angle = fmod(angle, cycle);
            while(angle < 0) { //pretty sure this comparison is valid for doubles and floats
                angle += cycle;
            }

            return angle;
        }

        static std::string exec_popen(const std::string &cmd);


        template <typename I>
        /**
         * @brief random_element Extract a random element from a collection
         * @param begin
         * @param end
         * @return
         */
        static I random_element(I begin, I end);
    };
}

#endif // SMARTASTRO_UTILS_H
