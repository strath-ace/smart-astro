/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2017 University of Strathclyde and Authors ------
------------------- Author: Francesco Torre --------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
*/

#ifndef SMARTASTRO_SMART_MATRIX_H
#define SMARTASTRO_SMART_MATRIX_H

#include "../exception.h"
#include "LinearAlgebra/Eigen/Core"


namespace smartastro
{
    namespace astrocore
    {
        namespace matrix_types
        {
            enum MATRIX_TYPE
            {
                ROT_MATRIX_2D,
                ROT_MATRIX_3D,
                COV_MATRIX,
                MAX_MATRIX_TYPE
            };

            const std::vector<std::string> matrix_type_str
            {
                "ROT_MATRIX_2D",
                "ROT_MATRIX_3D",
                "COV_MATRIX",
                "MAX_MATRIX_TYPE"
            };
        }
        
        template < class T >
        class smart_matrix
        {
        protected:
            // Matrix type
            matrix_types::MATRIX_TYPE m_type;

            // Time associated to the matrix
            double m_time;

            // Values of the matrix
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> m_value;


        public:

            /**
             *
             */
            smart_matrix();


            /**
             * 
             */
            smart_matrix(const matrix_types::MATRIX_TYPE &type);


            /**
             *
             */
            smart_matrix(const matrix_types::MATRIX_TYPE &type,
                         const double &time,
                         const unsigned int &rows,
                         const unsigned int &cols);


            /**
             *
             */
            smart_matrix(const matrix_types::MATRIX_TYPE &type,
                         const double &time,
                         const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix);


            /**
             *
             */
            ~smart_matrix();


            /**
             *
             */
            matrix_types::MATRIX_TYPE get_type() const;

            
            /**
             * 
             */
            std::string get_type_str();


            /**
             *
             */
            double get_time() const;


            /**
             *
             */
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> get_value() const;


            /**
             * 
             */
            T get_value(const unsigned int &x, const unsigned int &y);


            /**
             *
             */
            void set_new_value(const double &time, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix);


            /**
             * 
             */
            unsigned int rows();


            /**
             * 
             */
            unsigned int cols();

            
            /**
             * 
             */
            // T& operator()(const unsigned int &x, const unsigned int& y);
        };
    }
}

#endif // SMARTASTRO_SMART_MATRIX_H
