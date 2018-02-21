/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2016 University of Strathclyde and Authors ------
------------------- Author: Francesco Torre --------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
------------------- Author: Victor Rodriguez -------------------------
-------------- e-mail: victor.rodriguez@strath.ac.uk -----------------
*/

#ifndef SMARTASTRO_SMART_VECTOR_H
#define SMARTASTRO_SMART_VECTOR_H

#include <vector>
#include <limits>
#include <string>

#include "../config.h"

// Enables the use of polynomials
#ifdef ENABLE_SMARTUQ
    #include "smartuq.h"
#endif

#if defined(__GNUC__)
    #if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 9))
    #include "json/json.hpp"
    using json = nlohmann::json;
    #endif
#endif

namespace smartastro
{
    namespace astrocore
    {
        namespace vector_types
        {
            enum VECTOR_TYPE
            {
                // Scalar types
                MASS,
                GRAV_PARAMETER,
                ORBITAL_PERIOD,
                KEP_ELEM_SMA,
                KEP_ELEM_INCL,
                KEP_ELEM_ECC,
                KEP_ELEM_RAAN,
                KEP_ELEM_ARGP,
                // Vector types
                STATE_POS_VEL,
                KEPLERIAN_ELEMENTS,
                STATE_LIN_DERIV,
                STATE_ATT_VEL,
                STATE_ANG_DERIV,
                SENSOR_DATA,
                FILTER_OUTPUT,
                CONTROL,
                MAX_VECTOR_TYPE
            };

            const std::vector<std::string> vector_type_str
            {
                // Scalar types
                "MASS",
                "GRAV_PARAMETER",
                "ORBITAL_PERIOD",
                "KEP_ELEM_SMA",
                "KEP_ELEM_INCL",
                "KEP_ELEM_ECC",
                "KEP_ELEM_RAAN",
                "KEP_ELEM_ARGP",
                // Vector types
                "STATE_POS_VEL",
                "KEPLERIAN_ELEMENTS",
                "STATE_LIN_DERIV",
                "STATE_ATT_VEL",
                "STATE_ANG_DERIV",
                "SENSOR_DATA",
                "FILTER_OUTPUT",
                "CONTROL",
                "MAX_VECTOR_TYPE"
            };
        }
        
        template < class T >
        class smart_vector
        {
        protected:
            // State type
            vector_types::VECTOR_TYPE m_type;

            // Time associated to the vector
            double m_time;

            // Values of the vector
            std::vector<T> m_value;

        public:

            /**
             *
             */
            smart_vector();


            /**
             *
             */
            smart_vector(const vector_types::VECTOR_TYPE &type, const double &time, const int &length);


            /**
             *
             */
            smart_vector(const vector_types::VECTOR_TYPE &type, const double &time, const std::vector<T> &state);


            /**
             *
             */
            ~smart_vector();


            /**
             *
             */
            vector_types::VECTOR_TYPE get_type() const;


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
            std::vector<T> get_value() const;


            /**
             * 
             */
            unsigned int get_length() const;


            /**
             *
             */
            const std::vector<T>* get_pointer_to_value() const;


            /**
             *
             */
            T get_value(const unsigned int &index) const;


            /**
             *
             */
            void set_new_value(const double &time, const std::vector<T> &vector);
        };

        //JSON methods
        #if defined(__GNUC__)
            #if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 9))
                void to_json(json& j, const smart_vector<double>& val);
            #endif
        #endif
    }
}


#endif //SMARTASTRO_SMART_VECTOR_H
