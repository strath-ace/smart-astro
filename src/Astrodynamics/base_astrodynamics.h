/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------ Copyright (C) 2015 University of Strathclyde and Authors ------
-------------------- Author: Annalisa Riccardi -----------------------
-------------- e-mail: annalisa.riccardi@strath.ac.uk ----------------
-------------------- Author: Romain Serra ----------------------------
-------------- e-mail: romain.serra@strath.ac.uk ---------------------
-------------------- Author: Francesco Torre -------------------------
-------------- e-mail: francesco.torre@strath.ac.uk ------------------
*/

#ifndef SMARTASTRO_BASE_ASTRODYNAMICS_H
#define SMARTASTRO_BASE_ASTRODYNAMICS_H

// include dynamics from smartmath
#include "Dynamics/base_dynamics.h"
#include "smartmath.h"

#include "Astro-Core/smartastro_astrocore.h"
#include "config.h"
#include "exception.h"


#ifdef ENABLE_SMARTUQ
    #include "smartuq.h"
#endif

namespace smartastro
{
    namespace astrodynamics {

        /**
         * @brief The base_dynamics class is an abstract class. Any dynamics added to the toolbox needs to inherit from it and implement the method evaluate()
         *
         * The base_dynamics class is an abstract class. Any dynamical system added to the toolbox need to extend this class and implement the method evaluate.
         * The class has been designed to allow external forces and process noise in the dynamics. These are defined through the classes controller and noise,
         * providable in the contructor. It is possible to define for each dynamics also a scaling factor for lengths and time.
         */
	    template < class T >
        class base_astrodynamics : public smartmath::dynamics::base_dynamics<T>
        {
        protected:
            double m_L_scale;
            double m_T_scale;

            /**
             * @brief sc_number number of spacecraft.
             *
             * sc_number counts the number of spacecraft. It is used to apply the dynamics to multiple spacecraft.
             * The resulting state vector is composed by the stack of the states of the spacecraft and its length
             * must be consistent with the dynamics.
             */
            unsigned int m_sc_number;

            /**
             * @brief state_length indicates the length of the state.
             *
             * state_length stores the characteristic length of the state for each spacecraft.
             * It is meant to be used together with sc_number to propagate multiple spacecraft with the same dynamics.
             * Its value must be specified in the constructor of each dynamics.
             */
            unsigned int m_state_length;

            bool m_comments;


        public:

            /**
             * @brief base_dynamics constructor.
             *
             * In the constructor the name of the dynamics, the thrust and control profiles are initialized as class members.
             * @param name dynamical system name
             * @param L_scale length scaling factor. Default value = 1, meters (no scaling).
             * @param T_scale time scaling factor. Default value = 1, seconds (no scaling)
             * @param noise_gen process noise generator
             * @param controller controller
             */
            base_astrodynamics(const std::string &name,
                               const double &L_scale = 1.0,
                               const double &T_scale = 1.0
                              );
            

            virtual ~base_astrodynamics();


            /**
             * @brief Returns the number of spacecraft simultaneously propagated.
             * @param[out] sc_number number of spacecraft
             * @return
             */
            unsigned int get_sc_number() const;


            /**
             * @brief Defines the number of spacecraft simultaneously propagated.
             * @param[in] sc_number number of spacecraft
             * @return
             */
            void set_sc_number(const unsigned int &sc_number);


            /**
             * @brief Returns the length of the state of each spacecraft.
             * @param[out] state_length length of the state of each spacecraft.
             * @return
             */
            unsigned int get_state_length() const;

            /**
             * 
             */
            bool get_comments();


            /**
             * 
             */
            void enable_comments();


            /**
             * 
             */
            void disable_comments();
        };

    }
}

#endif // SMARTASTRO_BASE_ASTRODYNAMICS_H
