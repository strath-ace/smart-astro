/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/


#ifndef SMARTASTRO_BASE_PROPAGATOR_H
#define SMARTASTRO_BASE_PROPAGATOR_H

#include <vector>
#include <math.h>
#include <limits>
#include "../exception.h"

namespace smartastro
{
    namespace propagator {
        /**
         * @brief The base_propagator class is an abstract class. Any propagator added to the toolbox needs to inherit from it and implement the method propagate()
         *
         * @author Annalisa Riccardi
         */
        template < class T >
        class base_propagator
        {

        public:
            /**
             * @brief base_propagator constructor
             * @param name name of the propagator
             */
            base_propagator(const std::string &name)
            {
                m_name = name;
            }

            /**
             * @brief ~base_propagator deconstructor
             */
            virtual ~base_propagator()
            {
                // NOTHING YET
            }

            /**
             * @brief propagate virtual method. Propagation routine.
             * @param[in] t0 initial time
             * @param[in] tend final time
             * @param[in] initial_state initial state vector
             * @param[out] final_state final state vector
             * @return exit flag: 0 = success
             */
            virtual int propagate(const double &t0, const double&tend, const std::vector<T> &initial_state, std::vector<T> &final_state) const = 0;
            
            /**
             * @brief propagate propagate the trajectory from t0 to tend saving intermediate results.
             * @param[in] t0 initial time
             * @param[in] tend final time
             * @param[in] time_step constant step to save results
             * @param[in] initial_state initial state vector
             * @param[out] final_state final state matrix. The state vector at time t0, t0+time_step, t0+2*time_tep,...,tend are saved as rows
             * @return exit flag: 0 = success, n otherwise (propagator dependant)
             */
            int propagate_saving_intermediate_results(const double &t0,
                                                      const double&tend,
                                                      const double &time_step,
                                                      const std::vector<T> &initial_state,
                                                      std::vector<std::vector<T> > &final_state) const
            {
                //Sanity checks
                if(time_step<=0)
                    smartastro_throw("Base Propagator: the timestep for propagation must be a positive number");

                double ti = t0;
                double tf = ti + time_step;
                std::vector<T> xi = initial_state;
                std::vector<T> xf(final_state.size());
                int exit_flag = 0;
                bool last_step = false;

                if(tf>tend)
                    smartastro_throw("Base Propagator: the timestep for propagation is larger than the final propagation time");

                //clear final states matrix
                final_state.clear();

                while(tf<=tend){
                    std::cout << "base_propagator: propagating from ti=" << ti << "to tf=" << tf << std::endl << "[t0=" << t0 << "], tend = " << tend << "]";
                    exit_flag = propagate(ti,tf,xi,xf);
                    if(exit_flag == 0){ // check no error in the propagations
                        ti = tf;
                        tf = ti + time_step;
                        if(tf>tend && !last_step) {
                            tf = tend;
                            last_step = true;
                        }
                        xi = xf;
                        final_state.push_back(xf);
                    }
                    else break;
                }

                return exit_flag;
            }

            std::string get_name() const
            {
                return m_name;
            }

        protected:
            std::string m_name;
        };

    }
}

#endif // SMARTASTRO_BASE_PROPAGATOR_H
