/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Propagators/perturbation_propagator.h"

using namespace smartastro;
using namespace propagator;

perturbation_propagator::perturbation_propagator(const std::vector<bool> &flags, const int &zonal): base_propagator("perturbation propagator"), m_flags(flags), m_zonal(zonal)
{
    /* Sanity checks */
    if(zonal < 0)
        smartastro_throw("PERTURBATION_PROPAGATOR: degree of zonal harmonics must be positive");
    if(zonal > 5)
        smartastro_throw("PERTURBATION_PROPAGATOR: degree of zonal harmonics cannot be higher than 5");
    if(m_flags.size() < 1)
        smartastro_throw("PERTURBATION_PROPAGATOR: there must be at least one flag for perturbations");
}

perturbation_propagator::~perturbation_propagator(){
}


int perturbation_propagator::propagate(const double &L0, const double &Lf, const std::vector<double> &initial_state, std::vector<double> &final_state) const{
    
    /* Sanity checks */
    if(initial_state.size() != 6)
        smartastro_throw("PROPAGATE: initial state needs to have 6 dimensions"); 
    if(initial_state[5] != L0)
        smartastro_throw("PROPAGATE: initial true longitude does not equal value in initial state");     
    if(initial_state[1] * initial_state[1] + initial_state[2] * initial_state[2] >= 1.0)
        smartastro_throw("PROPAGATE: cannot handle non-closed orbits (ecc > 1)"); 
    if(m_timing.size() != 2)
        smartastro_throw("PROPAGATE: propagator must be initialized by associating a Julian date to the initial true longitude");          
    if(m_timing[0] != L0)
        smartastro_throw("PROPAGATE: initial true longitude does not equal value stored in vector for initial condition on the date");                

    /* Keplerian motion */
    final_state = initial_state;
    final_state[5] = Lf;

    /* Computing orbital perturbations */
    std::vector<double> delta(6), delta1(6, 0.0), delta2(6, 0.0);

    if(m_zonal >= 2) // zonal harmonics
        zonal_harmonics(L0, Lf, initial_state, delta1);

    if(m_flags[0])
        lunisolar(L0, Lf, initial_state, delta2);

    for(int i = 0; i < 6; i++)
        delta[i] = delta1[i] + delta2[i];

    /* Computing final state */
    final_state[1] += delta[2];
    final_state[2] += delta[1];
    final_state[3] += delta[4];
    final_state[4] += delta[3];
    double af = initial_state[0] / (1.0 - (initial_state[1] * initial_state[1] + initial_state[2] * initial_state[2])) + delta[0]; // final semi-major axis
    final_state[0] = af * (1.0 - final_state[1] * final_state[1] - final_state[2] * final_state[2]);  

    return 0;
}

int perturbation_propagator::set_timing(const double &L0, const double &jd0){

    m_timing.clear();
    m_timing.push_back(L0);
    m_timing.push_back(jd0);

    return 0;
}

int perturbation_propagator::get_timing(double &L0, double &jd0) const{

    L0 = m_timing[0];
    jd0 = m_timing[1];

    return 0;
}

int perturbation_propagator::zonal_harmonics(const double &L0, const double&Lf, const std::vector<double> &initial_state, std::vector<double> &delta) const{

    /* Gravitational constants */
    double R = 1.0;        
    // double muadim = 1.0;

    /* getting initial conditions */
    double p0  = initial_state[0];
    double P10 = initial_state[2];
    double P20 = initial_state[1];
    double Q10 = initial_state[4];
    double Q20 = initial_state[3];
    
    // Parameter B (equation [3] in Reference Zuiani)
    double B2 = 1.0 - (P10 * P10 + P20 * P20);
    // Semi-major axis (initial value)
    double a0 = p0 / B2;
    
    // Angular momentum (initial value)
    // double h0 = sqrt(muadim * p0);
    
     // Computation of the initial value of Omom = (Omega + omega)
     // Omom is the sum of RAAN and perigee argument.
     // If eccentricity is zero P10 and P20 are zero so Omom can not be
     // computed
    //double Omom = 0.0;
    //if(e_flag == true)
    //    Omom = atan2(P10,P20);

    double G = 1.0 + Q10 * Q10 + Q20 * Q20;

    /* Coefficients for J2 */
    double k1 = 3.0 * m_J2 * R * R / (pow(B2, 3) * G * G * a0);
    double k2 = k1 * B2 / a0;
    double k3 = k2 * G;

    /* Coefficients for higher zonal harmonics */
    double kJ3 = m_J3 * pow(R / a0, 3);
    double kJ4 = m_J4 * pow(R / a0, 4);
    double kJ5 = m_J5 * pow(R / a0, 5);

    // The following line compute e cos(theta0) where theta = L - (Omega + omega)
    //double e0_c_th0 = P10 * s_L0 + P20 * c_L0;
    // The following line compute e sin(theta0) where theta = L - (Omega + omega)
    //double e0_s_th0 = P20 * s_L0 - P10 * c_L0;
    
    // Initial true anomaly theta0
    // double th0 = atan2(e0_s_th0, e0_c_th0);   

    // Since (Omega + omega) remains always constant
    // L - theta = L0 - theta0
    // and therefore, the true anomaly theta is:    
    // double th = th0 + (Lf - L0); 

    /* Initializing variations of elements */
    double a1J = 0.0, P11J =0.0, P21J = 0.0, Q11J = 0.0, Q21J = 0.0;   

    /* Computing J2-related integrals */
    std::vector<double> I = integralsJ2(L0, Lf);

    /* Adding J2 contribution */
    a1J += k1 * ( (8.0*Q20*pow(P20, 3)*Q10 - 12.0*P10*pow(P20, 2)*pow(Q10, 2)) * I[17] +                                                      
        (- 24.0*pow(P10, 2)*P20*pow(Q10, 2) + 48.0*P10*pow(P20, 2)*Q10*Q20 + 20.0*pow(P20, 3)*pow(Q10, 2) - 8.0*pow(P20, 3)*pow(Q20, 2)) * I[18] +                       
        (- 24.0*P10*P20*pow(Q10, 2) + 24.0*pow(P20, 2)*Q10*Q20) * I[12] +
        (- 12.0*pow(P10, 3)*pow(Q10, 2) + 72.0*pow(P10, 2)*P20*Q10*Q20 + 48.0*P10*pow(P20, 2)*pow(Q10, 2) - 36.0*P10*pow(P20, 2)*pow(Q20, 2) - 32.0*pow(P20, 3)*Q10*Q20) * I[19] +
        (- 24.0*pow(P10, 2)*pow(Q10, 2) + 96.0*P10*P20*Q10*Q20 + 48.0*pow(P20, 2)*pow(Q10, 2) - 24.0*pow(P20, 2)*pow(Q20, 2)) * I[13] +
        (P10*pow(G, 2)*pow(P20, 2) + 24.0*Q20*P20*Q10 - 12.0*P10*pow(Q10, 2)) * I[8] +
        (32.0*pow(P10, 3)*Q10*Q20 + 36.0*pow(P10, 2)*P20*pow(Q10, 2) - 48.0*pow(P10, 2)*P20*pow(Q20, 2) - 72.0*P10*pow(P20, 2)*Q10*Q20 + 12.0*pow(P20, 3)*pow(Q20, 2)) * I[20] +
        (72.0*pow(P10, 2)*Q10*Q20 + 72.0*P10*P20*pow(Q10, 2) - 72.0*P10*P20*pow(Q20, 2) - 72.0*pow(P20, 2)*Q10*Q20) * I[14] +
        (2.0*pow(G, 2)*pow(P10, 2)*P20 - pow(G, 2)*pow(P20, 3) + 48.0*P10*Q10*Q20 + 36.0*P20*pow(Q10, 2) - 24.0*P20*pow(Q20, 2)) * I[9] +
        (2.0*P10*P20*pow(G, 2) + 8.0*Q10*Q20) * I[5] +
        (8.0*pow(P10, 3)*pow(Q10, 2) - 20.0*pow(P10, 3)*pow(Q20, 2) - 48.0*pow(P10, 2)*P20*Q10*Q20 + 24.0*P10*pow(P20, 2)*pow(Q20, 2)) * I[21] +
        (24.0*pow(P10, 2)*pow(Q10, 2) - 48.0*pow(P10, 2)*pow(Q20, 2) - 96.0*P10*P20*Q10*Q20 + 24.0*pow(P20, 2)*pow(Q20, 2)) * I[15] +
        (pow(G, 2)*pow(P10, 3) - 2.0*pow(G, 2)*P10*pow(P20, 2) + 24.0*P10*pow(Q10, 2) - 36.0*P10*pow(Q20, 2) - 48.0*P20*Q10*Q20) * I[10] +
        (2.0*pow(G, 2)*pow(P10, 2) - 2.0*pow(G, 2)*pow(P20, 2) + 8.0*pow(Q10, 2) - 8.0*pow(Q20, 2)) * I[6] +
        pow(G, 2)*P10 * I[3] +
        (- 8.0*Q10*pow(P10, 3)*Q20 + 12.0*P20*pow(P10, 2)*pow(Q20, 2)) * I[22] +
        (- 24.0*Q10*pow(P10, 2)*Q20 + 24.0*P20*P10*pow(Q20, 2)) * I[16] +
        (- P20*pow(G, 2)*pow(P10, 2) - 24.0*Q10*P10*Q20 + 12.0*P20*pow(Q20, 2)) * I[11] +
        (- 2.0*P10*P20*pow(G, 2) - 8.0*Q10*Q20) * I[7] +
        (-pow(G, 2)*P20) * I[4]);

    P11J += k2*( (-6.0*pow(P20, 2)*pow(Q10, 2)) * I[17] +
        (16.0*pow(P20, 2)*Q10*Q20 - 12.0*P10*P20*pow(Q10, 2)) * I[18] -
        12.0*P20*pow(Q10, 2) * I[12] +
        (- 6.0*pow(P10, 2)*pow(Q10, 2) + 32.0*P10*P20*Q10*Q20 + 4.0*pow(P20, 2)*pow(Q10, 2) - 10.0*pow(P20, 2)*pow(Q20, 2)) * I[19] +
        (- 12.0*P10*pow(Q10, 2) + 36.0*P20*Q20*Q10) * I[13] +
        ((pow(G, 2)*pow(P20, 2))/2.0 + 2.0*pow(P20, 2)*pow(Q10, 4) + 2.0*pow(P20, 2)*pow(Q10, 2)*pow(Q20, 2) - 2.0*pow(P20, 2)*pow(Q10, 2) + 4.0*P10*P20*Q10*Q20 - 6.0*pow(Q10, 2)) * I[8] +
        (16.0*pow(P10, 2)*Q10*Q20 + 8.0*P10*P20*pow(Q10, 2) - 20.0*P10*P20*pow(Q20, 2) - 4.0*pow(P20, 2)*Q10*Q20) * I[20] +
        (12.0*P20*pow(Q10, 2) + 36.0*P10*Q10*Q20 - 24.0*P20*pow(Q20, 2)) * I[14] +
        (pow(G, 2)*P10*P20 + 4.0*pow(P10, 2)*Q10*Q20 + 2.0*P10*P20*pow(Q10, 4) + 2.0*P10*P20*pow(Q10, 2)*pow(Q20, 2) + 2.0*P10*P20*pow(Q10, 2) - 4.0*P10*P20*pow(Q20, 2) - 4.0*pow(P20, 2)*pow(Q10, 3)*Q20 - 4.0*pow(P20, 2)*Q10*pow(Q20, 3) + 4.0*pow(P20, 2)*Q10*Q20 + 20.0*Q10*Q20) * I[9] +
        (P20*pow(G, 2) + 2.0*P20*pow(Q10, 4) + 2.0*P20*pow(Q10, 2)*pow(Q20, 2) - 2.0*P20*pow(Q10, 2) + 4.0*P10*Q10*Q20) * I[5] +
        (4.0*pow(P10, 2)*pow(Q10, 2) - 10.0*pow(P10, 2)*pow(Q20, 2) - 8.0*P20*P10*Q10*Q20) * I[21] +
        (12.0*P10*pow(Q10, 2) - 12.0*P20*Q10*Q20 - 24.0*P10*pow(Q20, 2)) * I[15] +
        ((pow(G, 2)*pow(P10, 2))/2.0 + 4.0*pow(P10, 2)*pow(Q10, 2) - 4.0*pow(P10, 2)*pow(Q20, 2) - 4.0*P10*P20*pow(Q10, 3)*Q20 - 4.0*P10*P20*Q10*pow(Q20, 3) + 2.0*pow(P20, 2)*pow(Q10, 2)*pow(Q20, 2) + 2.0*pow(P20, 2)*pow(Q20, 4) - 2.0*pow(P20, 2)*pow(Q20, 2) + 8.0*pow(Q10, 2) - 14.0*pow(Q20, 2)) * I[10] +
        (P10*pow(G, 2) - 4.0*P20*pow(Q10, 3)*Q20 + 4.0*P10*pow(Q10, 2) - 4.0*P20*Q10*pow(Q20, 3) + 4.0*P20*Q10*Q20 - 4.0*P10*pow(Q20, 2)) * I[6] +
         pow(G, 2)/2.0 * I[3] +
        (-4.0*pow(P10, 2)*Q10*Q20) * I[22] +
        (-12.0*P10*Q10*Q20) * I[16] +
        (- 4.0*pow(P10, 2)*Q10*Q20 + 2.0*P20*P10*pow(Q10, 2)*pow(Q20, 2) + 2.0*P20*P10*pow(Q20, 4) - 2.0*P20*P10*pow(Q20, 2) - 8.0*Q10*Q20) * I[11] +
        (2.0*P20*pow(Q10, 2)*pow(Q20, 2) - 4.0*P10*Q10*Q20 + 2.0*P20*pow(Q20, 4) - 2.0*P20*pow(Q20, 2)) * I[7]);
    
    P21J += k2*((4.0*pow(P20, 2)*Q10*Q20) * I[17] +
        (10.0*pow(P20, 2)*pow(Q10, 2) - 4.0*pow(P20, 2)*pow(Q20, 2) + 8.0*P10*P20*Q10*Q20) * I[18] +
        12.0*P20*Q10*Q20 * I[12] +
        (4.0*pow(P10, 2)*Q10*Q20 + 20.0*P10*P20*pow(Q10, 2) - 8.0*P10*P20*pow(Q20, 2) - 16.0*pow(P20, 2)*Q10*Q20) * I[19] +
        (24*P20*pow(Q10, 2) + 12*P10*Q10*Q20 - 12*P20*pow(Q20, 2)) * I[13] +
        (4.0*pow(P20, 2)*Q10*Q20 - 2.0*P10*P20*pow(Q10, 4) - 2.0*P10*P20*pow(Q10, 2)*pow(Q20, 2) + 2.0*P10*P20*pow(Q10, 2) + 8.0*Q10*Q20) * I[8] +
        (10.0*pow(P10, 2)*pow(Q10, 2) - 4.0*pow(P10, 2)*pow(Q20, 2) - 32.0*P10*P20*Q10*Q20 + 6.0*pow(P20, 2)*pow(Q20, 2)) * I[20] +
        (24.0*P10*pow(Q10, 2) - 36.0*P20*Q10*Q20 - 12.0*P10*pow(Q20, 2)) * I[14] +
        (- 2.0*pow(P10, 2)*pow(Q10, 4) - 2.0*pow(P10, 2)*pow(Q10, 2)*pow(Q20, 2) + 2.0*pow(P10, 2)*pow(Q10, 2) + 4.0*P10*P20*pow(Q10, 3)*Q20 + 4.0*P10*P20*Q10*pow(Q20, 3) + 4.0*pow(P20, 2)*pow(Q10, 2) - 4.0*pow(P20, 2)*pow(Q20, 2) - (pow(G, 2)*pow(P20, 2))/2.0 + 14.0*pow(Q10, 2) - 8.0*pow(Q20, 2)) * I[9] +
        (- 2.0*P10*pow(Q10, 4) - 2.0*P10*pow(Q10, 2)*pow(Q20, 2) + 2.0*P10*pow(Q10, 2) + 4.0*P20*Q10*Q20) * I[5] +
        (- 16.0*Q10*pow(P10, 2)*Q20 + 12.0*P20*P10*pow(Q20, 2)) * I[21] +
        (12.0*P20*pow(Q20, 2) - 36.0*P10*Q10*Q20) * I[15] +
        (- pow(G, 2)*P10*P20 + 4.0*pow(P10, 2)*pow(Q10, 3)*Q20 + 4.0*pow(P10, 2)*Q10*pow(Q20, 3) - 4.0*pow(P10, 2)*Q10*Q20 - 2.0*P10*P20*pow(Q10, 2)*pow(Q20, 2) + 4.0*P10*P20*pow(Q10, 2) - 2.0*P10*P20*pow(Q20, 4) - 2.0*P10*P20*pow(Q20, 2) - 4.0*pow(P20, 2)*Q10*Q20 - 20.0*Q10*Q20) * I[10] +
        (- P20*pow(G, 2) + 4.0*P10*pow(Q10, 3)*Q20 + 4.0*P20*pow(Q10, 2) + 4.0*P10*Q10*pow(Q20, 3) - 4.0*P10*Q10*Q20 - 4.0*P20*pow(Q20, 2)) * I[6] +
        (6.0*pow(P10, 2)*pow(Q20, 2)) * I[22] +
        (12.0*P10*pow(Q20, 2)) * I[16] +
        (- 2.0*pow(P10, 2)*pow(Q10, 2)*pow(Q20, 2) - 2.0*pow(P10, 2)*pow(Q20, 4) + 2.0*pow(P10, 2)*pow(Q20, 2) - (pow(G, 2)*pow(P10, 2))/2.0 - 4.0*P20*P10*Q10*Q20 + 6.0*pow(Q20, 2)) * I[11] +
        (- P10*pow(G, 2) - 2.0*P10*pow(Q10, 2)*pow(Q20, 2) - 4.0*P20*Q10*Q20 - 2.0*P10*pow(Q20, 4) + 2.0*P10*pow(Q20, 2)) * I[7] +
        (-pow(G, 2)/2.0) * I[4]);
  
    Q11J += k3*((- P20*pow(Q10, 3) - P20*Q10*pow(Q20, 2) + P20*Q10) * I[9] +
        (- P10*pow(Q10, 3) + P20*pow(Q10, 2)*Q20 - P10*Q10*pow(Q20, 2) + P10*Q10 + P20*pow(Q20, 3) - P20*Q20) * I[10] +
        (- pow(Q10, 3) - Q10*pow(Q20, 2) + Q10) * I[6] +
        (P10*pow(Q10, 2)*Q20 + P10*pow(Q20, 3) - P10*Q20) * I[11] +
        (pow(Q10, 2)*Q20 + pow(Q20, 3) - Q20) * I[7]);

    Q21J += k3*((- P20*pow(Q10, 3) - P20*Q10*pow(Q20, 2) + P20*Q10) * I[8] +
        (- P10*pow(Q10, 3) + P20*pow(Q20, 3) + P10*Q10 - P20*Q20 - P10*Q10*pow(Q20, 2) + P20*pow(Q10, 2)*Q20) * I[9] +
        (Q10 - Q10*pow(Q20, 2) - pow(Q10, 3) ) * I[5] +
        (P10*pow(Q10, 2)*Q20 + P10*pow(Q20, 3) - P10*Q20) * I[10] +
        (pow(Q10, 2)*Q20 + pow(Q20, 3) - Q20) * I[6]);

    if(m_zonal >= 3) // J3 contribution
    {
        /* Computing integrals */
        std::vector<double> IJ3 = integralsJ3(L0, Lf, P10, P20, Q10, Q20);

        // IJ3[0] = IJ3_a1, IJ3[1] = IJ3_a2, IJ3[2] = IJ3_P11, IJ3[3] = IJ3_P12, IJ3[4] = IJ3_P13, IJ3[5] = IJ3_P21, IJ3[6] = IJ3_Q1, IJ3[7] = IJ3_Q2
        a1J += kJ3 * a0 / (pow(B2, 4) * G) * (8.0 * IJ3[0] + 6.0 * IJ3[1]);
        P11J += kJ3 / (pow(B2, 3) * G) * (-4.0 * IJ3[2] + 3.0 * IJ3[3]) - 3.0 * kJ3 / (2.0 * pow(B2, 3) * G) * P20 * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ3[4];       
        P21J += kJ3 / (pow(B2, 3) * G) * IJ3[5] + 3.0 * kJ3 / (2.0 * pow(B2, 3) * G) * P10 * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ3[4];
        Q11J += 3.0 / 4.0 * kJ3 / (pow(B2, 3) * G) * (1.0 + Q10 * Q10 + Q20 * Q20) * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ3[6];
        Q21J += 3.0 / 4.0 * kJ3 / (pow(B2, 3) * G) * (1.0 + Q10 * Q10 + Q20 * Q20) * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ3[7];
    }

    if(m_zonal >= 4) // J4 contribution
    {
        /* Computing integrals */
        std::vector<double> IJ4 = integralsJ4(L0, Lf, P10, P20, Q10, Q20);

        a1J += kJ4 * a0 / pow(B2, 5) * (0.25 * IJ4[0] - 2.0 * IJ4[1] / pow(G, 2));
        P11J += kJ4 / pow(B2, 4) * (-IJ4[2]) - kJ4 / (pow(B2, 4) * pow(G, 2)) * P20 * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ4[3];
        P21J += kJ4 / pow(B2, 4) * (IJ4[4]) + kJ4 / (pow(B2, 4) * pow(G, 2)) * P10 * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ4[3];
        Q11J += -0.5 * kJ4 / (pow(B2, 4) * pow(G, 2)) * (1.0 + (Q10 * Q10 + Q20 * Q20)) * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ4[5];
        Q21J += -0.5 * kJ4 / (pow(B2, 4) * pow(G, 2)) * (1.0 + (Q10 * Q10 + Q20 * Q20)) * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ4[6];
    }

    if(m_zonal >= 5) // J5 contribution
    {
        /* Computing integrals */
        std::vector<double> IJ5 = integralsJ5(L0, Lf, P10, P20, Q10, Q20);

        a1J += kJ5 * a0 / (pow(B2, 6) * G) * IJ5[0];
        P11J += -kJ5 / (2.0 * G * pow(B2, 5)) * IJ5[1] + kJ5 / (8.0 * pow(B2, 5) * G) * P20 * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ5[2];
        P21J += -kJ5 / (2.0 * G * pow(B2, 5)) * IJ5[3] - kJ5 / (8.0 * pow(B2, 5) * G) * P10 * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ5[2];
        Q11J += -1.0 / 16.0 * kJ5 / (pow(B2, 5) * G) * (1.0 + (Q10 * Q10 + Q20 * Q20)) * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ5[4];
        Q21J += -1.0 / 16.0 * kJ5 / (pow(B2, 5) * G) * (1.0 + (Q10 * Q10 + Q20 * Q20)) * (1.0 - (Q10 * Q10 + Q20 * Q20)) * IJ5[5];
    }   
    
    /* Passing variations of elements */
    delta[0] = a1J;
    delta[1] = P11J;
    delta[2] = P21J;
    delta[3] = Q11J;     
    delta[4] = Q21J;
     
    return 0;   
}

int perturbation_propagator::lunisolar(const double &L0, const double &Lf, const std::vector<double> &initial_state, std::vector<double> &delta) const{

    double mu, R_3rd;
    std::vector<double> dir(3), f(3), g(3), w(3), pos(3);
    double mu_sun = constants::mu_sun * pow(constants::T_Earth, 2) / pow(constants::R_earth, 3);
    double mu_moon = constants::mu_moon * pow(constants::T_Earth, 2) / pow(constants::R_earth, 3);
    double L00, jd0, jd;
    get_timing(L00, jd0);
    jd = jd0;

    /* getting initial conditions */
    double p0  = initial_state[0];
    double P10 = initial_state[2];
    double P20 = initial_state[1];
    double Q10 = initial_state[4];
    double Q20 = initial_state[3];

    double e2 = P10 * P10 + P20 * P20;
    double B2 = 1.0 - e2;
    double a0 = p0 / B2;
    double G = 1.0 + Q10 * Q10 + Q20 * Q20;

    double I12 = Lf - L0;
    // if(sqrt(e2) >= 1.0e-5) // enigmatic computations from the Matlab code
    // {   
    //     double pi = constants::pi;
    //     double dL = Lf - L0;
    //     unsigned int n_periods = int(floor(fabs(dL) / (2.0 * pi)));
    //     double L_f = Lf - 2.0 * pi * floor(Lf / (2.0 * pi));
    //     if(L_f > pi)
    //         L_f -= 2.0 * pi;
    //     double L_0 = L0 - 2.0 * pi * floor(L0 / (2.0 * pi));
    //     if(L_0 > pi)
    //         L_0 -= 2.0 * pi;    

    //     if(fabs(P20) < 1.0e-13)
    //     {
    //         double H1 = sqrt(1.0 - P10 * P10);
    //         double H3 = pow(H1, 3);

    //         double M1 = 1.0 + P10 * sin(L_f);
    //         double M10 = 1.0 + P10 * sin(L_0);           
            
    //         double K = 2.0 * atan((P10 + tan(L_f / 2.0)) / H1); 
    //         if(K < 0.0)
    //             K += 2.0 * pi;
    //         double K0 = 2.0 * atan((P10 + tan(L_0 / 2.0)) / H1);   
    //         if(K0 < 0.0)
    //             K0 += 2.0 * pi;                  
    //         double K1 = K - K0;
    //         if(dL > 0.0) 
    //         {
    //             if(K1 < 0.0)
    //                 K1 += 2.0 * pi;
    //             K1 += 2.0 * pi * double(n_periods);
    //         }
    //         else
    //         {
    //             if(K1 > 0.0)
    //                 K1 -= 2.0 * pi;  
    //             K1 -= 2.0 * pi * double(n_periods);
    //         }                 
            
    //         I12 = (K1 + (P10 * H1 * (cos(L_f) / M1 - cos(L_0) / M10))) / H3;
    //     }
    //     else
    //     {
    //         double H0 = -B2;
    //         double H1 = sqrt(B2);
    //         double H3 = pow(H1, 3);
            
    //         double M1 = 1.0 + P20 * cos(L_f) + P10 * sin(L_f);
    //         double M10 = 1.0 + P20 * cos(L_0) + P10 * sin(L_0);
            
    //         double K = 2.0 * atan((P10 + (1.0 - P20) * tan(L_f / 2.0)) / H1);   
    //         if(K < 0.0)
    //             K += 2.0 * pi;                             
    //         double K0 = 2.0 * atan((P10 + (1.0 - P20) * tan(L_0 / 2.0)) / H1); 
    //         if(K0 < 0.0)
    //             K0 += 2.0 * pi;                              
    //         double K1 = K - K0;
    //         if(dL > 0.0) 
    //         {
    //             if(K1 < 0.0)
    //                 K1 += 2.0 * pi;
    //             K1 += 2.0 * pi * double(n_periods);
    //         }
    //         else
    //         {
    //             if(K1 > 0.0)
    //                 K1 -= 2.0 * pi;  
    //             K1 -= 2.0 * pi * double(n_periods);
    //         }  
            
    //         I12 = K1 / H3 + 1.0 / (P20 * H0) * ((P10 + e2 * sin(L_f)) / M1 - (P10 + e2 * sin(L_0)) / M10);
    //         std::cout << L0 << " " << Lf << " " << P10 << " " << P20 << std::endl;
    //         std::cout << I12 << std::endl;            
    //     }
    // } 

    f[0] = (1.0 - Q10 * Q10 + Q20 * Q20) / G;
    f[1] = 2.0 * Q10 * Q20 / G;
    f[2] = -2.0 * Q10 / G;
    g[0] = 2.0 * Q10 * Q20 / G;
    g[1] = (1.0 + Q10 * Q10 - Q20 * Q20) / G;
    g[2] = 2.0 * Q20 / G;
    w[0] = 2.0 * Q10 / G;
    w[1] = -2.0 * Q20 / G;
    w[2] = (1.0 - Q10 * Q10 - Q20 * Q20) / G;

    /** Sun's ephemerides **/
    double eps = constants::obliquity_ecliptic_2000;
    double r_sun, lambda_sun;
    astrodynamics::base_dearth<double>::Sun(jd, r_sun, lambda_sun);
    r_sun *= 1.0e3 / constants::R_earth;

    /** computing Sun's position in inertial frame **/
    std::vector<double> pos_sun(3);
    pos_sun[0] = r_sun * cos(lambda_sun);
    pos_sun[1] = r_sun * sin(lambda_sun) * cos(eps);
    pos_sun[2] = r_sun * sin(lambda_sun) * sin(eps);

    pos = pos_sun;
    mu = mu_sun;
    R_3rd = r_sun;
    
    dir[0] = (pos[0] * f[0] + pos[1] * f[1] + pos[2] * f[2]) / R_3rd;
    dir[1] = (pos[0] * g[0] + pos[1] * g[1] + pos[2] * g[2]) / R_3rd;
    dir[2] = (pos[0] * w[0] + pos[1] * w[1] + pos[2] * w[2]) / R_3rd;

    std::vector<double> I = integrals_3rd_body(a0, P10, P20, Q10, Q20, R_3rd, dir);

    double factor = B2 * (mu / 1.0) * I12 * pow(a0, 3) / (256.0 * pow(R_3rd, 7));
    delta[0] = 0.0;
    delta[1] = factor * I[0];
    delta[2] = -factor * I[1];
    delta[3] = factor * dir[2] * G * 0.5 * I[2];
    delta[4] = factor * dir[2] * G * 0.5 * I[3];   
    
    /** computing Moon's position from ephemerides in ecliptic plane **/
    astrodynamics::base_dearth<double>::moon_pos(jd, pos);

    /** computing Moon's scaled position in reference inertial frame **/
    std::vector<double> pos_moon(3);
    pos_moon[0] = pos[0] * 1.0e3 / constants::R_earth;
    pos_moon[1] = (cos(-eps)*pos[1]+sin(-eps)*pos[2]) * 1.0e3 / constants::R_earth;
    pos_moon[2] = (-sin(-eps)*pos[1]+cos(-eps)*pos[2]) * 1.0e3 / constants::R_earth;

    pos = pos_moon;
    mu = mu_moon;
    R_3rd = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

    dir[0] = (pos[0] * f[0] + pos[1] * f[1] + pos[2] * f[2]) / R_3rd;
    dir[1] = (pos[0] * g[0] + pos[1] * g[1] + pos[2] * g[2]) / R_3rd;
    dir[2] = (pos[0] * w[0] + pos[1] * w[1] + pos[2] * w[2]) / R_3rd; 
    I = integrals_3rd_body(a0, P10, P20, Q10, Q20, R_3rd, dir);
    // std::cout << a0 << " " << P10 << " " << P20 << " " << Q10 << " " << Q20 << " " << R_3rd << " " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl; 
    // std::cout << I[0] << " " << I[1] << " " << I[2] << " " << I[3] << std::endl; 
    // std::cout << std::endl;

    factor = B2 * (mu / 1.0) * I12 * pow(a0, 3) / (256.0 * pow(R_3rd, 7));
    delta[1] += factor * I[0];
    delta[2] += -factor * I[1];
    delta[3] += factor * dir[2] * G * 0.5 * I[2];
    delta[4] += factor * dir[2] * G * 0.5 * I[3];   

    return 0;   
}

std::vector<double>  perturbation_propagator::integralsJ2(const double &L0, const double &Lf) const{

    std::vector<double> I(23, 0.0);

    /* Useful precomputation */
    double pi = constants::pi;
    double L_0 = L0 - 2.0 * pi * floor(L0 / (2.0 * pi));
    double L_f = Lf - 2.0 * pi * floor(Lf / (2.0 * pi));
    double sL  = sin(L_f), sL0 = sin(L_0);
    double cL  = cos(L_f), cL0 = cos(L_0);
    double s2L  = sin(2.0 * L_f), s2L0 = sin(2.0 * L_0);
    double s3L  = sin(3.0 * L_f), s3L0 = sin(3.0 * L_0);
    double c3L  = cos(3.0 * L_f), c3L0 = cos(3.0 * L_0);
    double s4L  = sin(4.0 * L_f), s4L0 = sin(4.0 * L_0);
    double s5L  = sin(5.0 * L_f), s5L0 = sin(5.0 * L_0);
    double c5L  = cos(5.0 * L_f), c5L0 = cos(5.0 * L_0);

    /* Actual integrals */
    I[3] = sL - sL0; // IJc1
    I[4] = cL0 - cL; // IJs1
    I[5] = 0.5 * (Lf - L0) + 0.25 * (s2L - s2L0); // IJc2
    I[6] = -0.5 * (pow(cL, 2) - pow(cL0, 2)); // IJc1s1
    I[7] = 0.5 * (Lf - L0) - 0.25 * (s2L - s2L0); // IJs2
    I[8] = 0.75 * (sL - sL0) + 1.0 / 12.0 * (s3L - s3L0); // IJc3
    I[9] = -(pow(cL, 3) - pow(cL0, 3)) / 3.0; // IJc2s1
    I[10] = (pow(sL, 3) - pow(sL0, 3)) / 3.0; // IJc1s2
    I[11] = -0.75 * (cL - cL0) + 1.0 / 12.0* (c3L - c3L0); // IJs3
    I[12] = 3.0 / 8.0 * (Lf - L0) + 0.25 * (s2L - s2L0) + 1.0 / 32.0 * (s4L - s4L0); // IJc4
    I[13] = -(pow(cL, 4) - pow(cL0, 4)) / 4.0; // IJc3s1 
    I[14] = (Lf - L0) / 8.0 - (s4L - s4L0) / 32.0; // IJc2s2 
    I[15] = (pow(sL, 4) - pow(sL0, 4)) / 4.0; // IJc1s3
    I[16] = 3.0 / 8.0 * (Lf - L0) - 0.25 * (s2L - s2L0) + 1.0 / 32.0 * (s4L - s4L0); // IJs4 
    I[17] = 5.0 / 8.0 * (sL - sL0) + 5.0 / 48.0 * (s3L - s3L0) + 1.0 / 80.0 *(s5L - s5L0); // IJc5
    I[18] = -(pow(cL, 5) - pow(cL0, 5)) / 5.0; // IJc4s1 
    I[19] = 1.0 / 8.0 * (sL - sL0) - 1.0 / 48.0 *(s3L - s3L0) - 1.0 / 80.0 * (s5L - s5L0); // IJc3s2
    I[20] = -1.0 / 8.0 * (cL - cL0) - 1.0 / 48.0 * (c3L - c3L0) + 1.0 / 80.0 * (c5L - c5L0); // IJc2s3 
    I[21] = (pow(sL, 5) - pow(sL0, 5)) / 5.0; // IJc1s4 
    I[22] = -5.0 / 8.0 * (cL - cL0) + 5.0 / 48.0 * (c3L - c3L0) - 1.0 / 80.0 * (c5L - c5L0); // IJs5 

    return I;
}


std::vector<double> perturbation_propagator::integralsJ3(const double &L0, const double &Lf, const double &P1, const double &P2, const double &Q1, const double &Q2) const{

    std::vector<double> I(8, 0.0);

    /* Pre-computations */
    double G = 1.0 + Q1 * Q1 + Q2 * Q2;
    double pi = constants::pi;
    double L_0 = L0 - 2.0 * pi * floor(L0 / (2.0 * pi));
    double L_f = Lf - 2.0 * pi * floor(Lf / (2.0 * pi));
    double sL = sin(L_f), sL0 = sin(L_0);
    double cL = cos(L_f), cL0 = cos(L_0);
    double c2L = cos(2.0 * L_f), c2L0 = cos(2.0 * L_0);
    double s2L = sin(2.0 * L_f), s2L0 = sin(2.0 * L_0);
    double s3L = sin(3.0 * L_f), s3L0 = sin(3.0 * L_0);
    double c3L = cos(3.0 * L_f), c3L0 = cos(3.0 * L_0);
    double c4L = cos(4.0 * L_f), c4L0 = cos(4.0 * L_0);
    double s4L = sin(4.0 * L_f), s4L0 = sin(4.0 * L_0);
    double s5L = sin(5.0 * L_f), s5L0 = sin(5.0 * L_0);
    double c5L = cos(5.0 * L_f), c5L0 = cos(5.0 * L_0);
    double s6L = sin(6.0 * L_f), s6L0 = sin(6.0 * L_0);
    double c6L = cos(6.0 * L_f), c6L0 = cos(6.0 * L_0);
    double s7L = sin(7.0 * L_f), s7L0 = sin(7.0 * L_0);
    double c7L = cos(7.0 * L_f), c7L0 = cos(7.0 * L_0);

    /* Ia1 */
    I[0] = (3.0/2.0)*pow(G, -3)*(Lf - L0)*(P1*Q1+P2*Q2)*((-1.0)*pow(G, 2)*(4.0+3.0*pow(P1, 2)+ 
      3.0*pow(P2, 2))+10.0*((2+pow(P1, 2)+3.0*pow(P2, 2))*pow(Q1, 2)+(-4.0)*P1*P2*Q1*Q2+(2.0+ 
      3.0*pow(P1, 2)+pow(P2, 2))*pow(Q2, 2))) + 
    (cL - cL0) * ( 
    (1.0/4.0)*pow(G, -3)*(6.0*pow(G, 2)*(6.0+pow(P1, 2)+pow(P2, 2))*(pow(P1, 2)*Q1+(-1.0)* 
      pow(P2, 2)*Q1+2.0*P1*P2*Q2)+(-5.0)*(4.0*pow(P1, 3)*P2*Q2*((-3)*pow(Q1, 2)+5.0* 
      pow(Q2, 2))+6.0*pow(P1, 2)*Q1*((4.0+pow(P2, 2))*pow(Q1, 2)+(-3.0)*((-4.0)+pow(P2, 2))*pow(Q2, 2))+ 
      (-1.0)*pow(P2, 2)*Q1*((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+9.0*(8.0+pow(P2, 2))*pow(Q2, 2))+3.0* 
      pow(P1, 4)*(pow(Q1, 3)+5.0*Q1*pow(Q2, 2))+12.0*P1*(8.0*P2*pow(Q2, 3)+pow(P2, 3)*Q2*( 
      pow(Q1, 2)+pow(Q2, 2))))) 
    ) + 
    (c2L - c2L0) * ( 
    (-1.0/4.0)*pow(G, -3)*(6.0*pow(G, 2)*((2.0+(-3.0)*pow(P1, 2))*P2*Q1+3.0*pow(P2, 3)* 
      Q1+P1*(2.0+3.0*pow(P1, 2))*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+(-5.0)*((-3.0)*P1* 
      pow(P2, 2)*Q2*(3.0*pow(Q1, 2)+7.0*pow(Q2, 2))+P1*Q2*(3.0*(8.0+9.0*pow(P1, 2))*pow(Q1, 2)+(8.0+ 
      15.0*pow(P1, 2))*pow(Q2, 2))+3.0*pow(P2, 3)*(5.0*pow(Q1, 3)+9.0*Q1*pow(Q2, 2))+P2*((8.0+(-21.0) 
      *pow(P1, 2))*pow(Q1, 3)+3.0*(8.0+(-3.0)*pow(P1, 2))*Q1*pow(Q2, 2))))
     ) + 
    (c3L - c3L0) * ( 
    (1.0/4.0)*pow(G, -3.0)*(pow(P1, 4)*Q1+6.0*pow(P1, 2)*(2.0+pow(P2, 2))*Q1+(-3.0)* 
      pow(P2, 2)*(4.0+pow(P2, 2))*Q1+(-24.0)*P1*P2*Q2+(-8.0)*pow(P1, 3)*P2*Q2)*( 
      pow(G, 2)+(-5.0)*(pow(Q1, 2)+pow(Q2, 2)))
    ) + 
    (c4L - c4L0) * ( 
    (1.0/8.0)*pow(G, -3.0)*(9.0*pow(G, 2)*(3.0*pow(P1, 2)*P2*Q1+(-1.0)*pow(P2, 3)*Q1+ 
      pow(P1, 3)*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+(-20)*((-3.0)*pow(P2, 3)*pow(Q1, 3)+(-3.0)* 
      P1*pow(P2, 2)*Q2*(3.0*pow(Q1, 2)+2.0*pow(Q2, 2))+P1*Q2*((-3.0)*pow(Q1, 2)+(1.0+3.0* 
      pow(P1, 2))*pow(Q2, 2))+P2*Q1*(((-1.0)+6.0*pow(P1, 2))*pow(Q1, 2)+3.0*(1.0+3.0*pow(P1, 2))* 
      pow(Q2, 2)))) 
    ) + 
    (c5L - c5L0) * ( 
    (1.0/20)*pow(G, -3.0)*((-3.0)*pow(G, 2)*(pow(P1, 4)*Q1+(-6.0)*pow(P1, 2)*pow(P2, 2)* 
      Q1+pow(P2, 4)*Q1+(-4.0)*pow(P1, 3)*P2*Q2+4.0*P1*pow(P2, 3)*Q2)+5.0*((-16.0)* 
      pow(P1, 3)*P2*pow(Q2, 3)+8.0*P1*P2*Q2*(3.0*(3.0+pow(P2, 2))*pow(Q1, 2)+((-3.0)+pow(P2, 2)) 
      *pow(Q2, 2))+(-6.0)*pow(P1, 2)*Q1*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+3.0*((-2.0)+pow(P2, 2))* 
      pow(Q2, 2))+pow(P2, 2)*Q1*((12.0+5.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(12.0+pow(P2, 2))*pow(Q2, 2))+ 
      pow(P1, 4)*(pow(Q1, 3)+9.0*Q1*pow(Q2, 2)))) 
    ) + 
    (c6L - c6L0) * ( 
     (-5.0/4.0)*pow(G, -3.0)*(3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)* 
      pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3))) 
    ) + 
    (c7L - c7L0) * ( 
    (5.0/28.0)*pow(G, -3.0)*((-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+ 
      pow(P2, 4)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+4.0*pow(P1, 3)*P2*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+(-4.0)*P1*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 3)+(-3.0) 
      *Q1*pow(Q2, 2))) 
    ) + 
     (sL - sL0) * ( 
    (1.0/4.0)*pow(G, -3.0)*(6.0*pow(G, 2)*(6.0+pow(P1, 2)+pow(P2, 2))*((-2.0)*P1*P2*Q1+ 
      pow(P1, 2)*Q2+(-1.0)*pow(P2, 2)*Q2)+(-5.0)*((-12.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+ 
      pow(Q2, 2))+6.0*pow(P1, 2)*Q2*(3.0*(4.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*((-4.0)+pow(P2, 2))* 
      pow(Q2, 2))+(-3.0)*pow(P2, 2)*Q2*((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(8.0+pow(P2, 2))*pow(Q2, 2))+ 
      pow(P1, 4)*(9.0*pow(Q1, 2)*Q2+5.0*pow(Q2, 3))+(-4.0)*P1*(24.0*P2*pow(Q1, 3)+pow(P2, 3)*( 
      5.0*pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))))) 
    ) + 
    (s2L - s2L0) * ( 
     (1.0/4.0)*pow(G, -3.0)*(12.0*pow(G, 2)*((-1.0)*P1*(Q1+3.0*pow(P2, 2)*Q1)+P2*Q2+ 
      3.0*pow(P1, 2)*P2*Q2)+5.0*((-16.0)*P2*pow(Q2, 3)+(-3.0)*pow(P1, 2)*P2*Q2*(9.0* 
      pow(Q1, 2)+13.0*pow(Q2, 2))+3.0*pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))+P1*((16.0+39.0* 
      pow(P2, 2))*pow(Q1, 3)+27.0*pow(P2, 2)*Q1*pow(Q2, 2))+pow(P2, 3)*(9.0*pow(Q1, 2)*Q2+(-3.0)* 
      pow(Q2, 3)))) 
    ) + 
    (s3L - s3L0) * ( 
    (1.0/4.0)*pow(G, -3.0)*(8.0*P1*P2*(3.0+pow(P2, 2))*Q1+3.0*pow(P1, 4)*Q2+(-6.0)* 
      pow(P1, 2)*((-2.0)+pow(P2, 2))*Q2+(-1.0)*pow(P2, 2)*(12.0+pow(P2, 2))*Q2)*((-1.0)*pow(G, 2)+ 
      5.0*(pow(Q1, 2)+pow(Q2, 2))) 
    ) + 
    (s4L - s4L0) * ( 
    (1.0/8.0)*pow(G, -3.0)*(9.0*pow(G, 2)*(pow(P1, 3)*Q1+(-3.0)*P1*pow(P2, 2)*Q1+(-3.0)* 
      pow(P1, 2)*P2*Q2+pow(P2, 3)*Q2)+(-10)*((-3.0)*pow(P1, 2)*P2*Q2*(3.0*pow(Q1, 2)+ 
      5.0*pow(Q2, 2))+P2*Q2*((6.0+9.0*pow(P2, 2))*pow(Q1, 2)+((-2.0)+3.0*pow(P2, 2))*pow(Q2, 2))+( 
      -1.0)*P1*Q1*((2.0+15.0*pow(P2, 2))*pow(Q1, 2)+3.0*((-2.0)+3.0*pow(P2, 2))*pow(Q2, 2))+3.0* 
      pow(P1, 3)*(pow(Q1, 3)+3.0*Q1*pow(Q2, 2)))) 
    ) + 
    (s5L - s5L0) * ( 
    (1.0/20)*pow(G, -3.0)*(3.0*pow(G, 2)*(4.0*pow(P1, 3)*P2*Q1+(-4.0)*P1*pow(P2, 3)* 
      Q1+pow(P1, 4)*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q2+pow(P2, 4)*Q2)+5.0*(8.0*P1*P2*Q1* 
      ((3.0+2.0*pow(P2, 2))*pow(Q1, 2)+(-9.0)*pow(Q2, 2))+(-8.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+3.0* 
      pow(Q2, 2))+(-1.0)*pow(P2, 2)*Q2*(9.0*(4.0+pow(P2, 2))*pow(Q1, 2)+((-12.0)+pow(P2, 2))*pow(Q2, 2)) 
      +6.0*pow(P1, 2)*Q2*(3.0*(2.0+pow(P2, 2))*pow(Q1, 2)+((-2.0)+3.0*pow(P2, 2))*pow(Q2, 2))+ 
      pow(P1, 4)*(3.0*pow(Q1, 2)*Q2+(-5.0)*pow(Q2, 3))))  
    ) + 
    (s6L - s6L0) * ( 
     (-5.0/4.0)*pow(G, -3.0)*((-3.0)*P1*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0* 
      pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))) 
    ) + 
    (s7L - s7L0) * ( 
     (-5.0/28.0)*pow(G, -3.0)*(4.0*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-4.0)* 
      P1*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*Q2*((-3.0)* 
      pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 4)*(3.0* 
      pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3))) 
    );

    I[0] /= 4.0;

    /* Ia2 */
    I[1] = (Lf - L0) * ( 
      (-3.0/2.0)*pow(G, -3.0)*(P1*Q1+P2*Q2)*((-1.0)*pow(G, 2)*(4.0+3.0*pow(P1, 2)+3.0* 
      pow(P2, 2))+10.0*((2.0+pow(P1, 2)+3.0*pow(P2, 2))*pow(Q1, 2)+(-4.0)*P1*P2*Q1*Q2+(2.0+3.0* 
      pow(P1, 2)+pow(P2, 2))*pow(Q2, 2))) 
        ) + 
    (cL - cL0) * ( 
    (-3.0/16.0)*pow(G, -3.0)*(2.0*pow(G, 2)*((8.0+5.0*pow(P1, 4)+12.0*pow(P2, 2)+pow(P2, 4)+6.0* 
      pow(P1, 2)*(6.0+pow(P2, 2)))*Q1+4.0*P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))*Q2)+(-5.0)*((16.0+ 
      5.0*pow(P1, 4)+48.0*pow(P2, 2)+5.0*pow(P2, 4)+6.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2)))*pow(Q1, 3)+(-4.0)* 
      P1*P2*(24.0+7.0*pow(P1, 2)+pow(P2, 2))*pow(Q1, 2)*Q2+(16.0+25.0*pow(P1, 4)+(-48.0)*pow(P2, 2)+ 
      (-7.0)*pow(P2, 4)+(-6.0)*pow(P1, 2)*((-24.0)+pow(P2, 2)))*Q1*pow(Q2, 2)+4.0*P1*P2*(5.0* 
      pow(P1, 2)+3.0*(8.0+pow(P2, 2)))*pow(Q2, 3))) 
    ) + 
    (c2L - c2L0) * ( 
    (-3.0/4.0)*pow(G, -3.0)*(2.0*pow(G, 2)*((2.0+3.0*pow(P1, 2))*P2*Q1+pow(P2, 3)*Q1+P1* 
      (2.0+pow(P1, 2))*Q2+3.0*P1*pow(P2, 2)*Q2)+5.0*(3.0*P1*pow(P2, 2)*Q2*(pow(Q1, 2)+(-3.0) 
      *pow(Q2, 2))+P1*Q2*((8.0+7.0*pow(P1, 2))*pow(Q1, 2)+(-1.0)*(8.0+5.0*pow(P1, 2))*pow(Q2, 2))+ 
      pow(P2, 3)*((-5.0)*pow(Q1, 3)+7.0*Q1*pow(Q2, 2))+P2*((-1.0)*(8.0+9.0*pow(P1, 2))*pow(Q1, 3)+( 
      8.0+3.0*pow(P1, 2))*Q1*pow(Q2, 2)))) 
    ) + 
    (c3L - c3L0) * ( 
    (1.0/16.0)*pow(G, -3.0)*(pow(G, 2)*(5.0*pow(P1, 4)*Q1+(-6.0)*pow(P1, 2)*((-4.0)+pow(P2, 2)) 
      *Q1+(-3.0)*pow(P2, 2)*(8.0+pow(P2, 2))*Q1+(-4.0)*pow(P1, 3)*P2*Q2+(-12.0)*P1* 
      P2*(4.0+pow(P2, 2))*Q2)+5.0*((16.0+pow(P1, 4)+72.0*pow(P2, 2)+9.0*pow(P2, 4)+6.0*pow(P1, 2)*(4.0+ 
      3.0*pow(P2, 2)))*pow(Q1, 3)+4.0*P1*P2*(pow(P1, 2)+3.0*(4.0+pow(P2, 2)))*pow(Q1, 2)*Q2+(-1.0) 
      *(23.0*pow(P1, 4)+6.0*pow(P1, 2)*(28.0+5.0*pow(P2, 2))+3.0*(16.0+40.0*pow(P2, 2)+5.0*pow(P2, 4)))* 
      Q1*pow(Q2, 2)+4.0*P1*P2*(pow(P1, 2)+3.0*(4.0+pow(P2, 2)))*pow(Q2, 3))) 
    ) + 
    (c4L - c4L0) * ( 
    (3.0/8.0)*pow(G, -3.0)*(pow(G, 2)*(3.0*pow(P1, 2)*P2*Q1+(-1.0)*pow(P2, 3)*Q1+ 
      pow(P1, 3)*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+20.0*(3.0*P1*pow(P2, 2)*pow(Q1, 2)*Q2+P2* 
      Q1*(pow(Q1, 2)+(-3.0)*(1.0+pow(P1, 2))*pow(Q2, 2))+(-1.0)*P1*Q2*((-1.0)*(3.0+2.0* 
      pow(P1, 2))*pow(Q1, 2)+(1.0+pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(pow(Q1, 3)+(-2.0)*Q1*pow(Q2, 2)))) 
    ) +  
    (c5L - c5L0) * ( 
    (-3.0/80)*pow(G, -3.0)*(pow(G, 2)*(pow(P1, 4)*Q1+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1+ 
      pow(P2, 4)*Q1+(-4.0)*pow(P1, 3)*P2*Q2+4.0*P1*pow(P2, 3)*Q2)+5.0*(4.0*pow(P1, 3)*P2* 
      Q2*((-5.0)*pow(Q1, 2)+3.0*pow(Q2, 2))+6.0*pow(P1, 2)*Q1*((4.0+pow(P2, 2))*pow(Q1, 2)+((-12.0)+ 
      pow(P2, 2))*pow(Q2, 2))+4.0*P1*P2*Q2*((-1.0)*(36.0+7.0*pow(P2, 2))*pow(Q1, 2)+(12.0+ 
      pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)*Q1*((-1.0)*(24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(72.0+11.0* 
      pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(3.0*pow(Q1, 3)+(-13.0)*Q1*pow(Q2, 2)))) 
    ) + 
    (c6L - c6L0) * ( 
    (-5.0/4.0)*pow(G, -3.0)*(3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)* 
      pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3))) 
    ) + 
    (c7L - c7L0) * ( 
    (15.0/112.0)*pow(G, -3.0)*((-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2)) 
      +pow(P2, 4)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+4.0*pow(P1, 3)*P2*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+(-4.0)*P1*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 3)+(-3.0) 
      *Q1*pow(Q2, 2))) 
    ) + 
     (sL - sL0) * ( 
    (3.0/16.0)*pow(G, -3.0)*(2.0*pow(G, 2)*(4.0*pow(P1, 3)*P2*Q1+4.0*P1*P2*(6.0+ 
      pow(P2, 2))*Q1+pow(P1, 4)*Q2+6.0*pow(P1, 2)*(2.0+pow(P2, 2))*Q2+(8.0+36.0*pow(P2, 2)+5.0* 
      pow(P2, 4))*Q2)+5.0*(4.0*pow(P1, 3)*P2*Q1*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+6.0*pow(P1, 2)* 
      Q2*((8.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*(8.0+3.0*pow(P2, 2))*pow(Q2, 2))+(-4.0)*P1*P2* 
      Q1*((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(24.0+7.0*pow(P2, 2))*pow(Q2, 2))+(-1.0)*Q2*(( 
      16.0+144.0*pow(P2, 2)+25.0*pow(P2, 4))*pow(Q1, 2)+(16.0+48.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q2, 2))+ 
      pow(P1, 4)*(7.0*pow(Q1, 2)*Q2+(-5.0)*pow(Q2, 3)))) 
    ) + 
    (s2L - s2L0) * ( 
    (-3.0/4.0)*pow(G, -3.0)*(4.0*pow(G, 2)*(P1*Q1+pow(P1, 3)*Q1+(-1.0)*P2*(1.0+ 
      pow(P2, 2))*Q2)+(-5.0)*(3.0*pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)* 
      P2*Q2*((16.0+13.0*pow(P2, 2))*pow(Q1, 2)+pow(P2, 2)*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+13.0* 
      Q1*pow(Q2, 2))+P1*(16.0*Q1*pow(Q2, 2)+(-3.0)*pow(P2, 2)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)) 
      ))) 
    ) + 
    (s3L - s3L0) * ( 
    (-1.0/16.0)*pow(G, -3.0)*(pow(G, 2)*(12.0*pow(P1, 3)*P2*Q1+4.0*P1*P2*(12.0+ 
      pow(P2, 2))*Q1+3.0*pow(P1, 4)*Q2+6.0*pow(P1, 2)*(4.0+pow(P2, 2))*Q2+(-1.0)*pow(P2, 2)*(24.0+ 
      5.0*pow(P2, 2))*Q2)+5.0*(6.0*pow(P1, 2)*(4.0+pow(P2, 2))*Q2*(5.0*pow(Q1, 2)+(-3.0)*pow(Q2, 2)) 
      +(-12.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+pow(Q2, 2))+(-4.0)*P1*P2*(12.0+pow(P2, 2))* 
      Q1*(pow(Q1, 2)+pow(Q2, 2))+Q2*((48.0+168.0*pow(P2, 2)+23.0*pow(P2, 4))*pow(Q1, 2)+(-1.0)*(16.0+ 
      24.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+3.0*pow(P1, 4)*(5.0*pow(Q1, 2)*Q2+(-3.0)*pow(Q2, 3)))) 
    ) + 
    (s4L - s4L0) * ( 
    (3.0/8.0)*pow(G, -3.0)*(pow(G, 2)*(pow(P1, 3)*Q1+(-3.0)*P1*pow(P2, 2)*Q1+(-3.0)* 
      pow(P1, 2)*P2*Q2+pow(P2, 3)*Q2)+10.0*(3.0*pow(P1, 2)*P2*Q2*((-1.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+P1*Q1*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(2.0+pow(P2, 2))*pow(Q2, 2))+P2* 
      Q2*((-1.0)*(6.0+5.0*pow(P2, 2))*pow(Q1, 2)+(2.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+( 
      -5.0)*Q1*pow(Q2, 2)))) 
    ) + 
    (s5L - s5L0) * ( 
    (3.0/80)*pow(G, -3.0)*(pow(G, 2)*(4.0*pow(P1, 3)*P2*Q1+(-4.0)*P1*pow(P2, 3)*Q1+ 
      pow(P1, 4)*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q2+pow(P2, 4)*Q2)+5.0*(4.0*pow(P1, 3)*P2*Q1* 
      (pow(Q1, 2)+(-7.0)*pow(Q2, 2))+6.0*pow(P1, 2)*Q2*((12.0+pow(P2, 2))*pow(Q1, 2)+((-4.0)+pow(P2, 2)) 
      *pow(Q2, 2))+pow(P2, 2)*Q2*((-1.0)*(72.0+13.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+pow(P2, 2))* 
      pow(Q2, 2))+4.0*P1*P2*Q1*(3.0*(4.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*(36.0+5.0*pow(P2, 2))* 
      pow(Q2, 2))+pow(P1, 4)*(11.0*pow(Q1, 2)*Q2+(-5.0)*pow(Q2, 3)))) 
    ) + 
    (s6L - s6L0) * ( 
    (-5.0/4.0)*pow(G, -3.0)*((-3.0)*P1*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0* 
      pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))) 
    ) + 
    (s7L - s7L0) * ( 
    (-15.0/112.0)*pow(G, -3.0)*(4.0*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-4.0) 
      *P1*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*Q2*((-3.0)* 
      pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 4)*(3.0* 
      pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3))) 
    );

    I[1] /= 3.0;

    /* IP11 */
    I[2] = (Lf - L0) * ( 
    (3.0/8.0)*pow(G, -3.0)*(pow(G, 2)*((4.0+3.0*pow(P1, 2)+9.0*pow(P2, 2))*Q1+(-6.0)*P1*P2* 
      Q2)+(-10)*((2.0+pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 3)+(-6.0)*P1*P2*pow(Q1, 2)*Q2+(2.0+ 
      3.0*pow(P1, 2)+3.0*pow(P2, 2))*Q1*pow(Q2, 2)+(-2.0)*P1*P2*pow(Q2, 3)))    
        ) + 
    (cL - cL0) * ( 
    (-3.0/16.0)*pow(G, -3.0)*(2.0*pow(G, 2)*(pow(P1, 3)*Q1+3.0*P1*(2.0+pow(P2, 2))*Q1+(-3.0)* 
      pow(P1, 2)*P2*Q2+(-1.0)*P2*(6.0+pow(P2, 2))*Q2)+(-5.0)*((-1.0)*pow(P1, 2)*P2* 
      Q2*(9.0*pow(Q1, 2)+5.0*pow(Q2, 2))+(-1.0)*P2*Q2*((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(8.0+ 
      pow(P2, 2))*pow(Q2, 2))+P1*Q1*((8.0+5.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+3.0*pow(P2, 2))* 
      pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+5.0*Q1*pow(Q2, 2)))) 
    ) + 
    (c2L - c2L0) * ( 
    (1.0/16.0)*pow(G, -3.0)*(6.0*pow(G, 2)*((-6.0)*P1*P2*Q1+3.0*pow(P1, 2)*Q2+(2.0+3.0* 
      pow(P2, 2))*Q2)+(-5.0)*((-6.0)*P1*P2*Q1*(5.0*pow(Q1, 2)+9.0*pow(Q2, 2))+Q2*(3.0*( 
      8.0+15.0*pow(P2, 2))*pow(Q1, 2)+(8.0+9.0*pow(P2, 2))*pow(Q2, 2))+3.0*pow(P1, 2)*(9.0*pow(Q1, 2)*Q2+ 
      5.0*pow(Q2, 3)))) 
    ) + 
    (c3L - c3L0) * ( 
    (1.0/16.0)*pow(G, -3.0)*(pow(G, 2)*((-1.0)*pow(P1, 3)*Q1+(-3.0)*P1*(4.0+3.0*pow(P2, 2))* 
      Q1+3.0*pow(P1, 2)*P2*Q2+3.0*P2*(4.0+pow(P2, 2))*Q2)+5.0*(pow(P1, 3)*Q1*(pow(Q1, 2)+ 
      pow(Q2, 2))+3.0*P1*(4.0+3.0*pow(P2, 2))*Q1*(pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P1, 2)*P2* 
      Q2*(9.0*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*P2*(4.0+pow(P2, 2))*Q2*(9.0*pow(Q1, 2)+pow(Q2, 2)))) 
    ) + 
    (c4L - c4L0) * ( 
    (1.0/32.0)*pow(G, -3.0)*((-9.0)*pow(G, 2)*(2.0*P1*P2*Q1+pow(P1, 2)*Q2+(-1.0)* 
      pow(P2, 2)*Q2)+20.0*(6.0*P1*P2*pow(Q1, 3)+pow(Q2, 3)+3.0*pow(P1, 2)*pow(Q2, 3)+(-3.0)* 
      pow(Q1, 2)*(Q2+3.0*pow(P2, 2)*Q2))) 
    ) + 
    (c5L - c5L0) * ( 
    (1.0/80)*pow(G, -3.0)*(3.0*pow(G, 2)*(pow(P1, 3)*Q1+(-3.0)*P1*pow(P2, 2)*Q1+(-3.0)* 
      pow(P1, 2)*P2*Q2+pow(P2, 3)*Q2)+(-5.0)*((-9.0)*pow(P1, 2)*P2*Q2*(pow(Q1, 2)+pow(Q2, 2)) 
      +(-3.0)*P1*Q1*((4.0+5.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(4.0+pow(P2, 2))*pow(Q2, 2))+P2* 
      Q2*(3.0*(12.0+5.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(12.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*( 
      pow(Q1, 3)+9.0*Q1*pow(Q2, 2)))) 
    ) + 
    (c6L - c6L0) * ( 
    (5.0/16.0)*pow(G, -3.0)*(2.0*P1*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+pow(P2, 2)*Q2*(( 
      -3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3))) 
    ) + 
    (c7L - c7L0) * ( 
    (-5.0/112.0)*pow(G, -3.0)*((-3.0)*P1*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0* 
      pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))) 
    ) + 
     (sL - sL0) * ( 
    (1.0/16.0)*pow(G, -3.0)*(6.0*pow(G, 2)*(3.0*(6.0+pow(P1, 2))*P2*Q1+5.0*pow(P2, 3)*Q1+(-1.0) 
      *P1*(6.0+pow(P1, 2))*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+5.0*(9.0*P1*pow(P2, 2)*Q2*( 
      5.0*pow(Q1, 2)+pow(Q2, 2))+(-3.0)*P2*Q1*(5.0*(8.0+pow(P1, 2))*pow(Q1, 2)+3.0*(8.0+3.0*pow(P1, 2)) 
      *pow(Q2, 2))+P1*Q2*(9.0*(8.0+pow(P1, 2))*pow(Q1, 2)+(24.0+5.0*pow(P1, 2))*pow(Q2, 2))+(-5.0)* 
      pow(P2, 3)*(7.0*pow(Q1, 3)+3.0*Q1*pow(Q2, 2)))) 
    ) + 
    (s2L - s2L0) * ( 
    (1.0/16.0)*pow(G, -3.0)*(12.0*pow(G, 2)*(Q1+3.0*pow(P2, 2)*Q1)+(-5.0)*((16.0+3.0*pow(P1, 2)+ 
      45.0*pow(P2, 2))*pow(Q1, 3)+(-18.0)*P1*P2*pow(Q1, 2)*Q2+9.0*((-1.0)*pow(P1, 2)+pow(P2, 2)) 
      *Q1*pow(Q2, 2)+6.0*P1*P2*pow(Q2, 3))) 
    ) + 
    (s3L - s3L0) * ( 
    (1.0/16.0)*pow(G, -3.0)*(pow(G, 2)*((-3.0)*((-4.0)+pow(P1, 2))*P2*Q1+5.0*pow(P2, 3)*Q1+ 
      3.0*P1*(4.0+pow(P1, 2))*Q2+3.0*P1*pow(P2, 2)*Q2)+(-5.0)*(3.0*P1*(4.0+pow(P1, 2))* 
      Q2*(pow(Q1, 2)+pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*(pow(Q1, 2)+pow(Q2, 2))+(-1.0)*P2*Q1*(( 
      (-20)+pow(P1, 2))*pow(Q1, 2)+3.0*(4.0+3.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(7.0*pow(Q1, 3)+(-1.0) 
      *Q1*pow(Q2, 2)))) 
    ) + 
    (s4L - s4L0) * ( 
    (1.0/32.0)*pow(G, -3.0)*(pow(G, 2)*((-9.0)*pow(P1, 2)*Q1+9.0*pow(P2, 2)*Q1+18.0*P1*P2* 
      Q2)+10.0*(((-2.0)+3.0*pow(P1, 2)+(-9.0)*pow(P2, 2))*pow(Q1, 3)+(-18.0)*P1*P2*pow(Q1, 2)* 
      Q2+3.0*(2.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))*Q1*pow(Q2, 2)+(-6.0)*P1*P2*pow(Q2, 3))) 
    ) + 
    (s5L - s5L0) * ( 
    (1.0/80)*pow(G, -3.0)*(pow(G, 2)*((-9.0)*pow(P1, 2)*P2*Q1+3.0*pow(P2, 3)*Q1+(-3.0)* 
      pow(P1, 3)*Q2+9.0*P1*pow(P2, 2)*Q2)+(-5.0)*(3.0*P1*pow(P2, 2)*Q2*(9.0*pow(Q1, 2)+ 
      pow(Q2, 2))+(-3.0)*P2*Q1*(((-4.0)+3.0*pow(P1, 2))*pow(Q1, 2)+3.0*(4.0+pow(P1, 2))*pow(Q2, 2))+ 
      P1*Q2*(3.0*(12.0+pow(P1, 2))*pow(Q1, 2)+(-1.0)*(12.0+5.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*( 
      7.0*pow(Q1, 3)+(-9.0)*Q1*pow(Q2, 2)))) 
    ) + 
    (s6L - s6L0) * ( 
    (5.0/16.0)*pow(G, -3.0)*((-1.0)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+2.0*P1*P2* 
      Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))) 
    ) + 
    (s7L - s7L0) * ( 
    (5.0/112.0)*pow(G, -3.0)*(3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)* 
      pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3))) 
    );

    I[2] *= G;

    /* IP12 */
    I[3] = (Lf - L0) * ( 
      (1.0/2.0)*pow(G, -3.0)*(pow(G, 2)*((2.0+5.0*pow(P1, 2)+pow(P2, 2))*Q1+4.0*P1*P2*Q2)+( 
      -10)*((1.0+2.0*pow(P1, 2)+pow(P2, 2))*pow(Q1, 3)+(-1.0)*P1*P2*pow(Q1, 2)*Q2+(1.0+4.0* 
      pow(P1, 2)+(-1.0)*pow(P2, 2))*Q1*pow(Q2, 2)+3.0*P1*P2*pow(Q2, 3)))
        ) + 
    (cL - cL0) * ( 
    (1.0/16.0)*pow(G, -3.0)*((-2.0)*pow(G, 2)*(11.0*pow(P1, 3)*Q1+P1*(38.0+5.0*pow(P2, 2))* 
      Q1+7.0*pow(P1, 2)*P2*Q2+P2*(10+pow(P2, 2))*Q2)+5.0*(P1*Q1*((56.0+17.0* 
      pow(P2, 2))*pow(Q1, 2)+(136.0+(-11.0)*pow(P2, 2))*pow(Q2, 2))+P2*Q2*((-1.0)*(40+pow(P2, 2)) 
      *pow(Q1, 2)+(40+3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(13.0*pow(Q1, 3)+49.0*Q1*pow(Q2, 2))+ 
      pow(P1, 2)*P2*((-37.0)*pow(Q1, 2)*Q2+31.0*pow(Q2, 3))))
    ) + 
    (c2L - c2L0) * ( 
     (-1.0/4.0)*pow(G, -3.0)*(2.0*pow(G, 2)*(3.0*P1*P2*Q1+Q2+2.0*pow(P1, 2)*Q2+pow(P2, 2)* 
      Q2)+5.0*(Q2*((4.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*(4.0+3.0*pow(P2, 2))*pow(Q2, 2))+P1*P2*( 
      (-10)*pow(Q1, 3)+6.0*Q1*pow(Q2, 2))+pow(P1, 2)*(11.0*pow(Q1, 2)*Q2+(-9.0)*pow(Q2, 3))))
    ) + 
    (c3L - c3L0) * ( 
    (1.0/48.0)*pow(G, -3.0)*(pow(G, 2)*(9.0*pow(P1, 3)*Q1+P1*(20+(-7.0)*pow(P2, 2))*Q1+( 
      -11.0)*pow(P1, 2)*P2*Q2+(-1.0)*P2*(20+3.0*pow(P2, 2))*Q2)+5.0*(11.0*pow(P1, 2)* 
      P2*Q2*(pow(Q1, 2)+pow(Q2, 2))+P2*(20+3.0*pow(P2, 2))*Q2*(pow(Q1, 2)+pow(Q2, 2))+P1* 
      Q1*(3.0*(12.0+7.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(188.0+35.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)* 
      (5.0*pow(Q1, 3)+(-51.0)*Q1*pow(Q2, 2))))
    ) + 
    (c4L - c4L0) * ( 
     (1.0/8.0)*pow(G, -3.0)*(pow(G, 2)*(2.0*P1*P2*Q1+pow(P1, 2)*Q2+(-1.0)*pow(P2, 2)*Q2)+ 
      10.0*(P1*P2*Q1*(pow(Q1, 2)+(-7.0)*pow(Q2, 2))+Q2*((3.0+2.0*pow(P2, 2))*pow(Q1, 2)+(-1.0) 
      *pow(Q2, 2))+pow(P1, 2)*(7.0*pow(Q1, 2)*Q2+(-3.0)*pow(Q2, 3))))
    ) + 
    (c5L - c5L0) * ( 
    (-1.0/80)*pow(G, -3.0)*(pow(G, 2)*(pow(P1, 3)*Q1+(-3.0)*P1*pow(P2, 2)*Q1+(-3.0)* 
      pow(P1, 2)*P2*Q2+pow(P2, 3)*Q2)+5.0*(P2*Q2*((-1.0)*(60+7.0*pow(P2, 2))*pow(Q1, 2)+( 
      20+pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(7.0*pow(Q1, 3)+(-25.0)*Q1*pow(Q2, 2))+P1*((-1.0)*(( 
      -20)+pow(P2, 2))*pow(Q1, 3)+15.0*((-4.0)+pow(P2, 2))*Q1*pow(Q2, 2))+pow(P1, 2)*P2*((-39.0) 
      *pow(Q1, 2)*Q2+17.0*pow(Q2, 3)))) 
    ) + 
    (c6L - c6L0) * ( 
     (-5.0/12.0)*pow(G, -3.0)*(2.0*P1*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+pow(P2, 2)*Q2*( 
      (-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))
    ) + 
    (c7L - c7L0) * ( 
    (5.0/112.0)*pow(G, -3.0)*((-3.0)*P1*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0* 
      pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))) 
    ) + 
     (sL - sL0) * ( 
    (1.0/16.0)*pow(G, -3.0)*(2.0*pow(G, 2)*((10+7.0*pow(P1, 2))*P2*Q1+pow(P2, 3)*Q1+3.0* 
      P1*(6.0+pow(P1, 2))*Q2+9.0*P1*pow(P2, 2)*Q2)+5.0*((-1.0)*P1*pow(P2, 2)*Q2*(21.0* 
      pow(Q1, 2)+17.0*pow(Q2, 2))+P1*Q2*(3.0*(8.0+5.0*pow(P1, 2))*pow(Q1, 2)+(-1.0)*(56.0+13.0* 
      pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*((-5.0)*pow(Q1, 3)+7.0*Q1*pow(Q2, 2))+P2*((-5.0)*(8.0+5.0* 
      pow(P1, 2))*pow(Q1, 3)+(40+19.0*pow(P1, 2))*Q1*pow(Q2, 2)))) 
    ) + 
    (s2L - s2L0) * ( 
     (1.0/4.0)*pow(G, -3.0)*((-2.0)*pow(G, 2)*(Q1+3.0*pow(P1, 2)*Q1+(-1.0)*P1*P2*Q2)+ 
      5.0*(8.0*Q1*pow(Q2, 2)+2.0*P1*P2*Q2*((-7.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 2)*( 
      pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))+pow(P1, 2)*(pow(Q1, 3)+21.0*Q1*pow(Q2, 2))))
    ) + 
    (s3L - s3L0) * ( 
     (1.0/48.0)*pow(G, -3.0)*((-1.0)*pow(G, 2)*((20+17.0*pow(P1, 2))*P2*Q1+pow(P2, 3)*Q1+ 
      P1*(20+7.0*pow(P1, 2))*Q2+(-1.0)*P1*pow(P2, 2)*Q2)+5.0*((20+17.0*pow(P1, 2))* 
      P2*Q1*(pow(Q1, 2)+pow(Q2, 2))+pow(P2, 3)*Q1*(pow(Q1, 2)+pow(Q2, 2))+P1*Q2*((-1.0)*( 
      148.0+35.0*pow(P1, 2))*pow(Q1, 2)+(76.0+21.0*pow(P1, 2))*pow(Q2, 2))+P1*pow(P2, 2)*((-43.0)* 
      pow(Q1, 2)*Q2+13.0*pow(Q2, 3))))
    ) + 
    (s4L - s4L0) * ( 
    (1.0/8.0)*pow(G, -3.0)*(pow(G, 2)*(pow(P1, 2)*Q1+(-1.0)*pow(P2, 2)*Q1+(-2.0)*P1*P2* 
      Q2)+10.0*((1.0+2.0*pow(P1, 2)+pow(P2, 2))*pow(Q1, 3)+(-5.0)*P1*P2*pow(Q1, 2)*Q2+(-1.0)*( 
      3.0+8.0*pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 2)+3.0*P1*P2*pow(Q2, 3))) 
    ) + 
    (s5L - s5L0) * ( 
    (1.0/80)*pow(G, -3.0)*(pow(G, 2)*(3.0*pow(P1, 2)*P2*Q1+(-1.0)*pow(P2, 3)*Q1+pow(P1, 3)* 
      Q2+(-3.0)*P1*pow(P2, 2)*Q2)+5.0*(P2*Q1*((20+11.0*pow(P1, 2))*pow(Q1, 2)+(-15.0)* 
      (4.0+3.0*pow(P1, 2))*pow(Q2, 2))+P1*Q2*((60+23.0*pow(P1, 2))*pow(Q1, 2)+(-1.0)*(20+9.0* 
      pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(3.0*pow(Q1, 3)+(-5.0)*Q1*pow(Q2, 2))+P1*pow(P2, 2)*((-9.0) 
      *pow(Q1, 2)*Q2+7.0*pow(Q2, 3)))) 
    ) + 
    (s6L - s6L0) * ( 
     (-5.0/12.0)*pow(G, -3.0)*((-1.0)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+2.0*P1* 
      P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))
    ) + 
    (s7L - s7L0) * ( 
     (-5.0/112.0)*pow(G, -3.0)*(3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)* 
      pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))
    );

    I[3] *= G;

    /* IP13 */
    I[4] = (Lf - L0) * ( 
      pow(G, -2.0)*(P2*Q1+(-1.0)*P1*Q2)*(pow(G, 2)+(-15.0)*(pow(Q1, 2)+pow(Q2, 2))) 
        ) + 
    (cL - cL0) * ( 
    (1.0/4.0)*pow(G, -2.0)*(pow(G, 2)*((-2.0)*P1*P2*Q1+3.0*pow(P1, 2)*Q2+(4.0+ 
      pow(P2, 2))*Q2)+(-10)*((-2.0)*P1*P2*Q1*(pow(Q1, 2)+3.0*pow(Q2, 2))+Q2*(3.0*(2.0+ 
      pow(P2, 2))*pow(Q1, 2)+(6.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+5.0*pow(Q2, 3))))
    ) + 
    (c2L - c2L0) * ( 
    (1.0/2.0)*pow(G, -2.0)*(pow(G, 2)*((-1.0)*P1*Q1+P2*Q2)+10.0*((-1.0)*P2* 
      Q2*(3.0*pow(Q1, 2)+pow(Q2, 2))+P1*(pow(Q1, 3)+3.0*Q1*pow(Q2, 2)))) 
    ) + 
    (c3L - c3L0) * ( 
    (1.0/12.0)*pow(G, -2.0)*(pow(G, 2)*((-2.0)*P1*P2*Q1+(-1.0)*pow(P1, 2)*Q2+ 
      pow(P2, 2)*Q2)+(-5.0)*((-6.0)*P1*P2*Q1*(pow(Q1, 2)+pow(Q2, 2))+Q2*(3.0*(4.0+3.0* 
      pow(P2, 2))*pow(Q1, 2)+((-4.0)+pow(P2, 2))*pow(Q2, 2))+pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+(-5.0)* 
      pow(Q2, 3))))
    ) + 
    (c4L - c4L0) * ( 
     (-5.0/4.0)*pow(G, -2.0)*((-1.0)*P1*pow(Q1, 3)+3.0*P2*pow(Q1, 2)*Q2+3.0*P1*Q1* 
      pow(Q2, 2)+(-1.0)*P2*pow(Q2, 3))
    ) + 
    (c5L - c5L0) * ( 
     (-1.0/4.0)*pow(G, -2.0)*((-2.0)*P1*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)* 
      pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*((-3.0)*pow(Q1, 2)*Q2+pow(Q2, 3)))
    ) + 
     (sL - sL0) * ( 
     (1.0/4.0)*pow(G, -2.0)*(pow(G, 2)*((4.0+pow(P1, 2)+3.0*pow(P2, 2))*Q1+(-2.0)*P1*P2* 
      Q2)+(-10)*((6.0+pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 3)+(-6.0)*P1*P2*pow(Q1, 2)*Q2+3.0*( 
      2.0+pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 2)+(-2.0)*P1*P2*pow(Q2, 3)))
    ) + 
    (s2L - s2L0) * ( 
    (1.0/2.0)*pow(G, -2.0)*(pow(G, 2)*(P2*Q1+P1*Q2)+(-20)*(P2*pow(Q1, 3)+P1* 
      pow(Q2, 3)))
    ) + 
    (s3L - s3L0) * ( 
     (1.0/12.0)*pow(G, -2.0)*(pow(G, 2)*((-1.0)*pow(P1, 2)*Q1+pow(P2, 2)*Q1+2.0*P1*P2* 
      Q2)+5.0*(((-4.0)+pow(P1, 2)+(-5.0)*pow(P2, 2))*pow(Q1, 3)+(-6.0)*P1*P2*pow(Q1, 2)*Q2+ 
      3.0*(4.0+3.0*pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 2)+(-6.0)*P1*P2*pow(Q2, 3)))
    ) + 
    (s4L - s4L0) * ( 
    (-5.0/4.0)*pow(G, -2.0)*(P2*pow(Q1, 3)+3.0*P1*pow(Q1, 2)*Q2+(-3.0)*P2*Q1* 
      pow(Q2, 2)+(-1.0)*P1*pow(Q2, 3))
    ) + 
    (s5L - s5L0) * ( 
     (-1.0/4.0)*pow(G, -2.0)*(pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)*pow(P1, 2)* 
      (pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))+P1*P2*(6.0*pow(Q1, 2)*Q2+(-2.0)*pow(Q2, 3)))
    ); 

    /* IP2 */
    I[5] = (Lf - L0) * ( 
     (-3.0)*pow(G, -2.0)*((-5.0)*P1*P2*Q1+Q2+4.0*pow(P1, 2)*Q2+(-1.0)*pow(P2, 2)*Q2) 
      *(pow(G, 2)+(-5.0)*(pow(Q1, 2)+pow(Q2, 2))) 
        ) + 
    (cL - cL0) * ( 
    (-1.0/16.0)*pow(G, -2.0)*(6.0*pow(G, 2)*(21.0*(2.0+pow(P1, 2))*P2*Q1+7.0*pow(P2, 3)*Q1+( 
      -1.0)*P1*(62.0+19.0*pow(P1, 2))*Q2+(-5.0)*P1*pow(P2, 2)*Q2)+5.0*((-15.0)*P1* 
      pow(P2, 2)*Q2*((-11.0)*pow(Q1, 2)+pow(Q2, 2))+(-3.0)*P2*Q1*((88.0+29.0*pow(P1, 2))* 
      pow(Q1, 2)+9.0*(8.0+9.0*pow(P1, 2))*pow(Q2, 2))+P1*Q2*((408.0+81.0*pow(P1, 2))*pow(Q1, 2)+5.0*( 
      72.0+25.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*((-59.0)*pow(Q1, 3)+9.0*Q1*pow(Q2, 2))))
    ) + 
    (c2L - c2L0) * ( 
     (1.0/2.0)*pow(G, -2.0)*((-3.0)*pow(G, 2)*((3.0+4.0*pow(P1, 2)+5.0*pow(P2, 2))*Q1+(-3.0)*P1* 
      P2*Q2)+5.0*Q1*((10+9.0*pow(P1, 2)+21.0*pow(P2, 2))*pow(Q1, 2)+(-36.0)*P1*P2*Q1* 
      Q2+3.0*(2.0+7.0*pow(P1, 2)+(-1.0)*pow(P2, 2))*pow(Q2, 2)))
    ) + 
    (c3L - c3L0) * ( 
    (-1.0/16.0)*pow(G, -2.0)*(pow(G, 2)*((68.0+11.0*pow(P1, 2))*P2*Q1+19.0*pow(P2, 3)*Q1+ 
      P1*(68.0+21.0*pow(P1, 2))*Q2+5.0*P1*pow(P2, 2)*Q2)+5.0*(P1*pow(P2, 2)*Q2*(19.0* 
      pow(Q1, 2)+(-13.0)*pow(Q2, 2))+P1*Q2*((28.0+3.0*pow(P1, 2))*pow(Q1, 2)+(-1.0)*(100+29.0* 
      pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*((-33.0)*pow(Q1, 3)+23.0*Q1*pow(Q2, 2))+P2*((-1.0)*( 
      124.0+25.0*pow(P1, 2))*pow(Q1, 3)+(100+31.0*pow(P1, 2))*Q1*pow(Q2, 2))))
    ) + 
    (c4L - c4L0) * ( 
    (1.0/4.0)*pow(G, -2.0)*(6.0*pow(G, 2)*(pow(P1, 2)*Q1+(-1.0)*pow(P2, 2)*Q1+(-2.0)*P1* 
      P2*Q2)+(-5.0)*((-5.0)*(1.0+3.0*pow(P2, 2))*pow(Q1, 3)+(-21.0)*P1*P2*pow(Q1, 2)*Q2+ 
      3.0*(5.0+8.0*pow(P1, 2)+7.0*pow(P2, 2))*Q1*pow(Q2, 2)+(-9.0)*P1*P2*pow(Q2, 3))) 
    ) + 
    (c5L - c5L0) * ( 
    (1.0/16.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(3.0*pow(P1, 2)*P2*Q1+(-1.0)*pow(P2, 3)*Q1+ 
      pow(P1, 3)*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+3.0*P1*pow(P2, 2)*Q2*(57.0*pow(Q1, 2)+pow(Q2, 2)) 
      +(-3.0)*P2*Q1*(((-36.0)+11.0*pow(P1, 2))*pow(Q1, 2)+27.0*(4.0+pow(P1, 2))*pow(Q2, 2))+ 
      P1*Q2*((324.0+51.0*pow(P1, 2))*pow(Q1, 2)+(-1.0)*(108.0+37.0*pow(P1, 2))*pow(Q2, 2))+ 
      pow(P2, 3)*(47.0*pow(Q1, 3)+(-81.0)*Q1*pow(Q2, 2))) 
    ) + 
    (c6L - c6L0) * ( 
    (-5.0/2.0)*pow(G, -2.0)*((-1.0)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+2.0*P1*P2* 
      Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))
    ) + 
    (c7L - c7L0) * ( 
     (-5.0/16.0)*pow(G, -2.0)*(3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)* 
      pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))
    ) + 
     (sL - sL0) * ( 
     (3.0/16.0)*pow(G, -2.0)*(2.0*pow(G, 2)*(5.0*pow(P1, 3)*Q1+P1*(34.0+19.0*pow(P2, 2))*Q1+( 
      -7.0)*pow(P1, 2)*P2*Q2+7.0*P2*(2.0+pow(P2, 2))*Q2)+(-5.0)*((-1.0)*pow(P1, 2)*P2* 
      Q2*(47.0*pow(Q1, 2)+3.0*pow(Q2, 2))+P1*Q1*((72.0+51.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(( 
      -56.0)+pow(P2, 2))*pow(Q2, 2))+P2*Q2*((40+29.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+3.0*pow(P2, 2)) 
      *pow(Q2, 2))+pow(P1, 3)*(7.0*pow(Q1, 3)+19.0*Q1*pow(Q2, 2))))
    ) + 
    (s2L - s2L0) * ( 
     (1.0/2.0)*pow(G, -2.0)*(3.0*pow(G, 2)*((-1.0)*P1*P2*Q1+6.0*pow(P1, 2)*Q2+3.0*(1.0+ 
      pow(P2, 2))*Q2)+(-5.0)*(6.0*P1*P2*Q1*(pow(Q1, 2)+(-5.0)*pow(Q2, 2))+Q2*(3.0*(4.0+ 
      9.0*pow(P2, 2))*pow(Q1, 2)+(8.0+3.0*pow(P2, 2))*pow(Q2, 2))+3.0*pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+7.0* 
      pow(Q2, 3)))) 
    ) + 
    (s3L - s3L0) * ( 
      (1.0/16.0)*pow(G, -2.0)*(pow(G, 2)*((-15.0)*pow(P1, 3)*Q1+(-1.0)*P1*(68.0+23.0*pow(P2, 2)) 
      *Q1+29.0*pow(P1, 2)*P2*Q2+P2*(68.0+13.0*pow(P2, 2))*Q2)+5.0*((-1.0)*pow(P1, 2)* 
      P2*Q2*(71.0*pow(Q1, 2)+15.0*pow(Q2, 2))+P2*Q2*((-1.0)*(236.0+55.0*pow(P2, 2))* 
      pow(Q1, 2)+((-12.0)+pow(P2, 2))*pow(Q2, 2))+P1*Q1*(3.0*(12.0+5.0*pow(P2, 2))*pow(Q1, 2)+(164.0+ 
      47.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(7.0*pow(Q1, 3)+39.0*Q1*pow(Q2, 2))))
    ) + 
    (s4L - s4L0) * ( 
     (1.0/4.0)*pow(G, -2.0)*((-6.0)*pow(G, 2)*(2.0*P1*P2*Q1+pow(P1, 2)*Q2+(-1.0)* 
      pow(P2, 2)*Q2)+5.0*((-3.0)*(5.0+11.0*pow(P2, 2))*pow(Q1, 2)*Q2+(5.0+3.0*pow(P2, 2))*pow(Q2, 3)+ 
      12.0*pow(P1, 2)*Q2*((-1.0)*pow(Q1, 2)+pow(Q2, 2))+3.0*P1*P2*Q1*(5.0*pow(Q1, 2)+pow(Q2, 2)) 
      ))
    ) + 
    (s5L - s5L0) * ( 
    (1.0/16.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(pow(P1, 3)*Q1+(-3.0)*P1*pow(P2, 2)*Q1+(-3.0)* 
      pow(P1, 2)*P2*Q2+pow(P2, 3)*Q2)+3.0*P1*Q1*((36.0+29.0*pow(P2, 2))*pow(Q1, 2)+(-27.0)* 
      (4.0+pow(P2, 2))*pow(Q2, 2))+P2*Q2*((-3.0)*(108.0+37.0*pow(P2, 2))*pow(Q1, 2)+(108.0+17.0* 
      pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(7.0*pow(Q1, 3)+(-81.0)*Q1*pow(Q2, 2))+pow(P1, 2)*P2*(9.0* 
      pow(Q1, 2)*Q2+57.0*pow(Q2, 3))) 
    ) + 
    (s6L - s6L0) * ( 
    (5.0/2.0)*pow(G, -2.0)*(2.0*P1*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+pow(P2, 2)*Q2*(( 
      -3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3))) 
    ) + 
    (s7L - s7L0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*((-3.0)*P1*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0* 
      pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+ 
      pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))) 
    );

    /* IQ1 */
    I[6] = (Lf - L0) * ( 
      pow(G, -2.0)*(pow(G, 2)*P1+(-5.0)*P1*pow(Q1, 2)+10.0*P2*Q1*Q2+(-15.0)*P1* 
      pow(Q2, 2)) 
        ) + 
    (cL - cL0) * ( 
    (-1.0/4.0)*pow(G, -2.0)*(pow(G, 2)*(4.0+3.0*pow(P1, 2)+pow(P2, 2))+(-10)*((2.0+pow(P1, 2)+pow(P2, 2)) 
      *pow(Q1, 2)+(-4.0)*P1*P2*Q1*Q2+(6.0+5.0*pow(P1, 2)+pow(P2, 2))*pow(Q2, 2)))
    ) + 
    (c2L - c2L0) * ( 
     (-1.0/2.0)*pow(G, -2.0)*(pow(G, 2)*P2+(-10)*((-2.0)*P1*Q1*Q2+P2*(pow(Q1, 2)+ 
      pow(Q2, 2)))) 
    ) + 
    (c3L - c3L0) * ( 
    (1.0/12.0)*pow(G, -2.0)*(pow(G, 2)*(pow(P1, 2)+(-1.0)*pow(P2, 2))+5.0*((4.0+pow(P1, 2)+3.0*pow(P2, 2)) 
      *pow(Q1, 2)+(-4.0)*P1*P2*Q1*Q2+((-4.0)+(-5.0)*pow(P1, 2)+pow(P2, 2))*pow(Q2, 2)))
    ) + 
    (c4L - c4L0) * ( 
     (5.0/4.0)*pow(G, -2.0)*(2.0*P1*Q1*Q2+P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)))
    ) + 
    (c5L - c5L0) * ( 
     (1.0/4.0)*pow(G, -2.0)*(4.0*P1*P2*Q1*Q2+pow(P2, 2)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
      pow(P1, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))
    ) + 
     (sL - sL0) * ( 
     (1.0/2.0)*pow(G, -2.0)*(pow(G, 2)*P1*P2+10.0*(pow(P1, 2)*Q1*Q2+(2.0+pow(P2, 2))*Q1* 
      Q2+(-1.0)*P1*P2*(pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (s2L - s2L0) * ( 
    P1*((-1.0/2.0)+10.0*pow(G, -2.0)*pow(Q2, 2))
    ) + 
    (s3L - s3L0) * ( 
     (-1.0/6.0)*pow(G, -2.0)*(pow(G, 2)*P1*P2+5.0*(3.0*pow(P1, 2)*Q1*Q2+(4.0+pow(P2, 2))* 
      Q1*Q2+(-1.0)*P1*P2*(pow(Q1, 2)+3.0*pow(Q2, 2))))
    ) + 
    (s4L - s4L0) * ( 
    (-5.0/4.0)*pow(G, -2.0)*(2.0*P2*Q1*Q2+P1*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))
    ) + 
    (s5L - s5L0) * ( 
     (1.0/2.0)*pow(G, -2.0)*(pow(P1, 2)*Q1*Q2+(-1.0)*pow(P2, 2)*Q1*Q2+P1*P2*(pow(Q1, 2)+ 
      (-1.0)*pow(Q2, 2)))
    );

    /*IQ2 */
    I[7] = (Lf - L0) * ( 
      pow(G, -2.0)*(pow(G, 2)*P2+(-15.0)*P2*pow(Q1, 2)+10.0*P1*Q1*Q2+(-5.0)*P2* 
      pow(Q2, 2)) 
        ) + 
    (cL - cL0) * ( 
    (-1.0/2.0)*pow(G, -2.0)*(pow(G, 2)*P1*P2+10.0*(pow(P1, 2)*Q1*Q2+(2.0+pow(P2, 2))*Q1* 
      Q2+(-1.0)*P1*P2*(pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (c2L - c2L0) * ( 
     (-1.0/2.0)*pow(G, -2.0)*(pow(G, 2)*P1+(-10)*((-2.0)*P2*Q1*Q2+P1*(pow(Q1, 2)+ 
      pow(Q2, 2)))) 
    ) + 
    (c3L - c3L0) * ( 
    (-1.0/6.0)*pow(G, -2.0)*(pow(G, 2)*P1*P2+5.0*(pow(P1, 2)*Q1*Q2+(4.0+3.0*pow(P2, 2))* 
      Q1*Q2+(-1.0)*P1*P2*(3.0*pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (c4L - c4L0) * ( 
     (-5.0/4.0)*pow(G, -2.0)*(2.0*P2*Q1*Q2+P1*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))
    ) + 
    (c5L - c5L0) * ( 
    (1.0/2.0)*pow(G, -2.0)*(pow(P1, 2)*Q1*Q2+(-1.0)*pow(P2, 2)*Q1*Q2+P1*P2*(pow(Q1, 2)+ 
      (-1.0)*pow(Q2, 2)))
    ) + 
     (sL - sL0) * ( 
     (1.0/4.0)*pow(G, -2.0)*(pow(G, 2)*(4.0+pow(P1, 2)+3.0*pow(P2, 2))+(-10)*((6.0+pow(P1, 2)+5.0* 
      pow(P2, 2))*pow(Q1, 2)+(-4.0)*P1*P2*Q1*Q2+(2.0+pow(P1, 2)+pow(P2, 2))*pow(Q2, 2)))
    ) + 
    (s2L - s2L0) * ( 
    P2*((1.0/2.0)+(-10)*pow(G, -2.0)*pow(Q1, 2))
    ) + 
    (s3L - s3L0) * ( 
     (1.0/12.0)*pow(G, -2.0)*(pow(G, 2)*((-1.0)*pow(P1, 2)+pow(P2, 2))+5.0*(((-4.0)+pow(P1, 2)+(-5.0)* 
      pow(P2, 2))*pow(Q1, 2)+(-4.0)*P1*P2*Q1*Q2+(4.0+3.0*pow(P1, 2)+pow(P2, 2))*pow(Q2, 2)))
    ) + 
    (s4L - s4L0) * ( 
    (-5.0/4.0)*pow(G, -2.0)*(2.0*P1*Q1*Q2+P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)))
    ) + 
    (s5L - s5L0) * ( 
     (1.0/4.0)*pow(G, -2.0)*((-4.0)*P1*P2*Q1*Q2+pow(P1, 2)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
      pow(P2, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))
    );

    return I;
}


std::vector<double>  perturbation_propagator::integralsJ4(const double &L0, const double &Lf, const double &P1, const double &P2, const double &Q1, const double &Q2) const{

    std::vector<double> I(7, 0.0);

    /* Pre-computations */
    double G = 1.0 + Q1 * Q1 + Q2 * Q2;
    double pi = constants::pi;
    double L_0 = L0 - 2.0 * pi * floor(L0 / (2.0 * pi));
    double L_f = Lf - 2.0 * pi * floor(Lf / (2.0 * pi));
    double sL = sin(L_f), sL0 = sin(L_0);
    double cL = cos(L_f), cL0 = cos(L_0);
    double c2L = cos(2.0 * L_f), c2L0 = cos(2.0 * L_0);
    double s2L = sin(2.0 * L_f), s2L0 = sin(2.0 * L_0);
    double s3L = sin(3.0 * L_f), s3L0 = sin(3.0 * L_0);
    double c3L = cos(3.0 * L_f), c3L0 = cos(3.0 * L_0);
    double c4L = cos(4.0 * L_f), c4L0 = cos(4.0 * L_0);
    double s4L = sin(4.0 * L_f), s4L0 = sin(4.0 * L_0);
    double s5L = sin(5.0 * L_f), s5L0 = sin(5.0 * L_0);
    double c5L = cos(5.0 * L_f), c5L0 = cos(5.0 * L_0);
    double s6L = sin(6.0 * L_f), s6L0 = sin(6.0 * L_0);
    double c6L = cos(6.0 * L_f), c6L0 = cos(6.0 * L_0);
    double s7L = sin(7.0 * L_f), s7L0 = sin(7.0 * L_0);
    double c7L = cos(7.0 * L_f), c7L0 = cos(7.0 * L_0);
    double s8L = sin(8.0 * L_f), s8L0 = sin(8.0 * L_0);
    double c8L = cos(8.0 * L_f), c8L0 = cos(8.0 * L_0);
    double s9L = sin(9.0 * L_f), s9L0 = sin(9.0 * L_0);
    double c9L = cos(9.0 * L_f), c9L0 = cos(9.0 * L_0);    

I[0] = (Lf - L0) * ( 
    50.0*pow(G, -4.0)*(pow(P1, 2)*Q1*Q2+(-1.0)*pow(P2, 2)*Q1*Q2+P1*P2*((-1.0)* 
    pow(Q1, 2)+pow(Q2, 2)))*((-6.0)*pow(G, 2)*(2.0+pow(P1, 2)+pow(P2, 2))+7.0*((8.0+3.0*pow(P1, 2)+5.0* 
    pow(P2, 2))*pow(Q1, 2)+(-4.0)*P1*P2*Q1*Q2+(8.0+5.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q2, 2)))    
    ) + 
    (cL - cL0) * ( 
    (-5.0/8.0)*pow(G, -4.0)*(3.0*pow(G, 4)*P2*(8.0+pow(P1, 4)+12.0*pow(P2, 2)+pow(P2, 4)+2.0*pow(P1, 2)* 
    (6.0+pow(P2, 2)))+(-15.0)*pow(G, 2)*(2.0*P1*(16.0+48.0*pow(P1, 2)+5.0*pow(P1, 4))*Q1*Q2+( 
    -4.0)*P1*(24.0+pow(P1, 2))*pow(P2, 2)*Q1*Q2+(-14.0)*P1*pow(P2, 4)*Q1*Q2+ 
    pow(P2, 5)*(5.0*pow(Q1, 2)+3.0*pow(Q2, 2))+pow(P2, 3)*((-2.0)*((-24.0)+pow(P1, 2))*pow(Q1, 2)+6.0*( 
    8.0+3.0*pow(P1, 2))*pow(Q2, 2))+P2*((16.0+(-48.0)*pow(P1, 2)+(-7.0)*pow(P1, 4))*pow(Q1, 2)+3.0*( 
    16.0+48.0*pow(P1, 2)+5.0*pow(P1, 4))*pow(Q2, 2)))+35.0*(8.0*P1*pow(P2, 2)*Q1*Q2*(3.0*(( 
    -2.0)+pow(P1, 2))*pow(Q1, 2)+(-1.0)*(42.0+5.0*pow(P1, 2))*pow(Q2, 2))+4.0*P1*Q1*Q2*((16.0+ 
    36.0*pow(P1, 2)+3.0*pow(P1, 4))*pow(Q1, 2)+(16.0+60.0*pow(P1, 2)+7.0*pow(P1, 4))*pow(Q2, 2))+(-4.0)* 
    P1*pow(P2, 4)*Q2*(5.0*pow(Q1, 3)+9.0*Q1*pow(Q2, 2))+pow(P2, 5)*(7.0*pow(Q1, 4)+18.0* 
    pow(Q1, 2)*pow(Q2, 2)+3.0*pow(Q2, 4))+pow(P2, 3)*((-10.0)*((-6.0)+pow(P1, 2))*pow(Q1, 4)+36.0*(6.0+ 
    pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+30.0*(2.0+pow(P1, 2))*pow(Q2, 4))+P2*((16.0+(-84.0)*pow(P1, 2)+( 
    -9.0)*pow(P1, 4))*pow(Q1, 4)+(-6.0)*((-16.0)+12.0*pow(P1, 2)+5.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+ 
    5.0*(16.0+60.0*pow(P1, 2)+7.0*pow(P1, 4))*pow(Q2, 4)))) 
    ) + 
    (c2L - c2L0) * ( 
    (5.0/2.0)*pow(G, -4.0)*(3.0*pow(G, 4)*(2.0*pow(P1, 2)+pow(P1, 4)+(-1.0)*pow(P2, 2)*(2.0+pow(P2, 2)))+ 
    (-15.0)*pow(G, 2)*((-8.0)*pow(P1, 3)*P2*Q1*Q2+8.0*P1*pow(P2, 3)*Q1*Q2+pow(P1, 4)* 
    (3.0*pow(Q1, 2)+5.0*pow(Q2, 2))+pow(P1, 2)*((8.0+6.0*pow(P2, 2))*pow(Q1, 2)+2.0*(4.0+(-3.0)*pow(P2, 2)) 
    *pow(Q2, 2))+(-1.0)*pow(P2, 2)*((8.0+5.0*pow(P2, 2))*pow(Q1, 2)+(8.0+3.0*pow(P2, 2))*pow(Q2, 2)))+ 
    35.0*((-8.0)*pow(P1, 3)*P2*Q1*Q2*(3.0*pow(Q1, 2)+pow(Q2, 2))+8.0*P1*P2*Q1*Q2* 
    (((-2.0)+pow(P2, 2))*pow(Q1, 2)+(2.0+3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(3.0*pow(Q1, 4)+18.0* 
    pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+2.0*pow(P1, 2)*((5.0+6.0*pow(P2, 2))*pow(Q1, 4)+18.0*pow(Q1, 2)* 
    pow(Q2, 2)+(5.0+(-6.0)*pow(P2, 2))*pow(Q2, 4))+(-1.0)*pow(P2, 2)*((10.0+7.0*pow(P2, 2))*pow(Q1, 4)+ 
    18.0*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(10.0+3.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c3L - c3L0) * ( 
    (5.0/48.0)*pow(G, -4.0)*(9.0*pow(G, 4)*P2*(3.0*pow(P1, 4)+(-1.0)*pow(P2, 2)*(8.0+pow(P2, 2))+ 
    2.0*pow(P1, 2)*(12.0+pow(P2, 2)))+30.0*pow(G, 2)*(2.0*P1*(16.0+24.0*pow(P1, 2)+pow(P1, 4))*Q1* 
    Q2+4.0*P1*(12.0+7.0*pow(P1, 2))*pow(P2, 2)*Q1*Q2+(-6.0)*P1*pow(P2, 4)*Q1*Q2+( 
    -6.0)*((-4.0)+pow(P1, 2))*pow(P2, 3)*(3.0*pow(Q1, 2)+pow(Q2, 2))+3.0*pow(P2, 5)*(3.0*pow(Q1, 2)+ 
    pow(Q2, 2))+(-1.0)*P2*(((-16.0)+120.0*pow(P1, 2)+11.0*pow(P1, 4))*pow(Q1, 2)+(16.0+168.0* 
    pow(P1, 2)+25.0*pow(P1, 4))*pow(Q2, 2)))+(-140.0)*((-4.0)*P1*pow(P2, 4)*Q1*Q2*( 
    pow(Q1, 2)+2.0*pow(Q2, 2))+8.0*P1*pow(P2, 2)*Q1*Q2*((9.0+4.0*pow(P1, 2))*pow(Q1, 2)+3.0*(1.0+ 
    pow(P1, 2))*pow(Q2, 2))+4.0*P1*Q1*Q2*((12.0+18.0*pow(P1, 2)+pow(P1, 4))*pow(Q1, 2)+2.0*(2.0+ 
    3.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 5)*(7.0*pow(Q1, 4)+12.0*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-2.0)* 
    pow(P2, 3)*(((-27.0)+8.0*pow(P1, 2))*pow(Q1, 4)+6.0*((-9.0)+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(( 
    -3.0)+2.0*pow(P1, 2))*pow(Q2, 4))+(-1.0)*P2*(((-12.0)+90.0*pow(P1, 2)+7.0*pow(P1, 4))* 
    pow(Q1, 4)+12.0*((-2.0)+15.0*pow(P1, 2)+2.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(20.0+138.0*pow(P1, 2)+ 
    21.0*pow(P1, 4))*pow(Q2, 4))))
    ) + 
    (c4L - c4L0) * ( 
    (-5.0/8.0)*pow(G, -4.0)*(3.0*pow(G, 4)*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+(-120.0) 
    *pow(G, 2)*(2.0*pow(P1, 3)*P2*Q1*Q2+2.0*P1*P2*(2.0+pow(P2, 2))*Q1*Q2+pow(P1, 4)* 
    pow(Q2, 2)+pow(P2, 2)*((1.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-1.0)*pow(P1, 2)*((1.0+3.0* 
    pow(P2, 2))*pow(Q1, 2)+((-1.0)+3.0*pow(P2, 2))*pow(Q2, 2)))+(-70.0)*(pow(Q1, 2)+pow(Q2, 2))*(( 
    -16.0)*pow(P1, 3)*P2*Q1*Q2+(-16.0)*P1*P2*(2.0+pow(P2, 2))*Q1*Q2+pow(P1, 4)*( 
    pow(Q1, 2)+(-7.0)*pow(Q2, 2))+pow(P2, 2)*((-1.0)*(8.0+7.0*pow(P2, 2))*pow(Q1, 2)+(8.0+pow(P2, 2))* 
    pow(Q2, 2))+2.0*pow(P1, 2)*((4.0+9.0*pow(P2, 2))*pow(Q1, 2)+((-4.0)+9.0*pow(P2, 2))*pow(Q2, 2))))
    ) + 
    (c5L - c5L0) * ( 
    (1.0/16.0)*pow(G, -4.0)*((-3.0)*pow(G, 4)*(5.0*pow(P1, 4)*P2+(-10.0)*pow(P1, 2)*pow(P2, 3)+ 
    pow(P2, 5))+(-30.0)*pow(G, 2)*(6.0*pow(P1, 5)*Q1*Q2+(-12.0)*pow(P1, 3)*((-4.0)+pow(P2, 2))* 
    Q1*Q2+(-18.0)*P1*pow(P2, 2)*(8.0+pow(P2, 2))*Q1*Q2+(-1.0)*pow(P1, 4)*P2*( 
    pow(Q1, 2)+19.0*pow(Q2, 2))+pow(P2, 3)*((-1.0)*(24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(24.0+pow(P2, 2))* 
    pow(Q2, 2))+2.0*pow(P1, 2)*P2*((36.0+13.0*pow(P2, 2))*pow(Q1, 2)+((-36.0)+7.0*pow(P2, 2))* 
    pow(Q2, 2)))+140.0*((-8.0)*P1*pow(P2, 2)*Q1*Q2*((21.0+2.0*pow(P1, 2))*pow(Q1, 2)+(15.0+ 
    pow(P1, 2))*pow(Q2, 2))+4.0*P1*Q1*Q2*(((-4.0)+6.0*pow(P1, 2)+pow(P1, 4))*pow(Q1, 2)+2.0*(2.0+ 
    9.0*pow(P1, 2)+pow(P1, 4))*pow(Q2, 2))+(-4.0)*P1*pow(P2, 4)*Q2*(5.0*pow(Q1, 3)+4.0*Q1* 
    pow(Q2, 2))+pow(P2, 5)*((-5.0)*pow(Q1, 4)+pow(Q2, 4))+2.0*pow(P2, 3)*(5.0*((-3.0)+2.0*pow(P1, 2))* 
    pow(Q1, 4)+18.0*(1.0+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(9.0+4.0*pow(P1, 2))*pow(Q2, 4))+P2*(((-4.0)+ 
    66.0*pow(P1, 2)+pow(P1, 4))*pow(Q1, 4)+(-12.0)*((-2.0)+(-3.0)*pow(P1, 2)+pow(P1, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+(-1.0)*(4.0+78.0*pow(P1, 2)+17.0*pow(P1, 4))*pow(Q2, 4))))
    ) + 
    (c6L - c6L0) * ( 
    (25.0/6.0)*pow(G, -4.0)*(3.0*pow(G, 2)*((-8.0)*pow(P1, 3)*P2*Q1*Q2+8.0*P1*pow(P2, 3)* 
    Q1*Q2+pow(P1, 4)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0* 
    pow(P1, 2)*pow(P2, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+(-7.0)*((-8.0)*pow(P1, 3)*P2*Q1* 
    Q2*(pow(Q1, 2)+3.0*pow(Q2, 2))+8.0*P1*P2*Q1*Q2*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+((-2.0) 
    +pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+6.0*pow(Q1, 2)*pow(Q2, 2)+(-3.0)*pow(Q2, 4))+(-2.0)* 
    pow(P1, 2)*((1.0+6.0*pow(P2, 2))*pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+(1.0+(-6.0)*pow(P2, 2))* 
    pow(Q2, 4))+pow(P2, 2)*((2.0+3.0*pow(P2, 2))*pow(Q1, 4)+(-6.0)*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    -1.0)*((-2.0)+pow(P2, 2))*pow(Q2, 4)))) 
    ) + 
    (c7L - c7L0) * ( 
    (25.0/112.0)*pow(G, -4.0)*(6.0*pow(G, 2)*(2.0*pow(P1, 5)*Q1*Q2+(-20.0)*pow(P1, 3)* 
    pow(P2, 2)*Q1*Q2+10.0*P1*pow(P2, 4)*Q1*Q2+5.0*pow(P1, 4)*P2*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+(-10.0)*pow(P1, 2)*pow(P2, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 5)*(pow(Q1, 2)+(-1.0) 
    *pow(Q2, 2)))+7.0*(4.0*pow(P1, 5)*Q1*Q2*(pow(Q1, 2)+(-5.0)*pow(Q2, 2))+(-4.0)*P1* 
    pow(P2, 2)*Q1*Q2*((72.0+19.0*pow(P2, 2))*pow(Q1, 2)+((-72.0)+pow(P2, 2))*pow(Q2, 2))+8.0* 
    pow(P1, 3)*Q1*Q2*((12.0+7.0*pow(P2, 2))*pow(Q1, 2)+((-12.0)+13.0*pow(P2, 2))*pow(Q2, 2))+ 
    pow(P1, 4)*P2*((-11.0)*pow(Q1, 4)+(-54.0)*pow(Q1, 2)*pow(Q2, 2)+29.0*pow(Q2, 4))+2.0*pow(P1, 2)* 
    P2*((36.0+23.0*pow(P2, 2))*pow(Q1, 4)+(-18.0)*(12.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(36.0+( 
    -17.0)*pow(P2, 2))*pow(Q2, 4))+pow(P2, 3)*((-1.0)*(24.0+7.0*pow(P2, 2))*pow(Q1, 4)+18.0*(8.0+ 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-24.0)+pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c8L - c8L0) * ( 
    (-175.0/16.0)*pow(G, -4.0)*((-16.0)*pow(P1, 3)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)) 
    +16.0*P1*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+pow(P2, 4)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (c9L - c9L0) * ( 
    (-175.0/144.0)*pow(G, -4.0)*(4.0*pow(P1, 5)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-40.0) 
    *pow(P1, 3)*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+20.0*P1*pow(P2, 4)*Q1* 
    Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+5.0*pow(P1, 4)*P2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-10.0)*pow(P1, 2)*pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+ 
    pow(P2, 5)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (sL - sL0) * ( 
    (-5.0/8.0)*pow(G, -4.0)*(3.0*pow(G, 4)*P1*(8.0+pow(P1, 4)+12.0*pow(P2, 2)+pow(P2, 4)+2.0*pow(P1, 2)* 
    (6.0+pow(P2, 2)))+(-15.0)*pow(G, 2)*((-14.0)*pow(P1, 4)*P2*Q1*Q2+(-4.0)*pow(P1, 2)* 
    P2*(24.0+pow(P2, 2))*Q1*Q2+2.0*P2*(16.0+48.0*pow(P2, 2)+5.0*pow(P2, 4))*Q1*Q2+ 
    pow(P1, 5)*(3.0*pow(Q1, 2)+5.0*pow(Q2, 2))+2.0*pow(P1, 3)*(3.0*(8.0+3.0*pow(P2, 2))*pow(Q1, 2)+(-1.0) 
    *((-24.0)+pow(P2, 2))*pow(Q2, 2))+P1*(3.0*(16.0+48.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q1, 2)+( 
    16.0+(-48.0)*pow(P2, 2)+(-7.0)*pow(P2, 4))*pow(Q2, 2)))+35.0*((-8.0)*pow(P1, 2)*P2*Q1* 
    Q2*((42.0+5.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*((-2.0)+pow(P2, 2))*pow(Q2, 2))+4.0*P2*Q1* 
    Q2*((16.0+60.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 2)+(16.0+36.0*pow(P2, 2)+3.0*pow(P2, 4))* 
    pow(Q2, 2))+(-4.0)*pow(P1, 4)*P2*Q2*(9.0*pow(Q1, 3)+5.0*Q1*pow(Q2, 2))+pow(P1, 5)*(3.0* 
    pow(Q1, 4)+18.0*pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+2.0*pow(P1, 3)*(15.0*(2.0+pow(P2, 2))*pow(Q1, 4)+ 
    18.0*(6.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*((-6.0)+pow(P2, 2))*pow(Q2, 4))+P1*(5.0*( 
    16.0+60.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 4)+(-6.0)*((-16.0)+12.0*pow(P2, 2)+5.0*pow(P2, 4))* 
    pow(Q1, 2)*pow(Q2, 2)+(16.0+(-84.0)*pow(P2, 2)+(-9.0)*pow(P2, 4))*pow(Q2, 4))))
    ) + 
    (s2L - s2L0) * ( 
    (-5.0)*pow(G, -4.0)*(3.0*pow(G, 4)*P1*P2*(2.0+pow(P1, 2)+pow(P2, 2))+(-15.0)*pow(G, 2)*( 
    pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*Q2+pow(P2, 4)*Q1*Q2+2.0*pow(P1, 3)* 
    P2*(pow(Q1, 2)+3.0*pow(Q2, 2))+2.0*P1*P2*((4.0+3.0*pow(P2, 2))*pow(Q1, 2)+(4.0+pow(P2, 2))* 
    pow(Q2, 2)))+35.0*(4.0*pow(P1, 4)*Q1*pow(Q2, 3)+4.0*pow(P2, 2)*Q1*Q2*((1.0+pow(P2, 2))* 
    pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-4.0)*pow(P1, 2)*Q1*Q2*((1.0+3.0*pow(P2, 2))*pow(Q1, 2)+((-1.0) 
    +3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*P2*(3.0*pow(Q1, 4)+6.0*pow(Q1, 2)*pow(Q2, 2)+11.0*pow(Q2, 4)) 
    +P1*P2*((14.0+11.0*pow(P2, 2))*pow(Q1, 4)+6.0*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(14.0+3.0* 
    pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s3L - s3L0) * ( 
    (5.0/48.0)*pow(G, -4.0)*(9.0*pow(G, 4)*P1*(pow(P1, 4)+(-2.0)*pow(P1, 2)*((-4.0)+pow(P2, 2))+( 
    -3.0)*pow(P2, 2)*(8.0+pow(P2, 2)))+(-30.0)*pow(G, 2)*((-6.0)*pow(P1, 4)*P2*Q1*Q2+4.0* 
    pow(P1, 2)*P2*(12.0+7.0*pow(P2, 2))*Q1*Q2+2.0*P2*(16.0+24.0*pow(P2, 2)+pow(P2, 4))*Q1* 
    Q2+3.0*pow(P1, 5)*(pow(Q1, 2)+3.0*pow(Q2, 2))+(-6.0)*pow(P1, 3)*((-4.0)+pow(P2, 2))*(pow(Q1, 2)+ 
    3.0*pow(Q2, 2))+(-1.0)*P1*((16.0+168.0*pow(P2, 2)+25.0*pow(P2, 4))*pow(Q1, 2)+((-16.0)+120.0* 
    pow(P2, 2)+11.0*pow(P2, 4))*pow(Q2, 2)))+140.0*((-4.0)*pow(P1, 4)*P2*Q1*Q2*(2.0* 
    pow(Q1, 2)+pow(Q2, 2))+8.0*pow(P1, 2)*P2*Q1*Q2*(3.0*(1.0+pow(P2, 2))*pow(Q1, 2)+(9.0+4.0* 
    pow(P2, 2))*pow(Q2, 2))+4.0*P2*Q1*Q2*((4.0+6.0*pow(P2, 2))*pow(Q1, 2)+(12.0+18.0*pow(P2, 2)+ 
    pow(P2, 4))*pow(Q2, 2))+pow(P1, 5)*(pow(Q1, 4)+12.0*pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+(-2.0)* 
    pow(P1, 3)*(((-3.0)+2.0*pow(P2, 2))*pow(Q1, 4)+6.0*((-9.0)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(( 
    -27.0)+8.0*pow(P2, 2))*pow(Q2, 4))+(-1.0)*P1*((20.0+138.0*pow(P2, 2)+21.0*pow(P2, 4))* 
    pow(Q1, 4)+12.0*((-2.0)+15.0*pow(P2, 2)+2.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+((-12.0)+90.0* 
    pow(P2, 2)+7.0*pow(P2, 4))*pow(Q2, 4))))
    ) + 
    (s4L - s4L0) * ( 
    (5.0/2.0)*pow(G, -4.0)*(3.0*pow(G, 4)*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2))+(-70.0)*( 
    pow(Q1, 2)+pow(Q2, 2))*(4.0*pow(P1, 2)*Q1*Q2+2.0*pow(P1, 4)*Q1*Q2+(-2.0)*pow(P2, 2)*(2.0+ 
    pow(P2, 2))*Q1*Q2+(-1.0)*pow(P1, 3)*P2*(pow(Q1, 2)+5.0*pow(Q2, 2))+P1*P2*((4.0+5.0* 
    pow(P2, 2))*pow(Q1, 2)+((-4.0)+pow(P2, 2))*pow(Q2, 2)))+30.0*pow(G, 2)*(2.0*pow(P1, 2)*Q1*Q2+ 
    pow(P1, 4)*Q1*Q2+(-1.0)*pow(P2, 2)*(2.0+pow(P2, 2))*Q1*Q2+(-1.0)*pow(P1, 3)*P2*( 
    pow(Q1, 2)+3.0*pow(Q2, 2))+P1*P2*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+((-2.0)+pow(P2, 2))*pow(Q2, 2))) 
    ) 
    ) + 
    (s5L - s5L0) * ( 
    (1.0/16.0)*pow(G, -4.0)*((-3.0)*pow(G, 4)*(pow(P1, 5)+(-10.0)*pow(P1, 3)*pow(P2, 2)+5.0*P1* 
    pow(P2, 4))+(-30.0)*pow(G, 2)*((-18.0)*pow(P1, 4)*P2*Q1*Q2+6.0*pow(P2, 3)*(8.0+pow(P2, 2)) 
    *Q1*Q2+(-12.0)*pow(P1, 2)*P2*(12.0+pow(P2, 2))*Q1*Q2+pow(P1, 5)*(pow(Q1, 2)+(-5.0)* 
    pow(Q2, 2))+(-1.0)*P1*pow(P2, 2)*((72.0+19.0*pow(P2, 2))*pow(Q1, 2)+((-72.0)+pow(P2, 2))* 
    pow(Q2, 2))+2.0*pow(P1, 3)*((12.0+7.0*pow(P2, 2))*pow(Q1, 2)+((-12.0)+13.0*pow(P2, 2))*pow(Q2, 2)))+ 
    140.0*((-8.0)*pow(P1, 2)*P2*Q1*Q2*((15.0+pow(P2, 2))*pow(Q1, 2)+(21.0+2.0*pow(P2, 2))* 
    pow(Q2, 2))+4.0*P2*Q1*Q2*(2.0*(2.0+9.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)+((-4.0)+6.0* 
    pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+(-4.0)*pow(P1, 4)*P2*Q2*(4.0*pow(Q1, 3)+5.0*Q1*pow(Q2, 2))+ 
    pow(P1, 5)*(pow(Q1, 4)+(-5.0)*pow(Q2, 4))+2.0*pow(P1, 3)*((9.0+4.0*pow(P2, 2))*pow(Q1, 4)+18.0*(1.0+ 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*((-3.0)+2.0*pow(P2, 2))*pow(Q2, 4))+P1*((-1.0)*(4.0+78.0* 
    pow(P2, 2)+17.0*pow(P2, 4))*pow(Q1, 4)+(-12.0)*((-2.0)+(-3.0)*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+((-4.0)+66.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4)))) 
    ) + 
    (s6L - s6L0) * ( 
    (-25.0/3.0)*pow(G, -4.0)*(3.0*pow(G, 2)*(pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)* 
    Q1*Q2+pow(P2, 4)*Q1*Q2+2.0*pow(P1, 3)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*P1* 
    pow(P2, 3)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+(-7.0)*(4.0*pow(P1, 4)*Q1*pow(Q2, 3)+4.0*pow(P2, 2)* 
    Q1*Q2*((1.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-4.0)*pow(P1, 2)*Q1*Q2*((1.0+ 
    3.0*pow(P2, 2))*pow(Q1, 2)+((-1.0)+3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*P2*(3.0*pow(Q1, 4)+6.0* 
    pow(Q1, 2)*pow(Q2, 2)+(-5.0)*pow(Q2, 4))+P1*P2*((-1.0)*(2.0+5.0*pow(P2, 2))*pow(Q1, 4)+6.0*( 
    2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-2.0)+3.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s7L - s7L0) * ( 
    (25.0/112.0)*pow(G, -4.0)*(6.0*pow(G, 2)*((-10.0)*pow(P1, 4)*P2*Q1*Q2+20.0*pow(P1, 2)* 
    pow(P2, 3)*Q1*Q2+(-2.0)*pow(P2, 5)*Q1*Q2+pow(P1, 5)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-10.0) 
    *pow(P1, 3)*pow(P2, 2)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+5.0*P1*pow(P2, 4)*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2)))+(-7.0)*((-4.0)*pow(P1, 4)*P2*Q1*Q2*(pow(Q1, 2)+19.0*pow(Q2, 2))+4.0* 
    pow(P2, 3)*Q1*Q2*((-1.0)*(24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(24.0+pow(P2, 2))*pow(Q2, 2))+8.0* 
    pow(P1, 2)*P2*Q1*Q2*((36.0+13.0*pow(P2, 2))*pow(Q1, 2)+((-36.0)+7.0*pow(P2, 2))*pow(Q2, 2)) 
    +pow(P1, 5)*(pow(Q1, 4)+18.0*pow(Q1, 2)*pow(Q2, 2)+(-7.0)*pow(Q2, 4))+(-2.0)*pow(P1, 3)*((12.0+ 
    17.0*pow(P2, 2))*pow(Q1, 4)+18.0*((-4.0)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(12.0+(-23.0)*pow(P2, 2)) 
    *pow(Q2, 4))+P1*pow(P2, 2)*((72.0+29.0*pow(P2, 2))*pow(Q1, 4)+(-54.0)*(8.0+pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(72.0+(-11.0)*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s8L - s8L0) * ( 
    (175.0/4.0)*pow(G, -4.0)*(pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*Q1* 
    Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*Q1*Q2*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*P2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*P1* 
    pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (s9L - s9L0) * ( 
    (-175.0/144.0)*pow(G, -4.0)*((-20.0)*pow(P1, 4)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+40.0*pow(P1, 2)*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+4.0*pow(P2, 5)* 
    Q1*Q2*((-1.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 5)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-10.0)*pow(P1, 3)*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+5.0* 
    P1*pow(P2, 4)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    );


    I[1] = (Lf - L0) * ( 
    (25.0/4.0)*pow(G, -2.0)*(pow(P1, 2)*Q1*Q2+(-1.0)*pow(P2, 2)*Q1*Q2+P1*P2*((-1.0) 
    *pow(Q1, 2)+pow(Q2, 2)))*((-6.0)*pow(G, 2)*(2.0+pow(P1, 2)+pow(P2, 2))+7.0*((8.0+3.0*pow(P1, 2)+5.0* 
    pow(P2, 2))*pow(Q1, 2)+(-4.0)*P1*P2*Q1*Q2+(8.0+5.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q2, 2))) 
    ) + 
    (cL - cL0) * ( 
    (-5.0/32.0)*pow(G, -2.0)*(15.0*pow(G, 2)*((-2.0)*P1*(16.0+32.0*pow(P1, 2)+3.0*pow(P1, 4))* 
    Q1*Q2+(-4.0)*pow(P1, 3)*pow(P2, 2)*Q1*Q2+2.0*P1*pow(P2, 4)*Q1*Q2+(16.0+48.0* 
    pow(P1, 2)+5.0*pow(P1, 4))*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*(8.0+3.0*pow(P1, 2))*pow(P2, 3)*( 
    pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 5)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)))+14.0*(10.0*P1* 
    pow(P2, 4)*Q1*Q2*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+20.0*P1*pow(P2, 2)*Q1*Q2*(3.0*(4.0+ 
    pow(P1, 2))*pow(Q1, 2)+(-1.0)*(12.0+pow(P1, 2))*pow(Q2, 2))+2.0*P1*Q1*Q2*((80.0+120.0* 
    pow(P1, 2)+9.0*pow(P1, 4))*pow(Q1, 2)+(80.0+200.0*pow(P1, 2)+21.0*pow(P1, 4))*pow(Q2, 2))+pow(P2, 5)*(( 
    -7.0)*pow(Q1, 4)+12.0*pow(Q1, 2)*pow(Q2, 2)+3.0*pow(Q2, 4))+(-10.0)*pow(P2, 3)*((10.0+3.0*pow(P1, 2)) 
    *pow(Q1, 4)+(-12.0)*pow(Q1, 2)*pow(Q2, 2)+(-3.0)*(2.0+pow(P1, 2))*pow(Q2, 4))+(-5.0)*P2*(( 
    16.0+36.0*pow(P1, 2)+3.0*pow(P1, 4))*pow(Q1, 4)+12.0*pow(P1, 2)*(6.0+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    -1.0)*(16.0+60.0*pow(P1, 2)+7.0*pow(P1, 4))*pow(Q2, 4)))) 
    ) + 
    (c2L - c2L0) * ( 
    (-5.0/32.0)*pow(G, -2.0)*(3.0*pow(G, 2)*((16.0+25.0*pow(P1, 4)+80.0*pow(P2, 2)+25.0*pow(P2, 4)+ 
    10.0*pow(P1, 2)*(8.0+9.0*pow(P2, 2)))*pow(Q1, 2)+(-40.0)*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2)) 
    *Q1*Q2+(-1.0)*(16.0+25.0*pow(P1, 4)+80.0*pow(P2, 2)+25.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+9.0* 
    pow(P2, 2)))*pow(Q2, 2))+(-14.0)*((16.0+15.0*pow(P1, 4)+100.0*pow(P2, 2)+35.0*pow(P2, 4)+30.0* 
    pow(P1, 2)*(2.0+3.0*pow(P2, 2)))*pow(Q1, 4)+(-40.0)*P1*P2*(4.0+3.0*pow(P1, 2)+pow(P2, 2))* 
    pow(Q1, 3)*Q2+60.0*(2.0*pow(P1, 2)+pow(P1, 4)+(-1.0)*pow(P2, 2)*(2.0+pow(P2, 2)))*pow(Q1, 2)* 
    pow(Q2, 2)+40.0*P1*P2*(4.0+pow(P1, 2)+3.0*pow(P2, 2))*Q1*pow(Q2, 3)+(-1.0)*(16.0+35.0* 
    pow(P1, 4)+60.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(10.0+9.0*pow(P2, 2)))*pow(Q2, 4)))
    ) + 
    (c3L - c3L0) * ( 
    (-5.0/96.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(2.0*P1*(80.0+120.0*pow(P1, 2)+11.0*pow(P1, 4))* 
    Q1*Q2+20.0*P1*(12.0+pow(P1, 2))*pow(P2, 2)*Q1*Q2+30.0*P1*pow(P2, 4)*Q1*Q2+5.0* 
    (16.0+24.0*pow(P1, 2)+pow(P1, 4))*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+30.0*(4.0+pow(P1, 2))* 
    pow(P2, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+9.0*pow(P2, 5)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)))+(-28.0)*( 
    5.0*P1*pow(P2, 4)*Q1*Q2*(pow(Q1, 2)+5.0*pow(Q2, 2))+pow(P2, 5)*(7.0*pow(Q1, 4)+(-15.0)* 
    pow(Q1, 2)*pow(Q2, 2)+(-2.0)*pow(Q2, 4))+10.0*pow(P2, 3)*((9.0+2.0*pow(P1, 2))*pow(Q1, 4)+(-3.0)*( 
    6.0+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(3.0+pow(P1, 2))*pow(Q2, 4))+5.0*P2*((12.0+18.0* 
    pow(P1, 2)+pow(P1, 4))*pow(Q1, 4)+(-3.0)*(8.0+12.0*pow(P1, 2)+pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-2.0)* 
    (2.0+3.0*pow(P1, 2))*pow(Q2, 4))+(-10.0)*P1*pow(P2, 2)*Q1*Q2*((-24.0)*pow(Q2, 2)+ 
    pow(P1, 2)*(pow(Q1, 2)+(-3.0)*pow(Q2, 2)))+P1*Q1*Q2*(160.0*pow(Q2, 2)+240.0*pow(P1, 2)* 
    pow(Q2, 2)+pow(P1, 4)*(pow(Q1, 2)+21.0*pow(Q2, 2)))))
    ) + 
    (c4L - c4L0) * ( 
    (5.0/32.0)*pow(G, -2.0)*((-7.0)*(((-8.0)+5.0*pow(P1, 4)+(-80.0)*pow(P2, 2)+(-30.0)* 
    pow(P1, 2)*pow(P2, 2)+(-35.0)*pow(P2, 4))*pow(Q1, 4)+(-80.0)*P1*P2*(2.0+pow(P1, 2)+pow(P2, 2))* 
    pow(Q1, 3)*Q2+6.0*(8.0+15.0*pow(P1, 4)+40.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(4.0+3.0* 
    pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+(-80.0)*P1*P2*(2.0+pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 3)+( 
    -1.0)*(8.0+35.0*pow(P1, 4)+(-5.0)*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2)))*pow(Q2, 4))+ 
    30.0*pow(G, 2)*(2.0+pow(P1, 2)+pow(P2, 2))*((-4.0)*P1*P2*Q1*Q2+pow(P1, 2)*(pow(Q1, 2)+( 
    -1.0)*pow(Q2, 2))+pow(P2, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (c5L - c5L0) * ( 
    (5.0/32.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(8.0+pow(P1, 2)+pow(P2, 2))*(2.0*pow(P1, 3)*Q1*Q2+(-6.0) 
    *P1*pow(P2, 2)*Q1*Q2+3.0*pow(P1, 2)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 3)*((-1.0) 
    *pow(Q1, 2)+pow(Q2, 2)))+28.0*(P1*pow(P2, 4)*Q1*Q2*(5.0*pow(Q1, 2)+pow(Q2, 2))+P1*Q1* 
    Q2*((16.0+16.0*pow(P1, 2)+pow(P1, 4))*pow(Q1, 2)+(-1.0)*(16.0+32.0*pow(P1, 2)+3.0*pow(P1, 4))* 
    pow(Q2, 2))+pow(P2, 5)*(pow(Q1, 4)+(-3.0)*pow(Q1, 2)*pow(Q2, 2))+2.0*pow(P2, 3)*(5.0*pow(Q1, 4)+(-3.0) 
    *(6.0+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(1.0+pow(P1, 2))*pow(Q2, 4))+(-1.0)*P2*(((-4.0)+6.0* 
    pow(P1, 2)+pow(P1, 4))*pow(Q1, 4)+3.0*(8.0+12.0*pow(P1, 2)+pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-2.0)*(2.0+ 
    9.0*pow(P1, 2)+pow(P1, 4))*pow(Q2, 4))+pow(P2, 2)*(48.0*P1*pow(Q1, 3)*Q2+pow(P1, 3)*(6.0* 
    pow(Q1, 3)*Q2+(-2.0)*Q1*pow(Q2, 3)))))
    ) + 
    (c6L - c6L0) * ( 
    (-25.0/96.0)*pow(G, -2.0)*(3.0*pow(G, 2)*((-8.0)*pow(P1, 3)*P2*Q1*Q2+8.0*P1* 
    pow(P2, 3)*Q1*Q2+pow(P1, 4)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+14.0*((-8.0)*pow(P1, 3)* 
    P2*Q1*Q2*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-8.0)*P1*P2*Q1*Q2*((4.0+3.0* 
    pow(P2, 2))*pow(Q1, 2)+(-1.0)*(4.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+(-12.0)*pow(Q1, 2)* 
    pow(Q2, 2)+3.0*pow(Q2, 4))+pow(P1, 2)*((4.0+6.0*pow(P2, 2))*pow(Q1, 4)+(-24.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    2.0*(2.0+(-3.0)*pow(P2, 2))*pow(Q2, 4))+(-1.0)*pow(P2, 2)*((4.0+3.0*pow(P2, 2))*pow(Q1, 4)+(-12.0) 
    *(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(4.0+pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c7L - c7L0) * ( 
    (-5.0/224.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(2.0*pow(P1, 5)*Q1*Q2+(-20.0)*pow(P1, 3)* 
    pow(P2, 2)*Q1*Q2+10.0*P1*pow(P2, 4)*Q1*Q2+5.0*pow(P1, 4)*P2*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+(-10.0)*pow(P1, 2)*pow(P2, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 5)*(pow(Q1, 2)+(-1.0) 
    *pow(Q2, 2)))+7.0*(80.0*pow(P1, 3)*Q1*Q2*(2.0*pow(Q1, 2)+((-2.0)+pow(P2, 2))*pow(Q2, 2))+( 
    -40.0)*P1*pow(P2, 2)*Q1*Q2*(2.0*(6.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*(12.0+pow(P2, 2))* 
    pow(Q2, 2))+8.0*pow(P1, 5)*(2.0*pow(Q1, 3)*Q2+(-3.0)*Q1*pow(Q2, 3))+5.0*pow(P1, 4)*P2*( 
    pow(Q1, 4)+(-18.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+10.0*pow(P1, 2)*P2*(3.0*(4.0+pow(P2, 2))* 
    pow(Q1, 4)+(-6.0)*(12.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*((-12.0)+pow(P2, 2))*pow(Q2, 4))+( 
    -1.0)*pow(P2, 3)*((40.0+7.0*pow(P2, 2))*pow(Q1, 4)+(-30.0)*(8.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    40.0+3.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c8L - c8L0) * ( 
    (175.0/128.0)*pow(G, -2.0)*((-16.0)*pow(P1, 3)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)) 
    +16.0*P1*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+pow(P2, 4)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (c9L - c9L0) * ( 
    (35.0/288.0)*pow(G, -2.0)*(4.0*pow(P1, 5)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-40.0)* 
    pow(P1, 3)*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+20.0*P1*pow(P2, 4)*Q1*Q2*( 
    pow(Q1, 2)+(-1.0)*pow(Q2, 2))+5.0*pow(P1, 4)*P2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+ 
    (-10.0)*pow(P1, 2)*pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+pow(P2, 5)*( 
    pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (sL - sL0) * ( 
    (5.0/32.0)*pow(G, -2.0)*(15.0*pow(G, 2)*((-2.0)*pow(P1, 4)*P2*Q1*Q2+4.0*pow(P1, 2)* 
    pow(P2, 3)*Q1*Q2+2.0*P2*(16.0+32.0*pow(P2, 2)+3.0*pow(P2, 4))*Q1*Q2+pow(P1, 5)*( 
    pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*pow(P1, 3)*(8.0+3.0*pow(P2, 2))*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    P1*(16.0+48.0*pow(P2, 2)+5.0*pow(P2, 4))*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)))+(-14.0)*(10.0* 
    pow(P1, 4)*P2*Q1*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-20.0)*pow(P1, 2)*P2*Q1*Q2*( 
    (12.0+pow(P2, 2))*pow(Q1, 2)+(-3.0)*(4.0+pow(P2, 2))*pow(Q2, 2))+2.0*P2*Q1*Q2*((80.0+ 
    200.0*pow(P2, 2)+21.0*pow(P2, 4))*pow(Q1, 2)+(80.0+120.0*pow(P2, 2)+9.0*pow(P2, 4))*pow(Q2, 2))+ 
    pow(P1, 5)*(3.0*pow(Q1, 4)+12.0*pow(Q1, 2)*pow(Q2, 2)+(-7.0)*pow(Q2, 4))+10.0*pow(P1, 3)*(3.0*(2.0+ 
    pow(P2, 2))*pow(Q1, 4)+12.0*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(10.0+3.0*pow(P2, 2))*pow(Q2, 4))+5.0*P1* 
    ((16.0+60.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 4)+(-12.0)*pow(P2, 2)*(6.0+pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+(-1.0)*(16.0+36.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q2, 4)))) 
    ) + 
    (s2L - s2L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(35.0*pow(P1, 4)*Q1*Q2+10.0*pow(P1, 2)*(8.0+3.0* 
    pow(P2, 2))*Q1*Q2+(16.0+80.0*pow(P2, 2)+35.0*pow(P2, 4))*Q1*Q2+(-10.0)*pow(P1, 3)*P2* 
    (pow(Q1, 2)+(-1.0)*pow(Q2, 2))+10.0*P1*pow(P2, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)))+(-14.0)*( 
    20.0*pow(P1, 3)*P2*pow(Q2, 2)*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+10.0*pow(P1, 2)*Q1*Q2*((4.0+ 
    3.0*pow(P2, 2))*pow(Q1, 2)+3.0*(4.0+pow(P2, 2))*pow(Q2, 2))+Q1*Q2*((16.0+120.0*pow(P2, 2)+55.0* 
    pow(P2, 4))*pow(Q1, 2)+(16.0+40.0*pow(P2, 2)+15.0*pow(P2, 4))*pow(Q2, 2))+5.0*pow(P1, 4)*(3.0* 
    pow(Q1, 3)*Q2+11.0*Q1*pow(Q2, 3))+20.0*P1*P2*((1.0+pow(P2, 2))*pow(Q1, 4)+(-3.0)*(2.0+ 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))) 
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/96.0)*pow(G, -2.0)*(3.0*pow(G, 2)*((-30.0)*pow(P1, 4)*P2*Q1*Q2+(-20.0)* 
    pow(P1, 2)*P2*(12.0+pow(P2, 2))*Q1*Q2+(-2.0)*P2*(80.0+120.0*pow(P2, 2)+11.0*pow(P2, 4)) 
    *Q1*Q2+9.0*pow(P1, 5)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+30.0*pow(P1, 3)*(4.0+pow(P2, 2))*( 
    pow(Q1, 2)+(-1.0)*pow(Q2, 2))+5.0*P1*(16.0+24.0*pow(P2, 2)+pow(P2, 4))*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2)))+(-28.0)*((-5.0)*pow(P1, 4)*P2*Q1*Q2*(5.0*pow(Q1, 2)+pow(Q2, 2))+(-10.0)* 
    pow(P1, 2)*P2*Q1*Q2*(3.0*(8.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*pow(P2, 2)*pow(Q2, 2))+(-1.0)* 
    P2*Q1*Q2*((160.0+240.0*pow(P2, 2)+21.0*pow(P2, 4))*pow(Q1, 2)+pow(P2, 4)*pow(Q2, 2))+ 
    pow(P1, 5)*(2.0*pow(Q1, 4)+15.0*pow(Q1, 2)*pow(Q2, 2)+(-7.0)*pow(Q2, 4))+10.0*pow(P1, 3)*((3.0+ 
    pow(P2, 2))*pow(Q1, 4)+3.0*(6.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(9.0+2.0*pow(P2, 2))* 
    pow(Q2, 4))+5.0*P1*((4.0+6.0*pow(P2, 2))*pow(Q1, 4)+3.0*(8.0+12.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+(-1.0)*(12.0+18.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))))
    ) + 
    (s4L - s4L0) * ( 
    (-5.0/8.0)*pow(G, -2.0)*(15.0*pow(G, 2)*(2.0+pow(P1, 2)+pow(P2, 2))*(pow(P1, 2)*Q1*Q2+(-1.0)* 
    pow(P2, 2)*Q1*Q2+P1*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)))+7.0*(5.0*pow(P1, 4)*Q1*Q2*( 
    pow(Q1, 2)+(-5.0)*pow(Q2, 2))+10.0*pow(P1, 2)*Q1*Q2*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*( 
    2.0+pow(P2, 2))*pow(Q2, 2))+Q1*Q2*((8.0+60.0*pow(P2, 2)+25.0*pow(P2, 4))*pow(Q1, 2)+(-1.0)*(8.0+ 
    20.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q2, 2))+(-10.0)*pow(P1, 3)*P2*(pow(Q1, 4)+(-1.0)*pow(Q2, 4))+ 
    (-10.0)*P1*P2*(2.0+pow(P2, 2))*(pow(Q1, 4)+(-1.0)*pow(Q2, 4))))
    ) + 
    (s5L - s5L0) * ( 
    (5.0/32.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(8.0+pow(P1, 2)+pow(P2, 2))*((-6.0)*pow(P1, 2)*P2*Q1* 
    Q2+2.0*pow(P2, 3)*Q1*Q2+pow(P1, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*((-1.0) 
    *pow(Q1, 2)+pow(Q2, 2)))+(-28.0)*((-1.0)*pow(P1, 4)*P2*Q1*Q2*(pow(Q1, 2)+5.0*pow(Q2, 2))+ 
    P2*Q1*Q2*((16.0+32.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q1, 2)+(-1.0)*(16.0+16.0*pow(P2, 2)+ 
    pow(P2, 4))*pow(Q2, 2))+pow(P1, 5)*(3.0*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*pow(Q2, 4))+(-2.0)*pow(P1, 3)*( 
    (1.0+pow(P2, 2))*pow(Q1, 4)+(-3.0)*(6.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+P1*((-2.0) 
    *(2.0+9.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 4)+3.0*(8.0+12.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+( 
    (-4.0)+6.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+2.0*pow(P1, 2)*P2*Q1*Q2*((-24.0)*pow(Q2, 2)+ 
    pow(P2, 2)*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))))) 
    ) + 
    (s6L - s6L0) * ( 
    (25.0/48.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)* 
    Q1*Q2+pow(P2, 4)*Q1*Q2+2.0*pow(P1, 3)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*P1* 
    pow(P2, 3)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+14.0*(4.0*pow(P1, 3)*P2*pow(Q2, 2)*((-3.0)* 
    pow(Q1, 2)+pow(Q2, 2))+2.0*pow(P1, 2)*Q1*Q2*((4.0+3.0*pow(P2, 2))*pow(Q1, 2)+((-4.0)+3.0* 
    pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)*Q1*Q2*((-1.0)*(8.0+5.0*pow(P2, 2))*pow(Q1, 2)+(8.0+3.0* 
    pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(3.0*pow(Q1, 3)*Q2+(-5.0)*Q1*pow(Q2, 3))+4.0*P1*P2*(( 
    1.0+pow(P2, 2))*pow(Q1, 4)+(-3.0)*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))) 
    ) + 
    (s7L - s7L0) * ( 
    (-5.0/224.0)*pow(G, -2.0)*(3.0*pow(G, 2)*((-10.0)*pow(P1, 4)*P2*Q1*Q2+20.0*pow(P1, 2)* 
    pow(P2, 3)*Q1*Q2+(-2.0)*pow(P2, 5)*Q1*Q2+pow(P1, 5)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-10.0) 
    *pow(P1, 3)*pow(P2, 2)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+5.0*P1*pow(P2, 4)*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2)))+7.0*((-80.0)*pow(P1, 2)*P2*Q1*Q2*((6.0+pow(P2, 2))*pow(Q1, 2)+(-6.0)* 
    pow(Q2, 2))+(-40.0)*pow(P1, 4)*P2*Q1*Q2*(pow(Q1, 2)+(-2.0)*pow(Q2, 2))+8.0*pow(P2, 3)* 
    Q1*Q2*((20.0+3.0*pow(P2, 2))*pow(Q1, 2)+(-2.0)*(10.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 5)*(3.0* 
    pow(Q1, 4)+(-30.0)*pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+10.0*pow(P1, 3)*((4.0+pow(P2, 2))*pow(Q1, 4)+ 
    6.0*((-4.0)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(4.0+(-3.0)*pow(P2, 2))*pow(Q2, 4))+(-5.0)*P1* 
    pow(P2, 2)*((24.0+5.0*pow(P2, 2))*pow(Q1, 4)+(-18.0)*(8.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(24.0+ 
    pow(P2, 2))*pow(Q2, 4)))) 
    ) + 
    (s8L - s8L0) * ( 
    (-175.0/32.0)*pow(G, -2.0)*(pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)* 
    Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*Q1*Q2*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*P2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*P1* 
    pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))) 
    ) + 
    (s9L - s9L0) * ( 
    (35.0/288.0)*pow(G, -2.0)*((-20.0)*pow(P1, 4)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    40.0*pow(P1, 2)*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+4.0*pow(P2, 5)*Q1*Q2*(( 
    -1.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 5)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-10.0)* 
    pow(P1, 3)*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+5.0*P1*pow(P2, 4)*( 
    pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))) 
    );


    // P11
    I[2] = (Lf - L0) * ( 
    (15.0/16.0)*pow(G, -4.0)*(pow(G, 4)*P2*(4.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))+14.0*(pow(Q1, 2)+ 
    pow(Q2, 2))*(8.0*P1*((-1.0)+pow(P1, 2))*Q1*Q2+(-48.0)*P1*pow(P2, 2)*Q1*Q2+5.0* 
    pow(P2, 3)*(5.0*pow(Q1, 2)+pow(Q2, 2))+P2*((-3.0)*((-8.0)+pow(P1, 2))*pow(Q1, 2)+(16.0+33.0* 
    pow(P1, 2))*pow(Q2, 2)))+(-2.0)*pow(G, 2)*(12.0*P1*((-1.0)+pow(P1, 2))*Q1*Q2+(-72.0)* 
    P1*pow(P2, 2)*Q1*Q2+15.0*pow(P2, 3)*(3.0*pow(Q1, 2)+pow(Q2, 2))+P2*((46.0+3.0*pow(P1, 2))* 
    pow(Q1, 2)+(34.0+57.0*pow(P1, 2))*pow(Q2, 2))))  
    ) + 
    (cL - cL0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(6.0*pow(G, 4)*P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))+14.0*((-6.0)* 
    pow(P1, 2)*Q1*Q2*(((-22.0)+5.0*pow(P2, 2))*pow(Q1, 2)+((-26.0)+21.0*pow(P2, 2))*pow(Q2, 2))+ 
    (-1.0)*Q1*Q2*((16.0+228.0*pow(P2, 2)+33.0*pow(P2, 4))*pow(Q1, 2)+(16.0+252.0*pow(P2, 2)+ 
    21.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 4)*(27.0*pow(Q1, 3)*Q2+47.0*Q1*pow(Q2, 3))+3.0*pow(P1, 3)* 
    P2*((-5.0)*pow(Q1, 4)+2.0*pow(Q1, 2)*pow(Q2, 2)+31.0*pow(Q2, 4))+P1*P2*(((-6.0)+13.0* 
    pow(P2, 2))*pow(Q1, 4)+6.0*(54.0+17.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+3.0*(126.0+11.0*pow(P2, 2))* 
    pow(Q2, 4)))+(-3.0)*pow(G, 2)*(37.0*pow(P1, 4)*Q1*Q2+6.0*pow(P1, 2)*(24.0+(-13.0)*pow(P2, 2)) 
    *Q1*Q2+(-1.0)*(16.0+240.0*pow(P2, 2)+27.0*pow(P2, 4))*Q1*Q2+pow(P1, 3)*P2*((-14.0) 
    *pow(Q1, 2)+94.0*pow(Q2, 2))+P1*(10.0*pow(P2, 3)*(3.0*pow(Q1, 2)+5.0*pow(Q2, 2))+48.0*P2*( 
    pow(Q1, 2)+9.0*pow(Q2, 2))))) 
    ) + 
    (c2L - c2L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(6.0*pow(G, 4)*P1*(2.0+pow(P1, 2)+3.0*pow(P2, 2))+3.0*pow(G, 2)*(150.0* 
    pow(P1, 2)*P2*Q1*Q2+10.0*P2*(16.0+11.0*pow(P2, 2))*Q1*Q2+pow(P1, 3)*(19.0* 
    pow(Q1, 2)+(-99.0)*pow(Q2, 2))+(-1.0)*P1*((8.0+81.0*pow(P2, 2))*pow(Q1, 2)+(152.0+159.0* 
    pow(P2, 2))*pow(Q2, 2)))+(-28.0)*(P2*Q1*Q2*((72.0+65.0*pow(P2, 2))*pow(Q1, 2)+(88.0+ 
    45.0*pow(P2, 2))*pow(Q2, 2))+3.0*pow(P1, 2)*P2*Q2*(7.0*pow(Q1, 3)+43.0*Q1*pow(Q2, 2))+ 
    pow(P1, 3)*(9.0*pow(Q1, 4)+3.0*pow(Q1, 2)*pow(Q2, 2)+(-50.0)*pow(Q2, 4))+(-1.0)*P1*(((-4.0)+ 
    15.0*pow(P2, 2))*pow(Q1, 4)+3.0*(16.0+51.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+2.0*(34.0+27.0* 
    pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c3L - c3L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(pow(G, 4)*P1*P2*(pow(P1, 2)+3.0*(4.0+pow(P2, 2)))+pow(G, 2)*(51.0* 
    pow(P1, 4)*Q1*Q2+18.0*pow(P1, 2)*(24.0+7.0*pow(P2, 2))*Q1*Q2+3.0*(48.0+144.0*pow(P2, 2)+ 
    17.0*pow(P2, 4))*Q1*Q2+(-2.0)*pow(P1, 3)*P2*(7.0*pow(Q1, 2)+13.0*pow(Q2, 2))+(-6.0)* 
    P1*P2*((40.0+11.0*pow(P2, 2))*pow(Q1, 2)+(40.0+9.0*pow(P2, 2))*pow(Q2, 2)))+(-28.0)*(6.0* 
    pow(P1, 2)*Q1*Q2*((5.0+3.0*pow(P2, 2))*pow(Q1, 2)+(19.0+4.0*pow(P2, 2))*pow(Q2, 2))+Q1*Q2* 
    (2.0*(10.0+45.0*pow(P2, 2)+6.0*pow(P2, 4))*pow(Q1, 2)+(28.0+54.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q2, 2)) 
    +pow(P1, 4)*(2.0*pow(Q1, 3)*Q2+15.0*Q1*pow(Q2, 3))+pow(P1, 3)*P2*(pow(Q1, 4)+(-20.0)* 
    pow(Q1, 2)*pow(Q2, 2)+(-1.0)*pow(Q2, 4))+(-1.0)*P1*P2*(3.0*(5.0+2.0*pow(P2, 2))*pow(Q1, 4)+ 
    30.0*(5.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(15.0+4.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c4L - c4L0) * ( 
    (5.0/64.0)*pow(G, -4.0)*(3.0*pow(G, 4)*(pow(P1, 3)+(-3.0)*P1*pow(P2, 2))+12.0*pow(G, 2)*(( 
    -27.0)*pow(P1, 2)*P2*Q1*Q2+(-1.0)*P2*(34.0+25.0*pow(P2, 2))*Q1*Q2+2.0* 
    pow(P1, 3)*(4.0*pow(Q1, 2)+(-9.0)*pow(Q2, 2))+P1*((17.0+27.0*pow(P2, 2))*pow(Q1, 2)+((-17.0)+ 
    3.0*pow(P2, 2))*pow(Q2, 2)))+(-28.0)*((-4.0)*P2*Q1*Q2*((27.0+20.0*pow(P2, 2))* 
    pow(Q1, 2)+(7.0+5.0*pow(P2, 2))*pow(Q2, 2))+(-12.0)*pow(P1, 2)*P2*Q2*(7.0*pow(Q1, 3)+2.0*Q1* 
    pow(Q2, 2))+pow(P1, 3)*(5.0*pow(Q1, 4)+66.0*pow(Q1, 2)*pow(Q2, 2)+(-47.0)*pow(Q2, 4))+P1*((16.0+ 
    33.0*pow(P2, 2))*pow(Q1, 4)+18.0*(6.0+7.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(52.0+15.0* 
    pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c5L - c5L0) * ( 
    (1.0/32.0)*pow(G, -4.0)*(3.0*pow(G, 4)*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2))+3.0*pow(G, 2)*( 
    33.0*pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*((-32.0)+pow(P2, 2))*Q1*Q2+(-1.0)*pow(P2, 2)* 
    (192.0+31.0*pow(P2, 2))*Q1*Q2+2.0*pow(P1, 3)*P2*(7.0*pow(Q1, 2)+(-27.0)*pow(Q2, 2))+2.0* 
    P1*P2*((96.0+25.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(96.0+5.0*pow(P2, 2))*pow(Q2, 2)))+28.0*( 
    6.0*pow(P1, 2)*Q1*Q2*((15.0+7.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(47.0+6.0*pow(P2, 2))* 
    pow(Q2, 2))+Q1*Q2*((52.0+222.0*pow(P2, 2)+30.0*pow(P2, 4))*pow(Q1, 2)+((-52.0)+(-30.0)* 
    pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+pow(P1, 4)*(8.0*pow(Q1, 3)*Q2+(-41.0)*Q1*pow(Q2, 3))+ 
    pow(P1, 3)*P2*((-7.0)*pow(Q1, 4)+27.0*pow(Q2, 4))+(-1.0)*P1*P2*((81.0+20.0*pow(P2, 2))* 
    pow(Q1, 4)+30.0*(3.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(111.0+10.0*pow(P2, 2))*pow(Q2, 4)))) 
    ) + 
    (c6L - c6L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(15.0*pow(G, 2)*((-6.0)*pow(P1, 2)*P2*Q1*Q2+2.0*pow(P2, 3)* 
    Q1*Q2+pow(P1, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2)))+28.0*(P2*Q1*Q2*((-1.0)*(16.0+15.0*pow(P2, 2))*pow(Q1, 2)+(16.0+5.0* 
    pow(P2, 2))*pow(Q2, 2))+(-3.0)*pow(P1, 2)*P2*Q2*(pow(Q1, 3)+(-11.0)*Q1*pow(Q2, 2))+ 
    pow(P1, 3)*(pow(Q1, 4)+(-21.0)*pow(Q1, 2)*pow(Q2, 2)+6.0*pow(Q2, 4))+P1*((4.0+9.0*pow(P2, 2))* 
    pow(Q1, 4)+(-3.0)*(8.0+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+2.0*(2.0+(-3.0)*pow(P2, 2))*pow(Q2, 4))) 
    ) ) + 
    (c7L - c7L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(6.0*pow(P1, 2)*Q1*Q2*((44.0+7.0*pow(P2, 2))*pow(Q1, 2)+((-44.0)+ 
    21.0*pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)*Q1*Q2*((-3.0)*(88.0+17.0*pow(P2, 2))*pow(Q1, 2)+( 
    264.0+23.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(37.0*pow(Q1, 3)*Q2+(-65.0)*Q1*pow(Q2, 3))+ 
    pow(P1, 3)*P2*(pow(Q1, 4)+(-174.0)*pow(Q1, 2)*pow(Q2, 2)+57.0*pow(Q2, 4))+P1*P2*((132.0+ 
    43.0*pow(P2, 2))*pow(Q1, 4)+(-18.0)*(44.0+5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(132.0+(-13.0)* 
    pow(P2, 2))*pow(Q2, 4))+3.0*pow(G, 2)*(pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*Q2+ 
    pow(P2, 4)*Q1*Q2+2.0*pow(P1, 3)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*P1*pow(P2, 3)*(( 
    -1.0)*pow(Q1, 2)+pow(Q2, 2))))
    )  + 
    (c8L - c8L0) * ( 
    (175.0/64.0)*pow(G, -4.0)*((-12.0)*pow(P1, 2)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    4.0*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+(-3.0)*P1*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))  
    ) + 
    (c9L - c9L0) * ( 
    (35.0/32.0)*pow(G, -4.0)*(pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*Q1* 
    Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*Q1*Q2*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*P2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*P1* 
    pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) +  
    (sL - sL0) * ( 
    (5.0/64.0)*pow(G, -4.0)*(3.0*pow(G, 4)*(8.0+pow(P1, 4)+36.0*pow(P2, 2)+5.0*pow(P2, 4)+6.0*pow(P1, 2)*( 
    2.0+pow(P2, 2)))+3.0*pow(G, 2)*(11.0*((-16.0)+pow(P1, 4)+(-96.0)*pow(P2, 2)+(-6.0)*pow(P1, 2)* 
    pow(P2, 2)+(-15.0)*pow(P2, 4))*pow(Q1, 2)+8.0*P1*P2*(144.0+13.0*pow(P1, 2)+35.0*pow(P2, 2))* 
    Q1*Q2+(-1.0)*(144.0+51.0*pow(P1, 4)+384.0*pow(P2, 2)+35.0*pow(P2, 4)+6.0*pow(P1, 2)*(80.0+ 
    29.0*pow(P2, 2)))*pow(Q2, 2))+(-7.0)*(((-272.0)+21.0*pow(P1, 4)+(-1740.0)*pow(P2, 2)+(-287.0) 
    *pow(P2, 4)+(-18.0)*pow(P1, 2)*((-6.0)+pow(P2, 2)))*pow(Q1, 4)+16.0*P1*P2*(138.0+3.0* 
    pow(P1, 2)+43.0*pow(P2, 2))*pow(Q1, 3)*Q2+6.0*((-80.0)+pow(P1, 4)+(-372.0)*pow(P2, 2)+(-43.0)* 
    pow(P2, 4)+(-6.0)*pow(P1, 2)*(18.0+19.0*pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+16.0*P1*P2*(150.0+ 
    23.0*pow(P1, 2)+27.0*pow(P2, 2))*Q1*pow(Q2, 3)+(-1.0)*(208.0+103.0*pow(P1, 4)+396.0*pow(P2, 2)+ 
    27.0*pow(P2, 4)+6.0*pow(P1, 2)*(142.0+39.0*pow(P2, 2)))*pow(Q2, 4)))
    ) + 
    (s2L - s2L0) * ( 
    (5.0/32.0)*pow(G, -4.0)*(12.0*pow(G, 4)*(P2+pow(P2, 3))+(-3.0)*pow(G, 2)*((-2.0)*P1*(72.0+ 
    49.0*pow(P1, 2))*Q1*Q2+(-138.0)*P1*pow(P2, 2)*Q1*Q2+5.0*pow(P2, 3)*(29.0*pow(Q1, 2)+ 
    3.0*pow(Q2, 2))+5.0*P2*((32.0+9.0*pow(P1, 2))*pow(Q1, 2)+(-9.0)*pow(P1, 2)*pow(Q2, 2)))+(-28.0) 
    *(P1*Q1*Q2*((64.0+27.0*pow(P1, 2))*pow(Q1, 2)+(80.0+71.0*pow(P1, 2))*pow(Q2, 2))+3.0* 
    P1*pow(P2, 2)*Q2*(37.0*pow(Q1, 3)+9.0*Q1*pow(Q2, 2))+(-5.0)*pow(P2, 3)*(13.0*pow(Q1, 4)+ 
    9.0*pow(Q1, 2)*pow(Q2, 2))+(-1.0)*P2*((68.0+9.0*pow(P1, 2))*pow(Q1, 4)+9.0*(8.0+9.0*pow(P1, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(-12.0)*(1.0+3.0*pow(P1, 2))*pow(Q2, 4)))) 
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*pow(P1, 4)+6.0*pow(P1, 2)*(4.0+pow(P2, 2))+(-1.0)* 
    pow(P2, 2)*(24.0+5.0*pow(P2, 2)))+2.0*pow(G, 2)*((144.0+27.0*pow(P1, 4)+672.0*pow(P2, 2)+107.0* 
    pow(P2, 4)+6.0*pow(P1, 2)*(32.0+5.0*pow(P2, 2)))*pow(Q1, 2)+24.0*P1*P2*(pow(P1, 2)+(-1.0)* 
    pow(P2, 2))*Q1*Q2+(-1.0)*(144.0+87.0*pow(P1, 4)+192.0*pow(P2, 2)+7.0*pow(P2, 4)+6.0*pow(P1, 2)* 
    (112.0+25.0*pow(P2, 2)))*pow(Q2, 2))+(-28.0)*((44.0+5.0*pow(P1, 4)+222.0*pow(P2, 2)+35.0* 
    pow(P2, 4)+6.0*pow(P1, 2)*(7.0+2.0*pow(P2, 2)))*pow(Q1, 4)+(-16.0)*P1*P2*(2.0*pow(P1, 2)+3.0*( 
    5.0+pow(P2, 2)))*pow(Q1, 3)*Q2+4.0*(6.0+6.0*pow(P1, 4)+3.0*pow(P2, 2)+pow(P2, 4)+(-3.0)*pow(P1, 2)*(( 
    -11.0)+pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+16.0*P1*P2*(15.0+3.0*pow(P1, 2)+2.0*pow(P2, 2))* 
    Q1*pow(Q2, 3)+(-1.0)*(52.0+33.0*pow(P1, 4)+66.0*pow(P2, 2)+3.0*pow(P2, 4)+6.0*pow(P1, 2)*(41.0+8.0* 
    pow(P2, 2)))*pow(Q2, 4)))
    ) + 
    (s4L - s4L0) * ( 
    (-5.0/64.0)*pow(G, -4.0)*(pow(G, 4)*(9.0*pow(P1, 2)*P2+(-3.0)*pow(P2, 3))+6.0*pow(G, 2)*(4.0* 
    P1*(17.0+13.0*pow(P1, 2))*Q1*Q2+48.0*P1*pow(P2, 2)*Q1*Q2+5.0*pow(P2, 3)*(7.0* 
    pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)*P2*(((-34.0)+3.0*pow(P1, 2))*pow(Q1, 2)+(34.0+57.0* 
    pow(P1, 2))*pow(Q2, 2)))+56.0*(3.0*P1*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+(-17.0)*pow(Q2, 2))+ 
    P1*Q1*Q2*((2.0+pow(P1, 2))*pow(Q1, 2)+(-1.0)*(70.0+53.0*pow(P1, 2))*pow(Q2, 2))+5.0* 
    pow(P2, 3)*((-4.0)*pow(Q1, 4)+3.0*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+P2*((-2.0)*(11.0+3.0* 
    pow(P1, 2))*pow(Q1, 4)+15.0*(2.0+3.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+3.0*(4.0+7.0*pow(P1, 2))* 
    pow(Q2, 4))))
    ) + 
    (s5L - s5L0) * ( 
    (1.0/128.0)*pow(G, -4.0)*(3.0*pow(G, 4)*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+(-28.0) 
    *(((-52.0)+6.0*pow(P1, 2)+pow(P1, 4)+(-318.0)*pow(P2, 2)+(-53.0)*pow(P2, 4))*pow(Q1, 4)+(-16.0) 
    *P1*P2*(33.0+6.0*pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 3)*Q2+12.0*(26.0+11.0*pow(P1, 4)+63.0* 
    pow(P2, 2)+6.0*pow(P2, 4)+3.0*pow(P1, 2)*(31.0+9.0*pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+(-16.0)*P1* 
    P2*(63.0+11.0*pow(P1, 2)+10.0*pow(P2, 2))*Q1*pow(Q2, 3)+((-52.0)+(-65.0)*pow(P1, 4)+66.0* 
    pow(P2, 2)+9.0*pow(P2, 4)+6.0*pow(P1, 2)*((-63.0)+2.0*pow(P2, 2)))*pow(Q2, 4))+6.0*pow(G, 2)*(( 
    -136.0)*pow(P1, 3)*P2*Q1*Q2+(-24.0)*P1*P2*(32.0+5.0*pow(P2, 2))*Q1*Q2+ 
    pow(P1, 4)*(23.0*pow(Q1, 2)+(-43.0)*pow(Q2, 2))+pow(P2, 2)*((-1.0)*(192.0+41.0*pow(P2, 2))* 
    pow(Q1, 2)+3.0*(64.0+7.0*pow(P2, 2))*pow(Q2, 2))+6.0*pow(P1, 2)*((32.0+9.0*pow(P2, 2))*pow(Q1, 2)+(( 
    -32.0)+11.0*pow(P2, 2))*pow(Q2, 2)))) 
    ) + 
    (s6L - s6L0) * ( 
    (5.0/32.0)*pow(G, -4.0)*(15.0*pow(G, 2)*(2.0*pow(P1, 3)*Q1*Q2+(-6.0)*P1*pow(P2, 2)* 
    Q1*Q2+3.0*pow(P1, 2)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 3)*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2)))+28.0*(P1*Q1*Q2*((16.0+9.0*pow(P1, 2))*pow(Q1, 2)+(-1.0)*(16.0+19.0* 
    pow(P1, 2))*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*(7.0*pow(Q1, 3)+3.0*Q1*pow(Q2, 2))+5.0* 
    pow(P2, 3)*(pow(Q1, 4)+(-3.0)*pow(Q1, 2)*pow(Q2, 2))+P2*((4.0+(-3.0)*pow(P1, 2))*pow(Q1, 4)+(-3.0) 
    *(8.0+9.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+4.0*(1.0+3.0*pow(P1, 2))*pow(Q2, 4))))
    ) + 
    (s7L - s7L0) * ( 
    (-5.0/128.0)*pow(G, -4.0)*((-16.0)*P1*P2*Q1*Q2*((132.0+29.0*pow(P2, 2))* 
    pow(Q1, 2)+(-1.0)*(132.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*P2*((-240.0)*pow(Q1, 3)*Q2+ 
    688.0*Q1*pow(Q2, 3))+pow(P1, 4)*(23.0*pow(Q1, 4)+(-306.0)*pow(Q1, 2)*pow(Q2, 2)+79.0*pow(Q2, 4))+ 
    6.0*pow(P1, 2)*((44.0+21.0*pow(P2, 2))*pow(Q1, 4)+6.0*((-44.0)+7.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+(44.0+(-35.0)*pow(P2, 2))*pow(Q2, 4))+(-1.0)*pow(P2, 2)*((264.0+65.0*pow(P2, 2))* 
    pow(Q1, 4)+(-6.0)*(264.0+37.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+3.0*(88.0+3.0*pow(P2, 2))*pow(Q2, 4)) 
    +6.0*pow(G, 2)*((-8.0)*pow(P1, 3)*P2*Q1*Q2+8.0*P1*pow(P2, 3)*Q1*Q2+pow(P1, 4)*( 
    pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*(( 
    -1.0)*pow(Q1, 2)+pow(Q2, 2)))) 
    ) + 
    (s8L - s8L0) * ( 
    (-175.0/64.0)*pow(G, -4.0)*(4.0*pow(P1, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-12.0)* 
    P1*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*pow(P1, 2)*P2*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))) 
    ) + 
    (s9L - s9L0) * ( 
    (35.0/288.0)*pow(G, -4.0)*((-16.0)*pow(P1, 3)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    16.0*P1*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+pow(P2, 4)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    );


    // P12
    I[3] = (Lf - L0) * ( 
    (-15.0/8.0)*pow(G, -2.0)*((-14.0)*(pow(Q1, 2)+pow(Q2, 2))*((2.0+pow(P1, 2)+5.0*pow(P2, 2))* 
    pow(Q1, 2)+(-8.0)*P1*P2*Q1*Q2+(2.0+5.0*pow(P1, 2)+pow(P2, 2))*pow(Q2, 2))+pow(G, 2)*((4.0+ 
    3.0*pow(P1, 2)+9.0*pow(P2, 2))*pow(Q1, 2)+(-12.0)*P1*P2*Q1*Q2+(4.0+9.0*pow(P1, 2)+3.0* 
    pow(P2, 2))*pow(Q2, 2))) 
    ) + 
    (cL - cL0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(6.0*pow(G, 2)*(P1*(6.0+pow(P1, 2)+3.0*pow(P2, 2))*pow(Q1, 2)+(-2.0)* 
    P2*(6.0+3.0*pow(P1, 2)+pow(P2, 2))*Q1*Q2+P1*(5.0*pow(P1, 2)+3.0*(6.0+pow(P2, 2)))*pow(Q2, 2)) 
    +7.0*((-3.0)*P1*(8.0+pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 4)+4.0*P2*(24.0+9.0*pow(P1, 2)+5.0* 
    pow(P2, 2))*pow(Q1, 3)*Q2+(-6.0)*P1*(24.0+5.0*pow(P1, 2)+9.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+ 
    12.0*P2*(8.0+5.0*pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 3)+(-5.0)*P1*(7.0*pow(P1, 2)+3.0*(8.0+ 
    pow(P2, 2)))*pow(Q2, 4)))
    ) + 
    (c2L - c2L0) * ( 
    (1.0/8.0)*pow(G, -2.0)*(35.0*((-15.0)*P1*P2*pow(Q1, 4)+2.0*(8.0+9.0*pow(P1, 2)+15.0* 
    pow(P2, 2))*pow(Q1, 3)*Q2+(-54.0)*P1*P2*pow(Q1, 2)*pow(Q2, 2)+2.0*(8.0+15.0*pow(P1, 2)+9.0* 
    pow(P2, 2))*Q1*pow(Q2, 3)+(-15.0)*P1*P2*pow(Q2, 4))+(-30.0)*pow(G, 2)*(3.0*pow(P1, 2)* 
    Q1*Q2+(2.0+3.0*pow(P2, 2))*Q1*Q2+(-3.0)*P1*P2*(pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (c3L - c3L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(P1*(12.0+pow(P1, 2)+9.0*pow(P2, 2))*pow(Q1, 2)*(pow(G, 2)+(-7.0)* 
    pow(Q1, 2))+(-6.0)*P2*(4.0+pow(P1, 2)+pow(P2, 2))*Q1*(pow(G, 2)+(-14.0)*pow(Q1, 2))*Q2+(-1.0) 
    *P1*(pow(G, 2)*(12.0+5.0*pow(P1, 2)+(-3.0)*pow(P2, 2))+14.0*(12.0+pow(P1, 2)+9.0*pow(P2, 2))* 
    pow(Q1, 2))*pow(Q2, 2)+28.0*P2*(4.0+pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 3)+7.0*P1*(20.0+7.0* 
    pow(P1, 2)+(-1.0)*pow(P2, 2))*pow(Q2, 4))
    ) + 
    (c4L - c4L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(9.0*pow(G, 2)*(P2*Q1+P1*Q2)*(P1*Q1+(-1.0)*P2*Q2)+ 
    28.0*((-3.0)*P1*P2*pow(Q1, 4)+2.0*(1.0+3.0*pow(P2, 2))*pow(Q1, 3)*Q2+(-2.0)*(1.0+3.0* 
    pow(P1, 2))*Q1*pow(Q2, 3)+3.0*P1*P2*pow(Q2, 4)))
    ) + 
    (c5L - c5L0) * ( 
    (1.0/16.0)*pow(G, -2.0)*(7.0*(P1*(pow(P1, 2)+(-3.0)*(4.0+5.0*pow(P2, 2)))*pow(Q1, 4)+4.0* 
    P2*(12.0+(-3.0)*pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 3)*Q2+18.0*P1*(4.0+pow(P1, 2)+pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(-4.0)*P2*(12.0+9.0*pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 3)+P1*((-12.0)+( 
    -7.0)*pow(P1, 2)+9.0*pow(P2, 2))*pow(Q2, 4))+3.0*pow(G, 2)*(6.0*pow(P1, 2)*P2*Q1*Q2+(-2.0)* 
    pow(P2, 3)*Q1*Q2+3.0*P1*pow(P2, 2)*(Q1+(-1.0)*Q2)*(Q1+Q2)+pow(P1, 3)*((-1.0)* 
    pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (c6L - c6L0) * ( 
    (35.0/8.0)*pow(G, -2.0)*(2.0*P1*Q1*Q2+P2*(Q1+(-1.0)*Q2)*(Q1+Q2))*(2.0* 
    P2*Q1*Q2+P1*((-1.0)*pow(Q1, 2)+pow(Q2, 2))) 
    ) + 
    (c7L - c7L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(4.0*pow(P2, 3)*Q1*(Q1+(-1.0)*Q2)*Q2*(Q1+Q2)+12.0* 
    pow(P1, 2)*P2*Q1*Q2*((-1.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+(-3.0)*P1*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (sL - sL0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(6.0*pow(G, 2)*((-1.0)*P2*(18.0+3.0*pow(P1, 2)+5.0*pow(P2, 2))* 
    pow(Q1, 2)+2.0*P1*(6.0+pow(P1, 2)+3.0*pow(P2, 2))*Q1*Q2+(-1.0)*P2*(6.0+3.0*pow(P1, 2)+ 
    pow(P2, 2))*pow(Q2, 2))+7.0*(5.0*P2*(24.0+3.0*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q1, 4)+(-12.0)* 
    P1*(8.0+pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 3)*Q2+6.0*P2*(24.0+9.0*pow(P1, 2)+5.0*pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(-4.0)*P1*(24.0+5.0*pow(P1, 2)+9.0*pow(P2, 2))*Q1*pow(Q2, 3)+3.0*P2*( 
    8.0+5.0*pow(P1, 2)+pow(P2, 2))*pow(Q2, 4))) 
    ) + 
    (s2L - s2L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(12.0*pow(G, 2)*((-1.0)*(1.0+3.0*pow(P2, 2))*pow(Q1, 2)+(1.0+3.0* 
    pow(P1, 2))*pow(Q2, 2))+7.0*((16.0+3.0*pow(P1, 2)+45.0*pow(P2, 2))*pow(Q1, 4)+(-24.0)*P1*P2* 
    pow(Q1, 3)*Q2+18.0*((-1.0)*pow(P1, 2)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+24.0*P1*P2*Q1* 
    pow(Q2, 3)+(-1.0)*(16.0+45.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q2, 4)))
    ) + 
    (s3L - s3L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(pow(G, 2)*P2*((-12.0)+3.0*pow(P1, 2)+(-5.0)*pow(P2, 2))*pow(Q1, 2)+ 
    7.0*P2*(20.0+(-1.0)*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q1, 4)+2.0*P1*(4.0+pow(P1, 2)+pow(P2, 2))* 
    Q1*((-3.0)*pow(G, 2)+14.0*pow(Q1, 2))*Q2+P2*(12.0+9.0*pow(P1, 2)+pow(P2, 2))*(pow(G, 2)+( 
    -14.0)*pow(Q1, 2))*pow(Q2, 2)+84.0*P1*(4.0+pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 3)+(-7.0)*P2*( 
    12.0+9.0*pow(P1, 2)+pow(P2, 2))*pow(Q2, 4))
    ) + 
    (s4L - s4L0) * ( 
    (5.0/32.0)*pow(G, -2.0)*(14.0*((2.0+(-3.0)*pow(P1, 2)+9.0*pow(P2, 2))*pow(Q1, 4)+24.0*P1* 
    P2*pow(Q1, 3)*Q2+(-6.0)*(2.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+24.0*P1* 
    P2*Q1*pow(Q2, 3)+(2.0+9.0*pow(P1, 2)+(-3.0)*pow(P2, 2))*pow(Q2, 4))+9.0*pow(G, 2)*(P2*(Q1+( 
    -1.0)*Q2)+P1*(Q1+Q2))*(P1*(Q1+(-1.0)*Q2)+(-1.0)*P2*(Q1+Q2)))
    ) + 
    (s5L - s5L0) * ( 
    (1.0/16.0)*pow(G, -2.0)*(7.0*(P2*(12.0+(-9.0)*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q1, 4)+4.0*P1* 
    (12.0+pow(P1, 2)+9.0*pow(P2, 2))*pow(Q1, 3)*Q2+(-18.0)*P2*(4.0+pow(P1, 2)+pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+(-4.0)*P1*(12.0+5.0*pow(P1, 2)+(-3.0)*pow(P2, 2))*Q1*pow(Q2, 3)+P2*(12.0+15.0* 
    pow(P1, 2)+(-1.0)*pow(P2, 2))*pow(Q2, 4))+3.0*pow(G, 2)*(2.0*pow(P1, 3)*Q1*Q2+(-6.0)*P1* 
    pow(P2, 2)*Q1*Q2+3.0*pow(P1, 2)*P2*(Q1+(-1.0)*Q2)*(Q1+Q2)+pow(P2, 3)*((-1.0)* 
    pow(Q1, 2)+pow(Q2, 2)))) 
    ) + 
    (s6L - s6L0) * ( 
    (35.0/16.0)*pow(G, -2.0)*(8.0*P1*P2*Q1*(Q1+(-1.0)*Q2)*Q2*(Q1+Q2)+(-1.0) 
    *pow(P1, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+pow(P2, 2)*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (s7L - s7L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(12.0*P1*pow(P2, 2)*Q1*(Q1+(-1.0)*Q2)*Q2*(Q1+Q2)+ 
    4.0*pow(P1, 3)*Q1*Q2*((-1.0)*pow(Q1, 2)+pow(Q2, 2))+(-3.0)*pow(P1, 2)*P2*(pow(Q1, 4)+(-6.0) 
    *pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) ;


    // P21
    I[4] = (Lf - L0) * ( 
    (15.0/16.0)*pow(G, -4.0)*(pow(G, 4)*P1*(4.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))+14.0*(pow(Q1, 2)+ 
    pow(Q2, 2))*((-48.0)*pow(P1, 2)*P2*Q1*Q2+8.0*P2*((-1.0)+pow(P2, 2))*Q1*Q2+5.0* 
    pow(P1, 3)*(pow(Q1, 2)+5.0*pow(Q2, 2))+P1*((16.0+33.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*((-8.0)+ 
    pow(P2, 2))*pow(Q2, 2)))+(-2.0)*pow(G, 2)*((-72.0)*pow(P1, 2)*P2*Q1*Q2+12.0*P2*(( 
    -1.0)+pow(P2, 2))*Q1*Q2+15.0*pow(P1, 3)*(pow(Q1, 2)+3.0*pow(Q2, 2))+P1*((34.0+57.0*pow(P2, 2)) 
    *pow(Q1, 2)+(46.0+3.0*pow(P2, 2))*pow(Q2, 2)))) 
    ) + 
    (cL - cL0) * ( 
    (-5.0/64.0)*pow(G, -4.0)*(3.0*pow(G, 4)*(8.0+5.0*pow(P1, 4)+12.0*pow(P2, 2)+pow(P2, 4)+6.0*pow(P1, 2)* 
    (6.0+pow(P2, 2)))+(-3.0)*pow(G, 2)*((144.0+35.0*pow(P1, 4)+480.0*pow(P2, 2)+51.0*pow(P2, 4)+6.0* 
    pow(P1, 2)*(64.0+29.0*pow(P2, 2)))*pow(Q1, 2)+(-8.0)*P1*P2*(144.0+35.0*pow(P1, 2)+13.0* 
    pow(P2, 2))*Q1*Q2+11.0*(16.0+15.0*pow(P1, 4)+(-1.0)*pow(P2, 4)+6.0*pow(P1, 2)*(16.0+pow(P2, 2))) 
    *pow(Q2, 2))+7.0*((208.0+27.0*pow(P1, 4)+852.0*pow(P2, 2)+103.0*pow(P2, 4)+18.0*pow(P1, 2)*(22.0+ 
    13.0*pow(P2, 2)))*pow(Q1, 4)+(-16.0)*P1*P2*(150.0+27.0*pow(P1, 2)+23.0*pow(P2, 2))* 
    pow(Q1, 3)*Q2+6.0*(80.0+43.0*pow(P1, 4)+108.0*pow(P2, 2)+(-1.0)*pow(P2, 4)+6.0*pow(P1, 2)*(62.0+ 
    19.0*pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+(-16.0)*P1*P2*(43.0*pow(P1, 2)+3.0*(46.0+pow(P2, 2))) 
    *Q1*pow(Q2, 3)+(272.0+287.0*pow(P1, 4)+(-108.0)*pow(P2, 2)+(-21.0)*pow(P2, 4)+6.0*pow(P1, 2)*( 
    290.0+3.0*pow(P2, 2)))*pow(Q2, 4))) 
    ) + 
    (c2L - c2L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(6.0*pow(G, 4)*P2*(2.0+3.0*pow(P1, 2)+pow(P2, 2))+(-3.0)*pow(G, 2)*(( 
    -10.0)*P1*(16.0+11.0*pow(P1, 2))*Q1*Q2+(-150.0)*P1*pow(P2, 2)*Q1*Q2+pow(P2, 3)* 
    (99.0*pow(Q1, 2)+(-19.0)*pow(Q2, 2))+P2*((152.0+159.0*pow(P1, 2))*pow(Q1, 2)+(8.0+81.0* 
    pow(P1, 2))*pow(Q2, 2)))+(-28.0)*(P1*Q1*Q2*((88.0+45.0*pow(P1, 2))*pow(Q1, 2)+(72.0+ 
    65.0*pow(P1, 2))*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*(43.0*pow(Q1, 3)+7.0*Q1*pow(Q2, 2))+ 
    pow(P2, 3)*((-50.0)*pow(Q1, 4)+3.0*pow(Q1, 2)*pow(Q2, 2)+9.0*pow(Q2, 4))+(-1.0)*P2*((68.0+54.0* 
    pow(P1, 2))*pow(Q1, 4)+3.0*(16.0+51.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+((-4.0)+15.0*pow(P1, 2))* 
    pow(Q2, 4))))
    ) + 
    (c3L - c3L0) * ( 
    (5.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(5.0*pow(P1, 4)+(-6.0)*pow(P1, 2)*((-4.0)+pow(P2, 2))+(-3.0) 
    *pow(P2, 2)*(8.0+pow(P2, 2)))+2.0*pow(G, 2)*((144.0+7.0*pow(P1, 4)+672.0*pow(P2, 2)+87.0*pow(P2, 4)+ 
    6.0*pow(P1, 2)*(32.0+25.0*pow(P2, 2)))*pow(Q1, 2)+24.0*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2))* 
    Q1*Q2+(-1.0)*(107.0*pow(P1, 4)+6.0*pow(P1, 2)*(112.0+5.0*pow(P2, 2))+3.0*(48.0+64.0* 
    pow(P2, 2)+9.0*pow(P2, 4)))*pow(Q2, 2))+(-28.0)*((52.0+3.0*pow(P1, 4)+246.0*pow(P2, 2)+33.0* 
    pow(P2, 4)+6.0*pow(P1, 2)*(11.0+8.0*pow(P2, 2)))*pow(Q1, 4)+(-16.0)*P1*P2*(2.0*pow(P1, 2)+3.0* 
    (5.0+pow(P2, 2)))*pow(Q1, 3)*Q2+(-4.0)*(6.0+pow(P1, 4)+33.0*pow(P2, 2)+6.0*pow(P2, 4)+(-3.0)* 
    pow(P1, 2)*((-1.0)+pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+16.0*P1*P2*(15.0+3.0*pow(P1, 2)+2.0* 
    pow(P2, 2))*Q1*pow(Q2, 3)+(-1.0)*(44.0+35.0*pow(P1, 4)+42.0*pow(P2, 2)+5.0*pow(P2, 4)+6.0* 
    pow(P1, 2)*(37.0+2.0*pow(P2, 2)))*pow(Q2, 4)))
    ) + 
    (c4L - c4L0) * ( 
    (5.0/64.0)*pow(G, -4.0)*(pow(G, 4)*(9.0*pow(P1, 2)*P2+(-3.0)*pow(P2, 3))+12.0*pow(G, 2)*(P1* 
    (34.0+25.0*pow(P1, 2))*Q1*Q2+27.0*P1*pow(P2, 2)*Q1*Q2+2.0*pow(P2, 3)*(9.0*pow(Q1, 2)+( 
    -4.0)*pow(Q2, 2))+(-1.0)*P2*(((-17.0)+3.0*pow(P1, 2))*pow(Q1, 2)+(17.0+27.0*pow(P1, 2))* 
    pow(Q2, 2)))+(-28.0)*(4.0*P1*Q1*Q2*((7.0+5.0*pow(P1, 2))*pow(Q1, 2)+(27.0+20.0*pow(P1, 2)) 
    *pow(Q2, 2))+12.0*P1*pow(P2, 2)*Q2*(2.0*pow(Q1, 3)+7.0*Q1*pow(Q2, 2))+pow(P2, 3)*(47.0* 
    pow(Q1, 4)+(-66.0)*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*pow(Q2, 4))+P2*((52.0+15.0*pow(P1, 2))*pow(Q1, 4)+ 
    (-18.0)*(6.0+7.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(16.0+33.0*pow(P1, 2))*pow(Q2, 4))))
    ) + 
    (c5L - c5L0) * ( 
    (1.0/128.0)*pow(G, -4.0)*((-3.0)*pow(G, 4)*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+ 
    28.0*(((-52.0)+9.0*pow(P1, 4)+(-378.0)*pow(P2, 2)+(-65.0)*pow(P2, 4)+6.0*pow(P1, 2)*(11.0+2.0* 
    pow(P2, 2)))*pow(Q1, 4)+(-16.0)*P1*P2*(63.0+10.0*pow(P1, 2)+11.0*pow(P2, 2))*pow(Q1, 3)*Q2+ 
    12.0*(26.0+6.0*pow(P1, 4)+93.0*pow(P2, 2)+11.0*pow(P2, 4)+9.0*pow(P1, 2)*(7.0+3.0*pow(P2, 2)))* 
    pow(Q1, 2)*pow(Q2, 2)+(-16.0)*P1*P2*(33.0+5.0*pow(P1, 2)+6.0*pow(P2, 2))*Q1*pow(Q2, 3)+(( 
    -52.0)+(-318.0)*pow(P1, 2)+(-53.0)*pow(P1, 4)+6.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+(-6.0)* 
    pow(G, 2)*((-120.0)*pow(P1, 3)*P2*Q1*Q2+(-8.0)*P1*P2*(96.0+17.0*pow(P2, 2))* 
    Q1*Q2+pow(P1, 4)*(21.0*pow(Q1, 2)+(-41.0)*pow(Q2, 2))+6.0*pow(P1, 2)*((32.0+11.0*pow(P2, 2))* 
    pow(Q1, 2)+((-32.0)+9.0*pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)*((-1.0)*(192.0+43.0*pow(P2, 2))* 
    pow(Q1, 2)+(192.0+23.0*pow(P2, 2))*pow(Q2, 2)))) 
    ) + 
    (c6L - c6L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(15.0*pow(G, 2)*(2.0*pow(P1, 3)*Q1*Q2+(-6.0)*P1*pow(P2, 2)* 
    Q1*Q2+3.0*pow(P1, 2)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 3)*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2)))+28.0*(P1*Q1*Q2*((16.0+5.0*pow(P1, 2))*pow(Q1, 2)+(-1.0)*(16.0+15.0* 
    pow(P1, 2))*pow(Q2, 2))+P1*pow(P2, 2)*(33.0*pow(Q1, 3)*Q2+(-3.0)*Q1*pow(Q2, 3))+pow(P2, 3)*( 
    6.0*pow(Q1, 4)+(-21.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+P2*((4.0+(-6.0)*pow(P1, 2))*pow(Q1, 4)+( 
    -3.0)*(8.0+3.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(4.0+9.0*pow(P1, 2))*pow(Q2, 4)))) 
    ) + 
    (c7L - c7L0) * ( 
    (5.0/128.0)*pow(G, -4.0)*((-16.0)*pow(P1, 3)*P2*Q1*Q2*(pow(Q1, 2)+(-29.0)*pow(Q2, 2))+ 
    (-16.0)*P1*P2*Q1*Q2*((132.0+43.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(44.0+5.0*pow(P2, 2)) 
    *pow(Q2, 2))+pow(P1, 4)*(9.0*pow(Q1, 4)+(-222.0)*pow(Q1, 2)*pow(Q2, 2)+65.0*pow(Q2, 4))+6.0* 
    pow(P1, 2)*((44.0+35.0*pow(P2, 2))*pow(Q1, 4)+(-6.0)*(44.0+7.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    44.0+(-21.0)*pow(P2, 2))*pow(Q2, 4))+(-1.0)*pow(P2, 2)*((264.0+79.0*pow(P2, 2))*pow(Q1, 4)+( 
    -18.0)*(88.0+17.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(264.0+23.0*pow(P2, 2))*pow(Q2, 4))+6.0* 
    pow(G, 2)*((-8.0)*pow(P1, 3)*P2*Q1*Q2+8.0*P1*pow(P2, 3)*Q1*Q2+pow(P1, 4)*(pow(Q1, 2)+ 
    (-1.0)*pow(Q2, 2))+pow(P2, 4)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*((-1.0)* 
    pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (c8L - c8L0) * ( 
    (175.0/64.0)*pow(G, -4.0)*(4.0*pow(P1, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-12.0)* 
    P1*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*pow(P1, 2)*P2*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (c9L - c9L0) * ( 
    (-35.0/128.0)*pow(G, -4.0)*((-16.0)*pow(P1, 3)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)) 
    +16.0*P1*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+pow(P2, 4)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (sL - sL0) * ( 
    (5.0/32.0)*pow(G, -4.0)*(6.0*pow(G, 4)*P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))+3.0*pow(G, 2)*(27.0* 
    pow(P1, 4)*Q1*Q2+6.0*pow(P1, 2)*(40.0+13.0*pow(P2, 2))*Q1*Q2+(16.0+(-144.0)*pow(P2, 2)+( 
    -37.0)*pow(P2, 4))*Q1*Q2+(-10.0)*pow(P1, 3)*P2*(5.0*pow(Q1, 2)+3.0*pow(Q2, 2))+(-2.0)* 
    P1*P2*((216.0+47.0*pow(P2, 2))*pow(Q1, 2)+(24.0+(-7.0)*pow(P2, 2))*pow(Q2, 2)))+(-14.0)*( 
    6.0*pow(P1, 2)*Q1*Q2*(21.0*(2.0+pow(P2, 2))*pow(Q1, 2)+(38.0+5.0*pow(P2, 2))*pow(Q2, 2))+(-1.0) 
    *Q1*Q2*(((-16.0)+156.0*pow(P2, 2)+47.0*pow(P2, 4))*pow(Q1, 2)+((-16.0)+132.0*pow(P2, 2)+ 
    27.0*pow(P2, 4))*pow(Q2, 2))+3.0*pow(P1, 4)*(7.0*pow(Q1, 3)*Q2+11.0*Q1*pow(Q2, 3))+(-1.0)* 
    pow(P1, 3)*P2*(33.0*pow(Q1, 4)+102.0*pow(Q1, 2)*pow(Q2, 2)+13.0*pow(Q2, 4))+(-3.0)*P1*P2*( 
    (126.0+31.0*pow(P2, 2))*pow(Q1, 4)+2.0*(54.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(2.0+5.0* 
    pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s2L - s2L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(12.0*pow(G, 4)*(P1+pow(P1, 3))+28.0*((-1.0)*P2*Q1*Q2*(( 
    80.0+71.0*pow(P2, 2))*pow(Q1, 2)+(64.0+27.0*pow(P2, 2))*pow(Q2, 2))+(-3.0)*pow(P1, 2)*P2*Q2*( 
    9.0*pow(Q1, 3)+37.0*Q1*pow(Q2, 2))+5.0*pow(P1, 3)*(9.0*pow(Q1, 2)*pow(Q2, 2)+13.0*pow(Q2, 4))+P1* 
    ((-12.0)*(1.0+3.0*pow(P2, 2))*pow(Q1, 4)+9.0*(8.0+9.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(68.0+9.0* 
    pow(P2, 2))*pow(Q2, 4)))+(-3.0)*pow(G, 2)*((-138.0)*pow(P1, 2)*P2*Q1*Q2+(-2.0)*P2*( 
    72.0+49.0*pow(P2, 2))*Q1*Q2+5.0*pow(P1, 3)*(3.0*pow(Q1, 2)+29.0*pow(Q2, 2))+5.0*P1*(32.0* 
    pow(Q2, 2)+(-9.0)*pow(P2, 2)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))))) 
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(pow(G, 4)*P1*P2*(12.0+3.0*pow(P1, 2)+pow(P2, 2))+pow(G, 2)*(51.0* 
    pow(P1, 4)*Q1*Q2+18.0*pow(P1, 2)*(24.0+7.0*pow(P2, 2))*Q1*Q2+3.0*(48.0+144.0*pow(P2, 2)+ 
    17.0*pow(P2, 4))*Q1*Q2+(-6.0)*pow(P1, 3)*P2*(9.0*pow(Q1, 2)+11.0*pow(Q2, 2))+(-2.0)* 
    P1*P2*((120.0+13.0*pow(P2, 2))*pow(Q1, 2)+(120.0+7.0*pow(P2, 2))*pow(Q2, 2)))+(-28.0)*( 
    6.0*pow(P1, 2)*Q1*Q2*((9.0+4.0*pow(P2, 2))*pow(Q1, 2)+3.0*(5.0+pow(P2, 2))*pow(Q2, 2))+Q1* 
    Q2*((28.0+114.0*pow(P2, 2)+15.0*pow(P2, 4))*pow(Q1, 2)+2.0*(10.0+15.0*pow(P2, 2)+pow(P2, 4))* 
    pow(Q2, 2))+pow(P1, 4)*(5.0*pow(Q1, 3)*Q2+12.0*Q1*pow(Q2, 3))+(-2.0)*pow(P1, 3)*P2*(2.0* 
    pow(Q1, 4)+15.0*pow(Q1, 2)*pow(Q2, 2)+3.0*pow(Q2, 4))+(-1.0)*P1*P2*((15.0+pow(P2, 2))*pow(Q1, 4)+ 
    10.0*(15.0+2.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*((-15.0)+pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s4L - s4L0) * ( 
    (5.0/64.0)*pow(G, -4.0)*(3.0*pow(G, 4)*(pow(P1, 3)+(-3.0)*P1*pow(P2, 2))+6.0*pow(G, 2)*((-48.0) 
    *pow(P1, 2)*P2*Q1*Q2+(-4.0)*P2*(17.0+13.0*pow(P2, 2))*Q1*Q2+5.0*pow(P1, 3)*( 
    3.0*pow(Q1, 2)+(-7.0)*pow(Q2, 2))+P1*((34.0+57.0*pow(P2, 2))*pow(Q1, 2)+((-34.0)+3.0*pow(P2, 2)) 
    *pow(Q2, 2)))+(-56.0)*(3.0*pow(P1, 2)*P2*Q1*Q2*((-17.0)*pow(Q1, 2)+pow(Q2, 2))+P2* 
    Q1*Q2*((-1.0)*(70.0+53.0*pow(P2, 2))*pow(Q1, 2)+(2.0+pow(P2, 2))*pow(Q2, 2))+5.0*pow(P1, 3)*( 
    pow(Q1, 4)+3.0*pow(Q1, 2)*pow(Q2, 2)+(-4.0)*pow(Q2, 4))+P1*(3.0*(4.0+7.0*pow(P2, 2))*pow(Q1, 4)+ 
    15.0*(2.0+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-2.0)*(11.0+3.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s5L - s5L0) * ( 
    (1.0/32.0)*pow(G, -4.0)*(3.0*pow(G, 4)*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2))+3.0*pow(G, 2)*( 
    31.0*pow(P1, 4)*Q1*Q2+6.0*pow(P1, 2)*(32.0+pow(P2, 2))*Q1*Q2+(-3.0)*pow(P2, 2)*(64.0+ 
    11.0*pow(P2, 2))*Q1*Q2+10.0*pow(P1, 3)*P2*(pow(Q1, 2)+(-5.0)*pow(Q2, 2))+2.0*P1*P2*( 
    3.0*(32.0+9.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(96.0+7.0*pow(P2, 2))*pow(Q2, 2)))+(-28.0)*( 
    pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+30.0*pow(Q2, 2))+(-6.0)*pow(P1, 2)*Q1*Q2*((5.0+6.0* 
    pow(P2, 2))*pow(Q1, 2)+(-1.0)*(37.0+7.0*pow(P2, 2))*pow(Q2, 2))+Q1*Q2*((-1.0)*(52.0+282.0* 
    pow(P2, 2)+41.0*pow(P2, 4))*pow(Q1, 2)+2.0*(26.0+45.0*pow(P2, 2)+4.0*pow(P2, 4))*pow(Q2, 2))+10.0* 
    pow(P1, 3)*P2*(pow(Q1, 4)+(-3.0)*pow(Q1, 2)*pow(Q2, 2)+(-2.0)*pow(Q2, 4))+P1*P2*(3.0*(37.0+ 
    9.0*pow(P2, 2))*pow(Q1, 4)+(-90.0)*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(81.0+7.0*pow(P2, 2))*pow(Q2, 4)))) 
    ) + 
    (s6L - s6L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(15.0*pow(G, 2)*((-6.0)*pow(P1, 2)*P2*Q1*Q2+2.0*pow(P2, 3)* 
    Q1*Q2+pow(P1, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2)))+(-28.0)*(P2*Q1*Q2*((16.0+19.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(16.0+9.0* 
    pow(P2, 2))*pow(Q2, 2))+(-3.0)*pow(P1, 2)*P2*Q2*(3.0*pow(Q1, 3)+7.0*Q1*pow(Q2, 2))+5.0* 
    pow(P1, 3)*(3.0*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*pow(Q2, 4))+P1*((-4.0)*(1.0+3.0*pow(P2, 2))* 
    pow(Q1, 4)+3.0*(8.0+9.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-4.0)+3.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s7L - s7L0) * ( 
    (-5.0/32.0)*pow(G, -4.0)*(6.0*pow(P1, 2)*Q1*Q2*((44.0+21.0*pow(P2, 2))*pow(Q1, 2)+((-44.0) 
    +7.0*pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)*Q1*Q2*((-1.0)*(264.0+65.0*pow(P2, 2))*pow(Q1, 2)+( 
    264.0+37.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(23.0*pow(Q1, 3)*Q2+(-51.0)*Q1*pow(Q2, 3))+ 
    pow(P1, 3)*P2*((-13.0)*pow(Q1, 4)+(-90.0)*pow(Q1, 2)*pow(Q2, 2)+43.0*pow(Q2, 4))+P1*P2*( 
    3.0*(44.0+19.0*pow(P2, 2))*pow(Q1, 4)+(-6.0)*(132.0+29.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(132.0+ 
    pow(P2, 2))*pow(Q2, 4))+3.0*pow(G, 2)*(pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*Q2+ 
    pow(P2, 4)*Q1*Q2+2.0*pow(P1, 3)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*P1*pow(P2, 3)*(( 
    -1.0)*pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (s8L - s8L0) * ( 
    (175.0/64.0)*pow(G, -4.0)*((-12.0)*pow(P1, 2)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    4.0*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+(-3.0)*P1*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (s9L - s9L0) * ( 
    (35.0/32.0)*pow(G, -4.0)*(pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*Q1* 
    Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*Q1*Q2*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*P2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*P1* 
    pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))) 
    );


    // Q1
    I[5] = (Lf - L0) * ( 
    (-15.0/8.0)*pow(G, -2.0)*(pow(G, 2)*((-6.0)*P1*P2*Q1+9.0*pow(P1, 2)*Q2+(4.0+3.0* 
    pow(P2, 2))*Q2)+(-14.0)*((-2.0)*P1*P2*Q1*(pow(Q1, 2)+3.0*pow(Q2, 2))+Q2*((2.0+3.0* 
    pow(P2, 2))*pow(Q1, 2)+(2.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+5.0*pow(Q2, 3)))) 
    ) + 
    (cL - cL0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*(6.0*pow(G, 2)*(3.0*(2.0+pow(P1, 2))*P2*Q1+pow(P2, 3)*Q1+(-1.0)* 
    P1*(18.0+5.0*pow(P1, 2))*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+7.0*(3.0*P1*pow(P2, 2)*Q2*( 
    9.0*pow(Q1, 2)+5.0*pow(Q2, 2))+(-3.0)*P2*Q1*((8.0+3.0*pow(P1, 2))*pow(Q1, 2)+3.0*(8.0+5.0* 
    pow(P1, 2))*pow(Q2, 2))+P1*Q2*(3.0*(24.0+5.0*pow(P1, 2))*pow(Q1, 2)+5.0*(24.0+7.0*pow(P1, 2))* 
    pow(Q2, 2))+(-1.0)*pow(P2, 3)*(5.0*pow(Q1, 3)+9.0*Q1*pow(Q2, 2)))) 
    ) + 
    (c2L - c2L0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*(6.0*pow(G, 2)*((2.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))*Q1+(-6.0)*P1* 
    P2*Q2)+(-7.0)*((8.0+9.0*pow(P1, 2)+15.0*pow(P2, 2))*pow(Q1, 3)+(-54.0)*P1*P2* 
    pow(Q1, 2)*Q2+3.0*(8.0+15.0*pow(P1, 2)+9.0*pow(P2, 2))*Q1*pow(Q2, 2)+(-30.0)*P1*P2* 
    pow(Q2, 3)))
    ) + 
    (c3L - c3L0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*(pow(G, 2)*(3.0*(4.0+pow(P1, 2))*P2*Q1+3.0*pow(P2, 3)*Q1+P1*( 
    12.0+5.0*pow(P1, 2))*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+7.0*((-3.0)*(4.0+pow(P1, 2))*P2* 
    Q1*(pow(Q1, 2)+pow(Q2, 2))+(-3.0)*pow(P2, 3)*Q1*(pow(Q1, 2)+pow(Q2, 2))+P1*pow(P2, 2)*Q2*( 
    9.0*pow(Q1, 2)+pow(Q2, 2))+P1*Q2*((12.0+pow(P1, 2))*pow(Q1, 2)+(-1.0)*(20.0+7.0*pow(P1, 2))* 
    pow(Q2, 2))))
    ) + 
    (c4L - c4L0) * ( 
    (5.0/32.0)*pow(G, -2.0)*(9.0*pow(G, 2)*(pow(P1, 2)*Q1+(-1.0)*pow(P2, 2)*Q1+(-2.0)*P1* 
    P2*Q2)+28.0*((1.0+3.0*pow(P2, 2))*pow(Q1, 3)+(-3.0)*(1.0+3.0*pow(P1, 2))*Q1*pow(Q2, 2)+6.0* 
    P1*P2*pow(Q2, 3))) 
    ) + 
    (c5L - c5L0) * ( 
    (1.0/16.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(3.0*pow(P1, 2)*P2*Q1+(-1.0)*pow(P2, 3)*Q1+ 
    pow(P1, 3)*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+7.0*(9.0*P1*pow(P2, 2)*Q2*(pow(Q1, 2)+pow(Q2, 2)) 
    +(-3.0)*P2*Q1*(((-4.0)+pow(P1, 2))*pow(Q1, 2)+3.0*(4.0+3.0*pow(P1, 2))*pow(Q2, 2))+P1* 
    Q2*(9.0*(4.0+pow(P1, 2))*pow(Q1, 2)+(-1.0)*(12.0+7.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(5.0* 
    pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))))
    ) + 
    (c6L - c6L0) * ( 
    (-35.0/16.0)*pow(G, -2.0)*((-1.0)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+2.0*P1* 
    P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))) 
    ) + 
    (c7L - c7L0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*(3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)* 
    pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))
    ) + 
    (sL - sL0) * ( 
    (15.0/16.0)*pow(G, -2.0)*(2.0*pow(G, 2)*(pow(P1, 3)*Q1+3.0*P1*(2.0+pow(P2, 2))*Q1+(-3.0)* 
    pow(P1, 2)*P2*Q2+(-1.0)*P2*(6.0+pow(P2, 2))*Q2)+(-7.0)*((-1.0)*pow(P1, 2)*P2* 
    Q2*(9.0*pow(Q1, 2)+5.0*pow(Q2, 2))+(-1.0)*P2*Q2*((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(8.0+ 
    pow(P2, 2))*pow(Q2, 2))+P1*Q1*((8.0+5.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+3.0*pow(P2, 2))* 
    pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+5.0*Q1*pow(Q2, 2))))
    ) + 
    (s2L - s2L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*((-42.0)*P1*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+9.0* 
    pow(P1, 2)*Q2*(4.0*pow(G, 2)+(-7.0)*(pow(Q1, 2)+5.0*pow(Q2, 2)))+Q2*(12.0*pow(G, 2)+7.0*(( 
    -16.0)*pow(Q2, 2)+pow(P2, 2)*(9.0*pow(Q1, 2)+(-3.0)*pow(Q2, 2))))) 
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*(pow(G, 2)*(3.0*pow(P1, 3)*Q1+3.0*P1*(4.0+pow(P2, 2))*Q1+(-9.0)* 
    pow(P1, 2)*P2*Q2+(-1.0)*P2*(12.0+pow(P2, 2))*Q2)+(-7.0)*((-9.0)*pow(P1, 2)*P2* 
    Q2*(pow(Q1, 2)+pow(Q2, 2))+(-1.0)*P2*(12.0+pow(P2, 2))*Q2*(pow(Q1, 2)+pow(Q2, 2))+P1*(4.0+ 
    pow(P2, 2))*Q1*(pow(Q1, 2)+9.0*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+9.0*Q1*pow(Q2, 2))))
    ) + 
    (s4L - s4L0) * ( 
    (-5.0/32.0)*pow(G, -2.0)*(9.0*pow(G, 2)*(2.0*P1*P2*Q1+pow(P1, 2)*Q2+(-1.0)*pow(P2, 2)* 
    Q2)+14.0*(9.0*pow(P1, 2)*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-6.0)*P1*P2*Q1*( 
    pow(Q1, 2)+3.0*pow(Q2, 2))+Q2*((6.0+9.0*pow(P2, 2))*pow(Q1, 2)+((-2.0)+3.0*pow(P2, 2))*pow(Q2, 2)))) 
    ) + 
    (s5L - s5L0) * ( 
    (1.0/16.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(pow(P1, 3)*Q1+(-3.0)*P1*pow(P2, 2)*Q1+(-3.0)* 
    pow(P1, 2)*P2*Q2+pow(P2, 3)*Q2)+7.0*((-1.0)*P2*Q2*(9.0*(4.0+pow(P2, 2))*pow(Q1, 2)+(( 
    -12.0)+pow(P2, 2))*pow(Q2, 2))+3.0*P1*Q1*((4.0+3.0*pow(P2, 2))*pow(Q1, 2)+3.0*((-4.0)+ 
    pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-15.0)*Q1*pow(Q2, 2))+pow(P1, 2)*P2*((-9.0)* 
    pow(Q1, 2)*Q2+15.0*pow(Q2, 3)))) 
    ) + 
    (s6L - s6L0) * ( 
    (35.0/16.0)*pow(G, -2.0)*(2.0*P1*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+pow(P2, 2)*Q2*( 
    (-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))
    ) + 
    (s7L - s7L0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*((-3.0)*P1*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0* 
    pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))
    );


    I[6] = (Lf - L0) * ( 
    (15.0/8.0)*pow(G, -2.0)*(pow(G, 2)*((4.0+3.0*pow(P1, 2)+9.0*pow(P2, 2))*Q1+(-6.0)*P1*P2* 
    Q2)+(-14.0)*((2.0+pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 3)+(-6.0)*P1*P2*pow(Q1, 2)*Q2+(2.0+ 
    3.0*pow(P1, 2)+3.0*pow(P2, 2))*Q1*pow(Q2, 2)+(-2.0)*P1*P2*pow(Q2, 3))) 
    ) + 
    (cL - cL0) * ( 
    (-15.0/16.0)*pow(G, -2.0)*(2.0*pow(G, 2)*(pow(P1, 3)*Q1+3.0*P1*(2.0+pow(P2, 2))*Q1+(-3.0) 
    *pow(P1, 2)*P2*Q2+(-1.0)*P2*(6.0+pow(P2, 2))*Q2)+(-7.0)*((-1.0)*pow(P1, 2)*P2* 
    Q2*(9.0*pow(Q1, 2)+5.0*pow(Q2, 2))+(-1.0)*P2*Q2*((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(8.0+ 
    pow(P2, 2))*pow(Q2, 2))+P1*Q1*((8.0+5.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+3.0*pow(P2, 2))* 
    pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+5.0*Q1*pow(Q2, 2))))
    ) + 
    (c2L - c2L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(6.0*pow(G, 2)*((-6.0)*P1*P2*Q1+3.0*pow(P1, 2)*Q2+(2.0+3.0* 
    pow(P2, 2))*Q2)+(-7.0)*((-6.0)*P1*P2*Q1*(5.0*pow(Q1, 2)+9.0*pow(Q2, 2))+Q2*(3.0*( 
    8.0+15.0*pow(P2, 2))*pow(Q1, 2)+(8.0+9.0*pow(P2, 2))*pow(Q2, 2))+3.0*pow(P1, 2)*(9.0*pow(Q1, 2)*Q2+ 
    5.0*pow(Q2, 3))))
    ) + 
    (c3L - c3L0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*(pow(G, 2)*(pow(P1, 3)*Q1+3.0*P1*(4.0+3.0*pow(P2, 2))*Q1+(-3.0)* 
    pow(P1, 2)*P2*Q2+(-3.0)*P2*(4.0+pow(P2, 2))*Q2)+(-7.0)*(pow(P1, 3)*Q1*(pow(Q1, 2)+ 
    pow(Q2, 2))+3.0*P1*(4.0+3.0*pow(P2, 2))*Q1*(pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P1, 2)*P2* 
    Q2*(9.0*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*P2*(4.0+pow(P2, 2))*Q2*(9.0*pow(Q1, 2)+pow(Q2, 2))))
    ) + 
    (c4L - c4L0) * ( 
    (-5.0/32.0)*pow(G, -2.0)*(9.0*pow(G, 2)*(2.0*P1*P2*Q1+pow(P1, 2)*Q2+(-1.0)*pow(P2, 2)* 
    Q2)+(-28.0)*(6.0*P1*P2*pow(Q1, 3)+pow(Q2, 3)+3.0*pow(P1, 2)*pow(Q2, 3)+(-3.0)*pow(Q1, 2)*( 
    Q2+3.0*pow(P2, 2)*Q2))) 
    ) + 
    (c5L - c5L0) * ( 
    (1.0/16.0)*pow(G, -2.0)*(3.0*pow(G, 2)*(pow(P1, 3)*Q1+(-3.0)*P1*pow(P2, 2)*Q1+(-3.0)* 
    pow(P1, 2)*P2*Q2+pow(P2, 3)*Q2)+(-7.0)*((-9.0)*pow(P1, 2)*P2*Q2*(pow(Q1, 2)+pow(Q2, 2)) 
    +(-3.0)*P1*Q1*((4.0+5.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(4.0+pow(P2, 2))*pow(Q2, 2))+P2* 
    Q2*(3.0*(12.0+5.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(12.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*( 
    pow(Q1, 3)+9.0*Q1*pow(Q2, 2)))) 
    ) + 
    (c6L - c6L0) * ( 
    (35.0/16.0)*pow(G, -2.0)*(2.0*P1*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+pow(P2, 2)*Q2*( 
    (-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3))) 
    ) + 
    (c7L - c7L0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*((-3.0)*P1*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0* 
    pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))
    ) + 
    (sL - sL0) * ( 
    (-5.0/16.0)*pow(G, -2.0)*((-6.0)*pow(G, 2)*(3.0*(6.0+pow(P1, 2))*P2*Q1+5.0*pow(P2, 3)*Q1+ 
    (-1.0)*P1*(6.0+pow(P1, 2))*Q2+(-3.0)*P1*pow(P2, 2)*Q2)+(-7.0)*(9.0*P1*pow(P2, 2)* 
    Q2*(5.0*pow(Q1, 2)+pow(Q2, 2))+(-3.0)*P2*Q1*(5.0*(8.0+pow(P1, 2))*pow(Q1, 2)+3.0*(8.0+3.0* 
    pow(P1, 2))*pow(Q2, 2))+P1*Q2*(9.0*(8.0+pow(P1, 2))*pow(Q1, 2)+(24.0+5.0*pow(P1, 2))*pow(Q2, 2))+ 
    (-5.0)*pow(P2, 3)*(7.0*pow(Q1, 3)+3.0*Q1*pow(Q2, 2))))
    ) + 
    (s2L - s2L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(12.0*pow(G, 2)*(Q1+3.0*pow(P2, 2)*Q1)+(-7.0)*((16.0+3.0*pow(P1, 2)+ 
    45.0*pow(P2, 2))*pow(Q1, 3)+(-18.0)*P1*P2*pow(Q1, 2)*Q2+9.0*((-1.0)*pow(P1, 2)+pow(P2, 2)) 
    *Q1*pow(Q2, 2)+6.0*P1*P2*pow(Q2, 3))) 
    ) + 
    (s3L - s3L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(pow(G, 2)*((-3.0)*((-4.0)+pow(P1, 2))*P2*Q1+5.0*pow(P2, 3)*Q1+ 
    3.0*P1*(4.0+pow(P1, 2))*Q2+3.0*P1*pow(P2, 2)*Q2)+(-7.0)*(3.0*P1*(4.0+pow(P1, 2))* 
    Q2*(pow(Q1, 2)+pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*(pow(Q1, 2)+pow(Q2, 2))+(-1.0)*P2*Q1*(( 
    (-20.0)+pow(P1, 2))*pow(Q1, 2)+3.0*(4.0+3.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(7.0*pow(Q1, 3)+(-1.0) 
    *Q1*pow(Q2, 2))))
    ) + 
    (s4L - s4L0) * ( 
    (-5.0/32.0)*pow(G, -2.0)*(9.0*pow(G, 2)*(pow(P1, 2)*Q1+(-1.0)*pow(P2, 2)*Q1+(-2.0)*P1* 
    P2*Q2)+(-14.0)*(((-2.0)+3.0*pow(P1, 2)+(-9.0)*pow(P2, 2))*pow(Q1, 3)+(-18.0)*P1*P2* 
    pow(Q1, 2)*Q2+3.0*(2.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))*Q1*pow(Q2, 2)+(-6.0)*P1*P2*pow(Q2, 3)) 
    )
    ) + 
    (s5L - s5L0) * ( 
    (1.0/16.0)*pow(G, -2.0)*(pow(G, 2)*((-9.0)*pow(P1, 2)*P2*Q1+3.0*pow(P2, 3)*Q1+(-3.0)* 
    pow(P1, 3)*Q2+9.0*P1*pow(P2, 2)*Q2)+(-7.0)*(3.0*P1*pow(P2, 2)*Q2*(9.0*pow(Q1, 2)+ 
    pow(Q2, 2))+(-3.0)*P2*Q1*(((-4.0)+3.0*pow(P1, 2))*pow(Q1, 2)+3.0*(4.0+pow(P1, 2))*pow(Q2, 2))+ 
    P1*Q2*(3.0*(12.0+pow(P1, 2))*pow(Q1, 2)+(-1.0)*(12.0+5.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*( 
    7.0*pow(Q1, 3)+(-9.0)*Q1*pow(Q2, 2)))) 
    ) + 
    (s6L - s6L0) * ( 
    (35.0/16.0)*pow(G, -2.0)*((-1.0)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+2.0*P1* 
    P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 2)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))
    ) + 
    (s7L - s7L0) * ( 
    (5.0/16.0)*pow(G, -2.0)*(3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)* 
    pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))
    ) ;

    return I;
}

std::vector<double>  perturbation_propagator::integralsJ5(const double &L0, const double &Lf, const double &P1, const double &P2, const double &Q1, const double &Q2) const{

    std::vector<double> I(8, 0.0);

    /* Pre-computations */
    double G = 1.0 + Q1 * Q1 + Q2 * Q2;
    double pi = constants::pi;
    double L_0 = L0 - 2.0 * pi * floor(L0 / (2.0 * pi));
    double L_f = Lf - 2.0 * pi * floor(Lf / (2.0 * pi));
    double sL = sin(L_f), sL0 = sin(L_0);
    double cL = cos(L_f), cL0 = cos(L_0);
    double c2L = cos(2.0 * L_f), c2L0 = cos(2.0 * L_0);
    double s2L = sin(2.0 * L_f), s2L0 = sin(2.0 * L_0);
    double s3L = sin(3.0 * L_f), s3L0 = sin(3.0 * L_0);
    double c3L = cos(3.0 * L_f), c3L0 = cos(3.0 * L_0);
    double s4L = sin(4.0 * L_f), s4L0 = sin(4.0 * L_0);
    double c4L = cos(4.0 * L_f), c4L0 = cos(4.0 * L_0);
    double s5L = sin(5.0 * L_f), s5L0 = sin(5.0 * L_0);
    double c5L = cos(5.0 * L_f), c5L0 = cos(5.0 * L_0);
    double s6L = sin(6.0 * L_f), s6L0 = sin(6.0 * L_0);
    double c6L = cos(6.0 * L_f), c6L0 = cos(6.0 * L_0);
    double s7L = sin(7.0 * L_f), s7L0 = sin(7.0 * L_0);
    double c7L = cos(7.0 * L_f), c7L0 = cos(7.0 * L_0);
    double s8L = sin(8.0 * L_f), s8L0 = sin(8.0 * L_0);
    double c8L = cos(8.0 * L_f), c8L0 = cos(8.0 * L_0);
    double s9L = sin(9.0 * L_f), s9L0 = sin(9.0 * L_0);
    double c9L = cos(9.0 * L_f), c9L0 = cos(9.0 * L_0);     
    double s10L = sin(10.0 * L_f), s10L0 = sin(10.0 * L_0);
    double c10L = cos(10.0 * L_f), c10L0 = cos(10.0 * L_0);
    double s11L = sin(11.0 * L_f), s11L0 = sin(11.0 * L_0);
    double c11L = cos(11.0 * L_f), c11L0 = cos(11.0 * L_0);          

    // a
    I[0]= (cL - cL0) * ( 
    (3.0/128.0)*pow(G, -4.0)*(5.0*pow(G, 4)*((64.0+5.0*pow(P1, 6)+720.0*pow(P2, 2)+600.0*pow(P2, 4)+ 
    35.0*pow(P2, 6)+15.0*pow(P1, 4)*(8.0+3.0*pow(P2, 2))+15.0*pow(P1, 2)*(16.0+48.0*pow(P2, 2)+5.0* 
    pow(P2, 4)))*Q1+(-30.0)*P1*P2*(16.0+pow(P1, 4)+16.0*pow(P2, 2)+pow(P2, 4)+2.0*pow(P1, 2)*(8.0+ 
    pow(P2, 2)))*Q2)+(-140.0)*pow(G, 2)*((32.0+pow(P1, 6)+400.0*pow(P2, 2)+350.0*pow(P2, 4)+21.0* 
    pow(P2, 6)+15.0*pow(P1, 4)*(2.0+pow(P2, 2))+5.0*pow(P1, 2)*(16.0+60.0*pow(P2, 2)+7.0*pow(P2, 4)))* 
    pow(Q1, 3)+(-6.0)*P1*P2*(80.0+3.0*pow(P1, 4)+100.0*pow(P2, 2)+7.0*pow(P2, 4)+10.0*pow(P1, 2)*( 
    6.0+pow(P2, 2)))*pow(Q1, 2)*Q2+(32.0+7.0*pow(P1, 6)+240.0*pow(P2, 2)+150.0*pow(P2, 4)+7.0*pow(P2, 6)+ 
    15.0*pow(P1, 4)*(10.0+3.0*pow(P2, 2))+15.0*pow(P1, 2)*(16.0+36.0*pow(P2, 2)+3.0*pow(P2, 4)))*Q1* 
    pow(Q2, 2)+(-2.0)*P1*P2*(80.0+7.0*pow(P1, 4)+60.0*pow(P2, 2)+3.0*pow(P2, 4)+10.0*pow(P1, 2)*( 
    10.0+pow(P2, 2)))*pow(Q2, 3))+42.0*((320.0+5.0*pow(P1, 6)+4200.0*pow(P2, 2)+3780.0*pow(P2, 4)+ 
    231.0*pow(P2, 6)+15.0*pow(P1, 4)*(12.0+7.0*pow(P2, 2))+15.0*pow(P1, 2)*(40.0+168.0*pow(P2, 2)+21.0* 
    pow(P2, 4)))*pow(Q1, 5)+(-10.0)*P1*P2*(600.0+15.0*pow(P1, 4)+840.0*pow(P2, 2)+63.0*pow(P2, 4)+ 
    10.0*pow(P1, 2)*(36.0+7.0*pow(P2, 2)))*pow(Q1, 4)*Q2+10.0*(64.0+7.0*pow(P1, 6)+600.0*pow(P2, 2)+ 
    420.0*pow(P2, 4)+21.0*pow(P2, 6)+15.0*pow(P1, 4)*(12.0+5.0*pow(P2, 2))+15.0*pow(P1, 2)*(24.0+72.0* 
    pow(P2, 2)+7.0*pow(P2, 4)))*pow(Q1, 3)*pow(Q2, 2)+(-20.0)*P1*P2*(360.0+21.0*pow(P1, 4)+360.0* 
    pow(P2, 2)+21.0*pow(P2, 4)+10.0*pow(P1, 2)*(36.0+5.0*pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 3)+5.0*(64.0+ 
    21.0*pow(P1, 6)+360.0*pow(P2, 2)+180.0*pow(P2, 4)+7.0*pow(P2, 6)+105.0*pow(P1, 4)*(4.0+pow(P2, 2))+ 
    15.0*pow(P1, 2)*(40.0+72.0*pow(P2, 2)+5.0*pow(P2, 4)))*Q1*pow(Q2, 4)+(-2.0)*P1*P2*(63.0* 
    pow(P1, 4)+70.0*pow(P1, 2)*(12.0+pow(P2, 2))+15.0*(40.0+24.0*pow(P2, 2)+pow(P2, 4)))*pow(Q2, 5)))
    ) + 
    (c2L - c2L0) * ( 
    (15.0/32.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*(16.0+(-5.0)*pow(P1, 4))*P2*Q1+10.0*(16.0+3.0* 
    pow(P1, 2))*pow(P2, 3)*Q1+45.0*pow(P2, 5)*Q1+P1*(48.0+160.0*pow(P1, 2)+45.0*pow(P1, 4))*Q2+ 
    30.0*pow(P1, 3)*pow(P2, 2)*Q2+(-15.0)*P1*pow(P2, 4)*Q2)+42.0*(3.0*P1*pow(P2, 4)*Q2* 
    ((-35.0)*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+3.0*pow(P2, 5)*Q1*(21.0*pow(Q1, 4)+ 
    14.0*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+6.0*P1*pow(P2, 2)*Q2*((-5.0)*(8.0+pow(P1, 2))*pow(Q1, 4)+ 
    10.0*pow(P1, 2)*pow(Q1, 2)*pow(Q2, 2)+(8.0+7.0*pow(P1, 2))*pow(Q2, 4))+3.0*P2*Q1*((24.0+16.0* 
    pow(P1, 2)+pow(P1, 4))*pow(Q1, 4)+(-2.0)*((-8.0)+5.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(8.0+ 
    80.0*pow(P1, 2)+35.0*pow(P1, 4))*pow(Q2, 4))+P1*Q2*(3.0*((-8.0)+pow(P1, 4))*pow(Q1, 4)+2.0*( 
    24.0+80.0*pow(P1, 2)+21.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(72.0+224.0*pow(P1, 2)+63.0*pow(P1, 4))* 
    pow(Q2, 4))+2.0*pow(P2, 3)*(7.0*(16.0+3.0*pow(P1, 2))*pow(Q1, 5)+10.0*(8.0+3.0*pow(P1, 2))* 
    pow(Q1, 3)*pow(Q2, 2)+(-15.0)*pow(P1, 2)*Q1*pow(Q2, 4)))+(-56.0)*pow(G, 2)*((-15.0)*P1* 
    pow(P2, 4)*pow(Q1, 2)*Q2+10.0*pow(P2, 3)*Q1*((5.0+pow(P1, 2))*pow(Q1, 2)+pow(Q2, 2))+10.0*P1* 
    pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+(1.0+pow(P1, 2))*pow(Q2, 2))+P2*Q1*(2.0*(8.0+5.0*pow(P1, 2)) 
    *pow(Q1, 2)+(-15.0)*pow(P1, 2)*(2.0+pow(P1, 2))*pow(Q2, 2))+pow(P2, 5)*(14.0*pow(Q1, 3)+3.0*Q1* 
    pow(Q2, 2))+P1*Q2*(16.0*pow(Q2, 2)+10.0*pow(P1, 2)*(pow(Q1, 2)+5.0*pow(Q2, 2))+pow(P1, 4)*(3.0* 
    pow(Q1, 2)+14.0*pow(Q2, 2)))))
    ) + 
    (c3L - c3L0) * ( 
    (-5.0/128.0)*pow(G, -4.0)*(9.0*pow(G, 4)*(3.0*pow(P1, 6)*Q1+15.0*pow(P1, 4)*(4.0+pow(P2, 2))* 
    Q1+5.0*pow(P1, 2)*(16.0+24.0*pow(P2, 2)+pow(P2, 4))*Q1+(-1.0)*pow(P2, 2)*(80.0+100.0*pow(P2, 2)+ 
    7.0*pow(P2, 4))*Q1+(-18.0)*pow(P1, 5)*P2*Q2+(-20.0)*pow(P1, 3)*P2*(12.0+pow(P2, 2))* 
    Q2+(-2.0)*P1*P2*(80.0+40.0*pow(P2, 2)+pow(P2, 4))*Q2)+(-56.0)*pow(G, 2)*(((-16.0)+ 
    2.0*pow(P1, 6)+(-300.0)*pow(P2, 2)+(-315.0)*pow(P2, 4)+(-21.0)*pow(P2, 6)+15.0*pow(P1, 4)*(3.0+ 
    pow(P2, 2))+30.0*pow(P1, 2)*(2.0+3.0*pow(P2, 2)))*pow(Q1, 3)+(-12.0)*P1*P2*(3.0*pow(P1, 4)+ 
    15.0*(2.0+pow(P2, 2))+5.0*pow(P1, 2)*(9.0+pow(P2, 2)))*pow(Q1, 2)*Q2+3.0*(16.0+7.0*pow(P1, 6)+60.0* 
    pow(P2, 2)+15.0*pow(P2, 4)+15.0*pow(P1, 4)*(9.0+2.0*pow(P2, 2))+15.0*pow(P1, 2)*(12.0+18.0*pow(P2, 2)+ 
    pow(P2, 4)))*Q1*pow(Q2, 2)+(-2.0)*P1*P2*(21.0*pow(P1, 4)+10.0*pow(P1, 2)*(27.0+2.0* 
    pow(P2, 2))+3.0*(60.0+30.0*pow(P2, 2)+pow(P2, 4)))*pow(Q2, 3))+126.0*(((-32.0)+pow(P1, 6)+(-504.0) 
    *pow(P2, 2)+(-504.0)*pow(P2, 4)+(-33.0)*pow(P2, 6)+3.0*pow(P1, 4)*(8.0+3.0*pow(P2, 2))+pow(P1, 2)*( 
    24.0+(-9.0)*pow(P2, 4)))*pow(Q1, 5)+(-6.0)*P1*P2*(40.0+5.0*pow(P1, 4)+(-3.0)*pow(P2, 4)+ 
    10.0*pow(P1, 2)*(8.0+pow(P2, 2)))*pow(Q1, 4)*Q2+2.0*(32.0+11.0*pow(P1, 6)+120.0*pow(P2, 2)+(-3.0) 
    *pow(P2, 6)+15.0*pow(P1, 4)*(16.0+5.0*pow(P2, 2))+45.0*pow(P1, 2)*(8.0+16.0*pow(P2, 2)+pow(P2, 4)))* 
    pow(Q1, 3)*pow(Q2, 2)+(-4.0)*P1*P2*(360.0+33.0*pow(P1, 4)+240.0*pow(P2, 2)+9.0*pow(P2, 4)+10.0* 
    pow(P1, 2)*(48.0+5.0*pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 3)+3.0*(32.0+15.0*pow(P1, 6)+120.0*pow(P2, 2)+ 
    40.0*pow(P2, 4)+pow(P2, 6)+5.0*pow(P1, 4)*(56.0+11.0*pow(P2, 2))+5.0*pow(P1, 2)*(72.0+96.0*pow(P2, 2)+ 
    5.0*pow(P2, 4)))*Q1*pow(Q2, 4)+(-2.0)*P1*P2*(27.0*pow(P1, 4)+pow(P1, 2)*(336.0+22.0* 
    pow(P2, 2))+3.0*(72.0+32.0*pow(P2, 2)+pow(P2, 4)))*pow(Q2, 5)))
    ) + 
    (c4L - c4L0) * ( 
    (-3.0/16.0)*pow(G, -4.0)*(5.0*pow(G, 4)*(15.0*pow(P1, 4)*P2*Q1+30.0*pow(P1, 2)*P2*(2.0+ 
    pow(P2, 2))*Q1+(-1.0)*pow(P2, 3)*(20.0+9.0*pow(P2, 2))*Q1+9.0*pow(P1, 5)*Q2+10.0*pow(P1, 3)* 
    (2.0+(-3.0)*pow(P2, 2))*Q2+(-15.0)*P1*pow(P2, 2)*(4.0+pow(P2, 2))*Q2)+140.0*pow(G, 2)*( 
    10.0*P1*(2.0+pow(P1, 2))*pow(P2, 2)*Q2*(3.0*pow(Q1, 2)+pow(Q2, 2))+5.0*P1*pow(P2, 4)*Q2*( 
    3.0*pow(Q1, 2)+pow(Q2, 2))+(-10.0)*pow(P2, 3)*Q1*(((-2.0)+pow(P1, 2))*pow(Q1, 2)+(2.0+3.0* 
    pow(P1, 2))*pow(Q2, 2))+(-1.0)*P2*Q1*(((-4.0)+20.0*pow(P1, 2)+5.0*pow(P1, 4))*pow(Q1, 2)+3.0* 
    (4.0+20.0*pow(P1, 2)+5.0*pow(P1, 4))*pow(Q2, 2))+P1*Q2*((12.0+20.0*pow(P1, 2)+3.0*pow(P1, 4))* 
    pow(Q1, 2)+(-1.0)*(4.0+20.0*pow(P1, 2)+7.0*pow(P1, 4))*pow(Q2, 2))+pow(P2, 5)*(7.0*pow(Q1, 3)+(-3.0)* 
    Q1*pow(Q2, 2)))+(-84.0)*(30.0*P1*(2.0+pow(P1, 2))*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+10.0* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+15.0*P1*pow(P2, 4)*Q2*(5.0*pow(Q1, 4)+10.0*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-10.0)*pow(P2, 3)*Q1*(((-14.0)+3.0*pow(P1, 2))*pow(Q1, 4)+10.0*(2.0+3.0*pow(P1, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+5.0*(2.0+3.0*pow(P1, 2))*pow(Q2, 4))+(-3.0)*P2*Q1*(((-12.0)+20.0* 
    pow(P1, 2)+5.0*pow(P1, 4))*pow(Q1, 4)+10.0*(4.0+20.0*pow(P1, 2)+5.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+ 
    5.0*(4.0+20.0*pow(P1, 2)+5.0*pow(P1, 4))*pow(Q2, 4))+P1*Q2*(5.0*(12.0+20.0*pow(P1, 2)+3.0* 
    pow(P1, 4))*pow(Q1, 4)+10.0*(12.0+20.0*pow(P1, 2)+3.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(36.0+ 
    140.0*pow(P1, 2)+45.0*pow(P1, 4))*pow(Q2, 4))+15.0*pow(P2, 5)*(3.0*pow(Q1, 5)+(-2.0)*pow(Q1, 3)* 
    pow(Q2, 2)+(-1.0)*Q1*pow(Q2, 4))))
    ) + 
    (c5L - c5L0) * ( 
    (3.0/128.0)*pow(G, -4.0)*(5.0*pow(G, 4)*(5.0*pow(P1, 6)*Q1+(-15.0)*pow(P1, 4)*((-4.0)+ 
    pow(P2, 2))*Q1+(-45.0)*pow(P1, 2)*pow(P2, 2)*(8.0+pow(P2, 2))*Q1+pow(P2, 4)*(60.0+7.0*pow(P2, 2)) 
    *Q1+(-30.0)*pow(P1, 5)*P2*Q2+20.0*pow(P1, 3)*P2*((-12.0)+pow(P2, 2))*Q2+6.0*P1* 
    pow(P2, 3)*(40.0+3.0*pow(P2, 2))*Q2)+(-21.0)*(((-64.0)+5.0*pow(P1, 6)+(-1680.0)*pow(P2, 2)+( 
    -2160.0)*pow(P2, 4)+(-165.0)*pow(P2, 6)+15.0*pow(P1, 4)*(16.0+11.0*pow(P2, 2))+45.0*pow(P1, 2)*( 
    16.0+64.0*pow(P2, 2)+7.0*pow(P2, 4)))*pow(Q1, 5)+(-10.0)*P1*P2*(720.0+15.0*pow(P1, 4)+960.0* 
    pow(P2, 2)+63.0*pow(P2, 4)+10.0*pow(P1, 2)*(48.0+11.0*pow(P2, 2)))*pow(Q1, 4)*Q2+10.0*(64.0+(-5.0) 
    *pow(P1, 6)+720.0*pow(P2, 2)+75.0*pow(P1, 4)*pow(P2, 2)+480.0*pow(P2, 4)+21.0*pow(P2, 6)+15.0* 
    pow(P1, 2)*(16.0+96.0*pow(P2, 2)+11.0*pow(P2, 4)))*pow(Q1, 3)*pow(Q2, 2)+20.0*P1*P2*(15.0* 
    pow(P1, 4)+(-50.0)*pow(P1, 2)*pow(P2, 2)+(-3.0)*(80.0+160.0*pow(P2, 2)+11.0*pow(P2, 4)))* 
    pow(Q1, 2)*pow(Q2, 3)+(-5.0)*(64.0+75.0*pow(P1, 6)+(-240.0)*pow(P2, 2)+(-240.0)*pow(P2, 4)+(-11.0) 
    *pow(P2, 6)+75.0*pow(P1, 4)*(16.0+pow(P2, 2))+(-75.0)*pow(P1, 2)*((-16.0)+pow(P2, 4)))*Q1* 
    pow(Q2, 4)+10.0*P1*P2*(45.0*pow(P1, 4)+10.0*pow(P1, 2)*(48.0+pow(P2, 2))+(-3.0)*((-80.0)+ 
    pow(P2, 4)))*pow(Q2, 5))+(-280.0)*pow(G, 2)*(5.0*pow(P1, 6)*Q1*pow(Q2, 2)+(-10.0)*pow(P1, 5)* 
    P2*pow(Q2, 3)+20.0*pow(P1, 3)*P2*Q2*((3.0+pow(P2, 2))*pow(Q1, 2)+(-5.0)*pow(Q2, 2))+(-5.0)* 
    pow(P1, 2)*Q1*(2.0*(2.0+9.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)+3.0*((-4.0)+6.0*pow(P2, 2)+pow(P2, 4)) 
    *pow(Q2, 2))+2.0*P1*P2*Q2*(6.0*(10.0+15.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)+((-20.0)+ 
    10.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+pow(P2, 2)*Q1*((20.0+35.0*pow(P2, 2)+3.0*pow(P2, 4))* 
    pow(Q1, 2)+(-1.0)*(60.0+45.0*pow(P2, 2)+2.0*pow(P2, 4))*pow(Q2, 2))+(-5.0)*pow(P1, 4)*((1.0+ 
    pow(P2, 2))*pow(Q1, 3)+(-15.0)*Q1*pow(Q2, 2))))
    ) + 
    (c6L - c6L0) * ( 
    (1.0/32.0)*pow(G, -4.0)*(45.0*pow(G, 4)*(5.0*pow(P1, 4)*P2*Q1+(-10.0)*pow(P1, 2)* 
    pow(P2, 3)*Q1+pow(P2, 5)*Q1+pow(P1, 5)*Q2+(-10.0)*pow(P1, 3)*pow(P2, 2)*Q2+5.0*P1* 
    pow(P2, 4)*Q2)+280.0*pow(G, 2)*((-45.0)*pow(P1, 4)*P2*Q1*pow(Q2, 2)+30.0*pow(P1, 2)*P2* 
    Q1*((1.0+pow(P2, 2))*pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-15.0)*P1*pow(P2, 2)*Q2*(3.0*(2.0+ 
    pow(P2, 2))*pow(Q1, 2)+(-2.0)*pow(Q2, 2))+pow(P2, 3)*Q1*((-2.0)*(5.0+3.0*pow(P2, 2))*pow(Q1, 2)+ 
    3.0*(10.0+3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 5)*(9.0*pow(Q1, 2)*Q2+(-6.0)*pow(Q2, 3))+10.0* 
    pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+((-1.0)+3.0*pow(P2, 2))*pow(Q2, 3)))+(-63.0)*(15.0*P1* 
    pow(P2, 4)*Q2*((-65.0)*pow(Q1, 4)+(-30.0)*pow(Q1, 2)*pow(Q2, 2)+3.0*pow(Q2, 4))+30.0*P1* 
    pow(P2, 2)*Q2*((-5.0)*(16.0+3.0*pow(P1, 2))*pow(Q1, 4)+30.0*pow(P1, 2)*pow(Q1, 2)*pow(Q2, 2)+( 
    16.0+13.0*pow(P1, 2))*pow(Q2, 4))+3.0*P2*Q1*(((-16.0)+160.0*pow(P1, 2)+15.0*pow(P1, 4))* 
    pow(Q1, 4)+10.0*(16.0+(-15.0)*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*(16.0+160.0*pow(P1, 2)+ 
    65.0*pow(P1, 4))*pow(Q2, 4))+P1*Q2*(15.0*((-16.0)+3.0*pow(P1, 4))*pow(Q1, 4)+10.0*(48.0+ 
    160.0*pow(P1, 2)+39.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(48.0+320.0*pow(P1, 2)+135.0* 
    pow(P1, 4))*pow(Q2, 4))+(-15.0)*pow(P2, 5)*(9.0*pow(Q1, 5)+(-26.0)*pow(Q1, 3)*pow(Q2, 2)+(-3.0)* 
    Q1*pow(Q2, 4))+10.0*pow(P2, 3)*(((-32.0)+39.0*pow(P1, 2))*pow(Q1, 5)+10.0*(16.0+9.0*pow(P1, 2)) 
    *pow(Q1, 3)*pow(Q2, 2)+(-45.0)*pow(P1, 2)*Q1*pow(Q2, 4))))
    ) + 
    (c7L - c7L0) * ( 
    (-15.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 6)*Q1+(-15.0)*pow(P1, 4)*pow(P2, 2)*Q1+15.0* 
    pow(P1, 2)*pow(P2, 4)*Q1+(-1.0)*pow(P2, 6)*Q1+(-6.0)*pow(P1, 5)*P2*Q2+20.0*pow(P1, 3)* 
    pow(P2, 3)*Q2+(-6.0)*P1*pow(P2, 5)*Q2)+14.0*pow(G, 2)*(2.0*pow(P1, 5)*P2*Q2*((-9.0) 
    *pow(Q1, 2)+7.0*pow(Q2, 2))+(-20.0)*pow(P1, 3)*P2*Q2*((12.0+pow(P2, 2))*pow(Q1, 2)+((-4.0)+ 
    pow(P2, 2))*pow(Q2, 2))+pow(P2, 4)*Q1*((20.0+3.0*pow(P2, 2))*pow(Q1, 2)+(-5.0)*(12.0+pow(P2, 2))* 
    pow(Q2, 2))+(-5.0)*pow(P1, 2)*pow(P2, 2)*Q1*((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(24.0+ 
    pow(P2, 2))*pow(Q2, 2))+(-2.0)*P1*pow(P2, 3)*Q2*((-15.0)*(8.0+pow(P2, 2))*pow(Q1, 2)+(40.0+ 
    pow(P2, 2))*pow(Q2, 2))+5.0*pow(P1, 4)*Q1*((4.0+pow(P2, 2))*pow(Q1, 2)+3.0*((-4.0)+3.0*pow(P2, 2)) 
    *pow(Q2, 2))+pow(P1, 6)*(pow(Q1, 3)+(-7.0)*Q1*pow(Q2, 2)))+(-21.0)*(6.0*pow(P1, 5)*P2* 
    Q2*((-5.0)*pow(Q1, 4)+(-14.0)*pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+(-4.0)*pow(P1, 3)*P2* 
    Q2*((-5.0)*((-12.0)+pow(P2, 2))*pow(Q1, 4)+10.0*(36.0+5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+ 
    7.0*((-12.0)+pow(P2, 2))*pow(Q2, 4))+(-3.0)*pow(P1, 4)*Q1*(((-4.0)+pow(P2, 2))*pow(Q1, 4)+( 
    -10.0)*(12.0+5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-35.0)*((-4.0)+pow(P2, 2))*pow(Q2, 4))+ 
    pow(P2, 2)*Q1*((48.0+108.0*pow(P2, 2)+11.0*pow(P2, 4))*pow(Q1, 4)+(-2.0)*(240.0+300.0* 
    pow(P2, 2)+19.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(240.0+60.0*pow(P2, 2)+(-1.0)*pow(P2, 4))*pow(Q2, 4)) 
    +6.0*P1*P2*Q2*((80.0+200.0*pow(P2, 2)+19.0*pow(P2, 4))*pow(Q1, 4)+2.0*((-80.0)+(-40.0) 
    *pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*((-16.0)+24.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4)) 
    +(-3.0)*pow(P1, 2)*Q1*((16.0+120.0*pow(P2, 2)+19.0*pow(P2, 4))*pow(Q1, 4)+10.0*((-16.0)+( 
    -24.0)*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*((-16.0)+72.0*pow(P2, 2)+5.0*pow(P2, 4)) 
    *pow(Q2, 4))+pow(P1, 6)*(pow(Q1, 5)+14.0*pow(Q1, 3)*pow(Q2, 2)+(-35.0)*Q1*pow(Q2, 4))))
    ) + 
    (c8L - c8L0) * ( 
    (-105.0/16.0)*pow(G, -4.0)*(pow(G, 2)*(5.0*pow(P1, 4)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+ 
    (-10.0)*pow(P1, 2)*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+pow(P2, 5)*Q1*(pow(Q1, 2)+( 
    -3.0)*pow(Q2, 2))+10.0*pow(P1, 3)*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-5.0)*P1* 
    pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 5)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))+ 
    3.0*((-3.0)*pow(P1, 4)*P2*Q1*(pow(Q1, 4)+10.0*pow(Q1, 2)*pow(Q2, 2)+(-15.0)*pow(Q2, 4))+3.0* 
    pow(P1, 5)*Q2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+2.0*pow(P1, 3)*Q2*(5.0*(2.0+ 
    3.0*pow(P2, 2))*pow(Q1, 4)+10.0*((-2.0)+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(2.0+(-9.0)*pow(P2, 2)) 
    *pow(Q2, 4))+3.0*P1*pow(P2, 2)*Q2*((-5.0)*(4.0+3.0*pow(P2, 2))*pow(Q1, 4)+10.0*(4.0+ 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-4.0)+pow(P2, 2))*pow(Q2, 4))+6.0*pow(P1, 2)*P2*Q1*((2.0+ 
    3.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*((-2.0)+pow(P2, 2)) 
    *pow(Q2, 4))+(-1.0)*pow(P2, 3)*Q1*((4.0+3.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*(20.0+9.0*pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+(20.0+3.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c9L - c9L0) * ( 
    (7.0/128.0)*pow(G, -4.0)*(10.0*pow(G, 2)*((-15.0)*pow(P1, 4)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0) 
    *pow(Q2, 2))+15.0*pow(P1, 2)*pow(P2, 4)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-1.0)*pow(P2, 6)* 
    Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+6.0*pow(P1, 5)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+( 
    -20.0)*pow(P1, 3)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+6.0*P1*pow(P2, 5)*Q2*(( 
    -3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 6)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))+9.0*((-6.0)* 
    pow(P1, 5)*P2*Q2*(5.0*pow(Q1, 4)+(-50.0)*pow(Q1, 2)*pow(Q2, 2)+9.0*pow(Q2, 4))+15.0*pow(P1, 4)* 
    Q1*((4.0+3.0*pow(P2, 2))*pow(Q1, 4)+10.0*((-4.0)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(4.0+(-5.0) 
    *pow(P2, 2))*pow(Q2, 4))+(-20.0)*pow(P1, 3)*P2*Q2*(15.0*(4.0+pow(P2, 2))*pow(Q1, 4)+10.0*( 
    (-12.0)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(12.0+(-5.0)*pow(P2, 2))*pow(Q2, 4))+(-6.0)*P1* 
    pow(P2, 3)*Q2*((-5.0)*(40.0+7.0*pow(P2, 2))*pow(Q1, 4)+10.0*(40.0+3.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+((-40.0)+pow(P2, 2))*pow(Q2, 4))+(-15.0)*pow(P1, 2)*pow(P2, 2)*Q1*((24.0+7.0* 
    pow(P2, 2))*pow(Q1, 4)+(-30.0)*(8.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*((-24.0)+pow(P2, 2))* 
    pow(Q2, 4))+pow(P2, 4)*Q1*((60.0+11.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(60.0+7.0*pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+15.0*(20.0+pow(P2, 2))*pow(Q2, 4))+pow(P1, 6)*(pow(Q1, 5)+(-50.0)*pow(Q1, 3)* 
    pow(Q2, 2)+45.0*Q1*pow(Q2, 4))))
    ) + 
    (c10L - c10L0) * ( 
    (189.0/32.0)*pow(G, -4.0)*((-10.0)*pow(P1, 3)*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+(-10.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+5.0*P1*pow(P2, 4)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+5.0*pow(P1, 4)*P2*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+(-10.0)*pow(P1, 2)*pow(P2, 3)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+pow(P2, 5)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 5)*(5.0* 
    pow(Q1, 4)*Q2+(-10.0)*pow(Q1, 2)*pow(Q2, 3)+pow(Q2, 5)))
    ) + 
    (c11L - c11L0) * ( 
    (-63.0/128.0)*pow(G, -4.0)*((-6.0)*pow(P1, 5)*P2*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+20.0*pow(P1, 3)*pow(P2, 3)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-6.0)*P1*pow(P2, 5)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+( 
    -15.0)*pow(P1, 4)*pow(P2, 2)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+15.0* 
    pow(P1, 2)*pow(P2, 4)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+(-1.0)* 
    pow(P2, 6)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 6)*(pow(Q1, 5)+( 
    -10.0)*pow(Q1, 3)*pow(Q2, 2)+5.0*Q1*pow(Q2, 4)))
    ) + 
    (sL - sL0) * ( 
    (-3.0/128.0)*pow(G, -4.0)*(5.0*pow(G, 4)*((-30.0)*pow(P1, 5)*P2*Q1+(-60.0)*pow(P1, 3)* 
    P2*(8.0+pow(P2, 2))*Q1+(-30.0)*P1*P2*(16.0+16.0*pow(P2, 2)+pow(P2, 4))*Q1+35.0* 
    pow(P1, 6)*Q2+75.0*pow(P1, 4)*(8.0+pow(P2, 2))*Q2+45.0*pow(P1, 2)*(16.0+16.0*pow(P2, 2)+pow(P2, 4)) 
    *Q2+(64.0+240.0*pow(P2, 2)+120.0*pow(P2, 4)+5.0*pow(P2, 6))*Q2)+(-140.0)*pow(G, 2)*(7.0* 
    pow(P1, 6)*Q2*(pow(Q1, 2)+3.0*pow(Q2, 2))+(-6.0)*pow(P1, 5)*P2*Q1*(pow(Q1, 2)+7.0*pow(Q2, 2))+ 
    (-20.0)*pow(P1, 3)*P2*Q1*((6.0+pow(P2, 2))*pow(Q1, 2)+3.0*(10.0+pow(P2, 2))*pow(Q2, 2))+5.0* 
    pow(P1, 4)*Q2*((30.0+9.0*pow(P2, 2))*pow(Q1, 2)+7.0*(10.0+pow(P2, 2))*pow(Q2, 2))+5.0*pow(P1, 2)* 
    Q2*(3.0*(16.0+36.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q1, 2)+(80.0+60.0*pow(P2, 2)+3.0*pow(P2, 4))* 
    pow(Q2, 2))+(-2.0)*P1*P2*Q1*((80.0+100.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 2)+3.0*(80.0+ 
    60.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q2, 2))+Q2*((32.0+240.0*pow(P2, 2)+150.0*pow(P2, 4)+7.0* 
    pow(P2, 6))*pow(Q1, 2)+(32.0+80.0*pow(P2, 2)+30.0*pow(P2, 4)+pow(P2, 6))*pow(Q2, 2)))+42.0*((-30.0)* 
    pow(P1, 5)*P2*Q1*(pow(Q1, 4)+14.0*pow(Q1, 2)*pow(Q2, 2)+21.0*pow(Q2, 4))+15.0*pow(P1, 4)*Q2*( 
    5.0*(12.0+5.0*pow(P2, 2))*pow(Q1, 4)+70.0*(4.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+21.0*(12.0+pow(P2, 2)) 
    *pow(Q2, 4))+(-20.0)*pow(P1, 3)*P2*Q1*((36.0+7.0*pow(P2, 2))*pow(Q1, 4)+10.0*(36.0+5.0* 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+35.0*(12.0+pow(P2, 2))*pow(Q2, 4))+15.0*pow(P1, 2)*Q2*(5.0*( 
    24.0+72.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 4)+10.0*(40.0+72.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+7.0*(40.0+24.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+(-6.0)*P1*P2*Q1*((200.0+ 
    280.0*pow(P2, 2)+21.0*pow(P2, 4))*pow(Q1, 4)+10.0*(120.0+120.0*pow(P2, 2)+7.0*pow(P2, 4))* 
    pow(Q1, 2)*pow(Q2, 2)+25.0*(40.0+24.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+5.0*Q2*((64.0+600.0* 
    pow(P2, 2)+420.0*pow(P2, 4)+21.0*pow(P2, 6))*pow(Q1, 4)+2.0*(64.0+360.0*pow(P2, 2)+180.0*pow(P2, 4)+ 
    7.0*pow(P2, 6))*pow(Q1, 2)*pow(Q2, 2)+(64.0+120.0*pow(P2, 2)+36.0*pow(P2, 4)+pow(P2, 6))*pow(Q2, 4))+ 
    7.0*pow(P1, 6)*(5.0*pow(Q1, 4)*Q2+30.0*pow(Q1, 2)*pow(Q2, 3)+33.0*pow(Q2, 5))))
    ) + 
    (s2L - s2L0) * ( 
    (15.0/32.0)*pow(G, -4.0)*(pow(G, 4)*(15.0*pow(P1, 5)*Q1+10.0*pow(P1, 3)*(8.0+9.0*pow(P2, 2))* 
    Q1+3.0*P1*(16.0+80.0*pow(P2, 2)+25.0*pow(P2, 4))*Q1+(-75.0)*pow(P1, 4)*P2*Q2+(-30.0) 
    *pow(P1, 2)*P2*(8.0+3.0*pow(P2, 2))*Q2+(-1.0)*P2*(48.0+80.0*pow(P2, 2)+15.0*pow(P2, 4))* 
    Q2)+(-28.0)*pow(G, 2)*((-5.0)*pow(P1, 4)*P2*Q2*(9.0*pow(Q1, 2)+7.0*pow(Q2, 2))+(-10.0)* 
    pow(P1, 2)*P2*Q2*(9.0*(2.0+pow(P2, 2))*pow(Q1, 2)+(10.0+3.0*pow(P2, 2))*pow(Q2, 2))+10.0* 
    pow(P1, 3)*Q1*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+(10.0+9.0*pow(P2, 2))*pow(Q2, 2))+(-1.0)*P2* 
    Q2*((48.0+100.0*pow(P2, 2)+21.0*pow(P2, 4))*pow(Q1, 2)+(16.0+20.0*pow(P2, 2)+3.0*pow(P2, 4))* 
    pow(Q2, 2))+P1*Q1*((16.0+100.0*pow(P2, 2)+35.0*pow(P2, 4))*pow(Q1, 2)+3.0*(16.0+60.0*pow(P2, 2)+ 
    15.0*pow(P2, 4))*pow(Q2, 2))+3.0*pow(P1, 5)*(pow(Q1, 3)+7.0*Q1*pow(Q2, 2)))+42.0*((-3.0)* 
    pow(P1, 4)*P2*Q2*(25.0*pow(Q1, 4)+70.0*pow(Q1, 2)*pow(Q2, 2)+21.0*pow(Q2, 4))+(-6.0)* 
    pow(P1, 2)*P2*Q2*(5.0*(12.0+7.0*pow(P2, 2))*pow(Q1, 4)+10.0*(12.0+5.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+7.0*(4.0+pow(P2, 2))*pow(Q2, 4))+(-1.0)*P2*Q2*((120.0+280.0*pow(P2, 2)+63.0* 
    pow(P2, 4))*pow(Q1, 4)+6.0*(24.0+40.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+3.0*(8.0+8.0* 
    pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+3.0*P1*Q1*((8.0+56.0*pow(P2, 2)+21.0*pow(P2, 4))*pow(Q1, 4)+ 
    2.0*(24.0+120.0*pow(P2, 2)+35.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+5.0*(8.0+24.0*pow(P2, 2)+5.0* 
    pow(P2, 4))*pow(Q2, 4))+3.0*pow(P1, 5)*(pow(Q1, 5)+14.0*pow(Q1, 3)*pow(Q2, 2)+21.0*Q1*pow(Q2, 4))+ 
    pow(P1, 3)*(6.0*(4.0+7.0*pow(P2, 2))*pow(Q1, 5)+60.0*(4.0+5.0*pow(P2, 2))*pow(Q1, 3)*pow(Q2, 2)+70.0* 
    (4.0+3.0*pow(P2, 2))*Q1*pow(Q2, 4))))
    ) + 
    (s3L - s3L0) * ( 
    (5.0/128.0)*pow(G, -4.0)*(9.0*pow(G, 4)*(2.0*pow(P1, 5)*P2*Q1+20.0*pow(P1, 3)*P2*(4.0+ 
    pow(P2, 2))*Q1+2.0*P1*P2*(80.0+120.0*pow(P2, 2)+9.0*pow(P2, 4))*Q1+7.0*pow(P1, 6)*Q2+( 
    -5.0)*pow(P1, 4)*((-20.0)+pow(P2, 2))*Q2+(-5.0)*pow(P1, 2)*((-16.0)+24.0*pow(P2, 2)+3.0* 
    pow(P2, 4))*Q2+(-1.0)*pow(P2, 2)*(80.0+60.0*pow(P2, 2)+3.0*pow(P2, 4))*Q2)+(-56.0)*pow(G, 2)* 
    (6.0*pow(P1, 5)*P2*pow(Q1, 3)+21.0*pow(P1, 6)*pow(Q2, 3)+(-45.0)*pow(P1, 4)*Q2*((1.0+pow(P2, 2)) 
    *pow(Q1, 2)+(-7.0)*pow(Q2, 2))+20.0*pow(P1, 3)*P2*Q1*((9.0+2.0*pow(P2, 2))*pow(Q1, 2)+3.0*( 
    3.0+pow(P2, 2))*pow(Q2, 2))+(-15.0)*pow(P1, 2)*Q2*(6.0*(2.0+9.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)+( 
    (-20.0)+6.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+6.0*P1*P2*Q1*((60.0+90.0*pow(P2, 2)+7.0* 
    pow(P2, 4))*pow(Q1, 2)+6.0*(10.0+15.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+(-1.0)*Q2*(3.0*(16.0+ 
    180.0*pow(P2, 2)+135.0*pow(P2, 4)+7.0*pow(P2, 6))*pow(Q1, 2)+((-16.0)+60.0*pow(P2, 2)+45.0*pow(P2, 4)+ 
    2.0*pow(P2, 6))*pow(Q2, 2)))+(-126.0)*(3.0*pow(P1, 6)*Q2*(pow(Q1, 4)+(-2.0)*pow(Q1, 2)* 
    pow(Q2, 2)+(-11.0)*pow(Q2, 4))+(-6.0)*pow(P1, 5)*P2*Q1*(pow(Q1, 4)+6.0*pow(Q1, 2)*pow(Q2, 2)+( 
    -3.0)*pow(Q2, 4))+(-4.0)*pow(P1, 3)*P2*Q1*((48.0+11.0*pow(P2, 2))*pow(Q1, 4)+10.0*(24.0+ 
    5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+15.0*pow(P2, 2)*pow(Q2, 4))+3.0*pow(P1, 4)*Q2*(5.0*(8.0+5.0* 
    pow(P2, 2))*pow(Q1, 4)+30.0*pow(P2, 2)*pow(Q1, 2)*pow(Q2, 2)+(-3.0)*(56.0+pow(P2, 2))*pow(Q2, 4))+3.0* 
    pow(P1, 2)*Q2*(5.0*(24.0+96.0*pow(P2, 2)+11.0*pow(P2, 4))*pow(Q1, 4)+10.0*(8.0+48.0*pow(P2, 2)+ 
    5.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+3.0*((-56.0)+pow(P2, 4))*pow(Q2, 4))+(-6.0)*P1*P2* 
    Q1*((72.0+112.0*pow(P2, 2)+9.0*pow(P2, 4))*pow(Q1, 4)+2.0*(120.0+160.0*pow(P2, 2)+11.0*pow(P2, 4)) 
    *pow(Q1, 2)*pow(Q2, 2)+5.0*(8.0+16.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+Q2*(3.0*(32.0+360.0* 
    pow(P2, 2)+280.0*pow(P2, 4)+15.0*pow(P2, 6))*pow(Q1, 4)+2.0*(32.0+360.0*pow(P2, 2)+240.0*pow(P2, 4)+ 
    11.0*pow(P2, 6))*pow(Q1, 2)*pow(Q2, 2)+((-32.0)+24.0*pow(P2, 2)+24.0*pow(P2, 4)+pow(P2, 6))*pow(Q2, 4)) 
    ))
    ) + 
    (s4L - s4L0) * ( 
    (-3.0/8.0)*pow(G, -4.0)*(5.0*pow(G, 4)*(10.0*pow(P1, 3)*Q1+3.0*pow(P1, 5)*Q1+(-15.0)*P1* 
    pow(P2, 2)*(2.0+pow(P2, 2))*Q1+(-30.0)*pow(P1, 2)*P2*Q2+(-15.0)*pow(P1, 4)*P2*Q2+ 
    pow(P2, 3)*(10.0+3.0*pow(P2, 2))*Q2)+84.0*(15.0*pow(P1, 5)*Q1*pow(Q2, 2)*(pow(Q1, 2)+3.0* 
    pow(Q2, 2))+(-15.0)*pow(P1, 4)*P2*pow(Q2, 3)*(5.0*pow(Q1, 2)+3.0*pow(Q2, 2))+(-15.0)*pow(P1, 2)* 
    P2*Q2*((-5.0)*(1.0+pow(P2, 2))*pow(Q1, 4)+10.0*pow(Q1, 2)*pow(Q2, 2)+(7.0+pow(P2, 2))*pow(Q2, 4)) 
    +(-5.0)*pow(P1, 3)*Q1*((1.0+3.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*( 
    7.0+3.0*pow(P2, 2))*pow(Q2, 4))+(-3.0)*P1*Q1*((4.0+35.0*pow(P2, 2)+15.0*pow(P2, 4))*pow(Q1, 4)+ 
    25.0*pow(P2, 2)*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*(4.0+5.0*pow(P2, 2))*pow(Q2, 4))+P2* 
    Q2*(5.0*(12.0+35.0*pow(P2, 2)+9.0*pow(P2, 4))*pow(Q1, 4)+5.0*pow(P2, 2)*(10.0+3.0*pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(12.0+5.0*pow(P2, 2))*pow(Q2, 4)))+(-35.0)*pow(G, 2)*((-5.0)* 
    pow(P1, 4)*P2*Q2*(3.0*pow(Q1, 2)+7.0*pow(Q2, 2))+P2*Q2*((24.0+80.0*pow(P2, 2)+21.0* 
    pow(P2, 4))*pow(Q1, 2)+((-8.0)+pow(P2, 4))*pow(Q2, 2))+(-1.0)*P1*Q1*((8.0+80.0*pow(P2, 2)+ 
    35.0*pow(P2, 4))*pow(Q1, 2)+3.0*((-8.0)+5.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 5)*(pow(Q1, 3)+21.0* 
    Q1*pow(Q2, 2))+(-10.0)*pow(P1, 3)*Q1*((-8.0)*pow(Q2, 2)+pow(P2, 2)*(pow(Q1, 2)+(-3.0)* 
    pow(Q2, 2)))+(-10.0)*pow(P1, 2)*(8.0*P2*pow(Q2, 3)+pow(P2, 3)*((-3.0)*pow(Q1, 2)*Q2+pow(Q2, 3)) 
    ))) 
    ) + 
    (s5L - s5L0) * ( 
    (-3.0/128.0)*pow(G, -4.0)*(5.0*pow(G, 4)*(18.0*pow(P1, 5)*P2*Q1+(-30.0)*P1*pow(P2, 3)* 
    (8.0+pow(P2, 2))*Q1+20.0*pow(P1, 3)*P2*(12.0+pow(P2, 2))*Q1+7.0*pow(P1, 6)*Q2+15.0* 
    pow(P1, 4)*(4.0+(-3.0)*pow(P2, 2))*Q2+5.0*pow(P2, 4)*(12.0+pow(P2, 2))*Q2+(-15.0)*pow(P1, 2)* 
    pow(P2, 2)*(24.0+pow(P2, 2))*Q2)+280.0*pow(G, 2)*((-2.0)*pow(P1, 5)*P2*Q1*(pow(Q1, 2)+6.0* 
    pow(Q2, 2))+10.0*P1*P2*Q1*((4.0+10.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)+(-6.0)*(2.0+pow(P2, 2)) 
    *pow(Q2, 2))+(-5.0)*pow(P2, 2)*Q2*((12.0+15.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)+(-1.0)*(4.0+ 
    pow(P2, 2))*pow(Q2, 2))+(-20.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+(9.0+pow(P2, 2))*pow(Q2, 2))+5.0* 
    pow(P1, 4)*Q2*(3.0*(3.0+pow(P2, 2))*pow(Q1, 2)+((-7.0)+2.0*pow(P2, 2))*pow(Q2, 2))+5.0*pow(P1, 2)* 
    Q2*(6.0*(2.0+3.0*pow(P2, 2))*pow(Q1, 2)+((-4.0)+18.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+pow(P1, 6)* 
    (2.0*pow(Q1, 2)*Q2+(-3.0)*pow(Q2, 3)))+(-21.0)*((-30.0)*pow(P1, 5)*P2*Q1*(pow(Q1, 4)+ 
    22.0*pow(Q1, 2)*pow(Q2, 2)+21.0*pow(Q2, 4))+15.0*pow(P1, 4)*Q2*(5.0*(16.0+5.0*pow(P2, 2))* 
    pow(Q1, 4)+10.0*(32.0+11.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+3.0*((-48.0)+7.0*pow(P2, 2))*pow(Q2, 4)) 
    +150.0*P1*P2*Q1*((16.0+32.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q1, 4)+2.0*((-16.0)+pow(P2, 4)) 
    *pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(48.0+32.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+15.0*pow(P1, 2)*Q2*( 
    (80.0+(-25.0)*pow(P2, 4))*pow(Q1, 4)+10.0*(48.0+96.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+((-112.0)+192.0*pow(P2, 2)+11.0*pow(P2, 4))*pow(Q2, 4))+Q2*((-5.0)*(64.0+1200.0* 
    pow(P2, 2)+1200.0*pow(P2, 4)+75.0*pow(P2, 6))*pow(Q1, 4)+10.0*(64.0+240.0*pow(P2, 2)+(-5.0)* 
    pow(P2, 6))*pow(Q1, 2)*pow(Q2, 2)+((-64.0)+720.0*pow(P2, 2)+240.0*pow(P2, 4)+5.0*pow(P2, 6))* 
    pow(Q2, 4))+5.0*pow(P1, 6)*(11.0*pow(Q1, 4)*Q2+42.0*pow(Q1, 2)*pow(Q2, 3)+(-33.0)*pow(Q2, 5))+ 
    100.0*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+pow(Q2, 2))*((-96.0)*pow(Q2, 2)+pow(P2, 2)*(pow(Q1, 2)+( 
    -11.0)*pow(Q2, 2)))))  
    ) + 
    (s6L - s6L0) * ( 
    (1.0/32.0)*pow(G, -4.0)*(45.0*pow(G, 4)*(pow(P1, 5)*Q1+(-10.0)*pow(P1, 3)*pow(P2, 2)*Q1+5.0* 
    P1*pow(P2, 4)*Q1+(-5.0)*pow(P1, 4)*P2*Q2+10.0*pow(P1, 2)*pow(P2, 3)*Q2+(-1.0)* 
    pow(P2, 5)*Q2)+140.0*pow(G, 2)*(45.0*pow(P1, 4)*P2*Q2*((-1.0)*pow(Q1, 2)+pow(Q2, 2))+( 
    -30.0)*pow(P1, 2)*P2*Q2*(3.0*(2.0+pow(P2, 2))*pow(Q1, 2)+((-2.0)+pow(P2, 2))*pow(Q2, 2))+( 
    -15.0)*P1*pow(P2, 2)*Q1*((4.0+3.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(4.0+pow(P2, 2))*pow(Q2, 2))+ 
    10.0*pow(P1, 3)*Q1*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+3.0*((-2.0)+3.0*pow(P2, 2))*pow(Q2, 2))+ 
    pow(P2, 3)*Q2*(3.0*(20.0+9.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(20.0+3.0*pow(P2, 2))*pow(Q2, 2))+ 
    3.0*pow(P1, 5)*(pow(Q1, 3)+(-9.0)*Q1*pow(Q2, 2)))+(-63.0)*(15.0*pow(P1, 4)*P2*Q2*(( 
    -25.0)*pow(Q1, 4)+(-30.0)*pow(Q1, 2)*pow(Q2, 2)+27.0*pow(Q2, 4))+(-30.0)*pow(P1, 2)*P2*Q2*( 
    5.0*(8.0+3.0*pow(P2, 2))*pow(Q1, 4)+10.0*(8.0+5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+3.0*((-8.0)+ 
    pow(P2, 2))*pow(Q2, 4))+10.0*pow(P1, 3)*Q1*((8.0+9.0*pow(P2, 2))*pow(Q1, 4)+10.0*(8.0+15.0* 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+15.0*((-8.0)+3.0*pow(P2, 2))*pow(Q2, 4))+P2*Q2*(15.0*(16.0+ 
    80.0*pow(P2, 2)+27.0*pow(P2, 4))*pow(Q1, 4)+(-10.0)*(48.0+80.0*pow(P2, 2)+9.0*pow(P2, 4))* 
    pow(Q1, 2)*pow(Q2, 2)+(48.0+(-80.0)*pow(P2, 2)+(-15.0)*pow(P2, 4))*pow(Q2, 4))+(-3.0)*P1*Q1* 
    ((16.0+240.0*pow(P2, 2)+135.0*pow(P2, 4))*pow(Q1, 4)+(-10.0)*(16.0+80.0*pow(P2, 2)+15.0*pow(P2, 4)) 
    *pow(Q1, 2)*pow(Q2, 2)+(-5.0)*((-16.0)+80.0*pow(P2, 2)+25.0*pow(P2, 4))*pow(Q2, 4))+15.0* 
    pow(P1, 5)*(pow(Q1, 5)+6.0*pow(Q1, 3)*pow(Q2, 2)+(-27.0)*Q1*pow(Q2, 4)))) 
    ) + 
    (s7L - s7L0) * ( 
    (15.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(6.0*pow(P1, 5)*P2*Q1+(-20.0)*pow(P1, 3)*pow(P2, 3)* 
    Q1+6.0*P1*pow(P2, 5)*Q1+pow(P1, 6)*Q2+(-15.0)*pow(P1, 4)*pow(P2, 2)*Q2+15.0*pow(P1, 2)* 
    pow(P2, 4)*Q2+(-1.0)*pow(P2, 6)*Q2)+14.0*pow(G, 2)*(2.0*pow(P1, 5)*P2*Q1*(pow(Q1, 2)+( 
    -15.0)*pow(Q2, 2))+(-5.0)*pow(P1, 2)*pow(P2, 2)*Q2*(9.0*(8.0+pow(P2, 2))*pow(Q1, 2)+((-24.0)+ 
    pow(P2, 2))*pow(Q2, 2))+20.0*pow(P1, 3)*P2*Q1*((4.0+pow(P2, 2))*pow(Q1, 2)+((-12.0)+pow(P2, 2)) 
    *pow(Q2, 2))+pow(P2, 4)*Q2*((60.0+7.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(20.0+pow(P2, 2))*pow(Q2, 2)) 
    +(-2.0)*P1*pow(P2, 3)*Q1*((40.0+7.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(40.0+3.0*pow(P2, 2))* 
    pow(Q2, 2))+5.0*pow(P1, 4)*Q2*((-3.0)*((-4.0)+pow(P2, 2))*pow(Q1, 2)+((-4.0)+5.0*pow(P2, 2))* 
    pow(Q2, 2))+pow(P1, 6)*(5.0*pow(Q1, 2)*Q2+(-3.0)*pow(Q2, 3)))+(-21.0)*(6.0*pow(P1, 5)*P2* 
    Q1*(pow(Q1, 4)+(-2.0)*pow(Q1, 2)*pow(Q2, 2)+(-19.0)*pow(Q2, 4))+pow(P1, 6)*Q2*(pow(Q1, 4)+38.0* 
    pow(Q1, 2)*pow(Q2, 2)+(-11.0)*pow(Q2, 4))+4.0*pow(P1, 3)*P2*Q1*((36.0+7.0*pow(P2, 2))* 
    pow(Q1, 4)+10.0*(12.0+5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*(60.0+pow(P2, 2))*pow(Q2, 4))+( 
    -6.0)*P1*P2*Q1*((16.0+56.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 4)+(-2.0)*(80.0+120.0* 
    pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*((-16.0)+8.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4)) 
    +(-1.0)*pow(P2, 2)*Q2*((-5.0)*(48.0+84.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 4)+2.0*(240.0+ 
    180.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+((-48.0)+12.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4)) 
    +(-3.0)*pow(P1, 2)*Q2*(5.0*(16.0+72.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 4)+10.0*((-16.0)+ 
    24.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*((-16.0)+120.0*pow(P2, 2)+pow(P2, 4))* 
    pow(Q2, 4))+pow(P1, 4)*((-15.0)*(4.0+5.0*pow(P2, 2))*pow(Q1, 4)*Q2+30.0*(20.0+pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 3)+3.0*((-36.0)+19.0*pow(P2, 2))*pow(Q2, 5))))
    ) + 
    (s8L - s8L0) * ( 
    (-105.0/16.0)*pow(G, -4.0)*(pow(G, 2)*((-10.0)*pow(P1, 3)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)* 
    pow(Q2, 2))+5.0*P1*pow(P2, 4)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+5.0*pow(P1, 4)*P2*Q2*(( 
    -3.0)*pow(Q1, 2)+pow(Q2, 2))+(-10.0)*pow(P1, 2)*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+ 
    pow(P2, 5)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 5)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))+( 
    -12.0)*(3.0*pow(P1, 5)*Q1*pow(Q2, 2)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*pow(P1, 4)*P2* 
    pow(Q2, 3)*((-5.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P2, 3)*Q2*((-1.0)*(5.0+3.0*pow(P2, 2))*pow(Q1, 4)+( 
    10.0+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*pow(Q2, 4))+3.0*P1*pow(P2, 2)*Q1*((1.0+ 
    pow(P2, 2))*pow(Q1, 4)+(-5.0)*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+(-1.0)* 
    pow(P1, 3)*Q1*((1.0+3.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*(1.0+(-3.0)* 
    pow(P2, 2))*pow(Q2, 4))+3.0*pow(P1, 2)*P2*Q2*(5.0*(1.0+pow(P2, 2))*pow(Q1, 4)+(-10.0)* 
    pow(Q1, 2)*pow(Q2, 2)+(-1.0)*((-1.0)+pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s9L - s9L0) * ( 
    (-7.0/128.0)*pow(G, -4.0)*(10.0*pow(G, 2)*(6.0*pow(P1, 5)*P2*Q1*(pow(Q1, 2)+(-3.0)* 
    pow(Q2, 2))+(-20.0)*pow(P1, 3)*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+6.0*P1*pow(P2, 5)* 
    Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+15.0*pow(P1, 4)*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+ 
    (-15.0)*pow(P1, 2)*pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P2, 6)*Q2*((-3.0)* 
    pow(Q1, 2)+pow(Q2, 2))+pow(P1, 6)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))+9.0*((-6.0)*pow(P1, 5)* 
    P2*Q1*(pow(Q1, 4)+30.0*pow(Q1, 2)*pow(Q2, 2)+(-35.0)*pow(Q2, 4))+15.0*pow(P1, 4)*Q2*(5.0*( 
    4.0+pow(P2, 2))*pow(Q1, 4)+10.0*((-4.0)+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(4.0+(-7.0)*pow(P2, 2)) 
    *pow(Q2, 4))+15.0*pow(P1, 2)*pow(P2, 2)*Q2*((-5.0)*(24.0+5.0*pow(P2, 2))*pow(Q1, 4)+10.0*( 
    24.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+3.0*((-8.0)+pow(P2, 2))*pow(Q2, 4))+20.0*pow(P1, 3)*P2* 
    Q1*((12.0+5.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(12.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-15.0)*( 
    (-4.0)+pow(P2, 2))*pow(Q2, 4))+(-6.0)*P1*pow(P2, 3)*Q1*((40.0+9.0*pow(P2, 2))*pow(Q1, 4)+( 
    -50.0)*(8.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(40.0+pow(P2, 2))*pow(Q2, 4))+pow(P2, 4)*Q2*( 
    15.0*(20.0+3.0*pow(P2, 2))*pow(Q1, 4)+(-50.0)*(12.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(60.0+ 
    pow(P2, 2))*pow(Q2, 4))+pow(P1, 6)*(15.0*pow(Q1, 4)*Q2+(-70.0)*pow(Q1, 2)*pow(Q2, 3)+11.0* 
    pow(Q2, 5)))) 
    ) + 
    (s10L - s10L0) * ( 
    (189.0/32.0)*pow(G, -4.0)*((-5.0)*pow(P1, 4)*P2*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+10.0*pow(P1, 2)*pow(P2, 3)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-1.0)*pow(P2, 5)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-10.0) 
    *pow(P1, 3)*pow(P2, 2)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+5.0*P1* 
    pow(P2, 4)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 5)*(pow(Q1, 5)+( 
    -10.0)*pow(Q1, 3)*pow(Q2, 2)+5.0*Q1*pow(Q2, 4)))
    ) + 
    (s11L - s11L0) * ( 
    (63.0/128.0)*pow(G, -4.0)*((-15.0)*pow(P1, 4)*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+(-10.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+15.0*pow(P1, 2)*pow(P2, 4)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+(-1.0)*pow(P2, 6)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)) 
    +6.0*pow(P1, 5)*P2*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+(-20.0)* 
    pow(P1, 3)*pow(P2, 3)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+6.0*P1* 
    pow(P2, 5)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 6)*(5.0*pow(Q1, 4)* 
    Q2+(-10.0)*pow(Q1, 2)*pow(Q2, 3)+pow(Q2, 5))) 
    );

    //P11
    I[1] = (Lf - L0) * ( 
    (15.0/4.0)*pow(G, -4.0)*(pow(G, 4)*(((-4.0)+3.0*pow(P1, 4)+(-41.0)*pow(P2, 2)+(-18.0)*pow(P2, 4)+ 
    pow(P1, 2)*(1.0+(-15.0)*pow(P2, 2)))*Q1+21.0*P1*P2*(2.0+pow(P1, 2)+pow(P2, 2))*Q2)+( 
    -28.0)*pow(G, 2)*(((-2.0)+pow(P1, 4)+(-21.0)*pow(P2, 2)+(-10.0)*pow(P2, 4)+pow(P1, 2)*(1.0+(-3.0) 
    *pow(P2, 2)))*pow(Q1, 3)+3.0*P1*P2*(8.0+pow(P1, 2)+7.0*pow(P2, 2))*pow(Q1, 2)*Q2+((-2.0)+ 
    3.0*pow(P1, 4)+(-19.0)*pow(P2, 2)+(-6.0)*pow(P2, 4)+(-1.0)*pow(P1, 2)*(1.0+21.0*pow(P2, 2)))* 
    Q1*pow(Q2, 2)+P1*P2*(20.0+13.0*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q2, 3))+21.0*(pow(Q1, 2)+ 
    pow(Q2, 2))*(((-8.0)+3.0*pow(P1, 4)+(-85.0)*pow(P2, 2)+(-42.0)*pow(P2, 4)+pow(P1, 2)*(5.0+(-3.0)* 
    pow(P2, 2)))*pow(Q1, 3)+(-3.0)*P1*P2*((-34.0)+pow(P1, 2)+(-35.0)*pow(P2, 2))*pow(Q1, 2)* 
    Q2+((-8.0)+15.0*pow(P1, 4)+(-73.0)*pow(P2, 2)+(-18.0)*pow(P2, 4)+(-1.0)*pow(P1, 2)*(7.0+111.0* 
    pow(P2, 2)))*Q1*pow(Q2, 2)+3.0*P1*P2*(26.0+19.0*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q2, 3)))  
    ) + 
    (cL - cL0) * ( 
    (-15.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(45.0*pow(P1, 5)*Q1+pow(P1, 3)*(368.0+(-82.0)* 
    pow(P2, 2))*Q1+P1*(16.0+(-1008.0)*pow(P2, 2)+(-127.0)*pow(P2, 4))*Q1+207.0*pow(P1, 4)* 
    P2*Q2+2.0*pow(P1, 2)*P2*(984.0+121.0*pow(P2, 2))*Q2+P2*(624.0+592.0*pow(P2, 2)+35.0* 
    pow(P2, 4))*Q2)+(-28.0)*pow(G, 2)*(2.0*pow(P1, 3)*Q1*((58.0+3.0*pow(P2, 2))*pow(Q1, 2)+( 
    194.0+(-91.0)*pow(P2, 2))*pow(Q2, 2))+2.0*pow(P1, 2)*P2*Q2*((150.0+71.0*pow(P2, 2))* 
    pow(Q1, 2)+3.0*(202.0+19.0*pow(P2, 2))*pow(Q2, 2))+P2*Q2*((336.0+572.0*pow(P2, 2)+43.0* 
    pow(P2, 4))*pow(Q1, 2)+(304.0+204.0*pow(P2, 2)+9.0*pow(P2, 4))*pow(Q2, 2))+(-1.0)*P1*Q1*((( 
    -16.0)+252.0*pow(P2, 2)+45.0*pow(P2, 4))*pow(Q1, 2)+(16.0+1260.0*pow(P2, 2)+119.0*pow(P2, 4))* 
    pow(Q2, 2))+pow(P1, 5)*(11.0*pow(Q1, 3)+57.0*Q1*pow(Q2, 2))+pow(P1, 4)*P2*((-21.0)*pow(Q1, 2)* 
    Q2+145.0*pow(Q2, 3)))+42.0*(pow(P1, 4)*P2*Q2*((-109.0)*pow(Q1, 4)+50.0*pow(Q1, 2)* 
    pow(Q2, 2)+343.0*pow(Q2, 4))+2.0*pow(P1, 3)*Q1*((84.0+19.0*pow(P2, 2))*pow(Q1, 4)+2.0*(276.0+( 
    -59.0)*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(100.0+(-61.0)*pow(P2, 2))*pow(Q2, 4))+2.0* 
    pow(P1, 2)*P2*Q2*(3.0*((-4.0)+35.0*pow(P2, 2))*pow(Q1, 4)+2.0*(612.0+179.0*pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(1332.0+101.0*pow(P2, 2))*pow(Q2, 4))+P2*Q2*((696.0+1400.0*pow(P2, 2)+ 
    119.0*pow(P2, 4))*pow(Q1, 4)+2.0*(648.0+888.0*pow(P2, 2)+53.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+( 
    600.0+312.0*pow(P2, 2)+11.0*pow(P2, 4))*pow(Q2, 4))+(-1.0)*P1*Q1*(((-40.0)+264.0* 
    pow(P2, 2)+63.0*pow(P2, 4))*pow(Q1, 4)+2.0*(8.0+1704.0*pow(P2, 2)+225.0*pow(P2, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+(56.0+3336.0*pow(P2, 2)+251.0*pow(P2, 4))*pow(Q2, 4))+pow(P1, 5)*(13.0*pow(Q1, 5)+134.0* 
    pow(Q1, 3)*pow(Q2, 2)+161.0*Q1*pow(Q2, 4))))
    ) + 
    (c2L - c2L0) * ( 
    (-15.0/32.0)*pow(G, -4.0)*(pow(G, 4)*((-36.0)*pow(P1, 3)*P2*Q1+(-4.0)*P1*P2*(40.0+ 
    31.0*pow(P2, 2))*Q1+53.0*pow(P1, 4)*Q2+6.0*pow(P1, 2)*(28.0+31.0*pow(P2, 2))*Q2+(32.0+ 
    152.0*pow(P2, 2)+45.0*pow(P2, 4))*Q2)+56.0*pow(G, 2)*((-3.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+ 
    (-15.0)*pow(Q2, 2))+(-1.0)*pow(P1, 2)*Q2*((9.0+51.0*pow(P2, 2))*pow(Q1, 2)+(53.0+45.0* 
    pow(P2, 2))*pow(Q2, 2))+P1*P2*Q1*(3.0*(6.0+7.0*pow(P2, 2))*pow(Q1, 2)+(106.0+61.0*pow(P2, 2)) 
    *pow(Q2, 2))+(-1.0)*Q2*((8.0+71.0*pow(P2, 2)+27.0*pow(P2, 4))*pow(Q1, 2)+(8.0+27.0*pow(P2, 2)+ 
    6.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 4)*(4.0*pow(Q1, 2)*Q2+(-19.0)*pow(Q2, 3)))+(-42.0)*(( 
    -24.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 4)+(-4.0)*pow(Q1, 2)*pow(Q2, 2)+(-13.0)*pow(Q2, 4))+8.0* 
    P1*P2*Q1*((4.0+7.0*pow(P2, 2))*pow(Q1, 4)+4.0*(17.0+14.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+ 
    3.0*(24.0+11.0*pow(P2, 2))*pow(Q2, 4))+(-2.0)*pow(P1, 2)*Q2*(((-10.0)+69.0*pow(P2, 2))* 
    pow(Q1, 4)+2.0*(46.0+135.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(118.0+81.0*pow(P2, 2))*pow(Q2, 4))+( 
    -1.0)*Q2*((32.0+340.0*pow(P2, 2)+147.0*pow(P2, 4))*pow(Q1, 4)+2.0*(32.0+228.0*pow(P2, 2)+69.0* 
    pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(32.0+84.0*pow(P2, 2)+15.0*pow(P2, 4))*pow(Q2, 4))+pow(P1, 4)*(33.0* 
    pow(Q1, 4)*Q2+(-2.0)*pow(Q1, 2)*pow(Q2, 3)+(-91.0)*pow(Q2, 5))))
    ) + 
    (c3L - c3L0) * ( 
    (5.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(47.0*pow(P1, 5)*Q1+2.0*pow(P1, 3)*(308.0+73.0*pow(P2, 2)) 
    *Q1+3.0*P1*(208.0+632.0*pow(P2, 2)+81.0*pow(P2, 4))*Q1+(-51.0)*pow(P1, 4)*P2*Q2+( 
    -6.0)*pow(P1, 2)*P2*(180.0+43.0*pow(P2, 2))*Q2+(-3.0)*P2*(208.0+296.0*pow(P2, 2)+21.0* 
    pow(P2, 4))*Q2)+(-56.0)*pow(G, 2)*((-1.0)*pow(P1, 4)*P2*Q2*(39.0*pow(Q1, 2)+4.0* 
    pow(Q2, 2))+2.0*pow(P1, 3)*Q1*((11.0+6.0*pow(P2, 2))*pow(Q1, 2)+55.0*(5.0+pow(P2, 2))*pow(Q2, 2))+( 
    -6.0)*pow(P1, 2)*P2*Q2*((105.0+22.0*pow(P2, 2))*pow(Q1, 2)+(25.0+7.0*pow(P2, 2))*pow(Q2, 2))+ 
    P1*Q1*((68.0+342.0*pow(P2, 2)+51.0*pow(P2, 4))*pow(Q1, 2)+30.0*(14.0+29.0*pow(P2, 2)+3.0* 
    pow(P2, 4))*pow(Q2, 2))+(-1.0)*P2*Q2*((396.0+582.0*pow(P2, 2)+45.0*pow(P2, 4))*pow(Q1, 2)+ 
    2.0*(38.0+51.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 5)*(pow(Q1, 3)+44.0*Q1*pow(Q2, 2)))+( 
    -42.0)*(3.0*pow(P1, 4)*P2*Q2*(17.0*pow(Q1, 4)+174.0*pow(Q1, 2)*pow(Q2, 2)+(-11.0)* 
    pow(Q2, 4))+2.0*pow(P1, 2)*P2*Q2*((960.0+269.0*pow(P2, 2))*pow(Q1, 4)+2.0*(1560.0+259.0* 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(48.0+49.0*pow(P2, 2))*pow(Q2, 4))+2.0*pow(P1, 3)*Q1*((32.0+7.0* 
    pow(P2, 2))*pow(Q1, 4)+(-2.0)*(424.0+179.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-3.0)*(592.0+87.0* 
    pow(P2, 2))*pow(Q2, 4))+P2*Q2*((1944.0+3248.0*pow(P2, 2)+271.0*pow(P2, 4))*pow(Q1, 4)+2.0*( 
    1224.0+1408.0*pow(P2, 2)+89.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(120.0+208.0*pow(P2, 2)+11.0* 
    pow(P2, 4))*pow(Q2, 4))+(-3.0)*P1*Q1*((40.0+304.0*pow(P2, 2)+53.0*pow(P2, 4))*pow(Q1, 4)+2.0* 
    (344.0+1216.0*pow(P2, 2)+143.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(776.0+1104.0*pow(P2, 2)+97.0* 
    pow(P2, 4))*pow(Q2, 4))+pow(P1, 5)*(5.0*pow(Q1, 5)+(-98.0)*pow(Q1, 3)*pow(Q2, 2)+(-303.0)*Q1* 
    pow(Q2, 4))))
    ) + 
    (c4L - c4L0) * ( 
    (3.0/16.0)*pow(G, -4.0)*(5.0*pow(G, 4)*(9.0*pow(P1, 3)*P2*Q1+P1*P2*(38.0+29.0* 
    pow(P2, 2))*Q1+10.0*pow(P1, 4)*Q2+pow(P1, 2)*(19.0+(-3.0)*pow(P2, 2))*Q2+(-1.0)*pow(P2, 2)* 
    (19.0+9.0*pow(P2, 2))*Q2)+140.0*pow(G, 2)*((-2.0)*pow(P1, 3)*P2*Q1*(2.0*pow(Q1, 2)+3.0* 
    pow(Q2, 2))+(-2.0)*P1*P2*Q1*((9.0+7.0*pow(P2, 2))*pow(Q1, 2)+(11.0+8.0*pow(P2, 2))* 
    pow(Q2, 2))+Q2*((6.0+38.0*pow(P2, 2)+15.0*pow(P2, 4))*pow(Q1, 2)+((-2.0)+pow(P2, 4))*pow(Q2, 2))+ 
    pow(P1, 4)*(7.0*pow(Q1, 2)*Q2+(-9.0)*pow(Q2, 3))+pow(P1, 2)*(2.0*(11.0+12.0*pow(P2, 2))* 
    pow(Q1, 2)*Q2+(-2.0)*(10.0+3.0*pow(P2, 2))*pow(Q2, 3)))+(-42.0)*((-15.0)*pow(P1, 3)*P2* 
    Q1*(pow(Q1, 4)+22.0*pow(Q1, 2)*pow(Q2, 2)+(-3.0)*pow(Q2, 4))+5.0*pow(P1, 2)*Q2*((43.0+63.0* 
    pow(P2, 2))*pow(Q1, 4)+6.0*(15.0+11.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-3.0)*(19.0+7.0*pow(P2, 2)) 
    *pow(Q2, 4))+(-5.0)*P1*P2*Q1*((26.0+23.0*pow(P2, 2))*pow(Q1, 4)+2.0*(86.0+53.0* 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(2.0+11.0*pow(P2, 2))*pow(Q2, 4))+Q2*(5.0*(16.0+117.0*pow(P2, 2)+ 
    48.0*pow(P2, 4))*pow(Q1, 4)+10.0*(8.0+35.0*pow(P2, 2)+12.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)* 
    (32.0+35.0*pow(P2, 2))*pow(Q2, 4))+5.0*pow(P1, 4)*(11.0*pow(Q1, 4)*Q2+34.0*pow(Q1, 2)*pow(Q2, 3)+( 
    -25.0)*pow(Q2, 5))))
    ) + 
    (c5L - c5L0) * ( 
    (-3.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(29.0*pow(P1, 5)*Q1+2.0*pow(P1, 3)*(148.0+3.0*pow(P2, 2)) 
    *Q1+(-1.0)*P1*pow(P2, 2)*(888.0+151.0*pow(P2, 2))*Q1+(-121.0)*pow(P1, 4)*P2*Q2+( 
    -6.0)*pow(P1, 2)*P2*(148.0+9.0*pow(P2, 2))*Q2+pow(P2, 3)*(296.0+35.0*pow(P2, 2))*Q2)+ 
    56.0*pow(G, 2)*(P2*Q2*((-1.0)*(228.0+398.0*pow(P2, 2)+35.0*pow(P2, 4))*pow(Q1, 2)+2.0*( 
    38.0+17.0*pow(P2, 2))*pow(Q2, 2))+2.0*pow(P1, 2)*P2*Q2*((-3.0)*(29.0+8.0*pow(P2, 2))* 
    pow(Q1, 2)+(177.0+17.0*pow(P2, 2))*pow(Q2, 2))+2.0*pow(P1, 3)*Q1*((23.0+8.0*pow(P2, 2))*pow(Q1, 2)+ 
    (-1.0)*(217.0+27.0*pow(P2, 2))*pow(Q2, 2))+P1*Q1*((76.0+318.0*pow(P2, 2)+45.0*pow(P2, 4))* 
    pow(Q1, 2)+2.0*((-114.0)+(-33.0)*pow(P2, 2)+8.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 5)*(3.0*pow(Q1, 3)+( 
    -38.0)*Q1*pow(Q2, 2))+pow(P1, 4)*P2*((-5.0)*pow(Q1, 2)*Q2+42.0*pow(Q2, 3)))+(-21.0)*( 
    pow(P1, 4)*P2*Q2*((-301.0)*pow(Q1, 4)+442.0*pow(Q1, 2)*pow(Q2, 2)+359.0*pow(Q2, 4))+6.0* 
    pow(P1, 2)*P2*Q2*((-1.0)*(736.0+145.0*pow(P2, 2))*pow(Q1, 4)+34.0*(16.0+pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(512.0+51.0*pow(P2, 2))*pow(Q2, 4))+2.0*pow(P1, 3)*Q1*(5.0*(32.0+11.0* 
    pow(P2, 2))*pow(Q1, 4)+2.0*(304.0+109.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(3776.0+541.0* 
    pow(P2, 2))*pow(Q2, 4))+P2*Q2*((-1.0)*(3696.0+5920.0*pow(P2, 2)+505.0*pow(P2, 4))* 
    pow(Q1, 4)+(-2.0)*((-48.0)+448.0*pow(P2, 2)+55.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(720.0+416.0* 
    pow(P2, 2)+11.0*pow(P2, 4))*pow(Q2, 4))+P1*Q1*(5.0*(112.0+480.0*pow(P2, 2)+69.0*pow(P2, 4))* 
    pow(Q1, 4)+2.0*(848.0+3264.0*pow(P2, 2)+435.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(4496.0+ 
    4320.0*pow(P2, 2)+179.0*pow(P2, 4))*pow(Q2, 4))+pow(P1, 5)*(21.0*pow(Q1, 5)+78.0*pow(Q1, 3)* 
    pow(Q2, 2)+(-647.0)*Q1*pow(Q2, 4)))) 
    ) + 
    (c6L - c6L0) * ( 
    (1.0/32.0)*pow(G, -4.0)*((-45.0)*pow(G, 4)*(4.0*pow(P1, 3)*P2*Q1+(-4.0)*P1*pow(P2, 3)* 
    Q1+pow(P1, 4)*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q2+pow(P2, 4)*Q2)+(-280.0)*pow(G, 2)*( 
    pow(P1, 3)*P2*Q1*(pow(Q1, 2)+(-39.0)*pow(Q2, 2))+pow(P2, 2)*Q2*((-3.0)*(9.0+5.0*pow(P2, 2)) 
    *pow(Q1, 2)+(9.0+2.0*pow(P2, 2))*pow(Q2, 2))+3.0*pow(P1, 2)*Q2*(3.0*(3.0+pow(P2, 2))*pow(Q1, 2)+(( 
    -3.0)+5.0*pow(P2, 2))*pow(Q2, 2))+P1*P2*Q1*((18.0+17.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*( 
    18.0+5.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(12.0*pow(Q1, 2)*Q2+(-7.0)*pow(Q2, 3)))+21.0*(60.0* 
    pow(P1, 3)*P2*Q1*(3.0*pow(Q1, 4)+(-22.0)*pow(Q1, 2)*pow(Q2, 2)+(-41.0)*pow(Q2, 4))+15.0* 
    pow(P1, 4)*Q2*(pow(Q1, 4)+126.0*pow(Q1, 2)*pow(Q2, 2)+(-35.0)*pow(Q2, 4))+60.0*P1*P2*Q1* 
    ((16.0+13.0*pow(P2, 2))*pow(Q1, 4)+2.0*((-8.0)+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(64.0+ 
    23.0*pow(P2, 2))*pow(Q2, 4))+10.0*pow(P1, 2)*Q2*((-1.0)*(28.0+93.0*pow(P2, 2))*pow(Q1, 4)+2.0* 
    (244.0+165.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-92.0)+39.0*pow(P2, 2))*pow(Q2, 4))+Q2*((-5.0) 
    *(64.0+584.0*pow(P2, 2)+261.0*pow(P2, 4))*pow(Q1, 4)+10.0*(64.0+152.0*pow(P2, 2)+21.0*pow(P2, 4)) 
    *pow(Q1, 2)*pow(Q2, 2)+((-64.0)+280.0*pow(P2, 2)+75.0*pow(P2, 4))*pow(Q2, 4)))) 
    ) + 
    (c7L - c7L0) * ( 
    (15.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 5)*Q1+(-10.0)*pow(P1, 3)*pow(P2, 2)*Q1+5.0* 
    P1*pow(P2, 4)*Q1+(-5.0)*pow(P1, 4)*P2*Q2+10.0*pow(P1, 2)*pow(P2, 3)*Q2+(-1.0)* 
    pow(P2, 5)*Q2)+2.0*pow(G, 2)*(pow(P2, 3)*Q2*((408.0+55.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*( 
    136.0+9.0*pow(P2, 2))*pow(Q2, 2))+(-2.0)*pow(P1, 2)*P2*Q2*((612.0+71.0*pow(P2, 2))* 
    pow(Q1, 2)+((-204.0)+23.0*pow(P2, 2))*pow(Q2, 2))+2.0*pow(P1, 3)*Q1*((68.0+13.0*pow(P2, 2))* 
    pow(Q1, 2)+((-204.0)+101.0*pow(P2, 2))*pow(Q2, 2))+P1*pow(P2, 2)*Q1*((-3.0)*(136.0+27.0* 
    pow(P2, 2))*pow(Q1, 2)+(1224.0+103.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 5)*(11.0*pow(Q1, 3)+(-61.0)* 
    Q1*pow(Q2, 2))+pow(P1, 4)*P2*((-133.0)*pow(Q1, 2)*Q2+91.0*pow(Q2, 3)))+(-3.0)*( 
    pow(P1, 4)*P2*Q2*((-179.0)*pow(Q1, 4)+(-706.0)*pow(Q1, 2)*pow(Q2, 2)+289.0*pow(Q2, 4))+2.0* 
    pow(P1, 2)*P2*Q2*(((-516.0)+7.0*pow(P2, 2))*pow(Q1, 4)+(-6.0)*(644.0+97.0*pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+3.0*(292.0+pow(P2, 2))*pow(Q2, 4))+P2*Q2*((1200.0+2744.0*pow(P2, 2)+ 
    273.0*pow(P2, 4))*pow(Q1, 4)+(-2.0)*(1200.0+1112.0*pow(P2, 2)+53.0*pow(P2, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+(240.0+(-104.0)*pow(P2, 2)+(-11.0)*pow(P2, 4))*pow(Q2, 4))+P1*Q1*((-1.0)*( 
    240.0+1416.0*pow(P2, 2)+217.0*pow(P2, 4))*pow(Q1, 4)+2.0*(1200.0+2184.0*pow(P2, 2)+113.0* 
    pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+((-1200.0)+2712.0*pow(P2, 2)+299.0*pow(P2, 4))*pow(Q2, 4))+ 
    pow(P1, 5)*(3.0*pow(Q1, 5)+234.0*pow(Q1, 3)*pow(Q2, 2)+(-361.0)*Q1*pow(Q2, 4))+pow(P1, 3)*((-2.0) 
    *(4.0+19.0*pow(P2, 2))*pow(Q1, 5)+4.0*(836.0+251.0*pow(P2, 2))*pow(Q1, 3)*pow(Q2, 2)+2.0*(( 
    -1652.0)+153.0*pow(P2, 2))*Q1*pow(Q2, 4))))
    ) + 
    (c8L - c8L0) * ( 
    (105.0/32.0)*pow(G, -4.0)*(2.0*pow(G, 2)*(4.0*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+(-3.0)* 
    pow(Q2, 2))+(-4.0)*P1*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)* 
    Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+ 
    pow(P1, 4)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))+3.0*((-3.0)*pow(P1, 3)*P2*Q1*( 
    pow(Q1, 4)+22.0*pow(Q1, 2)*pow(Q2, 2)+(-27.0)*pow(Q2, 4))+pow(P2, 2)*Q2*((-1.0)*(35.0+24.0* 
    pow(P2, 2))*pow(Q1, 4)+2.0*(35.0+12.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-7.0)*pow(Q2, 4))+P1* 
    P2*Q1*((14.0+17.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*(70.0+37.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+ 
    (70.0+(-11.0)*pow(P2, 2))*pow(Q2, 4))+pow(P1, 2)*Q2*((35.0+39.0*pow(P2, 2))*pow(Q1, 4)+2.0*(( 
    -35.0)+33.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+7.0*(1.0+(-3.0)*pow(P2, 2))*pow(Q2, 4))+pow(P1, 4)*( 
    11.0*pow(Q1, 4)*Q2+(-46.0)*pow(Q1, 2)*pow(Q2, 3)+7.0*pow(Q2, 5))))
    ) + 
    (c9L - c9L0) * ( 
    (-7.0/128.0)*pow(G, -4.0)*((-15.0)*pow(P1, 4)*P2*Q2*(35.0*pow(Q1, 4)+(-190.0)* 
    pow(Q1, 2)*pow(Q2, 2)+31.0*pow(Q2, 4))+10.0*pow(P1, 3)*Q1*((52.0+23.0*pow(P2, 2))*pow(Q1, 4)+ 
    130.0*((-4.0)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(52.0+(-49.0)*pow(P2, 2))*pow(Q2, 4))+(-10.0) 
    *pow(P1, 2)*P2*Q2*(5.0*(156.0+31.0*pow(P2, 2))*pow(Q1, 4)+10.0*((-156.0)+5.0*pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+(156.0+(-41.0)*pow(P2, 2))*pow(Q2, 4))+(-15.0)*P1*pow(P2, 2)*Q1*(( 
    104.0+25.0*pow(P2, 2))*pow(Q1, 4)+(-130.0)*(8.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(104.0+ 
    pow(P2, 2))*pow(Q2, 4))+pow(P2, 3)*Q2*(5.0*(520.0+83.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(520.0+ 
    47.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(520.0+11.0*pow(P2, 2))*pow(Q2, 4))+pow(P1, 5)*(29.0* 
    pow(Q1, 5)+(-650.0)*pow(Q1, 3)*pow(Q2, 2)+505.0*Q1*pow(Q2, 4))+10.0*pow(G, 2)*((-10.0)* 
    pow(P1, 3)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+5.0*P1*pow(P2, 4)*Q1*(pow(Q1, 2)+( 
    -3.0)*pow(Q2, 2))+5.0*pow(P1, 4)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-10.0)*pow(P1, 2)* 
    pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P2, 5)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+ 
    pow(P1, 5)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2))))
    ) + 
    (c10L - c10L0) * ( 
    (-189.0/32.0)*pow(G, -4.0)*((-6.0)*pow(P1, 2)*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+(-10.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+pow(P2, 4)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+4.0*pow(P1, 3)*P2*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+(-4.0) 
    *P1*pow(P2, 3)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 4)*(5.0* 
    pow(Q1, 4)*Q2+(-10.0)*pow(Q1, 2)*pow(Q2, 3)+pow(Q2, 5)))
    ) + 
    (c11L - c11L0) * ( 
    (63.0/128.0)*pow(G, -4.0)*((-5.0)*pow(P1, 4)*P2*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+10.0*pow(P1, 2)*pow(P2, 3)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-1.0)*pow(P2, 5)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-10.0) 
    *pow(P1, 3)*pow(P2, 2)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+5.0*P1* 
    pow(P2, 4)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 5)*(pow(Q1, 5)+( 
    -10.0)*pow(Q1, 3)*pow(Q2, 2)+5.0*Q1*pow(Q2, 4)))
    ) + 
    (sL - sL0) * ( 
    (3.0/128.0)*pow(G, -4.0)*((-5.0)*pow(G, 4)*(3.0*(432.0+304.0*pow(P1, 2)+11.0*pow(P1, 4))* 
    P2*Q1+2.0*(1144.0+119.0*pow(P1, 2))*pow(P2, 3)*Q1+205.0*pow(P2, 5)*Q1+(-43.0)*P1*( 
    16.0+16.0*pow(P1, 2)+pow(P1, 4))*Q2+(-258.0)*P1*(8.0+pow(P1, 2))*pow(P2, 2)*Q2+(-215.0)* 
    P1*pow(P2, 4)*Q2)+140.0*pow(G, 2)*((-5.0)*P1*pow(P2, 4)*Q2*(53.0*pow(Q1, 2)+11.0* 
    pow(Q2, 2))+(-6.0)*P1*pow(P2, 2)*Q2*((334.0+23.0*pow(P1, 2))*pow(Q1, 2)+(118.0+21.0* 
    pow(P1, 2))*pow(Q2, 2))+2.0*pow(P2, 3)*Q1*(35.0*(18.0+pow(P1, 2))*pow(Q1, 2)+(398.0+133.0* 
    pow(P1, 2))*pow(Q2, 2))+(-1.0)*P1*Q2*((368.0+68.0*pow(P1, 2)+(-7.0)*pow(P1, 4))*pow(Q1, 2)+( 
    336.0+436.0*pow(P1, 2)+31.0*pow(P1, 4))*pow(Q2, 2))+P2*Q1*((656.0+156.0*pow(P1, 2)+(-9.0)* 
    pow(P1, 4))*pow(Q1, 2)+3.0*(208.0+452.0*pow(P1, 2)+31.0*pow(P1, 4))*pow(Q2, 2))+pow(P2, 5)*(119.0* 
    pow(Q1, 3)+53.0*Q1*pow(Q2, 2)))+(-42.0)*((-5.0)*P1*pow(P2, 4)*Q2*(735.0*pow(Q1, 4)+ 
    650.0*pow(Q1, 2)*pow(Q2, 2)+67.0*pow(Q2, 4))+(-10.0)*P1*pow(P2, 2)*Q2*(5.0*(492.0+17.0* 
    pow(P1, 2))*pow(Q1, 4)+2.0*(1548.0+191.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(540.0+113.0*pow(P1, 2)) 
    *pow(Q2, 4))+10.0*pow(P2, 3)*Q1*(7.0*(188.0+5.0*pow(P1, 2))*pow(Q1, 4)+490.0*(4.0+pow(P1, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(612.0+287.0*pow(P1, 2))*pow(Q2, 4))+5.0*P2*Q1*((1320.0+24.0*pow(P1, 2)+ 
    (-31.0)*pow(P1, 4))*pow(Q1, 4)+2.0*(1272.0+1752.0*pow(P1, 2)+47.0*pow(P1, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+(1224.0+3672.0*pow(P1, 2)+325.0*pow(P1, 4))*pow(Q2, 4))+(-1.0)*P1*Q2*((-5.0)*( 
    (-760.0)+120.0*pow(P1, 2)+29.0*pow(P1, 4))*pow(Q1, 4)+10.0*(712.0+392.0*pow(P1, 2)+pow(P1, 4))* 
    pow(Q1, 2)*pow(Q2, 2)+(3320.0+4840.0*pow(P1, 2)+371.0*pow(P1, 4))*pow(Q2, 4))+pow(P2, 5)*(1281.0* 
    pow(Q1, 5)+1470.0*pow(Q1, 3)*pow(Q2, 2)+325.0*Q1*pow(Q2, 4))))
    ) + 
    (s2L - s2L0) * ( 
    (-15.0/32.0)*pow(G, -4.0)*(pow(G, 4)*((32.0+31.0*pow(P1, 4)+240.0*pow(P2, 2)+111.0*pow(P2, 4)+ 
    pow(P1, 2)*(80.0+54.0*pow(P2, 2)))*Q1+4.0*P1*P2*((-4.0)+9.0*pow(P1, 2)+(-13.0)*pow(P2, 2)) 
    *Q2)+(-28.0)*pow(G, 2)*((16.0+5.0*pow(P1, 4)+142.0*pow(P2, 2)+67.0*pow(P2, 4)+6.0*pow(P1, 2)*( 
    3.0+4.0*pow(P2, 2)))*pow(Q1, 3)+(-14.0)*P1*P2*(10.0+3.0*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q1, 2)* 
    Q2+(16.0+47.0*pow(P1, 4)+54.0*pow(P2, 2)+21.0*pow(P2, 4)+2.0*pow(P1, 2)*(53.0+18.0*pow(P2, 2)))* 
    Q1*pow(Q2, 2)+2.0*P1*P2*(18.0+19.0*pow(P1, 2)+(-1.0)*pow(P2, 2))*pow(Q2, 3))+42.0*((32.0+ 
    3.0*pow(P1, 4)+304.0*pow(P2, 2)+147.0*pow(P2, 4)+2.0*pow(P1, 2)*(8.0+15.0*pow(P2, 2)))*pow(Q1, 5)+( 
    -28.0)*P1*P2*(14.0+3.0*pow(P1, 2)+11.0*pow(P2, 2))*pow(Q1, 4)*Q2+2.0*(32.0+45.0*pow(P1, 4)+ 
    184.0*pow(P2, 2)+69.0*pow(P2, 4)+2.0*pow(P1, 2)*(68.0+69.0*pow(P2, 2)))*pow(Q1, 3)*pow(Q2, 2)+( 
    -168.0)*P1*P2*(2.0+pow(P1, 2)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 3)+(32.0+143.0*pow(P1, 4)+32.0* 
    pow(P2, 2)+15.0*pow(P2, 4)+6.0*pow(P1, 2)*(48.0+pow(P2, 2)))*Q1*pow(Q2, 4)+12.0*P1*P2*(10.0+ 
    9.0*pow(P1, 2)+pow(P2, 2))*pow(Q2, 5)))
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*(208.0+(-8.0)*pow(P1, 2)+pow(P1, 4))*P2*Q1+(-2.0) 
    *((-628.0)+7.0*pow(P1, 2))*pow(P2, 3)*Q1+127.0*pow(P2, 5)*Q1+3.0*P1*(208.0+328.0* 
    pow(P1, 2)+25.0*pow(P1, 4))*Q2+18.0*P1*(44.0+13.0*pow(P1, 2))*pow(P2, 2)*Q2+15.0*P1* 
    pow(P2, 4)*Q2)+(-56.0)*pow(G, 2)*(pow(P2, 5)*Q1*(42.0*pow(Q1, 2)+pow(Q2, 2))+6.0*P1* 
    pow(P2, 2)*Q2*((-33.0)*pow(Q1, 2)+(55.0+13.0*pow(P1, 2))*pow(Q2, 2))+2.0*pow(P2, 3)*Q1*(7.0* 
    (31.0+pow(P1, 2))*pow(Q1, 2)+(-1.0)*(23.0+28.0*pow(P1, 2))*pow(Q2, 2))+P2*Q1*(2.0*(118.0+ 
    57.0*pow(P1, 2)+6.0*pow(P1, 4))*pow(Q1, 2)+(-3.0)*(28.0+122.0*pow(P1, 2)+11.0*pow(P1, 4))*pow(Q2, 2)) 
    +P1*Q2*((-3.0)*(36.0+50.0*pow(P1, 2)+5.0*pow(P1, 4))*pow(Q1, 2)+2.0*(122.0+189.0*pow(P1, 2)+ 
    15.0*pow(P1, 4))*pow(Q2, 2))+P1*pow(P2, 4)*((-33.0)*pow(Q1, 2)*Q2+16.0*pow(Q2, 3)))+(-42.0) 
    *((-3.0)*P1*pow(P2, 4)*Q2*((-111.0)*pow(Q1, 4)+46.0*pow(Q1, 2)*pow(Q2, 2)+21.0*pow(Q2, 4)) 
    +(-2.0)*P1*pow(P2, 2)*Q2*((-1.0)*(1272.0+91.0*pow(P1, 2))*pow(Q1, 4)+2.0*(480.0+91.0* 
    pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(696.0+169.0*pow(P1, 2))*pow(Q2, 4))+(-2.0)*pow(P2, 3)*Q1*( 
    3.0*(336.0+11.0*pow(P1, 2))*pow(Q1, 4)+6.0*(56.0+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(352.0+ 
    227.0*pow(P1, 2))*pow(Q2, 4))+(-3.0)*P2*Q1*((360.0+144.0*pow(P1, 2)+13.0*pow(P1, 4))* 
    pow(Q1, 4)+2.0*(88.0+192.0*pow(P1, 2)+31.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(312.0+1168.0* 
    pow(P1, 2)+119.0*pow(P1, 4))*pow(Q2, 4))+P1*Q2*((840.0+832.0*pow(P1, 2)+65.0*pow(P1, 4))* 
    pow(Q1, 4)+2.0*(24.0+368.0*pow(P1, 2)+55.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(1176.0+ 
    1888.0*pow(P1, 2)+155.0*pow(P1, 4))*pow(Q2, 4))+pow(P2, 5)*((-195.0)*pow(Q1, 5)+(-66.0)* 
    pow(Q1, 3)*pow(Q2, 2)+25.0*Q1*pow(Q2, 4))))
    ) + 
    (s4L - s4L0) * ( 
    (3.0/16.0)*pow(G, -4.0)*(5.0*pow(G, 4)*(7.0*pow(P1, 4)*Q1+(-1.0)*pow(P2, 2)*(19.0+12.0* 
    pow(P2, 2))*Q1+pow(P1, 2)*(19.0+15.0*pow(P2, 2))*Q1+(-21.0)*pow(P1, 3)*P2*Q2+(-1.0)* 
    P1*P2*(38.0+17.0*pow(P2, 2))*Q2)+70.0*pow(G, 2)*((4.0+2.0*pow(P1, 2)+pow(P1, 4)+38.0* 
    pow(P2, 2)+19.0*pow(P2, 4))*pow(Q1, 3)+2.0*P1*P2*(16.0+9.0*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q1, 2)* 
    Q2+(-1.0)*(12.0+31.0*pow(P1, 4)+38.0*pow(P2, 2)+9.0*pow(P2, 4)+pow(P1, 2)*(82.0+60.0*pow(P2, 2)))* 
    Q1*pow(Q2, 2)+2.0*P1*P2*(20.0+11.0*pow(P1, 2)+9.0*pow(P2, 2))*pow(Q2, 3))+(-84.0)*((14.0+ 
    5.0*pow(P1, 4)+125.0*pow(P2, 2)+60.0*pow(P2, 4)+15.0*pow(P1, 2)*(1.0+pow(P2, 2)))*pow(Q1, 5)+(-5.0)* 
    P1*P2*(8.0+3.0*pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 4)*Q2+(-10.0)*(2.0+2.0*pow(P1, 4)+11.0* 
    pow(P2, 2)+3.0*pow(P2, 4)+3.0*pow(P1, 2)*(3.0+5.0*pow(P2, 2)))*pow(Q1, 3)*pow(Q2, 2)+10.0*P1*P2*( 
    40.0+21.0*pow(P1, 2)+19.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 3)+(-5.0)*(10.0+29.0*pow(P1, 4)+27.0* 
    pow(P2, 2)+6.0*pow(P2, 4)+pow(P1, 2)*(73.0+45.0*pow(P2, 2)))*Q1*pow(Q2, 4)+5.0*P1*P2*(16.0+ 
    9.0*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q2, 5))) 
    ) + 
    (s5L - s5L0) * ( 
    (3.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(71.0*pow(P1, 4)*P2*Q1+(-1.0)*pow(P2, 3)*(296.0+45.0* 
    pow(P2, 2))*Q1+2.0*pow(P1, 2)*P2*(444.0+77.0*pow(P2, 2))*Q1+39.0*pow(P1, 5)*Q2+pow(P1, 3)* 
    (296.0+(-94.0)*pow(P2, 2))*Q2+(-1.0)*P1*pow(P2, 2)*(888.0+101.0*pow(P2, 2))*Q2)+56.0* 
    pow(G, 2)*(P1*pow(P2, 4)*Q2*(47.0*pow(Q1, 2)+18.0*pow(Q2, 2))+2.0*P1*pow(P2, 2)*Q2*(( 
    255.0+38.0*pow(P1, 2))*pow(Q1, 2)+3.0*(21.0+pow(P1, 2))*pow(Q2, 2))+(-2.0)*pow(P2, 3)*Q1*((( 
    -91.0)+9.0*pow(P1, 2))*pow(Q1, 2)+25.0*(5.0+2.0*pow(P1, 2))*pow(Q2, 2))+P1*Q2*((228.0+286.0* 
    pow(P1, 2)+21.0*pow(P1, 4))*pow(Q1, 2)+(-2.0)*(38.0+97.0*pow(P1, 2)+10.0*pow(P1, 4))*pow(Q2, 2))+( 
    -1.0)*P2*Q1*(((-76.0)+90.0*pow(P1, 2)+6.0*pow(P1, 4))*pow(Q1, 2)+(228.0+618.0*pow(P1, 2)+ 
    53.0*pow(P1, 4))*pow(Q2, 2))+5.0*pow(P2, 5)*(4.0*pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))+(-21.0)*( 
    P1*pow(P2, 4)*Q2*(243.0*pow(Q1, 4)+1018.0*pow(Q1, 2)*pow(Q2, 2)+71.0*pow(Q2, 4))+2.0*P1* 
    pow(P2, 2)*Q2*(27.0*(48.0+7.0*pow(P1, 2))*pow(Q1, 4)+2.0*(2784.0+419.0*pow(P1, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+(48.0+(-55.0)*pow(P1, 2))*pow(Q2, 4))+(-2.0)*pow(P2, 3)*Q1*(((-1088.0)+7.0* 
    pow(P1, 2))*pow(Q1, 4)+2.0*(1072.0+397.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(928.0+403.0*pow(P1, 2)) 
    *pow(Q2, 4))+P1*Q2*((976.0+1088.0*pow(P1, 2)+71.0*pow(P1, 4))*pow(Q1, 4)+2.0*(2672.0+ 
    3488.0*pow(P1, 2)+265.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(1264.0+2560.0*pow(P1, 2)+ 
    245.0*pow(P1, 4))*pow(Q2, 4))+P2*Q1*((1104.0+96.0*pow(P1, 2)+23.0*pow(P1, 4))*pow(Q1, 4)+( 
    -2.0)*(1872.0+4800.0*pow(P1, 2)+403.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(1776.0+ 
    5088.0*pow(P1, 2)+445.0*pow(P1, 4))*pow(Q2, 4))+3.0*pow(P2, 5)*(73.0*pow(Q1, 5)+(-90.0)* 
    pow(Q1, 3)*pow(Q2, 2)+(-35.0)*Q1*pow(Q2, 4))))  
    ) + 
    (s6L - s6L0) * ( 
    (1.0/32.0)*pow(G, -4.0)*((-45.0)*pow(G, 4)*(pow(P1, 4)*Q1+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1+ 
    pow(P2, 4)*Q1+(-4.0)*pow(P1, 3)*P2*Q2+4.0*P1*pow(P2, 3)*Q2)+21.0*(((-64.0)+45.0* 
    pow(P1, 4)+(-800.0)*pow(P2, 2)+(-435.0)*pow(P2, 4)+10.0*pow(P1, 2)*(16.0+21.0*pow(P2, 2)))* 
    pow(Q1, 5)+(-60.0)*P1*P2*(44.0+19.0*pow(P1, 2)+25.0*pow(P2, 2))*pow(Q1, 4)*Q2+10.0*(64.0+ 
    75.0*pow(P1, 4)+368.0*pow(P2, 2)+123.0*pow(P2, 4)+pow(P1, 2)*(272.0+366.0*pow(P2, 2)))*pow(Q1, 3)* 
    pow(Q2, 2)+(-120.0)*P1*P2*(28.0+9.0*pow(P1, 2)+19.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 3)+(-5.0)* 
    (64.0+339.0*pow(P1, 4)+(-64.0)*pow(P2, 2)+(-45.0)*pow(P2, 4)+pow(P1, 2)*(704.0+78.0*pow(P2, 2)))* 
    Q1*pow(Q2, 4)+60.0*P1*P2*(20.0+17.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q2, 5))+(-140.0)* 
    pow(G, 2)*((-2.0)*P1*P2*Q2*((54.0+33.0*pow(P2, 2))*pow(Q1, 2)+((-18.0)+pow(P2, 2))* 
    pow(Q2, 2))+6.0*pow(P1, 2)*Q1*((3.0+4.0*pow(P2, 2))*pow(Q1, 2)+3.0*((-3.0)+2.0*pow(P2, 2))* 
    pow(Q2, 2))+pow(P2, 2)*Q1*((-1.0)*(18.0+13.0*pow(P2, 2))*pow(Q1, 2)+3.0*(18.0+7.0*pow(P2, 2))* 
    pow(Q2, 2))+pow(P1, 4)*(5.0*pow(Q1, 3)+(-33.0)*Q1*pow(Q2, 2))+pow(P1, 3)*P2*((-42.0)* 
    pow(Q1, 2)*Q2+38.0*pow(Q2, 3))))
    ) + 
    (s7L - s7L0) * ( 
    (-15.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(5.0*pow(P1, 4)*P2*Q1+(-10.0)*pow(P1, 2)*pow(P2, 3)* 
    Q1+pow(P2, 5)*Q1+pow(P1, 5)*Q2+(-10.0)*pow(P1, 3)*pow(P2, 2)*Q2+5.0*P1*pow(P2, 4)*Q2)+ 
    2.0*pow(G, 2)*(7.0*pow(P1, 4)*P2*Q1*(3.0*pow(Q1, 2)+(-29.0)*pow(Q2, 2))+2.0*pow(P1, 2)*P2* 
    Q1*((204.0+47.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(612.0+pow(P2, 2))*pow(Q2, 2))+P1*pow(P2, 2)* 
    Q2*((-1.0)*(1224.0+173.0*pow(P2, 2))*pow(Q1, 2)+(408.0+11.0*pow(P2, 2))*pow(Q2, 2))+ 
    pow(P2, 3)*Q1*((-1.0)*(136.0+23.0*pow(P2, 2))*pow(Q1, 2)+(408.0+41.0*pow(P2, 2))*pow(Q2, 2))+ 
    2.0*pow(P1, 3)*Q2*((204.0+(-31.0)*pow(P2, 2))*pow(Q1, 2)+((-68.0)+57.0*pow(P2, 2))*pow(Q2, 2)) 
    +pow(P1, 5)*(47.0*pow(Q1, 2)*Q2+(-25.0)*pow(Q2, 3)))+3.0*(P1*pow(P2, 4)*Q2*(529.0* 
    pow(Q1, 4)+326.0*pow(Q1, 2)*pow(Q2, 2)+(-59.0)*pow(Q2, 4))+2.0*P1*pow(P2, 2)*Q2*(3.0*(772.0+ 
    81.0*pow(P1, 2))*pow(Q1, 4)+2.0*(132.0+(-119.0)*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(516.0+ 
    113.0*pow(P1, 2))*pow(Q2, 4))+P1*Q2*((1200.0+856.0*pow(P1, 2)+37.0*pow(P1, 4))*pow(Q1, 4)+( 
    -2.0)*(1200.0+2488.0*pow(P1, 2)+225.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(240.0+824.0*pow(P1, 2)+ 
    105.0*pow(P1, 4))*pow(Q2, 4))+P2*Q1*((240.0+(-696.0)*pow(P1, 2)+(-47.0)*pow(P1, 4))* 
    pow(Q1, 4)+(-2.0)*(1200.0+1416.0*pow(P1, 2)+17.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(1200.0+ 
    6312.0*pow(P1, 2)+829.0*pow(P1, 4))*pow(Q2, 4))+pow(P2, 5)*(85.0*pow(Q1, 5)+(-298.0)*pow(Q1, 3)* 
    pow(Q2, 2)+(-15.0)*Q1*pow(Q2, 4))+pow(P2, 3)*((712.0+(-138.0)*pow(P1, 2))*pow(Q1, 5)+(-4.0)*( 
    964.0+219.0*pow(P1, 2))*pow(Q1, 3)*pow(Q2, 2)+2.0*(148.0+223.0*pow(P1, 2))*Q1*pow(Q2, 4))))
    ) + 
    (s8L - s8L0) * ( 
    (105.0/32.0)*pow(G, -4.0)*(2.0*pow(G, 2)*((-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0) 
    *pow(Q2, 2))+pow(P2, 4)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+4.0*pow(P1, 3)*P2*Q2*((-3.0)* 
    pow(Q1, 2)+pow(Q2, 2))+(-4.0)*P1*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 4)*( 
    pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))+3.0*((-3.0)*pow(P1, 3)*P2*Q2*(3.0*pow(Q1, 4)+(-38.0)* 
    pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+pow(P1, 2)*Q1*((7.0+15.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*(35.0+ 
    3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(35.0+(-69.0)*pow(P2, 2))*pow(Q2, 4))+P1*P2*Q2*(( 
    -1.0)*(70.0+61.0*pow(P2, 2))*pow(Q1, 4)+2.0*(70.0+13.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+7.0*(( 
    -2.0)+pow(P2, 2))*pow(Q2, 4))+(-1.0)*pow(P2, 2)*Q1*((7.0+6.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*( 
    35.0+18.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(35.0+6.0*pow(P2, 2))*pow(Q2, 4))+pow(P1, 4)*(pow(Q1, 5)+( 
    -34.0)*pow(Q1, 3)*pow(Q2, 2)+29.0*Q1*pow(Q2, 4)))) 
    ) + 
    (s9L - s9L0) * ( 
    (7.0/128.0)*pow(G, -4.0)*(15.0*pow(P1, 4)*P2*Q1*(pow(Q1, 4)+(-130.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    125.0*pow(Q2, 4))+10.0*pow(P1, 3)*Q2*(5.0*(52.0+5.0*pow(P2, 2))*pow(Q1, 4)+10.0*((-52.0)+ 
    31.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(52.0+(-67.0)*pow(P2, 2))*pow(Q2, 4))+10.0*pow(P1, 2)*P2* 
    Q1*((156.0+49.0*pow(P2, 2))*pow(Q1, 4)+(-130.0)*(12.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*( 
    156.0+(-23.0)*pow(P2, 2))*pow(Q2, 4))+15.0*P1*pow(P2, 2)*Q2*((-5.0)*(104.0+19.0* 
    pow(P2, 2))*pow(Q1, 4)+10.0*(104.0+7.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-104.0)+5.0*pow(P2, 2))* 
    pow(Q2, 4))+(-1.0)*pow(P2, 3)*Q1*((520.0+101.0*pow(P2, 2))*pow(Q1, 4)+(-650.0)*(8.0+pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+5.0*(520.0+29.0*pow(P2, 2))*pow(Q2, 4))+pow(P1, 5)*(235.0*pow(Q1, 4)*Q2+( 
    -830.0)*pow(Q1, 2)*pow(Q2, 3)+119.0*pow(Q2, 5))+10.0*pow(G, 2)*(5.0*pow(P1, 4)*P2*Q1*( 
    pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-10.0)*pow(P1, 2)*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+ 
    pow(P2, 5)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+10.0*pow(P1, 3)*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+(-5.0)*P1*pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 5)*(3.0*pow(Q1, 2)* 
    Q2+(-1.0)*pow(Q2, 3)))) 
    ) + 
    (s10L - s10L0) * ( 
    (-189.0/32.0)*pow(G, -4.0)*((-4.0)*pow(P1, 3)*P2*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+4.0*P1*pow(P2, 3)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+pow(P2, 4)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 4)*( 
    pow(Q1, 5)+(-10.0)*pow(Q1, 3)*pow(Q2, 2)+5.0*Q1*pow(Q2, 4))) 
    ) + 
    (s11L - s11L0) * ( 
    (-63.0/128.0)*pow(G, -4.0)*((-10.0)*pow(P1, 3)*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+(-10.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+5.0*P1*pow(P2, 4)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+5.0*pow(P1, 4)*P2*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+(-10.0)*pow(P1, 2)*pow(P2, 3)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+pow(P2, 5)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 5)*(5.0* 
    pow(Q1, 4)*Q2+(-10.0)*pow(Q1, 2)*pow(Q2, 3)+pow(Q2, 5))) 
    );


    // P12
    I[2] = (Lf - L0) * ( 
    (-15.0/2.0)*pow(G, -4.0)*(pow(G, 4)*(4.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))*((-1.0)*P2*Q1+P1* 
    Q2)+28.0*pow(G, 2)*((-3.0)*P1*pow(P2, 2)*Q2*(3.0*pow(Q1, 2)+pow(Q2, 2))+3.0*P2*Q1*(( 
    2.0+pow(P1, 2))*pow(Q1, 2)+(2.0+3.0*pow(P1, 2))*pow(Q2, 2))+(-1.0)*P1*Q2*(3.0*(2.0+pow(P1, 2))* 
    pow(Q1, 2)+(6.0+5.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(5.0*pow(Q1, 3)+3.0*Q1*pow(Q2, 2)))+(-105.0) 
    *(pow(Q1, 2)+pow(Q2, 2))*((-3.0)*P1*pow(P2, 2)*Q2*(5.0*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*P1* 
    Q2*((8.0+3.0*pow(P1, 2))*pow(Q1, 2)+(8.0+7.0*pow(P1, 2))*pow(Q2, 2))+P2*Q1*((8.0+3.0* 
    pow(P1, 2))*pow(Q1, 2)+(8.0+15.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(7.0*pow(Q1, 3)+3.0*Q1*pow(Q2, 2)) 
    )) 
    ) + 
    (cL - cL0) * ( 
    (15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*((-4.0)*pow(P1, 3)*P2*Q1+(-4.0)*P1*P2*(6.0+ 
    pow(P2, 2))*Q1+5.0*pow(P1, 4)*Q2+6.0*pow(P1, 2)*(6.0+pow(P2, 2))*Q2+(8.0+12.0*pow(P2, 2)+ 
    pow(P2, 4))*Q2)+(-7.0)*pow(G, 2)*((-12.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+5.0*pow(Q2, 2))+ 
    6.0*pow(P1, 2)*Q2*(3.0*(8.0+3.0*pow(P2, 2))*pow(Q1, 2)+5.0*(8.0+pow(P2, 2))*pow(Q2, 2))+(-4.0)* 
    P1*P2*Q1*((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+9.0*(8.0+pow(P2, 2))*pow(Q2, 2))+3.0*Q2*(( 
    16.0+48.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q1, 2)+(16.0+16.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+5.0* 
    pow(P1, 4)*(3.0*pow(Q1, 2)*Q2+7.0*pow(Q2, 3)))+21.0*((-4.0)*pow(P1, 3)*P2*Q1*(3.0* 
    pow(Q1, 4)+30.0*pow(Q1, 2)*pow(Q2, 2)+35.0*pow(Q2, 4))+6.0*pow(P1, 2)*Q2*(15.0*(2.0+pow(P2, 2))* 
    pow(Q1, 4)+10.0*(10.0+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+7.0*(10.0+pow(P2, 2))*pow(Q2, 4))+(-4.0)* 
    P1*P2*Q1*((30.0+7.0*pow(P2, 2))*pow(Q1, 4)+30.0*(6.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+15.0* 
    (10.0+pow(P2, 2))*pow(Q2, 4))+Q2*(5.0*(16.0+60.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 4)+10.0*(16.0+ 
    36.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(80.0+60.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q2, 4))+ 
    pow(P1, 4)*(15.0*pow(Q1, 4)*Q2+70.0*pow(Q1, 2)*pow(Q2, 3)+63.0*pow(Q2, 5))))
    ) + 
    (c2L - c2L0) * ( 
    (-15.0/2.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 3)*Q1+P1*(2.0+3.0*pow(P2, 2))*Q1+(-3.0)* 
    pow(P1, 2)*P2*Q2+(-1.0)*P2*(2.0+pow(P2, 2))*Q2)+(-7.0)*pow(G, 2)*((-3.0)*pow(P1, 2)* 
    P2*Q2*(9.0*pow(Q1, 2)+5.0*pow(Q2, 2))+(-1.0)*P2*Q2*(3.0*(8.0+5.0*pow(P2, 2))*pow(Q1, 2)+ 
    (8.0+3.0*pow(P2, 2))*pow(Q2, 2))+P1*Q1*((8.0+15.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+9.0*pow(P2, 2)) 
    *pow(Q2, 2))+3.0*pow(P1, 3)*(pow(Q1, 3)+5.0*Q1*pow(Q2, 2)))+21.0*((-3.0)*pow(P1, 2)*P2* 
    Q2*(15.0*pow(Q1, 4)+30.0*pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+(-1.0)*P2*Q2*(5.0*(10.0+ 
    7.0*pow(P2, 2))*pow(Q1, 4)+30.0*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(10.0+3.0*pow(P2, 2))*pow(Q2, 4)) 
    +P1*Q1*((10.0+21.0*pow(P2, 2))*pow(Q1, 4)+30.0*(2.0+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    (10.0+9.0*pow(P2, 2))*pow(Q2, 4))+pow(P1, 3)*(3.0*pow(Q1, 5)+30.0*pow(Q1, 3)*pow(Q2, 2)+35.0*Q1* 
    pow(Q2, 4))))
    ) + 
    (c3L - c3L0) * ( 
    (-5.0/16.0)*pow(G, -4.0)*(pow(G, 4)*(4.0*pow(P1, 3)*P2*Q1+12.0*P1*P2*(4.0+pow(P2, 2))* 
    Q1+5.0*pow(P1, 4)*Q2+(-6.0)*pow(P1, 2)*((-4.0)+pow(P2, 2))*Q2+(-3.0)*pow(P2, 2)*(8.0+ 
    pow(P2, 2))*Q2)+14.0*pow(G, 2)*(3.0*pow(P1, 4)*Q2*(pow(Q1, 2)+(-7.0)*pow(Q2, 2))+(-12.0)* 
    pow(P1, 3)*P2*Q1*(pow(Q1, 2)+pow(Q2, 2))+(-36.0)*P1*P2*(4.0+pow(P2, 2))*Q1*(pow(Q1, 2)+ 
    pow(Q2, 2))+6.0*pow(P1, 2)*Q2*(3.0*(4.0+3.0*pow(P2, 2))*pow(Q1, 2)+((-20.0)+pow(P2, 2))*pow(Q2, 2)) 
    +Q2*(3.0*(16.0+72.0*pow(P2, 2)+9.0*pow(P2, 4))*pow(Q1, 2)+((-16.0)+24.0*pow(P2, 2)+3.0*pow(P2, 4)) 
    *pow(Q2, 2)))+(-84.0)*((-8.0)*pow(P1, 3)*P2*pow(Q1, 3)*(pow(Q1, 2)+5.0*pow(Q2, 2))+6.0* 
    pow(P1, 2)*Q2*(5.0*(3.0+2.0*pow(P2, 2))*pow(Q1, 4)+10.0*(1.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    -21.0)*pow(Q2, 4))+(-4.0)*P1*P2*Q1*((27.0+7.0*pow(P2, 2))*pow(Q1, 4)+10.0*(9.0+2.0* 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(3.0+pow(P2, 2))*pow(Q2, 4))+Q2*(5.0*(12.0+54.0*pow(P2, 2)+ 
    7.0*pow(P2, 4))*pow(Q1, 4)+20.0*(2.0+9.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+((-20.0)+6.0* 
    pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+pow(P1, 4)*(5.0*pow(Q1, 4)*Q2+(-21.0)*pow(Q2, 5))))
    ) + 
    (c4L - c4L0) * ( 
    (15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 3)*Q1+(-3.0)*P1*pow(P2, 2)*Q1+(-3.0)* 
    pow(P1, 2)*P2*Q2+pow(P2, 3)*Q2)+(-56.0)*pow(G, 2)*(3.0*pow(P1, 3)*Q1*pow(Q2, 2)+(-3.0)* 
    pow(P1, 2)*P2*pow(Q2, 3)+P2*Q2*(3.0*(1.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-1.0)* 
    P1*((1.0+3.0*pow(P2, 2))*pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))+(-42.0)*(3.0*pow(P1, 2)*P2* 
    Q2*((-5.0)*pow(Q1, 4)+10.0*pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+P2*Q2*((-5.0)*(8.0+7.0* 
    pow(P2, 2))*pow(Q1, 4)+(-10.0)*pow(P2, 2)*pow(Q1, 2)*pow(Q2, 2)+(8.0+pow(P2, 2))*pow(Q2, 4))+P1* 
    Q1*((8.0+21.0*pow(P2, 2))*pow(Q1, 4)+30.0*pow(P2, 2)*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*(8.0+3.0* 
    pow(P2, 2))*pow(Q2, 4))+pow(P1, 3)*(pow(Q1, 5)+(-10.0)*pow(Q1, 3)*pow(Q2, 2)+(-35.0)*Q1*pow(Q2, 4)) 
    ))
    ) + 
    (c5L - c5L0) * ( 
    (3.0/16.0)*pow(G, -4.0)*(pow(G, 4)*(4.0*pow(P1, 3)*P2*Q1+(-4.0)*P1*pow(P2, 3)*Q1+ 
    pow(P1, 4)*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q2+pow(P2, 4)*Q2)+14.0*pow(G, 2)*((-4.0)* 
    pow(P1, 3)*P2*Q1*(pow(Q1, 2)+9.0*pow(Q2, 2))+4.0*P1*P2*Q1*((12.0+5.0*pow(P2, 2))* 
    pow(Q1, 2)+(-3.0)*(12.0+pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)*Q2*((-3.0)*(24.0+5.0*pow(P2, 2))* 
    pow(Q1, 2)+(24.0+pow(P2, 2))*pow(Q2, 2))+6.0*pow(P1, 2)*Q2*(3.0*(4.0+pow(P2, 2))*pow(Q1, 2)+((-4.0)+ 
    3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(9.0*pow(Q1, 2)*Q2+(-7.0)*pow(Q2, 3)))+(-84.0)*((-40.0) 
    *pow(P1, 3)*P2*Q1*pow(Q2, 2)*(pow(Q1, 2)+pow(Q2, 2))+20.0*P1*P2*Q1*(pow(Q1, 2)+ 
    pow(Q2, 2))*((3.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*(9.0+pow(P2, 2))*pow(Q2, 2))+Q2*((-5.0)*(4.0+ 
    30.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q1, 4)+20.0*(2.0+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-4.0)+ 
    18.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 4))+pow(P1, 4)*(5.0*pow(Q1, 4)*Q2+20.0*pow(Q1, 2)*pow(Q2, 3)+( 
    -9.0)*pow(Q2, 5))+6.0*pow(P1, 2)*(5.0*pow(Q1, 4)*Q2+10.0*(3.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 3)+(( 
    -7.0)+2.0*pow(P2, 2))*pow(Q2, 5)))) 
    ) + 
    (c6L - c6L0) * ( 
    (-35.0/2.0)*pow(G, -4.0)*(pow(G, 2)*((-3.0)*P1*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2)) 
    +3.0*pow(P1, 2)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 3)*Q2*((-3.0)* 
    pow(Q1, 2)+pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))+3.0*(3.0*pow(P1, 2)*P2* 
    Q2*(5.0*pow(Q1, 4)+10.0*pow(Q1, 2)*pow(Q2, 2)+(-3.0)*pow(Q2, 4))+P1*Q1*((2.0+9.0*pow(P2, 2)) 
    *pow(Q1, 4)+(-10.0)*(2.0+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(2.0+(-3.0)*pow(P2, 2))* 
    pow(Q2, 4))+P2*Q2*((-5.0)*(2.0+3.0*pow(P2, 2))*pow(Q1, 4)+10.0*(2.0+pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+((-2.0)+pow(P2, 2))*pow(Q2, 4))+(-1.0)*pow(P1, 3)*(pow(Q1, 5)+10.0*pow(Q1, 3)*pow(Q2, 2)+( 
    -15.0)*Q1*pow(Q2, 4)))) 
    ) + 
    (c7L - c7L0) * ( 
    (-15.0/16.0)*pow(G, -4.0)*(2.0*pow(G, 2)*(4.0*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+(-3.0)* 
    pow(Q2, 2))+(-4.0)*P1*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)* 
    Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+ 
    pow(P1, 4)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))+3.0*((-4.0)*pow(P1, 3)*P2*Q1*(3.0* 
    pow(Q1, 4)+10.0*pow(Q1, 2)*pow(Q2, 2)+(-25.0)*pow(Q2, 4))+6.0*pow(P1, 2)*Q2*(5.0*(4.0+3.0* 
    pow(P2, 2))*pow(Q1, 4)+10.0*((-4.0)+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(4.0+(-5.0)*pow(P2, 2))* 
    pow(Q2, 4))+pow(P2, 2)*Q2*((-5.0)*(24.0+7.0*pow(P2, 2))*pow(Q1, 4)+30.0*(8.0+pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+((-24.0)+pow(P2, 2))*pow(Q2, 4))+4.0*P1*P2*Q1*((12.0+7.0*pow(P2, 2))* 
    pow(Q1, 4)+(-30.0)*(4.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*((-12.0)+pow(P2, 2))*pow(Q2, 4))+ 
    pow(P1, 4)*(5.0*pow(Q1, 4)*Q2+(-50.0)*pow(Q1, 2)*pow(Q2, 3)+9.0*pow(Q2, 5))))
    ) + 
    (c8L - c8L0) * ( 
    (315.0/16.0)*pow(G, -4.0)*((-3.0)*pow(P1, 2)*P2*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+pow(P2, 3)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-3.0) 
    *P1*pow(P2, 2)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 3)*( 
    pow(Q1, 5)+(-10.0)*pow(Q1, 3)*pow(Q2, 2)+5.0*Q1*pow(Q2, 4)))
    ) + 
    (c9L - c9L0) * ( 
    (35.0/16.0)*pow(G, -4.0)*((-6.0)*pow(P1, 2)*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+(-10.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+pow(P2, 4)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+4.0*pow(P1, 3)*P2*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+(-4.0) 
    *P1*pow(P2, 3)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 4)*(5.0* 
    pow(Q1, 4)*Q2+(-10.0)*pow(Q1, 2)*pow(Q2, 3)+pow(Q2, 5)))
    ) + 
    (sL - sL0) * ( 
    (15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*((8.0+pow(P1, 4)+36.0*pow(P2, 2)+5.0*pow(P2, 4)+6.0*pow(P1, 2)*(2.0+ 
    pow(P2, 2)))*Q1+(-4.0)*P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))*Q2)+(-7.0)*pow(G, 2)*((48.0+ 
    3.0*pow(P1, 4)+240.0*pow(P2, 2)+35.0*pow(P2, 4)+6.0*pow(P1, 2)*(8.0+5.0*pow(P2, 2)))*pow(Q1, 3)+(-12.0) 
    *P1*P2*(24.0+3.0*pow(P1, 2)+5.0*pow(P2, 2))*pow(Q1, 2)*Q2+3.0*(16.0+5.0*pow(P1, 4)+48.0* 
    pow(P2, 2)+5.0*pow(P2, 4)+6.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2)))*Q1*pow(Q2, 2)+(-4.0)*P1*P2*( 
    5.0*pow(P1, 2)+3.0*(8.0+pow(P2, 2)))*pow(Q2, 3))+21.0*((80.0+3.0*pow(P1, 4)+420.0*pow(P2, 2)+63.0* 
    pow(P2, 4)+6.0*pow(P1, 2)*(10.0+7.0*pow(P2, 2)))*pow(Q1, 5)+(-20.0)*P1*P2*(30.0+3.0*pow(P1, 2)+ 
    7.0*pow(P2, 2))*pow(Q1, 4)*Q2+10.0*(16.0+3.0*pow(P1, 4)+60.0*pow(P2, 2)+7.0*pow(P2, 4)+18.0* 
    pow(P1, 2)*(2.0+pow(P2, 2)))*pow(Q1, 3)*pow(Q2, 2)+(-120.0)*P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 3)+5.0*(16.0+7.0*pow(P1, 4)+36.0*pow(P2, 2)+3.0*pow(P2, 4)+6.0*pow(P1, 2)*(10.0+3.0* 
    pow(P2, 2)))*Q1*pow(Q2, 4)+(-4.0)*P1*P2*(30.0+7.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q2, 5)))
    ) + 
    (s2L - s2L0) * ( 
    (15.0/2.0)*pow(G, -4.0)*(2.0*pow(G, 4)*(P2*Q1+pow(P2, 3)*Q1+P1*(1.0+pow(P1, 2))*Q2)+( 
    -7.0)*pow(G, 2)*(16.0*P1*pow(Q2, 3)+3.0*P1*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+ 
    3.0*pow(P2, 3)*Q1*(5.0*pow(Q1, 2)+pow(Q2, 2))+3.0*pow(P1, 3)*Q2*(pow(Q1, 2)+5.0*pow(Q2, 2))+P2* 
    ((16.0+3.0*pow(P1, 2))*pow(Q1, 3)+(-9.0)*pow(P1, 2)*Q1*pow(Q2, 2)))+42.0*(2.0*pow(P2, 3)*(7.0* 
    pow(Q1, 5)+5.0*pow(Q1, 3)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-5.0)*pow(Q1, 4)+pow(Q2, 4))+P2* 
    Q1*(3.0*(5.0+pow(P1, 2))*pow(Q1, 4)+10.0*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*(1.0+3.0*pow(P1, 2))* 
    pow(Q2, 4))+P1*Q2*((-5.0)*pow(Q1, 4)+10.0*(1.0+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(15.0+14.0* 
    pow(P1, 2))*pow(Q2, 4)))) 
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/16.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*pow(P1, 4)*Q1+6.0*pow(P1, 2)*(4.0+pow(P2, 2))*Q1+( 
    -1.0)*pow(P2, 2)*(24.0+5.0*pow(P2, 2))*Q1+(-12.0)*pow(P1, 3)*P2*Q2+(-4.0)*P1*P2*( 
    12.0+pow(P2, 2))*Q2)+(-14.0)*pow(G, 2)*(((-16.0)+3.0*pow(P1, 4)+(-120.0)*pow(P2, 2)+(-21.0)* 
    pow(P2, 4)+6.0*pow(P1, 2)*(4.0+pow(P2, 2)))*pow(Q1, 3)+(-12.0)*P1*P2*(12.0+3.0*pow(P1, 2)+ 
    pow(P2, 2))*pow(Q1, 2)*Q2+3.0*(16.0+9.0*pow(P1, 4)+24.0*pow(P2, 2)+pow(P2, 4)+18.0*pow(P1, 2)*(4.0+ 
    pow(P2, 2)))*Q1*pow(Q2, 2)+(-12.0)*P1*P2*(12.0+3.0*pow(P1, 2)+pow(P2, 2))*pow(Q2, 3))+84.0* 
    (((-20.0)+6.0*pow(P1, 2)+pow(P1, 4)+(-126.0)*pow(P2, 2)+(-21.0)*pow(P2, 4))*pow(Q1, 5)+(-20.0)* 
    P1*(3.0+pow(P1, 2))*P2*pow(Q1, 4)*Q2+20.0*(2.0+pow(P1, 4)+3.0*pow(P2, 2)+3.0*pow(P1, 2)*(3.0+ 
    pow(P2, 2)))*pow(Q1, 3)*pow(Q2, 2)+(-40.0)*P1*P2*(9.0+2.0*pow(P1, 2)+pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 3)+5.0*(12.0+7.0*pow(P1, 4)+18.0*pow(P2, 2)+pow(P2, 4)+6.0*pow(P1, 2)*(9.0+2.0*pow(P2, 2)))* 
    Q1*pow(Q2, 4)+(-4.0)*P1*P2*(27.0+7.0*pow(P1, 2)+2.0*pow(P2, 2))*pow(Q2, 5)))
    ) + 
    (s4L - s4L0) * ( 
    (-15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*pow(P1, 2)*P2*Q1+(-1.0)*pow(P2, 3)*Q1+pow(P1, 3)* 
    Q2+(-3.0)*P1*pow(P2, 2)*Q2)+28.0*pow(G, 2)*(3.0*P1*pow(P2, 2)*Q2*(3.0*pow(Q1, 2)+ 
    pow(Q2, 2))+P1*Q2*(3.0*(2.0+pow(P1, 2))*pow(Q1, 2)+(-1.0)*(2.0+3.0*pow(P1, 2))*pow(Q2, 2))+( 
    -1.0)*P2*Q1*(((-2.0)+3.0*pow(P1, 2))*pow(Q1, 2)+3.0*(2.0+3.0*pow(P1, 2))*pow(Q2, 2))+3.0* 
    pow(P2, 3)*(pow(Q1, 3)+(-1.0)*Q1*pow(Q2, 2)))+(-84.0)*(3.0*P1*pow(P2, 2)*Q2*(5.0* 
    pow(Q1, 4)+10.0*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*P2*Q1*(3.0*((-2.0)+pow(P1, 2))* 
    pow(Q1, 4)+10.0*(2.0+3.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(2.0+3.0*pow(P1, 2))*pow(Q2, 4))+P1* 
    Q2*(5.0*(2.0+pow(P1, 2))*pow(Q1, 4)+10.0*(2.0+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(6.0+7.0* 
    pow(P1, 2))*pow(Q2, 4))+pow(P2, 3)*(7.0*pow(Q1, 5)+(-10.0)*pow(Q1, 3)*pow(Q2, 2)+(-5.0)*Q1* 
    pow(Q2, 4)))) 
    ) + 
    (s5L - s5L0) * ( 
    (3.0/16.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 4)*Q1+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1+pow(P2, 4)* 
    Q1+(-4.0)*pow(P1, 3)*P2*Q2+4.0*P1*pow(P2, 3)*Q2)+(-84.0)*(((-4.0)+pow(P1, 4)+(-42.0) 
    *pow(P2, 2)+(-9.0)*pow(P2, 4)+6.0*pow(P1, 2)*(3.0+2.0*pow(P2, 2)))*pow(Q1, 5)+(-20.0)*P1*P2* 
    (9.0+pow(P1, 2)+2.0*pow(P2, 2))*pow(Q1, 4)*Q2+20.0*(2.0+9.0*pow(P2, 2)+pow(P2, 4)+3.0*pow(P1, 2)*(1.0+ 
    pow(P2, 2)))*pow(Q1, 3)*pow(Q2, 2)+(-40.0)*P1*P2*(3.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 3)+(-5.0) 
    *(4.0+30.0*pow(P1, 2)+5.0*pow(P1, 4)+(-6.0)*pow(P2, 2)+(-1.0)*pow(P2, 4))*Q1*pow(Q2, 4)+20.0* 
    P1*(3.0+pow(P1, 2))*P2*pow(Q2, 5))+14.0*pow(G, 2)*(4.0*pow(P1, 3)*P2*Q2*((-3.0)* 
    pow(Q1, 2)+5.0*pow(Q2, 2))+(-4.0)*P1*P2*Q2*(9.0*(4.0+pow(P2, 2))*pow(Q1, 2)+((-12.0)+ 
    pow(P2, 2))*pow(Q2, 2))+6.0*pow(P1, 2)*Q1*((4.0+3.0*pow(P2, 2))*pow(Q1, 2)+3.0*((-4.0)+pow(P2, 2)) 
    *pow(Q2, 2))+pow(P2, 2)*Q1*((-1.0)*(24.0+7.0*pow(P2, 2))*pow(Q1, 2)+9.0*(8.0+pow(P2, 2))* 
    pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 3)+(-15.0)*Q1*pow(Q2, 2))))  
    ) + 
    (s6L - s6L0) * ( 
    (35.0/2.0)*pow(G, -4.0)*(pow(G, 2)*(3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+( 
    -1.0)*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-3.0)* 
    pow(Q1, 2)+pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))+(-6.0)*((-2.0)* 
    pow(P2, 3)*(pow(Q1, 5)+(-5.0)*pow(Q1, 3)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*Q2*((-5.0)*pow(Q1, 4)+ 
    pow(Q2, 4))+(-1.0)*P1*Q2*(5.0*pow(Q1, 4)+(-10.0)*(1.0+pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(1.0+ 
    2.0*pow(P1, 2))*pow(Q2, 4))+P2*Q1*(((-1.0)+3.0*pow(P1, 2))*pow(Q1, 4)+10.0*pow(Q1, 2)* 
    pow(Q2, 2)+(-5.0)*(1.0+3.0*pow(P1, 2))*pow(Q2, 4))))
    ) + 
    (s7L - s7L0) * ( 
    (-15.0/16.0)*pow(G, -4.0)*(2.0*pow(G, 2)*((-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0) 
    *pow(Q2, 2))+pow(P2, 4)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+4.0*pow(P1, 3)*P2*Q2*((-3.0)* 
    pow(Q1, 2)+pow(Q2, 2))+(-4.0)*P1*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 4)*( 
    pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))+(-3.0)*(4.0*pow(P1, 3)*P2*Q2*((-5.0)*pow(Q1, 4)+( 
    -30.0)*pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+(-6.0)*pow(P1, 2)*Q1*((4.0+5.0*pow(P2, 2))* 
    pow(Q1, 4)+(-10.0)*(4.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(4.0+(-3.0)*pow(P2, 2))*pow(Q2, 4))+ 
    4.0*P1*P2*Q2*(5.0*(12.0+5.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(12.0+pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+(-3.0)*((-4.0)+pow(P2, 2))*pow(Q2, 4))+pow(P2, 2)*Q1*(3.0*(8.0+3.0*pow(P2, 2))* 
    pow(Q1, 4)+(-10.0)*(24.0+5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(24.0+pow(P2, 2))*pow(Q2, 4))+ 
    pow(P1, 4)*(pow(Q1, 5)+30.0*pow(Q1, 3)*pow(Q2, 2)+(-35.0)*Q1*pow(Q2, 4))))
    ) + 
    (s8L - s8L0) * ( 
    (-315.0/16.0)*pow(G, -4.0)*((-3.0)*P1*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+3.0*pow(P1, 2)*P2*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+(-1.0)*pow(P2, 3)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+ 
    pow(P1, 3)*(5.0*pow(Q1, 4)*Q2+(-10.0)*pow(Q1, 2)*pow(Q2, 3)+pow(Q2, 5)))
    ) + 
    (s9L - s9L0) * ( 
    (35.0/16.0)*pow(G, -4.0)*((-4.0)*pow(P1, 3)*P2*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+4.0*P1*pow(P2, 3)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+pow(P2, 4)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 4)*( 
    pow(Q1, 5)+(-10.0)*pow(Q1, 3)*pow(Q2, 2)+5.0*Q1*pow(Q2, 4)))
    );

    // P21
    I[3] = (Lf - L0) * ( 
    (-15.0/4.0)*pow(G, -4.0)*(pow(G, 4)*((-21.0)*pow(P1, 3)*P2*Q1+(-21.0)*P1*P2*(2.0+ 
    pow(P2, 2))*Q1+18.0*pow(P1, 4)*Q2+pow(P1, 2)*(41.0+15.0*pow(P2, 2))*Q2+(-1.0)*((-4.0)+ 
    pow(P2, 2)+3.0*pow(P2, 4))*Q2)+(-28.0)*pow(G, 2)*((-7.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+3.0* 
    pow(Q2, 2))+(-1.0)*((-1.0)+pow(P2, 2))*Q2*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+(2.0+pow(P2, 2))* 
    pow(Q2, 2))+pow(P1, 2)*Q2*((19.0+21.0*pow(P2, 2))*pow(Q1, 2)+3.0*(7.0+pow(P2, 2))*pow(Q2, 2))+(-1.0) 
    *P1*P2*Q1*((20.0+13.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+pow(P2, 2))*pow(Q2, 2))+2.0* 
    pow(P1, 4)*(3.0*pow(Q1, 2)*Q2+5.0*pow(Q2, 3)))+21.0*(pow(Q1, 2)+pow(Q2, 2))*((-21.0)*pow(P1, 3)* 
    P2*Q1*(pow(Q1, 2)+5.0*pow(Q2, 2))+(-3.0)*P1*P2*Q1*((26.0+19.0*pow(P2, 2))*pow(Q1, 2)+ 
    (-1.0)*((-34.0)+pow(P2, 2))*pow(Q2, 2))+(-1.0)*((-1.0)+pow(P2, 2))*Q2*((8.0+15.0*pow(P2, 2)) 
    *pow(Q1, 2)+(8.0+3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 2)*Q2*((73.0+111.0*pow(P2, 2))*pow(Q1, 2)+( 
    85.0+3.0*pow(P2, 2))*pow(Q2, 2))+6.0*pow(P1, 4)*(3.0*pow(Q1, 2)*Q2+7.0*pow(Q2, 3))))
    ) + 
    (cL - cL0) * ( 
    (3.0/128.0)*pow(G, -4.0)*((-5.0)*pow(G, 4)*(43.0*(16.0+48.0*pow(P1, 2)+5.0*pow(P1, 4))*P2* 
    Q1+86.0*(8.0+3.0*pow(P1, 2))*pow(P2, 3)*Q1+43.0*pow(P2, 5)*Q1+(-1.0)*P1*(1296.0+ 
    2288.0*pow(P1, 2)+205.0*pow(P1, 4))*Q2+(-2.0)*P1*(456.0+119.0*pow(P1, 2))*pow(P2, 2)*Q2+ 
    (-33.0)*P1*pow(P2, 4)*Q2)+140.0*pow(G, 2)*((-2.0)*P1*pow(P2, 2)*Q2*((678.0+133.0* 
    pow(P1, 2))*pow(Q1, 2)+(78.0+35.0*pow(P1, 2))*pow(Q2, 2))+2.0*pow(P2, 3)*Q1*((218.0+63.0* 
    pow(P1, 2))*pow(Q1, 2)+(34.0+69.0*pow(P1, 2))*pow(Q2, 2))+(-1.0)*P1*Q2*((624.0+796.0* 
    pow(P1, 2)+53.0*pow(P1, 4))*pow(Q1, 2)+(656.0+1260.0*pow(P1, 2)+119.0*pow(P1, 4))*pow(Q2, 2))+P2* 
    Q1*((336.0+708.0*pow(P1, 2)+55.0*pow(P1, 4))*pow(Q1, 2)+(368.0+2004.0*pow(P1, 2)+265.0* 
    pow(P1, 4))*pow(Q2, 2))+pow(P2, 5)*(31.0*pow(Q1, 3)+(-7.0)*Q1*pow(Q2, 2))+P1*pow(P2, 4)*(( 
    -93.0)*pow(Q1, 2)*Q2+9.0*pow(Q2, 3)))+42.0*(5.0*P1*pow(P2, 4)*Q2*(325.0*pow(Q1, 4)+94.0* 
    pow(Q1, 2)*pow(Q2, 2)+(-31.0)*pow(Q2, 4))+(-10.0)*pow(P2, 3)*Q1*((484.0+113.0*pow(P1, 2))* 
    pow(Q1, 4)+2.0*(196.0+191.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*((-12.0)+17.0*pow(P1, 2))* 
    pow(Q2, 4))+10.0*P1*pow(P2, 2)*Q2*((1836.0+287.0*pow(P1, 2))*pow(Q1, 4)+2.0*(876.0+245.0* 
    pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(12.0+35.0*pow(P1, 2))*pow(Q2, 4))+(-5.0)*P2*Q1*((664.0+ 
    1080.0*pow(P1, 2)+67.0*pow(P1, 4))*pow(Q1, 4)+2.0*(712.0+3096.0*pow(P1, 2)+325.0*pow(P1, 4))* 
    pow(Q1, 2)*pow(Q2, 2)+5.0*(152.0+984.0*pow(P1, 2)+147.0*pow(P1, 4))*pow(Q2, 4))+P1*Q2*(5.0*( 
    1224.0+1224.0*pow(P1, 2)+65.0*pow(P1, 4))*pow(Q1, 4)+10.0*(1272.0+1960.0*pow(P1, 2)+147.0* 
    pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(6600.0+13160.0*pow(P1, 2)+1281.0*pow(P1, 4))*pow(Q2, 4))+ 
    pow(P2, 5)*((-371.0)*pow(Q1, 5)+(-10.0)*pow(Q1, 3)*pow(Q2, 2)+145.0*Q1*pow(Q2, 4))))
    ) + 
    (c2L - c2L0) * ( 
    (-15.0/32.0)*pow(G, -4.0)*(pow(G, 4)*((32.0+45.0*pow(P1, 4)+168.0*pow(P2, 2)+53.0*pow(P2, 4)+2.0* 
    pow(P1, 2)*(76.0+93.0*pow(P2, 2)))*Q1+(-4.0)*P1*P2*(40.0+31.0*pow(P1, 2)+9.0*pow(P2, 2))* 
    Q2)+(-56.0)*pow(G, 2)*((8.0+6.0*pow(P1, 4)+53.0*pow(P2, 2)+19.0*pow(P2, 4)+9.0*pow(P1, 2)*(3.0+5.0* 
    pow(P2, 2)))*pow(Q1, 3)+(-1.0)*P1*P2*(106.0+61.0*pow(P1, 2)+45.0*pow(P2, 2))*pow(Q1, 2)*Q2+ 
    (8.0+27.0*pow(P1, 4)+9.0*pow(P2, 2)+(-4.0)*pow(P2, 4)+pow(P1, 2)*(71.0+51.0*pow(P2, 2)))*Q1* 
    pow(Q2, 2)+3.0*P1*P2*((-6.0)+(-7.0)*pow(P1, 2)+pow(P2, 2))*pow(Q2, 3))+42.0*((32.0+15.0* 
    pow(P1, 4)+236.0*pow(P2, 2)+91.0*pow(P2, 4)+6.0*pow(P1, 2)*(14.0+27.0*pow(P2, 2)))*pow(Q1, 5)+(-24.0) 
    *P1*P2*(24.0+11.0*pow(P1, 2)+13.0*pow(P2, 2))*pow(Q1, 4)*Q2+2.0*(32.0+69.0*pow(P1, 4)+ 
    92.0*pow(P2, 2)+pow(P2, 4)+6.0*pow(P1, 2)*(38.0+45.0*pow(P2, 2)))*pow(Q1, 3)*pow(Q2, 2)+(-32.0)* 
    P1*P2*(17.0+14.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 3)+(32.0+147.0*pow(P1, 4)+(-20.0) 
    *pow(P2, 2)+(-33.0)*pow(P2, 4)+2.0*pow(P1, 2)*(170.0+69.0*pow(P2, 2)))*Q1*pow(Q2, 4)+(-8.0)* 
    P1*P2*(4.0+7.0*pow(P1, 2)+(-3.0)*pow(P2, 2))*pow(Q2, 5)))
    ) + 
    (c3L - c3L0) * ( 
    (-5.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*(208.0+264.0*pow(P1, 2)+5.0*pow(P1, 4))*P2*Q1+ 
    6.0*(164.0+39.0*pow(P1, 2))*pow(P2, 3)*Q1+75.0*pow(P2, 5)*Q1+P1*(624.0+1256.0*pow(P1, 2)+ 
    127.0*pow(P1, 4))*Q2+(-2.0)*P1*(12.0+7.0*pow(P1, 2))*pow(P2, 2)*Q2+3.0*P1*pow(P2, 4)* 
    Q2)+(-56.0)*pow(G, 2)*(2.0*P1*pow(P2, 2)*Q2*((-1.0)*(183.0+28.0*pow(P1, 2))*pow(Q1, 2)+ 
    (57.0+7.0*pow(P1, 2))*pow(Q2, 2))+P2*Q1*(2.0*(122.0+165.0*pow(P1, 2)+8.0*pow(P1, 4))* 
    pow(Q1, 2)+(-3.0)*(36.0+66.0*pow(P1, 2)+11.0*pow(P1, 4))*pow(Q2, 2))+P1*Q2*(((-84.0)+(-46.0) 
    *pow(P1, 2)+pow(P1, 4))*pow(Q1, 2)+2.0*(118.0+217.0*pow(P1, 2)+21.0*pow(P1, 4))*pow(Q2, 2))+6.0* 
    pow(P2, 3)*((63.0+13.0*pow(P1, 2))*pow(Q1, 3)+(-25.0)*Q1*pow(Q2, 2))+15.0*pow(P2, 5)*(2.0* 
    pow(Q1, 3)+(-1.0)*Q1*pow(Q2, 2))+P1*pow(P2, 4)*((-33.0)*pow(Q1, 2)*Q2+12.0*pow(Q2, 3)))+( 
    -42.0)*((-3.0)*P1*pow(P2, 4)*Q2*((-119.0)*pow(Q1, 4)+62.0*pow(Q1, 2)*pow(Q2, 2)+13.0* 
    pow(Q2, 4))+(-2.0)*pow(P2, 3)*Q1*((944.0+169.0*pow(P1, 2))*pow(Q1, 4)+2.0*((-184.0)+91.0* 
    pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(-13.0)*(32.0+7.0*pow(P1, 2))*pow(Q2, 4))+(-2.0)*P1* 
    pow(P2, 2)*Q2*((-1.0)*(1752.0+227.0*pow(P1, 2))*pow(Q1, 4)+6.0*(96.0+pow(P1, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+3.0*(72.0+11.0*pow(P1, 2))*pow(Q2, 4))+(-1.0)*P1*Q2*((-1.0)*(936.0+704.0* 
    pow(P1, 2)+25.0*pow(P1, 4))*pow(Q1, 4)+6.0*(88.0+112.0*pow(P1, 2)+11.0*pow(P1, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+3.0*(360.0+672.0*pow(P1, 2)+65.0*pow(P1, 4))*pow(Q2, 4))+(-3.0)*P2*Q1*((392.0+ 
    464.0*pow(P1, 2)+21.0*pow(P1, 4))*pow(Q1, 4)+2.0*((-8.0)+320.0*pow(P1, 2)+23.0*pow(P1, 4))* 
    pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(280.0+848.0*pow(P1, 2)+111.0*pow(P1, 4))*pow(Q2, 4))+pow(P2, 5)*(( 
    -155.0)*pow(Q1, 5)+110.0*pow(Q1, 3)*pow(Q2, 2)+65.0*Q1*pow(Q2, 4))))
    ) + 
    (c4L - c4L0) * ( 
    (3.0/16.0)*pow(G, -4.0)*(5.0*pow(G, 4)*(9.0*pow(P1, 4)*Q1+pow(P1, 2)*(19.0+3.0*pow(P2, 2))*Q1+ 
    (-1.0)*pow(P2, 2)*(19.0+10.0*pow(P2, 2))*Q1+(-29.0)*pow(P1, 3)*P2*Q2+(-1.0)*P1* 
    P2*(38.0+9.0*pow(P2, 2))*Q2)+(-140.0)*pow(G, 2)*(((-2.0)+pow(P1, 4)+(-20.0)*pow(P2, 2)+( 
    -6.0)*pow(P1, 2)*pow(P2, 2)+(-9.0)*pow(P2, 4))*pow(Q1, 3)+(-2.0)*P1*P2*(11.0+8.0*pow(P1, 2)+ 
    3.0*pow(P2, 2))*pow(Q1, 2)*Q2+(6.0+15.0*pow(P1, 4)+22.0*pow(P2, 2)+7.0*pow(P2, 4)+pow(P1, 2)*(38.0+ 
    24.0*pow(P2, 2)))*Q1*pow(Q2, 2)+(-2.0)*P1*P2*(9.0+7.0*pow(P1, 2)+2.0*pow(P2, 2))*pow(Q2, 3)) 
    +42.0*((-1.0)*(32.0+285.0*pow(P2, 2)+125.0*pow(P2, 4)+35.0*pow(P1, 2)*(1.0+3.0*pow(P2, 2)))* 
    pow(Q1, 5)+(-5.0)*P1*P2*(2.0+11.0*pow(P1, 2)+(-9.0)*pow(P2, 2))*pow(Q1, 4)*Q2+10.0*(8.0+ 
    12.0*pow(P1, 4)+45.0*pow(P2, 2)+17.0*pow(P2, 4)+pow(P1, 2)*(35.0+33.0*pow(P2, 2)))*pow(Q1, 3)* 
    pow(Q2, 2)+(-10.0)*P1*P2*(86.0+53.0*pow(P1, 2)+33.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 3)+5.0*( 
    16.0+48.0*pow(P1, 4)+43.0*pow(P2, 2)+11.0*pow(P2, 4)+9.0*pow(P1, 2)*(13.0+7.0*pow(P2, 2)))*Q1* 
    pow(Q2, 4)+(-5.0)*P1*P2*(26.0+23.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q2, 5)))
    ) + 
    (c5L - c5L0) * ( 
    (3.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(101.0*pow(P1, 4)*P2*Q1+(-1.0)*pow(P2, 3)*(296.0+ 
    39.0*pow(P2, 2))*Q1+pow(P1, 2)*(888.0*P2*Q1+94.0*pow(P2, 3)*Q1)+45.0*pow(P1, 5)*Q2+ 
    2.0*pow(P1, 3)*(148.0+(-77.0)*pow(P2, 2))*Q2+(-1.0)*P1*pow(P2, 2)*(888.0+71.0*pow(P2, 2)) 
    *Q2)+56.0*pow(G, 2)*(P1*pow(P2, 4)*Q2*(53.0*pow(Q1, 2)+6.0*pow(Q2, 2))+2.0*P1* 
    pow(P2, 2)*Q2*((309.0+50.0*pow(P1, 2))*pow(Q1, 2)+9.0*(5.0+pow(P1, 2))*pow(Q2, 2))+(-2.0)* 
    pow(P2, 3)*Q1*(((-97.0)+3.0*pow(P1, 2))*pow(Q1, 2)+(143.0+38.0*pow(P1, 2))*pow(Q2, 2))+P1* 
    Q2*((228.0+250.0*pow(P1, 2)+15.0*pow(P1, 4))*pow(Q1, 2)+(-2.0)*(38.0+91.0*pow(P1, 2)+10.0* 
    pow(P1, 4))*pow(Q2, 2))+(-1.0)*P2*Q1*(2.0*((-38.0)+63.0*pow(P1, 2)+9.0*pow(P1, 4))* 
    pow(Q1, 2)+(228.0+510.0*pow(P1, 2)+47.0*pow(P1, 4))*pow(Q2, 2))+pow(P2, 5)*(20.0*pow(Q1, 3)+(-21.0)* 
    Q1*pow(Q2, 2)))+(-21.0)*(P1*pow(P2, 4)*Q2*(445.0*pow(Q1, 4)+806.0*pow(Q1, 2)*pow(Q2, 2)+( 
    -23.0)*pow(Q2, 4))+2.0*P1*pow(P2, 2)*Q2*((2544.0+403.0*pow(P1, 2))*pow(Q1, 4)+2.0*(2400.0+ 
    397.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+((-48.0)+7.0*pow(P1, 2))*pow(Q2, 4))+2.0*pow(P2, 3)*Q1*( 
    5.0*(256.0+11.0*pow(P1, 2))*pow(Q1, 4)+(-2.0)*(1744.0+419.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    -1.0)*(544.0+189.0*pow(P1, 2))*pow(Q2, 4))+P1*Q2*((1776.0+1856.0*pow(P1, 2)+105.0* 
    pow(P1, 4))*pow(Q1, 4)+2.0*(1872.0+2144.0*pow(P1, 2)+135.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0) 
    *(1104.0+2176.0*pow(P1, 2)+219.0*pow(P1, 4))*pow(Q2, 4))+(-1.0)*P2*Q1*(((-1264.0)+ 
    96.0*pow(P1, 2)+71.0*pow(P1, 4))*pow(Q1, 4)+2.0*(2672.0+5568.0*pow(P1, 2)+509.0*pow(P1, 4))* 
    pow(Q1, 2)*pow(Q2, 2)+(976.0+2592.0*pow(P1, 2)+243.0*pow(P1, 4))*pow(Q2, 4))+pow(P2, 5)*(245.0* 
    pow(Q1, 5)+(-530.0)*pow(Q1, 3)*pow(Q2, 2)+(-71.0)*Q1*pow(Q2, 4))))
    ) + 
    (c6L - c6L0) * ( 
    (1.0/32.0)*pow(G, -4.0)*((-45.0)*pow(G, 4)*(pow(P1, 4)*Q1+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1+ 
    pow(P2, 4)*Q1+(-4.0)*pow(P1, 3)*P2*Q2+4.0*P1*pow(P2, 3)*Q2)+21.0*(((-64.0)+75.0* 
    pow(P1, 4)+(-920.0)*pow(P2, 2)+(-525.0)*pow(P2, 4)+10.0*pow(P1, 2)*(28.0+39.0*pow(P2, 2)))* 
    pow(Q1, 5)+(-60.0)*P1*P2*(64.0+23.0*pow(P1, 2)+41.0*pow(P2, 2))*pow(Q1, 4)*Q2+10.0*(64.0+ 
    21.0*pow(P1, 4)+488.0*pow(P2, 2)+189.0*pow(P2, 4)+2.0*pow(P1, 2)*(76.0+165.0*pow(P2, 2)))* 
    pow(Q1, 3)*pow(Q2, 2)+120.0*P1*P2*((-8.0)+3.0*pow(P1, 2)+(-11.0)*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 3)+(-5.0)*(64.0+261.0*pow(P1, 4)+56.0*pow(P2, 2)+(-3.0)*pow(P2, 4)+2.0*pow(P1, 2)*(292.0+ 
    93.0*pow(P2, 2)))*Q1*pow(Q2, 4)+60.0*P1*P2*(16.0+13.0*pow(P1, 2)+3.0*pow(P2, 2))*pow(Q2, 5)) 
    +(-280.0)*pow(G, 2)*(3.0*pow(P1, 2)*Q1*((3.0+5.0*pow(P2, 2))*pow(Q1, 2)+3.0*((-3.0)+pow(P2, 2)) 
    *pow(Q2, 2))+P1*P2*Q2*((-3.0)*(18.0+13.0*pow(P2, 2))*pow(Q1, 2)+(18.0+pow(P2, 2))* 
    pow(Q2, 2))+pow(P2, 2)*Q1*((-1.0)*(9.0+7.0*pow(P2, 2))*pow(Q1, 2)+3.0*(9.0+4.0*pow(P2, 2))* 
    pow(Q2, 2))+pow(P1, 4)*(2.0*pow(Q1, 3)+(-15.0)*Q1*pow(Q2, 2))+pow(P1, 3)*P2*((-15.0)* 
    pow(Q1, 2)*Q2+17.0*pow(Q2, 3))))
    ) + 
    (c7L - c7L0) * ( 
    (-15.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(5.0*pow(P1, 4)*P2*Q1+(-10.0)*pow(P1, 2)*pow(P2, 3)* 
    Q1+pow(P2, 5)*Q1+pow(P1, 5)*Q2+(-10.0)*pow(P1, 3)*pow(P2, 2)*Q2+5.0*P1*pow(P2, 4)*Q2)+ 
    2.0*pow(G, 2)*(pow(P1, 4)*P2*Q1*(11.0*pow(Q1, 2)+(-173.0)*pow(Q2, 2))+(-2.0)*pow(P1, 3)* 
    Q2*(((-204.0)+pow(P2, 2))*pow(Q1, 2)+(68.0+(-47.0)*pow(P2, 2))*pow(Q2, 2))+P1*pow(P2, 2)* 
    Q2*((-1.0)*(1224.0+203.0*pow(P2, 2))*pow(Q1, 2)+3.0*(136.0+7.0*pow(P2, 2))*pow(Q2, 2))+2.0* 
    pow(P1, 2)*P2*Q1*(3.0*(68.0+19.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(612.0+31.0*pow(P2, 2))* 
    pow(Q2, 2))+pow(P2, 3)*Q1*((-1.0)*(136.0+25.0*pow(P2, 2))*pow(Q1, 2)+(408.0+47.0*pow(P2, 2))* 
    pow(Q2, 2))+pow(P1, 5)*(41.0*pow(Q1, 2)*Q2+(-23.0)*pow(Q2, 3)))+(-3.0)*(P1*pow(P2, 4)*Q2* 
    ((-829.0)*pow(Q1, 4)+34.0*pow(Q1, 2)*pow(Q2, 2)+47.0*pow(Q2, 4))+2.0*P1*pow(P2, 2)*Q2*((-1.0) 
    *(3156.0+223.0*pow(P1, 2))*pow(Q1, 4)+6.0*(236.0+73.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+3.0*( 
    116.0+23.0*pow(P1, 2))*pow(Q2, 4))+P1*Q2*(((-1200.0)+(-296.0)*pow(P1, 2)+15.0*pow(P1, 4)) 
    *pow(Q1, 4)+2.0*(1200.0+1928.0*pow(P1, 2)+149.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(240.0+ 
    712.0*pow(P1, 2)+85.0*pow(P1, 4))*pow(Q2, 4))+P2*Q1*(((-240.0)+1032.0*pow(P1, 2)+59.0* 
    pow(P1, 4))*pow(Q1, 4)+(-2.0)*((-1200.0)+264.0*pow(P1, 2)+163.0*pow(P1, 4))*pow(Q1, 2)*pow(Q2, 2)+ 
    (-1.0)*(1200.0+4632.0*pow(P1, 2)+529.0*pow(P1, 4))*pow(Q2, 4))+pow(P2, 5)*((-105.0)*pow(Q1, 5)+ 
    450.0*pow(Q1, 3)*pow(Q2, 2)+(-37.0)*Q1*pow(Q2, 4))+pow(P2, 3)*(((-824.0)+226.0*pow(P1, 2))* 
    pow(Q1, 5)+4.0*(1244.0+119.0*pow(P1, 2))*pow(Q1, 3)*pow(Q2, 2)+(-2.0)*(428.0+243.0*pow(P1, 2))* 
    Q1*pow(Q2, 4))))
    ) + 
    (c8L - c8L0) * ( 
    (105.0/32.0)*pow(G, -4.0)*(2.0*pow(G, 2)*((-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*(pow(Q1, 2)+(-3.0) 
    *pow(Q2, 2))+pow(P2, 4)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+4.0*pow(P1, 3)*P2*Q2*((-3.0)* 
    pow(Q1, 2)+pow(Q2, 2))+(-4.0)*P1*pow(P2, 3)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 4)*( 
    pow(Q1, 3)+(-3.0)*Q1*pow(Q2, 2)))+(-3.0)*(24.0*pow(P1, 4)*Q1*pow(Q2, 2)*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+pow(P1, 3)*P2*Q2*((-11.0)*pow(Q1, 4)+(-74.0)*pow(Q1, 2)*pow(Q2, 2)+17.0*pow(Q2, 4)) 
    +P1*P2*Q2*((70.0+81.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*(70.0+33.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+(14.0+(-3.0)*pow(P2, 2))*pow(Q2, 4))+pow(P2, 2)*Q1*(7.0*(1.0+pow(P2, 2))*pow(Q1, 4)+( 
    -2.0)*(35.0+23.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(35.0+11.0*pow(P2, 2))*pow(Q2, 4))+pow(P1, 2)*(( 
    -7.0)*(1.0+3.0*pow(P2, 2))*pow(Q1, 5)+2.0*(35.0+33.0*pow(P2, 2))*pow(Q1, 3)*pow(Q2, 2)+((-35.0)+ 
    39.0*pow(P2, 2))*Q1*pow(Q2, 4))))
    ) + 
    (c9L - c9L0) * ( 
    (7.0/128.0)*pow(G, -4.0)*((-75.0)*pow(P1, 4)*P2*Q1*(pow(Q1, 4)+14.0*pow(Q1, 2)*pow(Q2, 2)+( 
    -19.0)*pow(Q2, 4))+10.0*pow(P1, 3)*Q2*(5.0*(52.0+23.0*pow(P2, 2))*pow(Q1, 4)+130.0*((-4.0)+ 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(52.0+(-49.0)*pow(P2, 2))*pow(Q2, 4))+10.0*pow(P1, 2)*P2*Q1* 
    ((156.0+67.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(156.0+31.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*( 
    156.0+(-5.0)*pow(P2, 2))*pow(Q2, 4))+(-15.0)*P1*pow(P2, 2)*Q2*(5.0*(104.0+25.0*pow(P2, 2)) 
    *pow(Q1, 4)+(-130.0)*(8.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(104.0+pow(P2, 2))*pow(Q2, 4))+(-1.0)* 
    pow(P2, 3)*Q1*((520.0+119.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(520.0+83.0*pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+5.0*(520.0+47.0*pow(P2, 2))*pow(Q2, 4))+pow(P1, 5)*(145.0*pow(Q1, 4)*Q2+( 
    -650.0)*pow(Q1, 2)*pow(Q2, 3)+101.0*pow(Q2, 5))+10.0*pow(G, 2)*(5.0*pow(P1, 4)*P2*Q1*( 
    pow(Q1, 2)+(-3.0)*pow(Q2, 2))+(-10.0)*pow(P1, 2)*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+ 
    pow(P2, 5)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+10.0*pow(P1, 3)*pow(P2, 2)*Q2*((-3.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+(-5.0)*P1*pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 5)*(3.0*pow(Q1, 2)* 
    Q2+(-1.0)*pow(Q2, 3))))
    ) + 
    (c10L - c10L0) * ( 
    (-189.0/32.0)*pow(G, -4.0)*((-4.0)*pow(P1, 3)*P2*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+4.0*P1*pow(P2, 3)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+pow(P2, 4)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 4)*( 
    pow(Q1, 5)+(-10.0)*pow(Q1, 3)*pow(Q2, 2)+5.0*Q1*pow(Q2, 4)))
    ) + 
    (c11L - c11L0) * ( 
    (-63.0/128.0)*pow(G, -4.0)*((-10.0)*pow(P1, 3)*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+(-10.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+5.0*P1*pow(P2, 4)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+5.0*pow(P1, 4)*P2*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+(-10.0)*pow(P1, 2)*pow(P2, 3)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0* 
    pow(Q2, 4))+pow(P2, 5)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 5)*(5.0* 
    pow(Q1, 4)*Q2+(-10.0)*pow(Q1, 2)*pow(Q2, 3)+pow(Q2, 5)))
    ) + 
    (sL - sL0) * ( 
    (15.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(35.0*pow(P1, 5)*Q1+2.0*pow(P1, 3)*(296.0+121.0* 
    pow(P2, 2))*Q1+3.0*P1*(208.0+656.0*pow(P2, 2)+69.0*pow(P2, 4))*Q1+(-127.0)*pow(P1, 4)* 
    P2*Q2+(-2.0)*pow(P1, 2)*P2*(504.0+41.0*pow(P2, 2))*Q2+P2*(16.0+368.0*pow(P2, 2)+ 
    45.0*pow(P2, 4))*Q2)+(-28.0)*pow(G, 2)*((-1.0)*pow(P1, 4)*P2*Q2*(119.0*pow(Q1, 2)+ 
    45.0*pow(Q2, 2))+(-2.0)*pow(P1, 2)*P2*Q2*(7.0*(90.0+13.0*pow(P2, 2))*pow(Q1, 2)+(-3.0)*(( 
    -42.0)+pow(P2, 2))*pow(Q2, 2))+2.0*pow(P1, 3)*Q1*(3.0*(34.0+19.0*pow(P2, 2))*pow(Q1, 2)+(286.0+ 
    71.0*pow(P2, 2))*pow(Q2, 2))+P1*Q1*((304.0+1212.0*pow(P2, 2)+145.0*pow(P2, 4))*pow(Q1, 2)+ 
    3.0*(112.0+100.0*pow(P2, 2)+(-7.0)*pow(P2, 4))*pow(Q2, 2))+P2*Q2*(((-16.0)+388.0* 
    pow(P2, 2)+57.0*pow(P2, 4))*pow(Q1, 2)+(16.0+116.0*pow(P2, 2)+11.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 5)*( 
    9.0*pow(Q1, 3)+43.0*Q1*pow(Q2, 2)))+42.0*((-1.0)*pow(P1, 4)*P2*Q2*(251.0*pow(Q1, 4)+ 
    450.0*pow(Q1, 2)*pow(Q2, 2)+63.0*pow(Q2, 4))+(-2.0)*pow(P1, 2)*P2*Q2*((1668.0+305.0* 
    pow(P2, 2))*pow(Q1, 4)+2.0*(852.0+59.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(132.0+(-19.0)*pow(P2, 2)) 
    *pow(Q2, 4))+2.0*pow(P1, 3)*Q1*((156.0+101.0*pow(P2, 2))*pow(Q1, 4)+2.0*(444.0+179.0* 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+35.0*(20.0+3.0*pow(P2, 2))*pow(Q2, 4))+P1*Q1*((600.0+ 
    2664.0*pow(P2, 2)+343.0*pow(P2, 4))*pow(Q1, 4)+2.0*(648.0+1224.0*pow(P2, 2)+25.0*pow(P2, 4))* 
    pow(Q1, 2)*pow(Q2, 2)+(696.0+(-24.0)*pow(P2, 2)+(-109.0)*pow(P2, 4))*pow(Q2, 4))+P2*Q2*((( 
    -56.0)+1000.0*pow(P2, 2)+161.0*pow(P2, 4))*pow(Q1, 4)+2.0*((-8.0)+552.0*pow(P2, 2)+67.0*pow(P2, 4)) 
    *pow(Q1, 2)*pow(Q2, 2)+(40.0+168.0*pow(P2, 2)+13.0*pow(P2, 4))*pow(Q2, 4))+pow(P1, 5)*(11.0* 
    pow(Q1, 5)+106.0*pow(Q1, 3)*pow(Q2, 2)+119.0*Q1*pow(Q2, 4)))) 
    ) + 
    (s2L - s2L0) * ( 
    (15.0/32.0)*pow(G, -4.0)*(pow(G, 4)*((-52.0)*pow(P1, 3)*P2*Q1+4.0*P1*P2*((-4.0)+ 
    9.0*pow(P2, 2))*Q1+111.0*pow(P1, 4)*Q2+6.0*pow(P1, 2)*(40.0+9.0*pow(P2, 2))*Q2+(32.0+80.0* 
    pow(P2, 2)+31.0*pow(P2, 4))*Q2)+(-28.0)*pow(G, 2)*((-2.0)*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+ 
    49.0*pow(Q2, 2))+2.0*P1*P2*Q1*((18.0+19.0*pow(P2, 2))*pow(Q1, 2)+(-7.0)*(10.0+3.0* 
    pow(P2, 2))*pow(Q2, 2))+2.0*pow(P1, 2)*Q2*(9.0*(3.0+2.0*pow(P2, 2))*pow(Q1, 2)+(71.0+12.0* 
    pow(P2, 2))*pow(Q2, 2))+Q2*((16.0+106.0*pow(P2, 2)+47.0*pow(P2, 4))*pow(Q1, 2)+(16.0+18.0* 
    pow(P2, 2)+5.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 4)*(21.0*pow(Q1, 2)*Q2+67.0*pow(Q2, 3)))+42.0*(4.0* 
    pow(P1, 3)*P2*Q1*(3.0*pow(Q1, 4)+(-42.0)*pow(Q1, 2)*pow(Q2, 2)+(-77.0)*pow(Q2, 4))+4.0*P1* 
    P2*Q1*(3.0*(10.0+9.0*pow(P2, 2))*pow(Q1, 4)+(-42.0)*(2.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    -7.0)*(14.0+3.0*pow(P2, 2))*pow(Q2, 4))+Q2*((32.0+288.0*pow(P2, 2)+143.0*pow(P2, 4))*pow(Q1, 4)+ 
    2.0*(32.0+136.0*pow(P2, 2)+45.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(32.0+16.0*pow(P2, 2)+3.0*pow(P2, 4)) 
    *pow(Q2, 4))+3.0*pow(P1, 4)*(5.0*pow(Q1, 4)*Q2+46.0*pow(Q1, 2)*pow(Q2, 3)+49.0*pow(Q2, 5))+ 
    pow(P1, 2)*(2.0*(16.0+3.0*pow(P2, 2))*pow(Q1, 4)*Q2+92.0*(4.0+3.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 3)+2.0*(152.0+15.0*pow(P2, 2))*pow(Q2, 5)))) 
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(63.0*pow(P1, 5)*Q1+6.0*pow(P1, 3)*(148.0+43.0*pow(P2, 2)) 
    *Q1+3.0*P1*(208.0+360.0*pow(P2, 2)+17.0*pow(P2, 4))*Q1+(-243.0)*pow(P1, 4)*P2*Q2+( 
    -2.0)*pow(P1, 2)*P2*(948.0+73.0*pow(P2, 2))*Q2+(-1.0)*P2*(624.0+616.0*pow(P2, 2)+47.0* 
    pow(P2, 4))*Q2)+(-56.0)*pow(G, 2)*((-3.0)*pow(P1, 4)*P2*Q2*(30.0*pow(Q1, 2)+17.0* 
    pow(Q2, 2))+(-2.0)*pow(P1, 2)*P2*Q2*(5.0*(87.0+11.0*pow(P2, 2))*pow(Q1, 2)+3.0*(57.0+2.0* 
    pow(P2, 2))*pow(Q2, 2))+6.0*pow(P1, 3)*Q1*((17.0+7.0*pow(P2, 2))*pow(Q1, 2)+(97.0+22.0*pow(P2, 2)) 
    *pow(Q2, 2))+(-1.0)*P2*Q2*((420.0+550.0*pow(P2, 2)+44.0*pow(P2, 4))*pow(Q1, 2)+(68.0+22.0* 
    pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+P1*Q1*(2.0*(38.0+75.0*pow(P2, 2)+2.0*pow(P2, 4))*pow(Q1, 2)+ 
    3.0*(132.0+210.0*pow(P2, 2)+13.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 5)*(6.0*pow(Q1, 3)+45.0*Q1* 
    pow(Q2, 2)))+42.0*((-3.0)*pow(P1, 4)*P2*Q2*(97.0*pow(Q1, 4)+286.0*pow(Q1, 2)*pow(Q2, 2)+ 
    53.0*pow(Q2, 4))+(-2.0)*pow(P1, 2)*P2*Q2*(9.0*(184.0+29.0*pow(P2, 2))*pow(Q1, 4)+2.0*( 
    1824.0+179.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(456.0+(-7.0)*pow(P2, 2))*pow(Q2, 4))+2.0* 
    pow(P1, 3)*Q1*((104.0+49.0*pow(P2, 2))*pow(Q1, 4)+2.0*(704.0+259.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+(1624.0+269.0*pow(P2, 2))*pow(Q2, 4))+P2*Q2*((-3.0)*(776.0+1184.0*pow(P2, 2)+ 
    101.0*pow(P2, 4))*pow(Q1, 4)+(-2.0)*(1032.0+848.0*pow(P2, 2)+49.0*pow(P2, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+((-120.0)+64.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q2, 4))+3.0*P1*Q1*((40.0+32.0* 
    pow(P2, 2)+(-11.0)*pow(P2, 4))*pow(Q1, 4)+2.0*(408.0+1040.0*pow(P2, 2)+87.0*pow(P2, 4))*pow(Q1, 2)* 
    pow(Q2, 2)+(648.0+640.0*pow(P2, 2)+17.0*pow(P2, 4))*pow(Q2, 4))+pow(P1, 5)*(11.0*pow(Q1, 5)+178.0* 
    pow(Q1, 3)*pow(Q2, 2)+271.0*Q1*pow(Q2, 4))))
    ) + 
    (s4L - s4L0) * ( 
    (-3.0/16.0)*pow(G, -4.0)*(5.0*pow(G, 4)*(17.0*pow(P1, 3)*P2*Q1+P1*P2*(38.0+21.0* 
    pow(P2, 2))*Q1+12.0*pow(P1, 4)*Q2+pow(P1, 2)*(19.0+(-15.0)*pow(P2, 2))*Q2+(-1.0)* 
    pow(P2, 2)*(19.0+7.0*pow(P2, 2))*Q2)+70.0*pow(G, 2)*((-2.0)*pow(P1, 3)*P2*Q1*(9.0* 
    pow(Q1, 2)+7.0*pow(Q2, 2))+(-2.0)*P1*P2*Q1*((20.0+11.0*pow(P2, 2))*pow(Q1, 2)+(16.0+9.0* 
    pow(P2, 2))*pow(Q2, 2))+Q2*((12.0+82.0*pow(P2, 2)+31.0*pow(P2, 4))*pow(Q1, 2)+(-1.0)*(4.0+2.0* 
    pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+pow(P1, 2)*(2.0*(19.0+30.0*pow(P2, 2))*pow(Q1, 2)*Q2+(-38.0)* 
    pow(Q2, 3))+pow(P1, 4)*(9.0*pow(Q1, 2)*Q2+(-19.0)*pow(Q2, 3)))+(-84.0)*((-5.0)*pow(P1, 3)* 
    P2*Q1*(7.0*pow(Q1, 4)+38.0*pow(Q1, 2)*pow(Q2, 2)+(-5.0)*pow(Q2, 4))+30.0*pow(P1, 4)*Q2*( 
    pow(Q1, 4)+pow(Q1, 2)*pow(Q2, 2)+(-2.0)*pow(Q2, 4))+(-5.0)*P1*P2*Q1*((16.0+9.0*pow(P2, 2))* 
    pow(Q1, 4)+2.0*(40.0+21.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(8.0+3.0*pow(P2, 2))*pow(Q2, 4))+ 
    5.0*pow(P1, 2)*Q2*(9.0*(3.0+5.0*pow(P2, 2))*pow(Q1, 4)+2.0*(11.0+15.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+(-1.0)*(25.0+3.0*pow(P2, 2))*pow(Q2, 4))+Q2*(5.0*(10.0+73.0*pow(P2, 2)+29.0*pow(P2, 4)) 
    *pow(Q1, 4)+10.0*(2.0+9.0*pow(P2, 2)+2.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(14.0+15.0* 
    pow(P2, 2)+5.0*pow(P2, 4))*pow(Q2, 4))))
    ) + 
    (s5L - s5L0) * ( 
    (3.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(35.0*pow(P1, 5)*Q1+pow(P1, 3)*(296.0+(-54.0)*pow(P2, 2)) 
    *Q1+(-1.0)*P1*pow(P2, 2)*(888.0+121.0*pow(P2, 2))*Q1+(-151.0)*pow(P1, 4)*P2*Q2+ 
    6.0*pow(P1, 2)*P2*((-148.0)+pow(P2, 2))*Q2+pow(P2, 3)*(296.0+29.0*pow(P2, 2))*Q2)+(-56.0) 
    *pow(G, 2)*(35.0*pow(P1, 5)*Q1*pow(Q2, 2)+(-1.0)*pow(P1, 4)*P2*Q2*(16.0*pow(Q1, 2)+45.0* 
    pow(Q2, 2))+2.0*pow(P1, 2)*P2*Q2*(3.0*(11.0+9.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(159.0+8.0* 
    pow(P2, 2))*pow(Q2, 2))+P2*Q2*((228.0+434.0*pow(P2, 2)+38.0*pow(P2, 4))*pow(Q1, 2)+(-1.0)*( 
    76.0+46.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q2, 2))+P1*Q1*((-2.0)*(38.0+177.0*pow(P2, 2)+21.0* 
    pow(P2, 4))*pow(Q1, 2)+(228.0+174.0*pow(P2, 2)+5.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 3)*((-34.0)*(1.0+ 
    pow(P2, 2))*pow(Q1, 3)+2.0*(199.0+24.0*pow(P2, 2))*Q1*pow(Q2, 2)))+(-21.0)*(pow(P1, 4)*P2* 
    Q2*((-179.0)*pow(Q1, 4)+870.0*pow(Q1, 2)*pow(Q2, 2)+345.0*pow(Q2, 4))+2.0*pow(P1, 2)*P2* 
    Q2*((-1.0)*(2160.0+541.0*pow(P2, 2))*pow(Q1, 4)+2.0*(1632.0+109.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+5.0*(240.0+11.0*pow(P2, 2))*pow(Q2, 4))+P2*Q2*((-1.0)*(4496.0+7552.0*pow(P2, 2)+ 
    647.0*pow(P2, 4))*pow(Q1, 4)+2.0*(848.0+608.0*pow(P2, 2)+39.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+( 
    560.0+320.0*pow(P2, 2)+21.0*pow(P2, 4))*pow(Q2, 4))+P1*Q1*((720.0+3072.0*pow(P2, 2)+359.0* 
    pow(P2, 4))*pow(Q1, 4)+2.0*(48.0+1632.0*pow(P2, 2)+221.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*( 
    3696.0+4416.0*pow(P2, 2)+301.0*pow(P2, 4))*pow(Q2, 4))+pow(P1, 5)*(11.0*pow(Q1, 5)+(-110.0)* 
    pow(Q1, 3)*pow(Q2, 2)+(-505.0)*Q1*pow(Q2, 4))+pow(P1, 3)*((416.0+306.0*pow(P2, 2))*pow(Q1, 5)+ 
    4.0*((-224.0)+51.0*pow(P2, 2))*pow(Q1, 3)*pow(Q2, 2)+(-10.0)*(592.0+87.0*pow(P2, 2))*Q1* 
    pow(Q2, 4)))) 
    ) + 
    (s6L - s6L0) * ( 
    (1.0/32.0)*pow(G, -4.0)*(45.0*pow(G, 4)*(4.0*pow(P1, 3)*P2*Q1+(-4.0)*P1*pow(P2, 3)*Q1+ 
    pow(P1, 4)*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q2+pow(P2, 4)*Q2)+140.0*pow(G, 2)*((-2.0)* 
    pow(P1, 3)*P2*Q1*(pow(Q1, 2)+33.0*pow(Q2, 2))+6.0*pow(P1, 2)*Q2*((9.0+6.0*pow(P2, 2))* 
    pow(Q1, 2)+((-3.0)+4.0*pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)*Q2*((-3.0)*(18.0+11.0*pow(P2, 2))* 
    pow(Q1, 2)+(18.0+5.0*pow(P2, 2))*pow(Q2, 2))+2.0*P1*P2*Q1*((18.0+19.0*pow(P2, 2))*pow(Q1, 2)+ 
    (-3.0)*(18.0+7.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(21.0*pow(Q1, 2)*Q2+(-13.0)*pow(Q2, 3)))+( 
    -21.0)*(60.0*pow(P1, 3)*P2*Q1*(3.0*pow(Q1, 4)+(-38.0)*pow(Q1, 2)*pow(Q2, 2)+(-25.0)* 
    pow(Q2, 4))+60.0*P1*P2*Q1*((20.0+17.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*(28.0+9.0*pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(44.0+19.0*pow(P2, 2))*pow(Q2, 4))+10.0*pow(P1, 2)*Q2*((32.0+( 
    -39.0)*pow(P2, 2))*pow(Q1, 4)+2.0*(184.0+183.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-80.0)+21.0* 
    pow(P2, 2))*pow(Q2, 4))+Q2*((-5.0)*(64.0+704.0*pow(P2, 2)+339.0*pow(P2, 4))*pow(Q1, 4)+10.0*( 
    64.0+272.0*pow(P2, 2)+75.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+((-64.0)+160.0*pow(P2, 2)+45.0* 
    pow(P2, 4))*pow(Q2, 4))+15.0*pow(P1, 4)*(15.0*pow(Q1, 4)*Q2+82.0*pow(Q1, 2)*pow(Q2, 3)+(-29.0)* 
    pow(Q2, 5))))
    ) + 
    (s7L - s7L0) * ( 
    (-15.0/128.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 5)*Q1+(-10.0)*pow(P1, 3)*pow(P2, 2)*Q1+5.0* 
    P1*pow(P2, 4)*Q1+(-5.0)*pow(P1, 4)*P2*Q2+10.0*pow(P1, 2)*pow(P2, 3)*Q2+(-1.0)* 
    pow(P2, 5)*Q2)+2.0*pow(G, 2)*(pow(P2, 3)*Q2*((408.0+61.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*( 
    136.0+11.0*pow(P2, 2))*pow(Q2, 2))+(-2.0)*pow(P1, 2)*P2*Q2*((612.0+101.0*pow(P2, 2))* 
    pow(Q1, 2)+((-204.0)+13.0*pow(P2, 2))*pow(Q2, 2))+2.0*pow(P1, 3)*Q1*((68.0+23.0*pow(P2, 2))* 
    pow(Q1, 2)+((-204.0)+71.0*pow(P2, 2))*pow(Q2, 2))+P1*pow(P2, 2)*Q1*((-1.0)*(408.0+91.0* 
    pow(P2, 2))*pow(Q1, 2)+(1224.0+133.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 5)*(9.0*pow(Q1, 3)+(-55.0)* 
    Q1*pow(Q2, 2))+pow(P1, 4)*P2*((-103.0)*pow(Q1, 2)*Q2+81.0*pow(Q2, 3)))+(-3.0)*( 
    pow(P1, 4)*P2*Q2*((-299.0)*pow(Q1, 4)+(-226.0)*pow(Q1, 2)*pow(Q2, 2)+217.0*pow(Q2, 4))+( 
    -2.0)*pow(P1, 3)*Q1*(((-52.0)+3.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*(556.0+291.0*pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+7.0*(196.0+pow(P2, 2))*pow(Q2, 4))+(-2.0)*pow(P1, 2)*P2*Q2*(3.0*(452.0+ 
    51.0*pow(P2, 2))*pow(Q1, 4)+2.0*(1092.0+251.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(708.0+ 
    19.0*pow(P2, 2))*pow(Q2, 4))+P2*Q2*((1200.0+3304.0*pow(P2, 2)+361.0*pow(P2, 4))*pow(Q1, 4)+( 
    -2.0)*(1200.0+1672.0*pow(P2, 2)+117.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+(240.0+8.0*pow(P2, 2)+( 
    -3.0)*pow(P2, 4))*pow(Q2, 4))+P1*Q1*((-1.0)*(240.0+1752.0*pow(P2, 2)+289.0*pow(P2, 4))* 
    pow(Q1, 4)+2.0*(1200.0+3864.0*pow(P2, 2)+353.0*pow(P2, 4))*pow(Q1, 2)*pow(Q2, 2)+((-1200.0)+ 
    1032.0*pow(P2, 2)+179.0*pow(P2, 4))*pow(Q2, 4))+pow(P1, 5)*(11.0*pow(Q1, 5)+106.0*pow(Q1, 3)* 
    pow(Q2, 2)+(-273.0)*Q1*pow(Q2, 4))))
    ) + 
    (s8L - s8L0) * ( 
    (-105.0/32.0)*pow(G, -4.0)*(2.0*pow(G, 2)*(4.0*pow(P1, 3)*P2*Q1*(pow(Q1, 2)+(-3.0)* 
    pow(Q2, 2))+(-4.0)*P1*pow(P2, 3)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)* 
    Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-1.0)*pow(P2, 4)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+ 
    pow(P1, 4)*(3.0*pow(Q1, 2)*Q2+(-1.0)*pow(Q2, 3)))+3.0*(6.0*pow(P1, 4)*Q2*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+pow(P1, 3)*P2*Q1*((-7.0)*pow(Q1, 4)+(-26.0)*pow(Q1, 2)* 
    pow(Q2, 2)+61.0*pow(Q2, 4))+pow(P1, 2)*Q2*((35.0+69.0*pow(P2, 2))*pow(Q1, 4)+2.0*((-35.0)+3.0* 
    pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(7.0+(-15.0)*pow(P2, 2))*pow(Q2, 4))+(-1.0)*pow(P2, 2)*Q2*(( 
    35.0+29.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*(35.0+17.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(7.0+pow(P2, 2)) 
    *pow(Q2, 4))+P1*P2*Q1*(7.0*(2.0+3.0*pow(P2, 2))*pow(Q1, 4)+(-2.0)*(70.0+57.0*pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+(70.0+9.0*pow(P2, 2))*pow(Q2, 4)))) 
    ) + 
    (s9L - s9L0) * ( 
    (7.0/128.0)*pow(G, -4.0)*((-75.0)*pow(P1, 4)*P2*Q2*(pow(Q1, 4)+(-26.0)*pow(Q1, 2)* 
    pow(Q2, 2)+5.0*pow(Q2, 4))+10.0*pow(P1, 3)*Q1*((52.0+41.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(52.0+ 
    5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+5.0*(52.0+(-31.0)*pow(P2, 2))*pow(Q2, 4))+(-15.0)*P1* 
    pow(P2, 2)*Q1*((104.0+31.0*pow(P2, 2))*pow(Q1, 4)+(-10.0)*(104.0+19.0*pow(P2, 2))*pow(Q1, 2)* 
    pow(Q2, 2)+5.0*(104.0+7.0*pow(P2, 2))*pow(Q2, 4))+10.0*pow(P1, 2)*P2*Q2*((-5.0)*(156.0+ 
    49.0*pow(P2, 2))*pow(Q1, 4)+130.0*(12.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-156.0)+23.0*pow(P2, 2)) 
    *pow(Q2, 4))+pow(P2, 3)*Q2*(5.0*(520.0+101.0*pow(P2, 2))*pow(Q1, 4)+(-650.0)*(8.0+pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+(520.0+29.0*pow(P2, 2))*pow(Q2, 4))+pow(P1, 5)*(11.0*pow(Q1, 5)+(-470.0)* 
    pow(Q1, 3)*pow(Q2, 2)+415.0*Q1*pow(Q2, 4))+10.0*pow(G, 2)*((-10.0)*pow(P1, 3)*pow(P2, 2)*Q1*( 
    pow(Q1, 2)+(-3.0)*pow(Q2, 2))+5.0*P1*pow(P2, 4)*Q1*(pow(Q1, 2)+(-3.0)*pow(Q2, 2))+5.0* 
    pow(P1, 4)*P2*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+(-10.0)*pow(P1, 2)*pow(P2, 3)*Q2*((-3.0) 
    *pow(Q1, 2)+pow(Q2, 2))+pow(P2, 5)*Q2*((-3.0)*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 5)*(pow(Q1, 3)+(-3.0)* 
    Q1*pow(Q2, 2))))
    ) + 
    (s10L - s10L0) * ( 
    (189.0/32.0)*pow(G, -4.0)*((-6.0)*pow(P1, 2)*pow(P2, 2)*Q2*(5.0*pow(Q1, 4)+(-10.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+pow(P2, 4)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+4.0*pow(P1, 3)*P2*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+(-4.0) 
    *P1*pow(P2, 3)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 4)*(5.0* 
    pow(Q1, 4)*Q2+(-10.0)*pow(Q1, 2)*pow(Q2, 3)+pow(Q2, 5)))
    ) + 
    (s11L - s11L0) * ( 
    (-63.0/128.0)*pow(G, -4.0)*((-5.0)*pow(P1, 4)*P2*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+10.0*pow(P1, 2)*pow(P2, 3)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+(-1.0)*pow(P2, 5)*Q2*(5.0*pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-10.0) 
    *pow(P1, 3)*pow(P2, 2)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+5.0*P1* 
    pow(P2, 4)*Q1*(pow(Q1, 4)+(-10.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+pow(P1, 5)*(pow(Q1, 5)+( 
    -10.0)*pow(Q1, 3)*pow(Q2, 2)+5.0*Q1*pow(Q2, 4)))
    );


    // Q1
    I[4] = (Lf - L0) * ( 
    (15.0/2.0)*pow(G, -4.0)*(pow(G, 4)*P1*(4.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))+(-28.0)*pow(G, 2)*(( 
    -6.0)*pow(P1, 2)*P2*Q1*Q2+(-2.0)*P2*(2.0+pow(P2, 2))*Q1*Q2+pow(P1, 3)*(pow(Q1, 2)+ 
    5.0*pow(Q2, 2))+P1*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+3.0*(2.0+pow(P2, 2))*pow(Q2, 2)))+21.0*((-4.0) 
    *P2*Q1*Q2*((8.0+5.0*pow(P2, 2))*pow(Q1, 2)+(8.0+3.0*pow(P2, 2))*pow(Q2, 2))+(-12.0)* 
    pow(P1, 2)*P2*Q2*(3.0*pow(Q1, 3)+5.0*Q1*pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 4)+30.0* 
    pow(Q1, 2)*pow(Q2, 2)+35.0*pow(Q2, 4))+P1*((8.0+15.0*pow(P2, 2))*pow(Q1, 4)+6.0*(8.0+9.0*pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+5.0*(8.0+3.0*pow(P2, 2))*pow(Q2, 4))))  
    ) + 
    (cL - cL0) * ( 
    (-15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*(8.0+5.0*pow(P1, 4)+12.0*pow(P2, 2)+pow(P2, 4)+6.0*pow(P1, 2)*(6.0+ 
    pow(P2, 2)))+(-7.0)*pow(G, 2)*((16.0+5.0*pow(P1, 4)+48.0*pow(P2, 2)+5.0*pow(P2, 4)+6.0*pow(P1, 2)*(8.0+ 
    3.0*pow(P2, 2)))*pow(Q1, 2)+(-8.0)*P1*P2*(5.0*pow(P1, 2)+3.0*(8.0+pow(P2, 2)))*Q1*Q2+( 
    35.0*pow(P1, 4)+30.0*pow(P1, 2)*(8.0+pow(P2, 2))+3.0*(16.0+16.0*pow(P2, 2)+pow(P2, 4)))*pow(Q2, 2))+ 
    21.0*((16.0+3.0*pow(P1, 4)+60.0*pow(P2, 2)+7.0*pow(P2, 4)+18.0*pow(P1, 2)*(2.0+pow(P2, 2)))*pow(Q1, 4)+ 
    (-48.0)*P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))*pow(Q1, 3)*Q2+6.0*(16.0+7.0*pow(P1, 4)+36.0* 
    pow(P2, 2)+3.0*pow(P2, 4)+6.0*pow(P1, 2)*(10.0+3.0*pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+(-16.0)*P1* 
    P2*(30.0+7.0*pow(P1, 2)+3.0*pow(P2, 2))*Q1*pow(Q2, 3)+(80.0+63.0*pow(P1, 4)+60.0*pow(P2, 2)+3.0* 
    pow(P2, 4)+42.0*pow(P1, 2)*(10.0+pow(P2, 2)))*pow(Q2, 4)))
    ) + 
    (c2L - c2L0) * ( 
    (-15.0/2.0)*pow(G, -4.0)*(pow(G, 4)*P2*(2.0+3.0*pow(P1, 2)+pow(P2, 2))+(-7.0)*pow(G, 2)*((-2.0) 
    *P1*(8.0+5.0*pow(P1, 2))*Q1*Q2+(-18.0)*P1*pow(P2, 2)*Q1*Q2+pow(P2, 3)*(5.0* 
    pow(Q1, 2)+3.0*pow(Q2, 2))+P2*((8.0+9.0*pow(P1, 2))*pow(Q1, 2)+(8.0+15.0*pow(P1, 2))*pow(Q2, 2)))+ 
    21.0*((-36.0)*P1*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+pow(Q2, 2))+(-4.0)*P1*Q1*Q2*( 
    3.0*(2.0+pow(P1, 2))*pow(Q1, 2)+(10.0+7.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(7.0*pow(Q1, 4)+18.0* 
    pow(Q1, 2)*pow(Q2, 2)+3.0*pow(Q2, 4))+P2*((10.0+9.0*pow(P1, 2))*pow(Q1, 4)+18.0*(2.0+3.0*pow(P1, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+(10.0+21.0*pow(P1, 2))*pow(Q2, 4))))
    ) + 
    (c3L - c3L0) * ( 
    (5.0/16.0)*pow(G, -4.0)*(pow(G, 4)*(5.0*pow(P1, 4)+(-6.0)*pow(P1, 2)*((-4.0)+pow(P2, 2))+(-3.0)* 
    pow(P2, 2)*(8.0+pow(P2, 2)))+14.0*pow(G, 2)*((16.0+pow(P1, 4)+72.0*pow(P2, 2)+9.0*pow(P2, 4)+6.0* 
    pow(P1, 2)*(4.0+3.0*pow(P2, 2)))*pow(Q1, 2)+(-8.0)*P1*P2*(pow(P1, 2)+3.0*(4.0+pow(P2, 2)))* 
    Q1*Q2+((-16.0)+(-21.0)*pow(P1, 4)+24.0*pow(P2, 2)+3.0*pow(P2, 4)+6.0*pow(P1, 2)*((-20.0)+ 
    pow(P2, 2)))*pow(Q2, 2))+(-84.0)*((12.0+pow(P1, 4)+54.0*pow(P2, 2)+7.0*pow(P2, 4)+6.0*pow(P1, 2)*(3.0+ 
    2.0*pow(P2, 2)))*pow(Q1, 4)+(-16.0)*P1*P2*(9.0+pow(P1, 2)+2.0*pow(P2, 2))*pow(Q1, 3)*Q2+ 
    12.0*(2.0+9.0*pow(P2, 2)+pow(P2, 4)+3.0*pow(P1, 2)*(1.0+pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+(-16.0)* 
    P1*P2*(3.0+pow(P2, 2))*Q1*pow(Q2, 3)+((-20.0)+(-126.0)*pow(P1, 2)+(-21.0)*pow(P1, 4)+6.0* 
    pow(P2, 2)+pow(P2, 4))*pow(Q2, 4)))
    ) + 
    (c4L - c4L0) * ( 
    (15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*pow(P1, 2)*P2+(-1.0)*pow(P2, 3))+56.0*pow(G, 2)*( 
    pow(P2, 3)*pow(Q1, 2)+2.0*P1*(1.0+pow(P1, 2))*Q1*Q2+P2*(pow(Q1, 2)+(-1.0)*(1.0+3.0* 
    pow(P1, 2))*pow(Q2, 2)))+(-42.0)*((-12.0)*P1*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+P2*(pow(Q1, 2)+pow(Q2, 2))*((8.0+3.0*pow(P1, 2))*pow(Q1, 2)+(-1.0)*(8.0+21.0*pow(P1, 2)) 
    *pow(Q2, 2))+pow(P2, 3)*(7.0*pow(Q1, 4)+6.0*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*pow(Q2, 4))+4.0*P1*Q1* 
    Q2*(8.0*pow(Q2, 2)+pow(P1, 2)*(pow(Q1, 2)+7.0*pow(Q2, 2))))) 
    ) + 
    (c5L - c5L0) * ( 
    (-3.0/16.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+(-84.0)*( 
    ((-4.0)+6.0*pow(P1, 2)+pow(P1, 4)+(-30.0)*pow(P2, 2)+(-5.0)*pow(P2, 4))*pow(Q1, 4)+(-16.0)*P1*( 
    3.0+pow(P1, 2))*P2*pow(Q1, 3)*Q2+12.0*(2.0+pow(P1, 4)+3.0*pow(P2, 2)+3.0*pow(P1, 2)*(3.0+pow(P2, 2))) 
    *pow(Q1, 2)*pow(Q2, 2)+(-16.0)*P1*P2*(9.0+2.0*pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 3)+((-4.0)+ 
    (-9.0)*pow(P1, 4)+18.0*pow(P2, 2)+pow(P2, 4)+6.0*pow(P1, 2)*((-7.0)+2.0*pow(P2, 2)))*pow(Q2, 4))+ 
    14.0*pow(G, 2)*((-24.0)*pow(P1, 3)*P2*Q1*Q2+(-8.0)*P1*P2*(12.0+pow(P2, 2))*Q1* 
    Q2+pow(P1, 4)*(3.0*pow(Q1, 2)+(-7.0)*pow(Q2, 2))+pow(P2, 2)*((-1.0)*(24.0+5.0*pow(P2, 2))* 
    pow(Q1, 2)+(24.0+pow(P2, 2))*pow(Q2, 2))+6.0*pow(P1, 2)*((4.0+pow(P2, 2))*pow(Q1, 2)+((-4.0)+3.0* 
    pow(P2, 2))*pow(Q2, 2))))
    ) + 
    (c6L - c6L0) * ( 
    (-35.0/2.0)*pow(G, -4.0)*(pow(G, 2)*(2.0*pow(P1, 3)*Q1*Q2+(-6.0)*P1*pow(P2, 2)*Q1* 
    Q2+3.0*pow(P1, 2)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 3)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+ 
    3.0*(12.0*P1*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+pow(Q2, 2))+4.0*P1*Q1*Q2*((2.0+pow(P1, 2)) 
    *pow(Q1, 2)+(-1.0)*(2.0+3.0*pow(P1, 2))*pow(Q2, 2))+pow(P2, 3)*(3.0*pow(Q1, 4)+(-6.0)*pow(Q1, 2)* 
    pow(Q2, 2)+(-1.0)*pow(Q2, 4))+P2*((2.0+(-3.0)*pow(P1, 2))*pow(Q1, 4)+(-6.0)*(2.0+3.0*pow(P1, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+(2.0+9.0*pow(P1, 2))*pow(Q2, 4))))
    ) + 
    (c7L - c7L0) * ( 
    (15.0/16.0)*pow(G, -4.0)*(2.0*pow(G, 2)*((-8.0)*pow(P1, 3)*P2*Q1*Q2+8.0*P1* 
    pow(P2, 3)*Q1*Q2+pow(P1, 4)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+3.0*((-16.0)*pow(P1, 3)* 
    P2*Q1*Q2*(pow(Q1, 2)+(-5.0)*pow(Q2, 2))+(-16.0)*P1*P2*Q1*Q2*(3.0*(4.0+ 
    pow(P2, 2))*pow(Q1, 2)+((-12.0)+pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+(-30.0)*pow(Q1, 2)* 
    pow(Q2, 2)+9.0*pow(Q2, 4))+6.0*pow(P1, 2)*((4.0+3.0*pow(P2, 2))*pow(Q1, 4)+6.0*((-4.0)+pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(4.0+(-5.0)*pow(P2, 2))*pow(Q2, 4))+pow(P2, 2)*((-1.0)*(24.0+7.0*pow(P2, 2))* 
    pow(Q1, 4)+18.0*(8.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+((-24.0)+pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c8L - c8L0) * ( 
    (315.0/16.0)*pow(G, -4.0)*(4.0*pow(P1, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-12.0)* 
    P1*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*pow(P1, 2)*P2*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (c9L - c9L0) * ( 
    (-35.0/16.0)*pow(G, -4.0)*((-16.0)*pow(P1, 3)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    16.0*P1*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+pow(P2, 4)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (sL - sL0) * ( 
    (15.0/4.0)*pow(G, -4.0)*(2.0*pow(G, 4)*P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))+7.0*pow(G, 2)*(5.0* 
    pow(P1, 4)*Q1*Q2+6.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2))*Q1*Q2+(16.0+48.0*pow(P2, 2)+5.0* 
    pow(P2, 4))*Q1*Q2+(-2.0)*pow(P1, 3)*P2*(3.0*pow(Q1, 2)+5.0*pow(Q2, 2))+(-2.0)*P1*P2* 
    ((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+pow(P2, 2))*pow(Q2, 2)))+(-42.0)*(6.0*pow(P1, 2)*Q1* 
    Q2*(3.0*(2.0+pow(P2, 2))*pow(Q1, 2)+(10.0+3.0*pow(P2, 2))*pow(Q2, 2))+Q1*Q2*((16.0+60.0* 
    pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 2)+(16.0+36.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 4)*(3.0* 
    pow(Q1, 3)*Q2+7.0*Q1*pow(Q2, 3))+(-1.0)*pow(P1, 3)*P2*(3.0*pow(Q1, 4)+18.0*pow(Q1, 2)* 
    pow(Q2, 2)+7.0*pow(Q2, 4))+(-1.0)*P1*P2*((30.0+7.0*pow(P2, 2))*pow(Q1, 4)+18.0*(6.0+pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+3.0*(10.0+pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s2L - s2L0) * ( 
    (-15.0/2.0)*pow(G, -4.0)*(2.0*pow(G, 4)*(P1+pow(P1, 3))+42.0*((-12.0)*pow(P1, 2)*P2*Q1* 
    pow(Q2, 3)+4.0*P2*Q1*Q2*((1.0+pow(P2, 2))*pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*pow(P1, 3)*(3.0* 
    pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+P1*((-1.0)*(1.0+3.0*pow(P2, 2))*pow(Q1, 4)+6.0*pow(Q1, 2)* 
    pow(Q2, 2)+3.0*(5.0+pow(P2, 2))*pow(Q2, 4)))+(-7.0)*pow(G, 2)*((-6.0)*pow(P1, 2)*P2*Q1*Q2+ 
    2.0*pow(P2, 3)*Q1*Q2+pow(P1, 3)*(pow(Q1, 2)+15.0*pow(Q2, 2))+P1*(16.0*pow(Q2, 2)+(-3.0)* 
    pow(P2, 2)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2)))))
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/4.0)*pow(G, -4.0)*(pow(G, 4)*P1*P2*(12.0+3.0*pow(P1, 2)+pow(P2, 2))+7.0*pow(G, 2)*(9.0* 
    pow(P1, 4)*Q1*Q2+18.0*pow(P1, 2)*(4.0+pow(P2, 2))*Q1*Q2+(16.0+24.0*pow(P2, 2)+pow(P2, 4))* 
    Q1*Q2+(-6.0)*pow(P1, 3)*P2*(pow(Q1, 2)+3.0*pow(Q2, 2))+(-2.0)*P1*P2*(12.0+pow(P2, 2)) 
    *(pow(Q1, 2)+3.0*pow(Q2, 2)))+(-84.0)*(6.0*pow(P1, 2)*Q1*Q2*((3.0+pow(P2, 2))*pow(Q1, 2)+( 
    9.0+2.0*pow(P2, 2))*pow(Q2, 2))+Q1*Q2*((4.0+6.0*pow(P2, 2))*pow(Q1, 2)+(12.0+18.0*pow(P2, 2)+ 
    pow(P2, 4))*pow(Q2, 2))+pow(P1, 4)*(2.0*pow(Q1, 3)*Q2+7.0*Q1*pow(Q2, 3))+(-1.0)*pow(P1, 3)* 
    P2*(pow(Q1, 4)+12.0*pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+(-1.0)*P1*P2*(3.0*pow(Q1, 4)+6.0*( 
    9.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(27.0+2.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (s4L - s4L0) * ( 
    (15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 3)+(-3.0)*P1*pow(P2, 2))+28.0*pow(G, 2)*((-6.0)* 
    pow(P1, 2)*P2*Q1*Q2+(-2.0)*P2*(2.0+pow(P2, 2))*Q1*Q2+pow(P1, 3)*(pow(Q1, 2)+(-3.0)* 
    pow(Q2, 2))+P1*((2.0+3.0*pow(P2, 2))*pow(Q1, 2)+((-2.0)+3.0*pow(P2, 2))*pow(Q2, 2)))+(-84.0)*(( 
    -12.0)*pow(P1, 2)*P2*Q1*Q2*(pow(Q1, 2)+pow(Q2, 2))+(-4.0)*P2*(2.0+pow(P2, 2))*Q1* 
    Q2*(pow(Q1, 2)+pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 4)+6.0*pow(Q1, 2)*pow(Q2, 2)+(-7.0)*pow(Q2, 4))+P1*( 
    (2.0+3.0*pow(P2, 2))*pow(Q1, 4)+6.0*(2.0+3.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+3.0*((-2.0)+pow(P2, 2)) 
    *pow(Q2, 4))))
    ) + 
    (s5L - s5L0) * ( 
    (3.0/4.0)*pow(G, -4.0)*(pow(G, 4)*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2))+7.0*pow(G, 2)*(5.0* 
    pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*((-4.0)+pow(P2, 2))*Q1*Q2+(-3.0)*pow(P2, 2)*(8.0+ 
    pow(P2, 2))*Q1*Q2+2.0*pow(P1, 3)*P2*(pow(Q1, 2)+(-5.0)*pow(Q2, 2))+2.0*P1*P2*(3.0*( 
    4.0+pow(P2, 2))*pow(Q1, 2)+((-12.0)+pow(P2, 2))*pow(Q2, 2)))+(-84.0)*(5.0*pow(P1, 4)*Q1* 
    pow(Q2, 3)+(-6.0)*pow(P1, 2)*Q1*Q2*((1.0+pow(P2, 2))*pow(Q1, 2)+(-5.0)*pow(Q2, 2))+(-1.0)* 
    Q1*Q2*(2.0*(2.0+9.0*pow(P2, 2)+pow(P2, 4))*pow(Q1, 2)+((-4.0)+6.0*pow(P2, 2)+pow(P2, 4))* 
    pow(Q2, 2))+P1*P2*((9.0+2.0*pow(P2, 2))*pow(Q1, 4)+6.0*(3.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    -15.0)*pow(Q2, 4))+pow(P1, 3)*P2*(pow(Q1, 4)+(-5.0)*pow(Q2, 4)))) 
    ) + 
    (s6L - s6L0) * ( 
    (-35.0/2.0)*pow(G, -4.0)*(pow(G, 2)*((-6.0)*pow(P1, 2)*P2*Q1*Q2+2.0*pow(P2, 3)*Q1* 
    Q2+pow(P1, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+( 
    -6.0)*((-12.0)*pow(P1, 2)*P2*Q1*pow(Q2, 3)+4.0*P2*Q1*Q2*((1.0+pow(P2, 2))* 
    pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 3)*(6.0*pow(Q1, 2)*pow(Q2, 2)+(-2.0)*pow(Q2, 4))+(-1.0)*P1* 
    ((1.0+3.0*pow(P2, 2))*pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+(1.0+(-3.0)*pow(P2, 2))*pow(Q2, 4)))) 
    ) + 
    (s7L - s7L0) * ( 
    (-15.0/4.0)*pow(G, -4.0)*(pow(G, 2)*(pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1* 
    Q2+pow(P2, 4)*Q1*Q2+2.0*pow(P1, 3)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*P1*pow(P2, 3)*( 
    (-1.0)*pow(Q1, 2)+pow(Q2, 2)))+3.0*(pow(P2, 2)*Q1*Q2*((-1.0)*(24.0+5.0*pow(P2, 2))* 
    pow(Q1, 2)+(24.0+pow(P2, 2))*pow(Q2, 2))+6.0*pow(P1, 2)*Q1*Q2*((4.0+pow(P2, 2))*pow(Q1, 2)+((-4.0) 
    +3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(3.0*pow(Q1, 3)*Q2+(-7.0)*Q1*pow(Q2, 3))+(-1.0)* 
    pow(P1, 3)*P2*(pow(Q1, 4)+18.0*pow(Q1, 2)*pow(Q2, 2)+(-7.0)*pow(Q2, 4))+P1*P2*((12.0+5.0* 
    pow(P2, 2))*pow(Q1, 4)+(-6.0)*(12.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-3.0)*((-4.0)+pow(P2, 2))* 
    pow(Q2, 4)))) 
    ) + 
    (s8L - s8L0) * ( 
    (315.0/16.0)*pow(G, -4.0)*((-12.0)*pow(P1, 2)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    4.0*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+(-3.0)*P1*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))) 
    ) + 
    (s9L - s9L0) * ( 
    (35.0/4.0)*pow(G, -4.0)*(pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*Q1* 
    Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*Q1*Q2*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*P2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*P1* 
    pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) ;


    // Q2
    I[5] = (Lf - L0) * ( 
    (15.0/2.0)*pow(G, -4.0)*(pow(G, 4)*P2*(4.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))+(-28.0)*pow(G, 2)*(( 
    -2.0)*P1*(2.0+pow(P1, 2))*Q1*Q2+(-6.0)*P1*pow(P2, 2)*Q1*Q2+pow(P2, 3)*(5.0* 
    pow(Q1, 2)+pow(Q2, 2))+P2*(3.0*(2.0+pow(P1, 2))*pow(Q1, 2)+(2.0+3.0*pow(P1, 2))*pow(Q2, 2)))+(-21.0) 
    *(4.0*P1*Q1*Q2*((8.0+3.0*pow(P1, 2))*pow(Q1, 2)+(8.0+5.0*pow(P1, 2))*pow(Q2, 2))+12.0* 
    P1*pow(P2, 2)*Q2*(5.0*pow(Q1, 3)+3.0*Q1*pow(Q2, 2))+(-1.0)*pow(P2, 3)*(35.0*pow(Q1, 4)+ 
    30.0*pow(Q1, 2)*pow(Q2, 2)+3.0*pow(Q2, 4))+(-1.0)*P2*(5.0*(8.0+3.0*pow(P1, 2))*pow(Q1, 4)+6.0*( 
    8.0+9.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+(8.0+15.0*pow(P1, 2))*pow(Q2, 4)))) 
    ) + 
    (cL - cL0) * ( 
    (-15.0/4.0)*pow(G, -4.0)*(2.0*pow(G, 4)*P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))+7.0*pow(G, 2)*(5.0* 
    pow(P1, 4)*Q1*Q2+6.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2))*Q1*Q2+(16.0+48.0*pow(P2, 2)+5.0* 
    pow(P2, 4))*Q1*Q2+(-2.0)*pow(P1, 3)*P2*(3.0*pow(Q1, 2)+5.0*pow(Q2, 2))+(-2.0)*P1*P2* 
    ((24.0+5.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+pow(P2, 2))*pow(Q2, 2)))+(-42.0)*(6.0*pow(P1, 2)*Q1* 
    Q2*(3.0*(2.0+pow(P2, 2))*pow(Q1, 2)+(10.0+3.0*pow(P2, 2))*pow(Q2, 2))+Q1*Q2*((16.0+60.0* 
    pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 2)+(16.0+36.0*pow(P2, 2)+3.0*pow(P2, 4))*pow(Q2, 2))+pow(P1, 4)*(3.0* 
    pow(Q1, 3)*Q2+7.0*Q1*pow(Q2, 3))+(-1.0)*pow(P1, 3)*P2*(3.0*pow(Q1, 4)+18.0*pow(Q1, 2)* 
    pow(Q2, 2)+7.0*pow(Q2, 4))+(-1.0)*P1*P2*((30.0+7.0*pow(P2, 2))*pow(Q1, 4)+18.0*(6.0+pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+3.0*(10.0+pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c2L - c2L0) * ( 
    (-15.0/2.0)*pow(G, -4.0)*(pow(G, 4)*P1*(2.0+pow(P1, 2)+3.0*pow(P2, 2))+(-7.0)*pow(G, 2)*((-18.0) 
    *pow(P1, 2)*P2*Q1*Q2+(-2.0)*P2*(8.0+5.0*pow(P2, 2))*Q1*Q2+pow(P1, 3)*(3.0* 
    pow(Q1, 2)+5.0*pow(Q2, 2))+P1*((8.0+15.0*pow(P2, 2))*pow(Q1, 2)+(8.0+9.0*pow(P2, 2))*pow(Q2, 2)))+ 
    21.0*((-36.0)*pow(P1, 2)*P2*Q1*Q2*(pow(Q1, 2)+pow(Q2, 2))+(-4.0)*P2*Q1*Q2*(( 
    10.0+7.0*pow(P2, 2))*pow(Q1, 2)+3.0*(2.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(3.0*pow(Q1, 4)+18.0* 
    pow(Q1, 2)*pow(Q2, 2)+7.0*pow(Q2, 4))+P1*((10.0+21.0*pow(P2, 2))*pow(Q1, 4)+18.0*(2.0+3.0*pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+(10.0+9.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c3L - c3L0) * ( 
    (-5.0/4.0)*pow(G, -4.0)*(pow(G, 4)*P1*P2*(pow(P1, 2)+3.0*(4.0+pow(P2, 2)))+7.0*pow(G, 2)*( 
    pow(P1, 4)*Q1*Q2+6.0*pow(P1, 2)*(4.0+3.0*pow(P2, 2))*Q1*Q2+(16.0+72.0*pow(P2, 2)+9.0* 
    pow(P2, 4))*Q1*Q2+(-2.0)*pow(P1, 3)*P2*(3.0*pow(Q1, 2)+pow(Q2, 2))+(-6.0)*P1*P2*(4.0+ 
    pow(P2, 2))*(3.0*pow(Q1, 2)+pow(Q2, 2)))+(-84.0)*(pow(P1, 4)*pow(Q1, 3)*Q2+(-2.0)*pow(P1, 3)* 
    P2*pow(Q1, 2)*(pow(Q1, 2)+3.0*pow(Q2, 2))+6.0*pow(P1, 2)*Q1*Q2*((3.0+2.0*pow(P2, 2))* 
    pow(Q1, 2)+(1.0+pow(P2, 2))*pow(Q2, 2))+Q1*Q2*((12.0+54.0*pow(P2, 2)+7.0*pow(P2, 4))*pow(Q1, 2)+ 
    2.0*(2.0+9.0*pow(P2, 2)+pow(P2, 4))*pow(Q2, 2))+(-1.0)*P1*P2*((27.0+7.0*pow(P2, 2))*pow(Q1, 4)+ 
    6.0*(9.0+2.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(3.0+pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c4L - c4L0) * ( 
    (15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 3)+(-3.0)*P1*pow(P2, 2))+(-56.0)*pow(G, 2)*(2.0* 
    P2*(1.0+pow(P2, 2))*Q1*Q2+pow(P1, 3)*pow(Q2, 2)+P1*((-1.0)*(1.0+3.0*pow(P2, 2))*pow(Q1, 2)+ 
    pow(Q2, 2)))+(-42.0)*((-12.0)*pow(P1, 2)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-4.0) 
    *P2*Q1*Q2*((8.0+7.0*pow(P2, 2))*pow(Q1, 2)+pow(P2, 2)*pow(Q2, 2))+P1*(pow(Q1, 2)+pow(Q2, 2)) 
    *((8.0+21.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(8.0+3.0*pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 4)+( 
    -6.0)*pow(Q1, 2)*pow(Q2, 2)+(-7.0)*pow(Q2, 4)))) 
    ) + 
    (c5L - c5L0) * ( 
    (3.0/4.0)*pow(G, -4.0)*(pow(G, 4)*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2))+7.0*pow(G, 2)*(3.0* 
    pow(P1, 4)*Q1*Q2+6.0*pow(P1, 2)*(4.0+pow(P2, 2))*Q1*Q2+(-1.0)*pow(P2, 2)*(24.0+5.0* 
    pow(P2, 2))*Q1*Q2+(-2.0)*pow(P1, 3)*P2*(pow(Q1, 2)+3.0*pow(Q2, 2))+2.0*P1*P2*((12.0+ 
    5.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(12.0+pow(P2, 2))*pow(Q2, 2)))+(-84.0)*((-2.0)*pow(P1, 3)* 
    P2*pow(Q2, 2)*(3.0*pow(Q1, 2)+pow(Q2, 2))+pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+2.0*pow(Q2, 2))+6.0* 
    pow(P1, 2)*Q1*Q2*(pow(Q1, 2)+(3.0+pow(P2, 2))*pow(Q2, 2))+Q1*Q2*((-1.0)*(4.0+30.0* 
    pow(P2, 2)+5.0*pow(P2, 4))*pow(Q1, 2)+2.0*(2.0+3.0*pow(P2, 2))*pow(Q2, 2))+P1*P2*(5.0*(3.0+ 
    pow(P2, 2))*pow(Q1, 4)+(-18.0)*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(9.0+pow(P2, 2))*pow(Q2, 4)))) 
    ) + 
    (c6L - c6L0) * ( 
    (-35.0/2.0)*pow(G, -4.0)*(pow(G, 2)*((-6.0)*pow(P1, 2)*P2*Q1*Q2+2.0*pow(P2, 3)*Q1* 
    Q2+pow(P1, 3)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*P1*pow(P2, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+( 
    -3.0)*((-12.0)*pow(P1, 2)*P2*Q1*Q2*(pow(Q1, 2)+pow(Q2, 2))+4.0*P2*Q1*Q2*((2.0+ 
    3.0*pow(P2, 2))*pow(Q1, 2)+(-1.0)*(2.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 4)+6.0*pow(Q1, 2)* 
    pow(Q2, 2)+(-3.0)*pow(Q2, 4))+P1*((-1.0)*(2.0+9.0*pow(P2, 2))*pow(Q1, 4)+6.0*(2.0+3.0*pow(P2, 2)) 
    *pow(Q1, 2)*pow(Q2, 2)+((-2.0)+3.0*pow(P2, 2))*pow(Q2, 4))))
    ) + 
    (c7L - c7L0) * ( 
    (-15.0/4.0)*pow(G, -4.0)*(pow(G, 2)*(pow(P1, 4)*Q1*Q2+(-6.0)*pow(P1, 2)*pow(P2, 2)*Q1* 
    Q2+pow(P2, 4)*Q1*Q2+2.0*pow(P1, 3)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+2.0*P1*pow(P2, 3)*( 
    (-1.0)*pow(Q1, 2)+pow(Q2, 2)))+3.0*(pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+(-5.0)*pow(Q2, 2))+6.0* 
    pow(P1, 2)*Q1*Q2*((4.0+3.0*pow(P2, 2))*pow(Q1, 2)+((-4.0)+pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)* 
    Q1*Q2*((-1.0)*(24.0+7.0*pow(P2, 2))*pow(Q1, 2)+3.0*(8.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 3)* 
    P2*((-3.0)*pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+5.0*pow(Q2, 4))+P1*P2*((12.0+7.0* 
    pow(P2, 2))*pow(Q1, 4)+(-18.0)*(4.0+pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*((-12.0)+pow(P2, 2))* 
    pow(Q2, 4))))
    ) + 
    (c8L - c8L0) * ( 
    (315.0/16.0)*pow(G, -4.0)*((-12.0)*pow(P1, 2)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    4.0*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)* 
    pow(Q2, 2)+pow(Q2, 4))+(-3.0)*P1*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (c9L - c9L0) * ( 
    (35.0/4.0)*pow(G, -4.0)*(pow(P1, 4)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*Q1* 
    Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*Q1*Q2*((-1.0)*pow(Q1, 2)+ 
    pow(Q2, 2))+pow(P1, 3)*P2*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*P1* 
    pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (sL - sL0) * ( 
    (15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*(8.0+pow(P1, 4)+36.0*pow(P2, 2)+5.0*pow(P2, 4)+6.0*pow(P1, 2)*(2.0+ 
    pow(P2, 2)))+(-7.0)*pow(G, 2)*((48.0+3.0*pow(P1, 4)+240.0*pow(P2, 2)+35.0*pow(P2, 4)+6.0*pow(P1, 2)*( 
    8.0+5.0*pow(P2, 2)))*pow(Q1, 2)+(-8.0)*P1*P2*(24.0+3.0*pow(P1, 2)+5.0*pow(P2, 2))*Q1*Q2+( 
    16.0+5.0*pow(P1, 4)+48.0*pow(P2, 2)+5.0*pow(P2, 4)+6.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2)))*pow(Q2, 2))+ 
    21.0*((80.0+3.0*pow(P1, 4)+420.0*pow(P2, 2)+63.0*pow(P2, 4)+6.0*pow(P1, 2)*(10.0+7.0*pow(P2, 2)))* 
    pow(Q1, 4)+(-16.0)*P1*P2*(30.0+3.0*pow(P1, 2)+7.0*pow(P2, 2))*pow(Q1, 3)*Q2+6.0*(16.0+3.0* 
    pow(P1, 4)+60.0*pow(P2, 2)+7.0*pow(P2, 4)+18.0*pow(P1, 2)*(2.0+pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+(-48.0) 
    *P1*P2*(6.0+pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 3)+(16.0+7.0*pow(P1, 4)+36.0*pow(P2, 2)+3.0* 
    pow(P2, 4)+6.0*pow(P1, 2)*(10.0+3.0*pow(P2, 2)))*pow(Q2, 4))) 
    ) + 
    (s2L - s2L0) * ( 
    (15.0/2.0)*pow(G, -4.0)*(2.0*pow(G, 4)*(P2+pow(P2, 3))+(-7.0)*pow(G, 2)*(2.0*pow(P1, 3)*Q1* 
    Q2+(-6.0)*P1*pow(P2, 2)*Q1*Q2+pow(P2, 3)*(15.0*pow(Q1, 2)+pow(Q2, 2))+P2*((16.0+3.0* 
    pow(P1, 2))*pow(Q1, 2)+(-3.0)*pow(P1, 2)*pow(Q2, 2)))+42.0*((-12.0)*P1*pow(P2, 2)*pow(Q1, 3)* 
    Q2+4.0*P1*Q1*Q2*((-1.0)*pow(Q1, 2)+(1.0+pow(P1, 2))*pow(Q2, 2))+2.0*pow(P2, 3)*(7.0* 
    pow(Q1, 4)+3.0*pow(Q1, 2)*pow(Q2, 2))+P2*(3.0*(5.0+pow(P1, 2))*pow(Q1, 4)+6.0*pow(Q1, 2)*pow(Q2, 2)+( 
    -1.0)*(1.0+3.0*pow(P1, 2))*pow(Q2, 4)))) 
    ) + 
    (s3L - s3L0) * ( 
    (-5.0/16.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*pow(P1, 4)+6.0*pow(P1, 2)*(4.0+pow(P2, 2))+(-1.0)* 
    pow(P2, 2)*(24.0+5.0*pow(P2, 2)))+(-14.0)*pow(G, 2)*(((-16.0)+3.0*pow(P1, 4)+(-120.0)*pow(P2, 2)+ 
    (-21.0)*pow(P2, 4)+6.0*pow(P1, 2)*(4.0+pow(P2, 2)))*pow(Q1, 2)+(-8.0)*P1*P2*(12.0+3.0* 
    pow(P1, 2)+pow(P2, 2))*Q1*Q2+(16.0+9.0*pow(P1, 4)+24.0*pow(P2, 2)+pow(P2, 4)+18.0*pow(P1, 2)*(4.0+ 
    pow(P2, 2)))*pow(Q2, 2))+84.0*(((-20.0)+6.0*pow(P1, 2)+pow(P1, 4)+(-126.0)*pow(P2, 2)+(-21.0)* 
    pow(P2, 4))*pow(Q1, 4)+(-16.0)*P1*(3.0+pow(P1, 2))*P2*pow(Q1, 3)*Q2+12.0*(2.0+pow(P1, 4)+3.0* 
    pow(P2, 2)+3.0*pow(P1, 2)*(3.0+pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+(-16.0)*P1*P2*(9.0+2.0* 
    pow(P1, 2)+pow(P2, 2))*Q1*pow(Q2, 3)+(12.0+7.0*pow(P1, 4)+18.0*pow(P2, 2)+pow(P2, 4)+6.0*pow(P1, 2)*(9.0+ 
    2.0*pow(P2, 2)))*pow(Q2, 4)))
    ) + 
    (s4L - s4L0) * ( 
    (-15.0/8.0)*pow(G, -4.0)*(pow(G, 4)*(3.0*pow(P1, 2)*P2+(-1.0)*pow(P2, 3))+28.0*pow(G, 2)*(2.0* 
    P1*(2.0+pow(P1, 2))*Q1*Q2+6.0*P1*pow(P2, 2)*Q1*Q2+pow(P2, 3)*(3.0*pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+(-1.0)*P2*(((-2.0)+3.0*pow(P1, 2))*pow(Q1, 2)+(2.0+3.0*pow(P1, 2))*pow(Q2, 2)))+( 
    -84.0)*(4.0*P1*(2.0+pow(P1, 2))*Q1*Q2*(pow(Q1, 2)+pow(Q2, 2))+12.0*P1*pow(P2, 2)*Q1* 
    Q2*(pow(Q1, 2)+pow(Q2, 2))+pow(P2, 3)*(7.0*pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*pow(Q2, 4)) 
    +(-1.0)*P2*(3.0*((-2.0)+pow(P1, 2))*pow(Q1, 4)+6.0*(2.0+3.0*pow(P1, 2))*pow(Q1, 2)*pow(Q2, 2)+( 
    2.0+3.0*pow(P1, 2))*pow(Q2, 4)))) 
    ) + 
    (s5L - s5L0) * ( 
    (3.0/16.0)*pow(G, -4.0)*(pow(G, 4)*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+(-84.0)*(( 
    (-4.0)+pow(P1, 4)+(-42.0)*pow(P2, 2)+(-9.0)*pow(P2, 4)+6.0*pow(P1, 2)*(3.0+2.0*pow(P2, 2)))* 
    pow(Q1, 4)+(-16.0)*P1*P2*(9.0+pow(P1, 2)+2.0*pow(P2, 2))*pow(Q1, 3)*Q2+12.0*(2.0+9.0* 
    pow(P2, 2)+pow(P2, 4)+3.0*pow(P1, 2)*(1.0+pow(P2, 2)))*pow(Q1, 2)*pow(Q2, 2)+(-16.0)*P1*P2*(3.0+ 
    pow(P2, 2))*Q1*pow(Q2, 3)+((-4.0)+(-30.0)*pow(P1, 2)+(-5.0)*pow(P1, 4)+6.0*pow(P2, 2)+pow(P2, 4))* 
    pow(Q2, 4))+14.0*pow(G, 2)*((-8.0)*pow(P1, 3)*P2*Q1*Q2+(-24.0)*P1*P2*(4.0+pow(P2, 2)) 
    *Q1*Q2+pow(P1, 4)*(pow(Q1, 2)+(-5.0)*pow(Q2, 2))+6.0*pow(P1, 2)*((4.0+3.0*pow(P2, 2))* 
    pow(Q1, 2)+((-4.0)+pow(P2, 2))*pow(Q2, 2))+pow(P2, 2)*((-1.0)*(24.0+7.0*pow(P2, 2))*pow(Q1, 2)+3.0*( 
    8.0+pow(P2, 2))*pow(Q2, 2)))) 
    ) + 
    (s6L - s6L0) * ( 
    (35.0/2.0)*pow(G, -4.0)*(pow(G, 2)*(2.0*pow(P1, 3)*Q1*Q2+(-6.0)*P1*pow(P2, 2)*Q1*Q2+ 
    3.0*pow(P1, 2)*P2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 3)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+(-6.0) 
    *((-12.0)*P1*pow(P2, 2)*pow(Q1, 3)*Q2+4.0*P1*Q1*Q2*((-1.0)*pow(Q1, 2)+(1.0+ 
    pow(P1, 2))*pow(Q2, 2))+(-2.0)*pow(P2, 3)*(pow(Q1, 4)+(-3.0)*pow(Q1, 2)*pow(Q2, 2))+P2*(((-1.0)+ 
    3.0*pow(P1, 2))*pow(Q1, 4)+6.0*pow(Q1, 2)*pow(Q2, 2)+(-1.0)*(1.0+3.0*pow(P1, 2))*pow(Q2, 4))))
    ) + 
    (s7L - s7L0) * ( 
    (-15.0/16.0)*pow(G, -4.0)*(2.0*pow(G, 2)*((-8.0)*pow(P1, 3)*P2*Q1*Q2+8.0*P1* 
    pow(P2, 3)*Q1*Q2+pow(P1, 4)*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P2, 4)*(pow(Q1, 2)+(-1.0)* 
    pow(Q2, 2))+6.0*pow(P1, 2)*pow(P2, 2)*((-1.0)*pow(Q1, 2)+pow(Q2, 2)))+(-3.0)*((-16.0)*pow(P1, 3)* 
    P2*Q1*Q2*(pow(Q1, 2)+3.0*pow(Q2, 2))+16.0*P1*P2*Q1*Q2*((12.0+5.0*pow(P2, 2))* 
    pow(Q1, 2)+(-1.0)*(12.0+pow(P2, 2))*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+18.0*pow(Q1, 2)*pow(Q2, 2)+(-7.0) 
    *pow(Q2, 4))+(-6.0)*pow(P1, 2)*((4.0+5.0*pow(P2, 2))*pow(Q1, 4)+(-6.0)*(4.0+pow(P2, 2))* 
    pow(Q1, 2)*pow(Q2, 2)+(4.0+(-3.0)*pow(P2, 2))*pow(Q2, 4))+pow(P2, 2)*(3.0*(8.0+3.0*pow(P2, 2))* 
    pow(Q1, 4)+(-6.0)*(24.0+5.0*pow(P2, 2))*pow(Q1, 2)*pow(Q2, 2)+(24.0+pow(P2, 2))*pow(Q2, 4)))) 
    ) + 
    (s8L - s8L0) * ( 
    (-315.0/16.0)*pow(G, -4.0)*(4.0*pow(P1, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+(-12.0)* 
    P1*pow(P2, 2)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+3.0*pow(P1, 2)*P2*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-1.0)*pow(P2, 3)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) + 
    (s9L - s9L0) * ( 
    (35.0/16.0)*pow(G, -4.0)*((-16.0)*pow(P1, 3)*P2*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+ 
    16.0*P1*pow(P2, 3)*Q1*Q2*(pow(Q1, 2)+(-1.0)*pow(Q2, 2))+pow(P1, 4)*(pow(Q1, 4)+(-6.0)* 
    pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4))+(-6.0)*pow(P1, 2)*pow(P2, 2)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+ 
    pow(Q2, 4))+pow(P2, 4)*(pow(Q1, 4)+(-6.0)*pow(Q1, 2)*pow(Q2, 2)+pow(Q2, 4)))
    ) ;

    return I;
}

std::vector<double>  perturbation_propagator::integrals_3rd_body(const double &a, const double &P1, const double &P2, const double &Q1, const double &Q2, const double &R, const std::vector<double> &dir) const{

    std::vector<double> I(4, 0.0);
    double alpha = dir[0];
    double beta = dir[1];
    double gamma = dir[2];

    I[0] = ((1.0+(-1.0)*pow(P1, 2)+(-1.0)*pow(P2, 2))*((21.0/16.0) 
      *pow(a, 4)*(15.0*(pow(alpha, 2)+(-1.0)*pow(beta, 2))*(16.0+33.0*pow(alpha, 4)+(-48.0)* 
      pow(beta, 2)+33.0*pow(beta, 4)+6.0*pow(alpha, 2)*((-8.0)+11.0*pow(beta, 2)))*P2*(48.0+( 
      -15.0)*pow(P1, 4)+160.0*pow(P2, 2)+30.0*pow(P1, 2)*pow(P2, 2)+45.0*pow(P2, 4))+28314.0*alpha* 
      beta*(3.0*pow(alpha, 4)+(-10.0)*pow(alpha, 2)*pow(beta, 2)+3.0*pow(beta, 4))*(pow(P1, 5)+( 
      -10.0)*pow(P1, 3.0)*pow(P2, 2)+5.0*P1*pow(P2, 4))+14157.0*(pow(alpha, 2)+(-1.0)*pow(beta, 2))* 
      (pow(alpha, 4)+(-14.0)*pow(alpha, 2)*pow(beta, 2)+pow(beta, 4))*(5.0*pow(P1, 4)*P2+(-10.0)* 
      pow(P1, 2)*pow(P2, 3.0)+pow(P2, 5))+(-1584.0)*alpha*beta*(pow(alpha, 2)+(-1.0)*pow(beta, 2)) 
      *((-10.0)+11.0*pow(alpha, 2)+11.0*pow(beta, 2))*P1*(10.0*pow(P1, 2)+3.0*pow(P1, 4)+(-15.0) 
      *pow(P2, 2)*(2.0+pow(P2, 2)))+10.0*((-16.0)+168.0*(pow(alpha, 2)+pow(beta, 2))+(-378.0)*pow( 
      pow(alpha, 2)+pow(beta, 2), 2)+231.0*pow(pow(alpha, 2)+pow(beta, 2), 3))*P2*(8.0+20.0*( 
      pow(P1, 2)+pow(P2, 2))+5.0*pow(pow(P1, 2)+pow(P2, 2), 2))+30.0*alpha*beta*(16.0+33.0* 
      pow(alpha, 4)+(-48.0)*pow(beta, 2)+33.0*pow(beta, 4)+6.0*pow(alpha, 2)*((-8.0)+11.0* 
      pow(beta, 2)))*P1*(48.0+15.0*pow(P1, 4)+240.0*pow(P2, 2)+75.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+ 
      9.0*pow(P2, 2)))+198.0*((-10.0)+11.0*pow(alpha, 2)+11.0*pow(beta, 2))*(pow(alpha, 4)+(-6.0)* 
      pow(alpha, 2)*pow(beta, 2)+pow(beta, 4))*P2*((-15.0)*pow(P1, 4)+(-30.0)*pow(P1, 2)*(2.0+ 
      pow(P2, 2))+pow(P2, 2)*(20.0+9.0*pow(P2, 2))))+(-105.0/2.0)*pow(a, 3.0)*(21.0*pow(alpha, 4)* 
      beta*P1*P2*(236.0+(-101.0)*pow(P1, 2)+337.0*pow(P2, 2))+84.0*pow(alpha, 5)*(1.0+ 
      pow(P1, 4)+12.0*pow(P2, 2)+8.0*pow(P2, 4)+(-2.0)*pow(P1, 2)*(1.0+6.0*pow(P2, 2)))+beta*P1*P2* 
      (80.0*(2.0+pow(P1, 2)+pow(P2, 2))+(-21.0)*pow(beta, 4)*(52.0+41.0*pow(P1, 2)+11.0*pow(P2, 2))+56.0* 
      pow(beta, 2)*(14.0+pow(P1, 2)+13.0*pow(P2, 2)))+14.0*pow(alpha, 2)*beta*P1*P2*(3.0* 
      pow(beta, 2)*(92.0+193.0*pow(P1, 2)+(-101.0)*pow(P2, 2))+(-4.0)*(82.0+23.0*pow(P1, 2)+59.0* 
      pow(P2, 2)))+4.0*alpha*(8.0+5.0*pow(P1, 4)+60.0*pow(P2, 2)+25.0*pow(P2, 4)+10.0*pow(P1, 2)*(2.0+ 
      3.0*pow(P2, 2))+21.0*pow(beta, 4)*(1.0+16.0*pow(P1, 4)+(-6.0)*pow(P2, 2)+5.0*pow(P2, 4)+pow(P1, 2)*( 
      16.0+(-48.0)*pow(P2, 2)))+(-14.0)*pow(beta, 2)*(2.0+8.0*pow(P1, 4)+(-3.0)*pow(P2, 2)+(-5.0)* 
      pow(P2, 4)+pow(P1, 2)*(23.0+21.0*pow(P2, 2))))+(-56.0)*pow(alpha, 3.0)*(2.0+(-1.0)*pow(P1, 4)+21.0* 
      pow(P2, 2)+10.0*pow(P2, 4)+pow(P1, 2)*((-1.0)+3.0*pow(P2, 2))+3.0*pow(beta, 2)*((-1.0)+8.0*pow(P1, 4)+ 
      (-3.0)*pow(P2, 2)+10.0*pow(P2, 4)+(-1.0)*pow(P1, 2)*(7.0+69.0*pow(P2, 2)))))*R+240.0*pow(a, 2)* 
      ((-49.0)*pow(alpha, 3.0)*beta*P1*((-1.0)+pow(P1, 2)+(-6.0)*pow(P2, 2))+14.0* 
      pow(alpha, 4)*P2*(3.0+(-3.0)*pow(P1, 2)+4.0*pow(P2, 2))+(-1.0)*pow(alpha, 2)*P2*(41.0+ 
      15.0*pow(P1, 2)+36.0*pow(P2, 2)+(-7.0)*pow(beta, 2)*(5.0+51.0*pow(P1, 2)+(-12.0)*pow(P2, 2)))+7.0* 
      alpha*beta*P1*(7.0*pow(beta, 2)*(1.0+2.0*pow(P1, 2)+(-3.0)*pow(P2, 2))+(-3.0)*(2.0+ 
      pow(P1, 2)+3.0*pow(P2, 2)))+P2*(4.0+3.0*pow(P1, 2)+3.0*pow(P2, 2)+(-7.0)*pow(beta, 4)*(1.0+6.0* 
      pow(P1, 2)+(-1.0)*pow(P2, 2))+pow(beta, 2)*(1.0+(-15.0)*pow(P1, 2)+6.0*pow(P2, 2))))*pow(R, 2)+ 
      120.0*a*((-90.0)*pow(alpha, 2)*beta*P1*P2+2.0*beta*((-6.0)+25.0*pow(beta, 2)) 
      *P1*P2+5.0*pow(alpha, 3.0)*((-2.0)+5.0*pow(P1, 2)+pow(P2, 2))+(-1.0)*alpha*((-8.0)+6.0* 
      pow(P1, 2)+18.0*pow(P2, 2)+5.0*pow(beta, 2)*(2.0+9.0*pow(P1, 2)+(-15.0)*pow(P2, 2))))*pow(R, 3.0)+ 
      384.0*(5.0*alpha*beta*P1+4.0*pow(alpha, 2)*P2+(-1.0)*(1.0+pow(beta, 2))*P2)* 
      pow(R, 4))+(-1.0)*gamma*P2*(Q1*((21.0/16.0)*pow(a, 4)*(1584.0*beta*(30.0* 
      pow(alpha, 2)+(-55.0)*pow(alpha, 4)+(-10.0)*pow(beta, 2)+11.0*pow(beta, 4))*P1*P2*( 
      pow(P1, 2)+(-1.0)*pow(P2, 2))*(10.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))+(-14157.0)*(pow(alpha, 5)+( 
      -10.0)*pow(alpha, 3.0)*pow(beta, 2)+5.0*alpha*pow(beta, 4))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*( 
      pow(P1, 4)+(-14.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+198.0*alpha*(33.0*pow(alpha, 4)+60.0* 
      pow(beta, 2)+(-55.0)*pow(beta, 4)+(-10.0)*pow(alpha, 2)*(2.0+11.0*pow(beta, 2)))*(10.0+3.0* 
      pow(P1, 2)+3.0*pow(P2, 2))*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+28314.0*(5.0* 
      pow(alpha, 4)*beta+(-10.0)*pow(alpha, 2)*pow(beta, 3.0)+pow(beta, 5))*P1*P2*(3.0* 
      pow(P1, 4)+(-10.0)*pow(P1, 2)*pow(P2, 2)+3.0*pow(P2, 4))+10.0*alpha*(8.0+(-36.0)*( 
      pow(alpha, 2)+pow(beta, 2))+33.0*pow(pow(alpha, 2)+pow(beta, 2), 2))*(16.0+168.0*(pow(P1, 2)+ 
      pow(P2, 2))+210.0*pow(pow(P1, 2)+pow(P2, 2), 2)+35.0*pow(pow(P1, 2)+pow(P2, 2), 3))+30.0*beta*(16.0+ 
      165.0*pow(alpha, 4)+(-48.0)*pow(beta, 2)+33.0*pow(beta, 4)+18.0*pow(alpha, 2)*((-8.0)+11.0* 
      pow(beta, 2)))*P1*P2*(48.0+15.0*pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*( 
      8.0+3.0*pow(P2, 2)))+(-15.0)*alpha*(16.0+99.0*pow(alpha, 4)+(-33.0)*pow(beta, 4)+ 
      pow(alpha, 2)*((-96.0)+66.0*pow(beta, 2)))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*(48.0+15.0* 
      pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2))))+(-105.0/2.0)* 
      pow(a, 3.0)*(84.0*pow(alpha, 4)*P2*(5.0+5.0*pow(P1, 4)+20.0*pow(P2, 2)+8.0*pow(P2, 4)+(-10.0)* 
      pow(P1, 2)*(1.0+2.0*pow(P2, 2)))+(-21.0)*pow(alpha, 3.0)*beta*P1*((-16.0)+11.0*pow(P1, 4)+( 
      -472.0)*pow(P2, 2)+(-337.0)*pow(P2, 4)+2.0*pow(P1, 2)*(52.0+101.0*pow(P2, 2)))+4.0*P2*(8.0+ 
      5.0*pow(P1, 4)+20.0*pow(P2, 2)+5.0*pow(P2, 4)+10.0*pow(P1, 2)*(2.0+pow(P2, 2))+21.0*pow(beta, 4)*( 
      16.0*pow(P1, 4)+(-16.0)*pow(P1, 2)*((-1.0)+pow(P2, 2))+pow((-1.0)+pow(P2, 2), 2))+(-14.0)* 
      pow(beta, 2)*(2.0+8.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+(-1.0)*pow(P2, 4)+pow(P1, 2)*(23.0+7.0*pow(P2, 2))) 
      )+(-168.0)*pow(alpha, 2)*P2*(2.0+(-1.0)*pow(P1, 4)+7.0*pow(P2, 2)+2.0*pow(P2, 4)+pow(P1, 2)*(( 
      -1.0)+pow(P2, 2))+3.0*pow(beta, 2)*((-1.0)+8.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+2.0*pow(P2, 4)+(-1.0)* 
      pow(P1, 2)*(7.0+23.0*pow(P2, 2))))+(-7.0)*alpha*beta*P1*(3.0*pow(beta, 2)*((-16.0)+ 
      41.0*pow(P1, 4)+(-184.0)*pow(P2, 2)+101.0*pow(P2, 4)+pow(P1, 2)*(8.0+(-386.0)*pow(P2, 2)))+4.0*( 
      8.0+(-13.0)*pow(P1, 4)+164.0*pow(P2, 2)+59.0*pow(P2, 4)+pow(P1, 2)*((-28.0)+46.0*pow(P2, 2)))))* 
      R+240.0*pow(a, 2)*((-147.0)*pow(alpha, 2)*beta*P1*P2*((-1.0)+pow(P1, 2)+(-2.0)* 
      pow(P2, 2))+7.0*beta*P1*P2*(7.0*pow(beta, 2)*(1.0+2.0*pow(P1, 2)+(-1.0)*pow(P2, 2))+(-3.0) 
      *(2.0+pow(P1, 2)+pow(P2, 2)))+7.0*pow(alpha, 3.0)*(1.0+pow(P1, 4)+12.0*pow(P2, 2)+8.0*pow(P2, 4)+(-2.0)* 
      pow(P1, 2)*(1.0+6.0*pow(P2, 2)))+(-1.0)*alpha*(4.0+(-3.0)*pow(P1, 4)+41.0*pow(P2, 2)+18.0* 
      pow(P2, 4)+pow(P1, 2)*((-1.0)+15.0*pow(P2, 2))+7.0*pow(beta, 2)*((-1.0)+6.0*pow(P1, 4)+(-5.0)* 
      pow(P2, 2)+6.0*pow(P2, 4)+(-1.0)*pow(P1, 2)*(5.0+51.0*pow(P2, 2)))))*pow(R, 2)+120.0*a*(10.0* 
      alpha*beta*P1*((-2.0)+5.0*pow(P1, 2)+(-9.0)*pow(P2, 2))+5.0*pow(alpha, 2)*P2*(( 
      -6.0)+15.0*pow(P1, 2)+pow(P2, 2))+P2*(8.0+(-6.0)*pow(P1, 2)+(-6.0)*pow(P2, 2)+(-5.0)* 
      pow(beta, 2)*(2.0+9.0*pow(P1, 2)+(-5.0)*pow(P2, 2))))*pow(R, 3.0)+384.0*(alpha+(-1.0)* 
      alpha*pow(P1, 2)+5.0*beta*P1*P2+4.0*alpha*pow(P2, 2))*pow(R, 4))+(-1.0)*Q2*(( 
      21.0/16.0)*pow(a, 4)*((-1584.0)*alpha*((-10.0)*pow(alpha, 2)+11.0*pow(alpha, 4)+30.0* 
      pow(beta, 2)+(-55.0)*pow(beta, 4))*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2))*(10.0+3.0* 
      pow(P1, 2)+3.0*pow(P2, 2))+14157.0*(5.0*pow(alpha, 4)*beta+(-10.0)*pow(alpha, 2)* 
      pow(beta, 3.0)+pow(beta, 5))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*(pow(P1, 4)+(-14.0)*pow(P1, 2)*pow(P2, 2)+ 
      pow(P2, 4))+(-198.0)*beta*(55.0*pow(alpha, 4)+20.0*pow(beta, 2)+(-33.0)*pow(beta, 4)+10.0* 
      pow(alpha, 2)*((-6.0)+11.0*pow(beta, 2)))*(10.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))*(pow(P1, 4)+(-6.0) 
      *pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+28314.0*(pow(alpha, 5)+(-10.0)*pow(alpha, 3.0)*pow(beta, 2)+ 
      5.0*alpha*pow(beta, 4))*P1*P2*(3.0*pow(P1, 4)+(-10.0)*pow(P1, 2)*pow(P2, 2)+3.0* 
      pow(P2, 4))+10.0*beta*(8.0+(-36.0)*(pow(alpha, 2)+pow(beta, 2))+33.0*pow(pow(alpha, 2)+ 
      pow(beta, 2), 2))*(16.0+168.0*(pow(P1, 2)+pow(P2, 2))+210.0*pow(pow(P1, 2)+pow(P2, 2), 2)+35.0*pow( 
      pow(P1, 2)+pow(P2, 2), 3))+30.0*alpha*(16.0+33.0*pow(alpha, 4)+(-144.0)*pow(beta, 2)+165.0* 
      pow(beta, 4)+6.0*pow(alpha, 2)*((-8.0)+33.0*pow(beta, 2)))*P1*P2*(48.0+15.0*pow(P1, 4)+ 
      80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2)))+15.0*beta*(16.0+(-33.0)* 
      pow(alpha, 4)+(-96.0)*pow(beta, 2)+66.0*pow(alpha, 2)*pow(beta, 2)+99.0*pow(beta, 4))*( 
      pow(P1, 2)+(-1.0)*pow(P2, 2))*(48.0+15.0*pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*( 
      8.0+3.0*pow(P2, 2))))+(105.0/8.0)*pow(a, 3.0)*(((-80.0)+231.0*pow(alpha, 4)+1848.0*pow(beta, 2)+( 
      -3633.0)*pow(beta, 4)+14.0*pow(alpha, 2)*((-52.0)+369.0*pow(beta, 2)))*pow(P1, 5)+3584.0* 
      alpha*beta*(1.0+3.0*pow(alpha, 2)+(-6.0)*pow(beta, 2))*pow(P1, 4)*P2+448.0*alpha* 
      beta*P2*((-1.0)+pow(P2, 2))*((-2.0)+(-1.0)*pow(P2, 2)+(-3.0)*pow(beta, 2)*((-1.0)+ 
      pow(P2, 2))+pow(alpha, 2)*(3.0+6.0*pow(P2, 2)))+(-448.0)*alpha*beta*pow(P1, 2)*P2*(( 
      -23.0)+(-7.0)*pow(P2, 2)+(-48.0)*pow(beta, 2)*((-1.0)+pow(P2, 2))+3.0*pow(alpha, 2)*(7.0+23.0* 
      pow(P2, 2)))+2.0*pow(P1, 3.0)*((-168.0)*pow(beta, 2)*((-18.0)+pow(P2, 2))+(-80.0)*(2.0+pow(P2, 2)) 
      +105.0*pow(beta, 4)*((-44.0)+41.0*pow(P2, 2))+21.0*pow(alpha, 4)*(52.0+101.0*pow(P2, 2))+( 
      -14.0)*pow(alpha, 2)*(56.0+(-92.0)*pow(P2, 2)+9.0*pow(beta, 2)*((-4.0)+193.0*pow(P2, 2))))+( 
      -1.0)*P1*(16.0*(8.0+20.0*pow(P2, 2)+5.0*pow(P2, 4))+(-105.0)*pow(beta, 4)*((-16.0)+104.0* 
      pow(P2, 2)+11.0*pow(P2, 4))+168.0*pow(beta, 2)*((-8.0)+28.0*pow(P2, 2)+13.0*pow(P2, 4))+21.0* 
      pow(alpha, 4)*(16.0+472.0*pow(P2, 2)+337.0*pow(P2, 4))+(-14.0)*pow(alpha, 2)*(4.0*(8.0+164.0* 
      pow(P2, 2)+59.0*pow(P2, 4))+9.0*pow(beta, 2)*((-16.0)+(-184.0)*pow(P2, 2)+101.0*pow(P2, 4)))))* 
      R+240.0*pow(a, 2)*(147.0*alpha*pow(beta, 2)*P1*P2*(1.0+2.0*pow(P1, 2)+(-1.0)* 
      pow(P2, 2))+7.0*pow(beta, 3.0)*(8.0*pow(P1, 4)+(-12.0)*pow(P1, 2)*((-1.0)+pow(P2, 2))+pow((-1.0)+ 
      pow(P2, 2), 2))+(-7.0)*alpha*P1*P2*(7.0*pow(alpha, 2)*((-1.0)+pow(P1, 2)+(-2.0)* 
      pow(P2, 2))+3.0*(2.0+pow(P1, 2)+pow(P2, 2)))+(-1.0)*beta*(4.0+18.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+( 
      -3.0)*pow(P2, 4)+pow(P1, 2)*(41.0+15.0*pow(P2, 2))+7.0*pow(alpha, 2)*((-1.0)+6.0*pow(P1, 4)+(-5.0) 
      *pow(P2, 2)+6.0*pow(P2, 4)+(-1.0)*pow(P1, 2)*(5.0+51.0*pow(P2, 2)))))*pow(R, 2)+(-120.0)*a*(( 
      6.0+(-25.0)*pow(alpha, 2)+(-5.0)*pow(beta, 2))*pow(P1, 3.0)+90.0*alpha*beta*pow(P1, 2)* 
      P2+10.0*alpha*beta*P2*(2.0+(-5.0)*pow(P2, 2))+P1*((-8.0)+6.0*pow(P2, 2)+ 
      pow(beta, 2)*(30.0+(-75.0)*pow(P2, 2))+5.0*pow(alpha, 2)*(2.0+9.0*pow(P2, 2))))*pow(R, 3.0)+384.0* 
      (beta+4.0*beta*pow(P1, 2)+5.0*alpha*P1*P2+(-1.0)*beta*pow(P2, 2))*pow(R, 4)))); 


    I[1] = ((1.0+(-1.0)*pow(P1, 2)+(-1.0)*pow(P2, 2))*((21.0/16.0) 
      *pow(a, 4)*(1584.0*alpha*beta*(pow(alpha, 2)+(-1.0)*pow(beta, 2))*((-10.0)+11.0* 
      pow(alpha, 2)+11.0*pow(beta, 2))*P2*((-30.0)*pow(P1, 2)+(-15.0)*pow(P1, 4)+10.0*pow(P2, 2)+ 
      3.0*pow(P2, 4))+(-14157.0)*(pow(alpha, 2)+(-1.0)*pow(beta, 2))*(pow(alpha, 4)+(-14.0)* 
      pow(alpha, 2)*pow(beta, 2)+pow(beta, 4))*(pow(P1, 5)+(-10.0)*pow(P1, 3.0)*pow(P2, 2)+5.0*P1* 
      pow(P2, 4))+28314.0*alpha*beta*(3.0*pow(alpha, 4)+(-10.0)*pow(alpha, 2)*pow(beta, 2)+ 
      3.0*pow(beta, 4))*(5.0*pow(P1, 4)*P2+(-10.0)*pow(P1, 2)*pow(P2, 3.0)+pow(P2, 5))+198.0*((-10.0)+ 
      11.0*pow(alpha, 2)+11.0*pow(beta, 2))*(pow(alpha, 4)+(-6.0)*pow(alpha, 2)*pow(beta, 2)+ 
      pow(beta, 4))*P1*(9.0*pow(P1, 4)+pow(P1, 2)*(20.0+(-30.0)*pow(P2, 2))+(-15.0)*pow(P2, 2)*(4.0+ 
      pow(P2, 2)))+10.0*((-16.0)+168.0*(pow(alpha, 2)+pow(beta, 2))+(-378.0)*pow(pow(alpha, 2)+ 
      pow(beta, 2), 2)+231.0*pow(pow(alpha, 2)+pow(beta, 2), 3))*P1*(8.0+20.0*(pow(P1, 2)+pow(P2, 2))+ 
      5.0*pow(pow(P1, 2)+pow(P2, 2), 2))+30.0*alpha*beta*(16.0+33.0*pow(alpha, 4)+(-48.0)* 
      pow(beta, 2)+33.0*pow(beta, 4)+6.0*pow(alpha, 2)*((-8.0)+11.0*pow(beta, 2)))*P2*(48.0+75.0* 
      pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+30.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2)))+(-15.0)*( 
      pow(alpha, 2)+(-1.0)*pow(beta, 2))*(16.0+33.0*pow(alpha, 4)+(-48.0)*pow(beta, 2)+33.0* 
      pow(beta, 4)+6.0*pow(alpha, 2)*((-8.0)+11.0*pow(beta, 2)))*P1*(48.0+45.0*pow(P1, 4)+(-15.0) 
      *pow(P2, 4)+10.0*pow(P1, 2)*(16.0+3.0*pow(P2, 2))))+(-105.0/8.0)*pow(a, 3.0)*(10752.0*alpha* 
      pow(beta, 4)*P1*P2*(1.0+2.0*pow(P1, 2)+(-1.0)*pow(P2, 2))+(-448.0)*alpha*pow(beta, 2)* 
      P1*P2*(23.0+16.0*pow(P1, 2)+7.0*pow(P2, 2)+3.0*pow(alpha, 2)*((-7.0)+16.0*pow(P1, 2)+(-23.0) 
      *pow(P2, 2)))+64.0*alpha*P1*P2*(21.0*pow(alpha, 4)*((-1.0)+pow(P1, 2)+(-2.0)* 
      pow(P2, 2))+7.0*pow(alpha, 2)*(1.0+2.0*pow(P1, 2)+(-1.0)*pow(P2, 2))+5.0*(2.0+pow(P1, 2)+pow(P2, 2)))+ 
      21.0*pow(beta, 5)*(16.0+173.0*pow(P1, 4)+(-104.0)*pow(P2, 2)+(-11.0)*pow(P2, 4)+(-6.0)* 
      pow(P1, 2)*((-44.0)+41.0*pow(P2, 2)))+(-14.0)*pow(beta, 3.0)*(pow(alpha, 2)*((-48.0)+615.0* 
      pow(P1, 4)+(-552.0)*pow(P2, 2)+303.0*pow(P2, 4)+pow(P1, 2)*(72.0+(-3474.0)*pow(P2, 2)))+4.0*(8.0+ 
      55.0*pow(P1, 4)+(-28.0)*pow(P2, 2)+(-13.0)*pow(P2, 4)+(-6.0)*pow(P1, 2)*((-18.0)+pow(P2, 2))))+ 
      beta*(56.0*pow(alpha, 2)*((-8.0)+65.0*pow(P1, 4)+(-164.0)*pow(P2, 2)+(-59.0)*pow(P2, 4)+ 
      pow(P1, 2)*(84.0+(-138.0)*pow(P2, 2)))+16.0*(8.0+25.0*pow(P1, 4)+20.0*pow(P2, 2)+5.0*pow(P2, 4)+ 
      30.0*pow(P1, 2)*(2.0+pow(P2, 2)))+(-21.0)*pow(alpha, 4)*((-16.0)+55.0*pow(P1, 4)+(-472.0)* 
      pow(P2, 2)+(-337.0)*pow(P2, 4)+6.0*pow(P1, 2)*(52.0+101.0*pow(P2, 2)))))*R+30.0*pow(a, 2)*( 
      14.0*((-6.0)+7.0*(pow(alpha, 2)+pow(beta, 2)))*((-1.0)*pow(alpha, 2)*P1+pow(beta, 2)*P1+ 
      2.0*alpha*beta*P2)*(2.0+pow(P1, 2)+pow(P2, 2))+147.0*(pow(alpha, 4)+(-6.0)* 
      pow(alpha, 2)*pow(beta, 2)+pow(beta, 4))*(pow(P1, 3.0)+(-3.0)*P1*pow(P2, 2))+(-588.0)* 
      alpha*beta*(pow(alpha, 2)+(-1.0)*pow(beta, 2))*(3.0*pow(P1, 2)*P2+(-1.0)*pow(P2, 3.0))+ 
      14.0*((-6.0)+7.0*(pow(alpha, 2)+pow(beta, 2)))*P1*(4.0*alpha*beta*P1*P2+( 
      pow(alpha, 2)+(-1.0)*pow(beta, 2))*((-1.0)*pow(P1, 2)+pow(P2, 2)))+(8.0+(-40.0)*(pow(alpha, 2)+ 
      pow(beta, 2))+35.0*pow(pow(alpha, 2)+pow(beta, 2), 2))*P1*(4.0+3.0*(pow(P1, 2)+pow(P2, 2))))* 
      pow(R, 2)+120.0*a*(2.0*alpha*((-6.0)+25.0*pow(alpha, 2))*P1*P2+(-90.0)*alpha* 
      pow(beta, 2)*P1*P2+5.0*pow(beta, 3.0)*((-2.0)+pow(P1, 2)+5.0*pow(P2, 2))+beta*(8.0+(-18.0)* 
      pow(P1, 2)+(-6.0)*pow(P2, 2)+5.0*pow(alpha, 2)*((-2.0)+15.0*pow(P1, 2)+(-9.0)*pow(P2, 2))))* 
      pow(R, 3.0)+(-384.0)*((1.0+pow(alpha, 2)+(-4.0)*pow(beta, 2))*P1+(-5.0)*alpha*beta* 
      P2)*pow(R, 4))+(-1.0)*gamma*P1*(Q1*((21.0/16.0)*pow(a, 4)*(1584.0*beta*( 
      30.0*pow(alpha, 2)+(-55.0)*pow(alpha, 4)+(-10.0)*pow(beta, 2)+11.0*pow(beta, 4))*P1*P2* 
      (pow(P1, 2)+(-1.0)*pow(P2, 2))*(10.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))+(-14157.0)*(pow(alpha, 5)+( 
      -10.0)*pow(alpha, 3.0)*pow(beta, 2)+5.0*alpha*pow(beta, 4))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*( 
      pow(P1, 4)+(-14.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+198.0*alpha*(33.0*pow(alpha, 4)+60.0* 
      pow(beta, 2)+(-55.0)*pow(beta, 4)+(-10.0)*pow(alpha, 2)*(2.0+11.0*pow(beta, 2)))*(10.0+3.0* 
      pow(P1, 2)+3.0*pow(P2, 2))*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+28314.0*(5.0* 
      pow(alpha, 4)*beta+(-10.0)*pow(alpha, 2)*pow(beta, 3.0)+pow(beta, 5))*P1*P2*(3.0* 
      pow(P1, 4)+(-10.0)*pow(P1, 2)*pow(P2, 2)+3.0*pow(P2, 4))+10.0*alpha*(8.0+(-36.0)*( 
      pow(alpha, 2)+pow(beta, 2))+33.0*pow(pow(alpha, 2)+pow(beta, 2), 2))*(16.0+168.0*(pow(P1, 2)+ 
      pow(P2, 2))+210.0*pow(pow(P1, 2)+pow(P2, 2), 2)+35.0*pow(pow(P1, 2)+pow(P2, 2), 3))+30.0*beta*(16.0+ 
      165.0*pow(alpha, 4)+(-48.0)*pow(beta, 2)+33.0*pow(beta, 4)+18.0*pow(alpha, 2)*((-8.0)+11.0* 
      pow(beta, 2)))*P1*P2*(48.0+15.0*pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*( 
      8.0+3.0*pow(P2, 2)))+(-15.0)*alpha*(16.0+99.0*pow(alpha, 4)+(-33.0)*pow(beta, 4)+ 
      pow(alpha, 2)*((-96.0)+66.0*pow(beta, 2)))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*(48.0+15.0* 
      pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2))))+(-105.0/2.0)* 
      pow(a, 3.0)*(84.0*pow(alpha, 4)*P2*(5.0+5.0*pow(P1, 4)+20.0*pow(P2, 2)+8.0*pow(P2, 4)+(-10.0)* 
      pow(P1, 2)*(1.0+2.0*pow(P2, 2)))+(-21.0)*pow(alpha, 3.0)*beta*P1*((-16.0)+11.0*pow(P1, 4)+( 
      -472.0)*pow(P2, 2)+(-337.0)*pow(P2, 4)+2.0*pow(P1, 2)*(52.0+101.0*pow(P2, 2)))+4.0*P2*(8.0+ 
      5.0*pow(P1, 4)+20.0*pow(P2, 2)+5.0*pow(P2, 4)+10.0*pow(P1, 2)*(2.0+pow(P2, 2))+21.0*pow(beta, 4)*( 
      16.0*pow(P1, 4)+(-16.0)*pow(P1, 2)*((-1.0)+pow(P2, 2))+pow((-1.0)+pow(P2, 2), 2))+(-14.0)* 
      pow(beta, 2)*(2.0+8.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+(-1.0)*pow(P2, 4)+pow(P1, 2)*(23.0+7.0*pow(P2, 2))) 
      )+(-168.0)*pow(alpha, 2)*P2*(2.0+(-1.0)*pow(P1, 4)+7.0*pow(P2, 2)+2.0*pow(P2, 4)+pow(P1, 2)*(( 
      -1.0)+pow(P2, 2))+3.0*pow(beta, 2)*((-1.0)+8.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+2.0*pow(P2, 4)+(-1.0)* 
      pow(P1, 2)*(7.0+23.0*pow(P2, 2))))+(-7.0)*alpha*beta*P1*(3.0*pow(beta, 2)*((-16.0)+ 
      41.0*pow(P1, 4)+(-184.0)*pow(P2, 2)+101.0*pow(P2, 4)+pow(P1, 2)*(8.0+(-386.0)*pow(P2, 2)))+4.0*( 
      8.0+(-13.0)*pow(P1, 4)+164.0*pow(P2, 2)+59.0*pow(P2, 4)+pow(P1, 2)*((-28.0)+46.0*pow(P2, 2)))))* 
      R+240.0*pow(a, 2)*((-147.0)*pow(alpha, 2)*beta*P1*P2*((-1.0)+pow(P1, 2)+(-2.0)* 
      pow(P2, 2))+7.0*beta*P1*P2*(7.0*pow(beta, 2)*(1.0+2.0*pow(P1, 2)+(-1.0)*pow(P2, 2))+(-3.0) 
      *(2.0+pow(P1, 2)+pow(P2, 2)))+7.0*pow(alpha, 3)*(1.0+pow(P1, 4)+12.0*pow(P2, 2)+8.0*pow(P2, 4)+(-2.0)* 
      pow(P1, 2)*(1.0+6.0*pow(P2, 2)))+(-1.0)*alpha*(4.0+(-3.0)*pow(P1, 4)+41.0*pow(P2, 2)+18.0* 
      pow(P2, 4)+pow(P1, 2)*((-1.0)+15.0*pow(P2, 2))+7.0*pow(beta, 2)*((-1.0)+6.0*pow(P1, 4)+(-5.0)* 
      pow(P2, 2)+6.0*pow(P2, 4)+(-1.0)*pow(P1, 2)*(5.0+51.0*pow(P2, 2)))))*pow(R, 2)+120.0*a*(10.0* 
      alpha*beta*P1*((-2.0)+5.0*pow(P1, 2)+(-9.0)*pow(P2, 2))+5.0*pow(alpha, 2)*P2*(( 
      -6.0)+15.0*pow(P1, 2)+pow(P2, 2))+P2*(8.0+(-6.0)*pow(P1, 2)+(-6.0)*pow(P2, 2)+(-5.0)* 
      pow(beta, 2)*(2.0+9.0*pow(P1, 2)+(-5.0)*pow(P2, 2))))*pow(R, 3.0)+384.0*(alpha+(-1.0)* 
      alpha*pow(P1, 2)+5.0*beta*P1*P2+4.0*alpha*pow(P2, 2))*pow(R, 4))+(-1.0)*Q2*(( 
      21.0/16.0)*pow(a, 4)*((-1584.0)*alpha*((-10.0)*pow(alpha, 2)+11.0*pow(alpha, 4)+30.0* 
      pow(beta, 2)+(-55.0)*pow(beta, 4))*P1*P2*(pow(P1, 2)+(-1.0)*pow(P2, 2))*(10.0+3.0* 
      pow(P1, 2)+3.0*pow(P2, 2))+14157.0*(5.0*pow(alpha, 4)*beta+(-10.0)*pow(alpha, 2)* 
      pow(beta, 3.0)+pow(beta, 5))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*(pow(P1, 4)+(-14.0)*pow(P1, 2)*pow(P2, 2)+ 
      pow(P2, 4))+(-198.0)*beta*(55.0*pow(alpha, 4)+20.0*pow(beta, 2)+(-33.0)*pow(beta, 4)+10.0* 
      pow(alpha, 2)*((-6.0)+11.0*pow(beta, 2)))*(10.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))*(pow(P1, 4)+(-6.0) 
      *pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+28314.0*(pow(alpha, 5)+(-10.0)*pow(alpha, 3.0)*pow(beta, 2)+ 
      5.0*alpha*pow(beta, 4))*P1*P2*(3.0*pow(P1, 4)+(-10.0)*pow(P1, 2)*pow(P2, 2)+3.0* 
      pow(P2, 4))+10.0*beta*(8.0+(-36.0)*(pow(alpha, 2)+pow(beta, 2))+33.0*pow(pow(alpha, 2)+ 
      pow(beta, 2), 2))*(16.0+168.0*(pow(P1, 2)+pow(P2, 2))+210.0*pow(pow(P1, 2)+pow(P2, 2), 2)+35.0*pow( 
      pow(P1, 2)+pow(P2, 2), 3))+30.0*alpha*(16.0+33.0*pow(alpha, 4)+(-144.0)*pow(beta, 2)+165.0* 
      pow(beta, 4)+6.0*pow(alpha, 2)*((-8.0)+33.0*pow(beta, 2)))*P1*P2*(48.0+15.0*pow(P1, 4)+ 
      80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2)))+15.0*beta*(16.0+(-33.0)* 
      pow(alpha, 4)+(-96.0)*pow(beta, 2)+66.0*pow(alpha, 2)*pow(beta, 2)+99.0*pow(beta, 4))*( 
      pow(P1, 2)+(-1.0)*pow(P2, 2))*(48.0+15.0*pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*( 
      8.0+3.0*pow(P2, 2))))+(105.0/8.0)*pow(a, 3.0)*(((-80.0)+231.0*pow(alpha, 4)+1848.0*pow(beta, 2)+( 
      -3633.0)*pow(beta, 4)+14.0*pow(alpha, 2)*((-52.0)+369.0*pow(beta, 2)))*pow(P1, 5)+3584.0* 
      alpha*beta*(1.0+3.0*pow(alpha, 2)+(-6.0)*pow(beta, 2))*pow(P1, 4)*P2+448.0*alpha* 
      beta*P2*((-1.0)+pow(P2, 2))*((-2.0)+(-1.0)*pow(P2, 2)+(-3.0)*pow(beta, 2)*((-1.0)+ 
      pow(P2, 2))+pow(alpha, 2)*(3.0+6.0*pow(P2, 2)))+(-448.0)*alpha*beta*pow(P1, 2)*P2*(( 
      -23.0)+(-7.0)*pow(P2, 2)+(-48.0)*pow(beta, 2)*((-1.0)+pow(P2, 2))+3.0*pow(alpha, 2)*(7.0+23.0* 
      pow(P2, 2)))+2.0*pow(P1, 3.0)*((-168.0)*pow(beta, 2)*((-18.0)+pow(P2, 2))+(-80.0)*(2.0+pow(P2, 2)) 
      +105.0*pow(beta, 4)*((-44.0)+41.0*pow(P2, 2))+21.0*pow(alpha, 4)*(52.0+101.0*pow(P2, 2))+( 
      -14.0)*pow(alpha, 2)*(56.0+(-92.0)*pow(P2, 2)+9.0*pow(beta, 2)*((-4.0)+193.0*pow(P2, 2))))+( 
      -1.0)*P1*(16.0*(8.0+20.0*pow(P2, 2)+5.0*pow(P2, 4))+(-105.0)*pow(beta, 4)*((-16.0)+104.0* 
      pow(P2, 2)+11.0*pow(P2, 4))+168.0*pow(beta, 2)*((-8.0)+28.0*pow(P2, 2)+13.0*pow(P2, 4))+21.0* 
      pow(alpha, 4)*(16.0+472.0*pow(P2, 2)+337.0*pow(P2, 4))+(-14.0)*pow(alpha, 2)*(4.0*(8.0+164.0* 
      pow(P2, 2)+59.0*pow(P2, 4))+9.0*pow(beta, 2)*((-16.0)+(-184.0)*pow(P2, 2)+101.0*pow(P2, 4)))))* 
      R+240.0*pow(a, 2)*(147.0*alpha*pow(beta, 2)*P1*P2*(1.0+2.0*pow(P1, 2)+(-1.0)* 
      pow(P2, 2))+7.0*pow(beta, 3.0)*(8.0*pow(P1, 4)+(-12.0)*pow(P1, 2)*((-1.0)+pow(P2, 2))+pow((-1.0)+ 
      pow(P2, 2), 2))+(-7.0)*alpha*P1*P2*(7.0*pow(alpha, 2)*((-1.0)+pow(P1, 2)+(-2.0)* 
      pow(P2, 2))+3.0*(2.0+pow(P1, 2)+pow(P2, 2)))+(-1.0)*beta*(4.0+18.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+( 
      -3.0)*pow(P2, 4)+pow(P1, 2)*(41.0+15.0*pow(P2, 2))+7.0*pow(alpha, 2)*((-1.0)+6.0*pow(P1, 4)+(-5.0) 
      *pow(P2, 2)+6.0*pow(P2, 4)+(-1.0)*pow(P1, 2)*(5.0+51.0*pow(P2, 2)))))*pow(R, 2)+(-120.0)*a*(( 
      6.0+(-25.0)*pow(alpha, 2)+(-5.0)*pow(beta, 2))*pow(P1, 3.0)+90.0*alpha*beta*pow(P1, 2)* 
      P2+10.0*alpha*beta*P2*(2.0+(-5.0)*pow(P2, 2))+P1*((-8.0)+6.0*pow(P2, 2)+ 
      pow(beta, 2)*(30.0+(-75.0)*pow(P2, 2))+5.0*pow(alpha, 2)*(2.0+9.0*pow(P2, 2))))*pow(R, 3.0)+384.0* 
      (beta+4.0*beta*pow(P1, 2)+5.0*alpha*P1*P2+(-1.0)*beta*pow(P2, 2))*pow(R, 4)))); 


    I[2] = ((21.0/16.0)*pow(a, 4)*((-1584.0)*alpha*(( 
      -10.0)*pow(alpha, 2)+11.0*pow(alpha, 4)+30.0*pow(beta, 2)+(-55.0)*pow(beta, 4))*P1*P2*( 
      pow(P1, 2)+(-1.0)*pow(P2, 2))*(10.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))+14157.0*(5.0*pow(alpha, 4)* 
      beta+(-10.0)*pow(alpha, 2)*pow(beta, 3.0)+pow(beta, 5))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*( 
      pow(P1, 4)+(-14.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+(-198.0)*beta*(55.0*pow(alpha, 4)+20.0* 
      pow(beta, 2)+(-33.0)*pow(beta, 4)+10.0*pow(alpha, 2)*((-6.0)+11.0*pow(beta, 2)))*(10.0+3.0* 
      pow(P1, 2)+3.0*pow(P2, 2))*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+28314.0*( 
      pow(alpha, 5)+(-10.0)*pow(alpha, 3.0)*pow(beta, 2)+5.0*alpha*pow(beta, 4))*P1*P2*(3.0* 
      pow(P1, 4)+(-10.0)*pow(P1, 2)*pow(P2, 2)+3.0*pow(P2, 4))+10.0*beta*(8.0+(-36.0)*(pow(alpha, 2)+ 
      pow(beta, 2))+33.0*pow(pow(alpha, 2)+pow(beta, 2), 2))*(16.0+168.0*(pow(P1, 2)+pow(P2, 2))+210.0*pow( 
      pow(P1, 2)+pow(P2, 2), 2)+35.0*pow(pow(P1, 2)+pow(P2, 2), 3))+30.0*alpha*(16.0+33.0*pow(alpha, 4)+ 
      (-144.0)*pow(beta, 2)+165.0*pow(beta, 4)+6.0*pow(alpha, 2)*((-8.0)+33.0*pow(beta, 2)))* 
      P1*P2*(48.0+15.0*pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2))) 
      +15.0*beta*(16.0+(-33.0)*pow(alpha, 4)+(-96.0)*pow(beta, 2)+66.0*pow(alpha, 2)* 
      pow(beta, 2)+99.0*pow(beta, 4))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*(48.0+15.0*pow(P1, 4)+80.0* 
      pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2))))+(105.0/8.0)*pow(a, 3.0)*(((-80.0)+ 
      231.0*pow(alpha, 4)+1848.0*pow(beta, 2)+(-3633.0)*pow(beta, 4)+14.0*pow(alpha, 2)*((-52.0)+ 
      369.0*pow(beta, 2)))*pow(P1, 5)+3584.0*alpha*beta*(1.0+3.0*pow(alpha, 2)+(-6.0)* 
      pow(beta, 2))*pow(P1, 4)*P2+448.0*alpha*beta*P2*((-1.0)+pow(P2, 2))*((-2.0)+(-1.0) 
      *pow(P2, 2)+(-3.0)*pow(beta, 2)*((-1.0)+pow(P2, 2))+pow(alpha, 2)*(3.0+6.0*pow(P2, 2)))+(-448.0) 
      *alpha*beta*pow(P1, 2)*P2*((-23.0)+(-7.0)*pow(P2, 2)+(-48.0)*pow(beta, 2)*((-1.0) 
      +pow(P2, 2))+3.0*pow(alpha, 2)*(7.0+23.0*pow(P2, 2)))+2.0*pow(P1, 3.0)*((-168.0)*pow(beta, 2)*(( 
      -18.0)+pow(P2, 2))+(-80.0)*(2.0+pow(P2, 2))+105.0*pow(beta, 4)*((-44.0)+41.0*pow(P2, 2))+21.0* 
      pow(alpha, 4)*(52.0+101.0*pow(P2, 2))+(-14.0)*pow(alpha, 2)*(56.0+(-92.0)*pow(P2, 2)+9.0* 
      pow(beta, 2)*((-4.0)+193.0*pow(P2, 2))))+(-1.0)*P1*(16.0*(8.0+20.0*pow(P2, 2)+5.0*pow(P2, 4)) 
      +(-105.0)*pow(beta, 4)*((-16.0)+104.0*pow(P2, 2)+11.0*pow(P2, 4))+168.0*pow(beta, 2)*((-8.0) 
      +28.0*pow(P2, 2)+13.0*pow(P2, 4))+21.0*pow(alpha, 4)*(16.0+472.0*pow(P2, 2)+337.0*pow(P2, 4))+( 
      -14.0)*pow(alpha, 2)*(4.0*(8.0+164.0*pow(P2, 2)+59.0*pow(P2, 4))+9.0*pow(beta, 2)*((-16.0)+( 
      -184.0)*pow(P2, 2)+101.0*pow(P2, 4)))))*R+240.0*pow(a, 2)*(147.0*alpha*pow(beta, 2)* 
      P1*P2*(1.0+2.0*pow(P1, 2)+(-1.0)*pow(P2, 2))+7.0*pow(beta, 3.0)*(8.0*pow(P1, 4)+(-12.0)* 
      pow(P1, 2)*((-1.0)+pow(P2, 2))+pow((-1.0)+pow(P2, 2), 2))+(-7.0)*alpha*P1*P2*(7.0* 
      pow(alpha, 2)*((-1.0)+pow(P1, 2)+(-2.0)*pow(P2, 2))+3.0*(2.0+pow(P1, 2)+pow(P2, 2)))+(-1.0)* 
      beta*(4.0+18.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+(-3.0)*pow(P2, 4)+pow(P1, 2)*(41.0+15.0*pow(P2, 2))+ 
      7.0*pow(alpha, 2)*((-1.0)+6.0*pow(P1, 4)+(-5.0)*pow(P2, 2)+6.0*pow(P2, 4)+(-1.0)*pow(P1, 2)*(5.0+ 
      51.0*pow(P2, 2)))))*pow(R, 2)+(-120.0)*a*((6.0+(-25.0)*pow(alpha, 2)+(-5.0)*pow(beta, 2)) 
      *pow(P1, 3.0)+90.0*alpha*beta*pow(P1, 2)*P2+10.0*alpha*beta*P2*(2.0+(-5.0)* 
      pow(P2, 2))+P1*((-8.0)+6.0*pow(P2, 2)+pow(beta, 2)*(30.0+(-75.0)*pow(P2, 2))+5.0*pow(alpha, 2)* 
      (2.0+9.0*pow(P2, 2))))*pow(R, 3.0)+384.0*(beta+4.0*beta*pow(P1, 2)+5.0*alpha*P1*P2+( 
      -1.0)*beta*pow(P2, 2))*pow(R, 4));


    I[3] = ((21.0/16.0)*pow(a, 4)*(1584.0*beta*(30.0* 
      pow(alpha, 2)+(-55.0)*pow(alpha, 4)+(-10.0)*pow(beta, 2)+11.0*pow(beta, 4))*P1*P2*( 
      pow(P1, 2)+(-1.0)*pow(P2, 2))*(10.0+3.0*pow(P1, 2)+3.0*pow(P2, 2))+(-14157.0)*(pow(alpha, 5)+( 
      -10.0)*pow(alpha, 3.0)*pow(beta, 2)+5.0*alpha*pow(beta, 4))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*( 
      pow(P1, 4)+(-14.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+198.0*alpha*(33.0*pow(alpha, 4)+60.0* 
      pow(beta, 2)+(-55.0)*pow(beta, 4)+(-10.0)*pow(alpha, 2)*(2.0+11.0*pow(beta, 2)))*(10.0+3.0* 
      pow(P1, 2)+3.0*pow(P2, 2))*(pow(P1, 4)+(-6.0)*pow(P1, 2)*pow(P2, 2)+pow(P2, 4))+28314.0*(5.0* 
      pow(alpha, 4)*beta+(-10.0)*pow(alpha, 2)*pow(beta, 3.0)+pow(beta, 5))*P1*P2*(3.0* 
      pow(P1, 4)+(-10.0)*pow(P1, 2)*pow(P2, 2)+3.0*pow(P2, 4))+10.0*alpha*(8.0+(-36.0)*( 
      pow(alpha, 2)+pow(beta, 2))+33.0*pow(pow(alpha, 2)+pow(beta, 2), 2))*(16.0+168.0*(pow(P1, 2)+ 
      pow(P2, 2))+210.0*pow(pow(P1, 2)+pow(P2, 2), 2)+35.0*pow(pow(P1, 2)+pow(P2, 2), 3))+30.0*beta*(16.0+ 
      165.0*pow(alpha, 4)+(-48.0)*pow(beta, 2)+33.0*pow(beta, 4)+18.0*pow(alpha, 2)*((-8.0)+11.0* 
      pow(beta, 2)))*P1*P2*(48.0+15.0*pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*( 
      8.0+3.0*pow(P2, 2)))+(-15.0)*alpha*(16.0+99.0*pow(alpha, 4)+(-33.0)*pow(beta, 4)+ 
      pow(alpha, 2)*((-96.0)+66.0*pow(beta, 2)))*(pow(P1, 2)+(-1.0)*pow(P2, 2))*(48.0+15.0* 
      pow(P1, 4)+80.0*pow(P2, 2)+15.0*pow(P2, 4)+10.0*pow(P1, 2)*(8.0+3.0*pow(P2, 2))))+(-105.0/2.0)* 
      pow(a, 3.0)*(84.0*pow(alpha, 4)*P2*(5.0+5.0*pow(P1, 4)+20.0*pow(P2, 2)+8.0*pow(P2, 4)+(-10.0)* 
      pow(P1, 2)*(1.0+2.0*pow(P2, 2)))+(-21.0)*pow(alpha, 3.0)*beta*P1*((-16.0)+11.0*pow(P1, 4)+( 
      -472.0)*pow(P2, 2)+(-337.0)*pow(P2, 4)+2.0*pow(P1, 2)*(52.0+101.0*pow(P2, 2)))+4.0*P2*(8.0+ 
      5.0*pow(P1, 4)+20.0*pow(P2, 2)+5.0*pow(P2, 4)+10.0*pow(P1, 2)*(2.0+pow(P2, 2))+21.0*pow(beta, 4)*( 
      16.0*pow(P1, 4)+(-16.0)*pow(P1, 2)*((-1.0)+pow(P2, 2))+pow((-1.0)+pow(P2, 2), 2))+(-14.0)* 
      pow(beta, 2)*(2.0+8.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+(-1.0)*pow(P2, 4)+pow(P1, 2)*(23.0+7.0*pow(P2, 2))) 
      )+(-168.0)*pow(alpha, 2)*P2*(2.0+(-1.0)*pow(P1, 4)+7.0*pow(P2, 2)+2.0*pow(P2, 4)+pow(P1, 2)*(( 
      -1.0)+pow(P2, 2))+3.0*pow(beta, 2)*((-1.0)+8.0*pow(P1, 4)+(-1.0)*pow(P2, 2)+2.0*pow(P2, 4)+(-1.0)* 
      pow(P1, 2)*(7.0+23.0*pow(P2, 2))))+(-7.0)*alpha*beta*P1*(3.0*pow(beta, 2)*((-16.0)+ 
      41.0*pow(P1, 4)+(-184.0)*pow(P2, 2)+101.0*pow(P2, 4)+pow(P1, 2)*(8.0+(-386.0)*pow(P2, 2)))+4.0*( 
      8.0+(-13.0)*pow(P1, 4)+164.0*pow(P2, 2)+59.0*pow(P2, 4)+pow(P1, 2)*((-28.0)+46.0*pow(P2, 2)))))* 
      R+240.0*pow(a, 2)*((-147.0)*pow(alpha, 2)*beta*P1*P2*((-1.0)+pow(P1, 2)+(-2.0)* 
      pow(P2, 2))+7.0*beta*P1*P2*(7.0*pow(beta, 2)*(1.0+2.0*pow(P1, 2)+(-1.0)*pow(P2, 2))+(-3.0) 
      *(2.0+pow(P1, 2)+pow(P2, 2)))+7.0*pow(alpha, 3.0)*(1.0+pow(P1, 4)+12.0*pow(P2, 2)+8.0*pow(P2, 4)+(-2.0)* 
      pow(P1, 2)*(1.0+6.0*pow(P2, 2)))+(-1.0)*alpha*(4.0+(-3.0)*pow(P1, 4)+41.0*pow(P2, 2)+18.0* 
      pow(P2, 4)+pow(P1, 2)*((-1.0)+15.0*pow(P2, 2))+7.0*pow(beta, 2)*((-1.0)+6.0*pow(P1, 4)+(-5.0)* 
      pow(P2, 2)+6.0*pow(P2, 4)+(-1.0)*pow(P1, 2)*(5.0+51.0*pow(P2, 2)))))*pow(R, 2)+120.0*a*(10.0* 
      alpha*beta*P1*((-2.0)+5.0*pow(P1, 2)+(-9.0)*pow(P2, 2))+5.0*pow(alpha, 2)*P2*(( 
      -6.0)+15.0*pow(P1, 2)+pow(P2, 2))+P2*(8.0+(-6.0)*pow(P1, 2)+(-6.0)*pow(P2, 2)+(-5.0)* 
      pow(beta, 2)*(2.0+9.0*pow(P1, 2)+(-5.0)*pow(P2, 2))))*pow(R, 3.0)+384.0*(alpha+(-1.0)* 
      alpha*pow(P1, 2)+5.0*beta*P1*P2+4.0*alpha*pow(P2, 2))*pow(R, 4));

    return I;
}

