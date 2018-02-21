/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/base_dearth.h"
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

using namespace smartastro;
using namespace smartastro::astrodynamics;


template<class T>
base_dearth<T>::base_dearth(const std::string &name, const double &L_scale, const double &t_scale, const int &n, const std::vector<int> &flags, const std::vector<T> &p, const std::vector<T> &p2, const std::vector< std::vector<double> > &F10dot7): base_astrodynamics<T>(name, L_scale, t_scale), m_max_degree_Earth_gravity(n), m_flags(flags), m_drag(p), m_other(p2), m_F10dot7(F10dot7){

    if(m_max_degree_Earth_gravity>constants::n_Earth_gravity)
        smartastro_throw("BASE_DEARTH: Earth harmonics exceeding available data");
    if(m_max_degree_Earth_gravity<0)
        smartastro_throw("BASE_DEARTH: Earth harmonics order must be a natural integer");      
    if((m_flags.size()<4)||(m_flags.size()>5))
        smartastro_throw("BASE_DEARTH: there must be at least 4 or 5 flags for perturbations");
    if(m_drag.size() >= 4)
    {
        if((-m_drag[3] > 0.0) || (m_drag[3] > 9.0))
            smartastro_throw("BASE_DEARTH: geomagnetic index must be between 0 and 9");
    }
    if(m_F10dot7.size() != 2)
        smartastro_throw("BASE_DEARTH: parameters for solar flux variations needs to be a vector of 2 vectors");  
    if(m_F10dot7[0].size() > 0)
    {
        double ans = m_F10dot7[0][0];
        for(unsigned int k = 1; k < m_F10dot7[0].size(); k++)
        {
            if(m_F10dot7[0][k] < ans)
                smartastro_throw("BASE_DEARTH: dates for variable solar flux must be ordered");
            ans = m_F10dot7[0][k];
        }
    }

    std::vector<double> C((m_max_degree_Earth_gravity+1)*(m_max_degree_Earth_gravity+1),0.0),S((m_max_degree_Earth_gravity+1)*(m_max_degree_Earth_gravity+1),0.0);
    m_C=C; m_S=S;

    /** computing non-normalized coefficients for spherical harmonics **/
    double prod;
    for(int m=0; m<m_max_degree_Earth_gravity+1; m++)
    {
        for(int n=0; n<m_max_degree_Earth_gravity+1; n++)
        {
            prod=1.0;
            if(n>=m)
            {
                for(int k=n-m+1; k<n+m+1; k++)
                    prod*=double(k);   
            }
            if(m==0)
            {
                m_C[(m_max_degree_Earth_gravity+1)*m+n]=constants::C_Earth_norm[(constants::n_Earth_gravity+1)*m+n]*sqrt((2.0*double(n)+1.0)/prod);
                m_S[(m_max_degree_Earth_gravity+1)*m+n]=constants::S_Earth_norm[(constants::n_Earth_gravity+1)*m+n]*sqrt((2.0*double(n)+1.0)/prod);
            }
            else if(n>=m)
            {
                m_C[(m_max_degree_Earth_gravity+1)*m+n]=constants::C_Earth_norm[(constants::n_Earth_gravity+1)*m+n]*sqrt(2.0*(2.0*double(n)+1.0)/prod);
                m_S[(m_max_degree_Earth_gravity+1)*m+n]=constants::S_Earth_norm[(constants::n_Earth_gravity+1)*m+n]*sqrt(2.0*(2.0*double(n)+1.0)/prod);
            }
        }
    }

    m_max_degree_Earth_magnetic=constants::n_Earth_magnetic;

};


template < class T >
base_dearth<T>::~base_dearth()
{
      
}


template<class T>
int base_dearth<T>::harmonics(const std::vector<T> &pos, const int &n_max, const double &Re, const std::vector<double> &C, const std::vector<double> &S, std::vector<T> &grad){

  /** sanity checks **/
    if(pos.size()!=3)
        smartastro_throw("HARMONICS: position must be a 3 dimensional vector");
    if(grad.size()!=3)
        smartastro_throw("HARMONICS: output for potential's gradient must be a 3 dimensional vector");    
    if(C.size()!=(n_max+1)*(n_max+1))
        smartastro_throw("HARMONICS: number of C coefficients must match asked expansion order");    
    if(S.size()!=(n_max+1)*(n_max+1))
        smartastro_throw("HARMONICS: number of S coefficients must match asked expansion order");    

    /** declarations **/
    std::vector<T> V, W;

    T zero=0.0*pos[0];
    for(int k=0; k<(n_max+2)*(n_max+2); k++)
    {
       V.push_back(zero);
       W.push_back(zero);
    }

    T r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]); // distance to Earth center
    T x0=pos[0]; // 1st component of position
    T y0=pos[1]; // 2nd component of position
    T z0=pos[2]; // 3rd component of position
    T K=Re/(r*r); // multiplicative factor

    /* initializing harmonics terms */
    /** n=0, m=0 **/
    V[(n_max+2)*0+0] = Re/r;
    W[(n_max+2)*0+0] = 0.0;
    /** n=1, m=0 **/
    V[(n_max+2)*0+1] = z0*K*V[(n_max+2)*0+0];
    W[(n_max+2)*0+1] = 0.0;

    /* computing harmonics terms for n=m */
    for(int n=1; n<=n_max+1; n++)
    {
        V[(n_max+2)*n+n]=(2.0*double(n)-1.0)*K*(x0*V[(n_max+2)*(n-1)+n-1]-y0*W[(n_max+2)*(n-1)+n-1]);
        W[(n_max+2)*n+n]=(2.0*double(n)-1.0)*K*(x0*W[(n_max+2)*(n-1)+n-1]+y0*V[(n_max+2)*(n-1)+n-1]);
    }

    /* computing the rest of the harmonics */
    for(int n=2; n<=n_max+1; n++)
    {
        for(int m=0; m<=n-1; m++)
        {
            V[(n_max+2)*m+n]=((2.0*double(n)-1.0)*z0*V[(n_max+2)*m+n-1]-(-1.0+double(n)+double(m))*Re*V[(n_max+2)*m+n-2])*K/(double(n)-double(m));
            W[(n_max+2)*m+n]=((2.0*double(n)-1.0)*z0*W[(n_max+2)*m+n-1]-(-1.0+double(n)+double(m))*Re*W[(n_max+2)*m+n-2])*K/(double(n)-double(m));
        }
    }

    /* computing gradient of potential */
    grad[0]=0.0;grad[1]=0.0;grad[2]=0.0; // initializing
    for(int n=0; n<=n_max; n++)
    {
        /* m=0 */
        grad[0]+=-C[(n_max+1)*0+n]*V[(n_max+2)*1+n+1];
        grad[1]+=-C[(n_max+1)*0+n]*W[(n_max+2)*1+n+1];
        grad[2]+=-(1.0+double(n))*C[(n_max+1)*0+n]*V[(n_max+2)*0+n+1];
        /* m>0 */
        for(int m=1; m<=n; m++)
        {
            grad[0]+=(-C[(n_max+1)*m+n]*V[(n_max+2)*(m+1)+n+1]-S[(n_max+1)*m+n]*W[(n_max+2)*(m+1)+n+1]+(1.0+double(n)-double(m))*(2.0+double(n)-double(m))*(C[(n_max+1)*m+n]*V[(n_max+2)*(m-1)+n+1]+S[(n_max+1)*m+n]*W[(n_max+2)*(m-1)+n+1]))/2.0;
            grad[1]+=(-C[(n_max+1)*m+n]*W[(n_max+2)*(m+1)+n+1]+S[(n_max+1)*m+n]*V[(n_max+2)*(m+1)+n+1]+(1.0+double(n)-double(m))*(2.0+double(n)-double(m))*(-C[(n_max+1)*m+n]*W[(n_max+2)*(m-1)+n+1]+S[(n_max+1)*m+n]*V[(n_max+2)*(m-1)+n+1]))/2.0;
            grad[2]+=-(1.0+double(n)-double(m))*(C[(n_max+1)*m+n]*V[(n_max+2)*m+n+1]+S[(n_max+1)*m+n]*W[(n_max+2)*m+n+1]);
        }
    }

    return 0;
}


template<class T> 
int base_dearth<T>::tidal_corrections(const double &jd, std::vector<double> &C, std::vector<double> &S) const{
    // NOT WORKING PROPERLY, DO NOT USE AS IT IS
    // TODO: fix this
    if(m_max_degree_Earth_gravity >= 1)
    {
        /** declarations **/
        std::vector<double> pos_sun(3), pos_moon(3), pos(3), M(9);
        double ratio_sun, ratio_moon, eps, r_sun, lambda_sun, r_moon, lambda_moon, phi_moon, phi_sun, re, factor_sun, factor_moon, factor;
        double prod, kn, k2m_plus;
        int n_max;

        /** constants **/
        ratio_sun=constants::mu_sun/constants::mu_earth;
        ratio_moon=constants::mu_moon/constants::mu_earth;
        eps=constants::obliquity_ecliptic_2000;
        re = constants::R_earth / this->m_L_scale;

        /** computing rotation matrix **/
        astrocore::conversion_frames::inertial_to_bf(jd,M);

        /** Sun's ephemerides **/
        Sun(jd, r_sun, lambda_sun);
        r_sun *= 1.0e3 / m_L_scale;

        /** computing Sun's position in inertial frame **/
        pos[0]=r_sun*cos(lambda_sun);
        pos[1]=r_sun*sin(lambda_sun)*cos(eps);
        pos[2]=r_sun*sin(lambda_sun)*sin(eps);

        /** computing Sun's position in Earth-fixed frame **/
        pos_sun[0]=M[0]*pos[0]+M[3]*pos[1]+M[6]*pos[2];
        pos_sun[1]=M[1]*pos[0]+M[4]*pos[1]+M[7]*pos[2];
        pos_sun[2]=M[2]*pos[0]+M[5]*pos[1]+M[8]*pos[2];

        /** computing Sun's geocentric coordinates **/
        lambda_sun = atan2(pos_sun[1], pos_sun[0]);
        phi_sun = atan2(pos_sun[2], sqrt(pos_sun[0]*pos_sun[0]+pos_sun[1]*pos_sun[1]));       

        /** computing Moon's position from ephemerides **/
        Moon(jd, pos);

        /** computing Moon's position in Earth-fixed frame **/
        pos_moon[0]=M[0]*pos[0]+M[3]*pos[1]+M[6]*pos[2];
        pos_moon[1]=M[1]*pos[0]+M[4]*pos[1]+M[7]*pos[2];
        pos_moon[2]=M[2]*pos[0]+M[5]*pos[1]+M[8]*pos[2];

        /** computing Moon's geocentric coordinates **/
        r_moon = sqrt(pos_moon[0]*pos_moon[0]+pos_moon[1]*pos_moon[1]+pos_moon[2]*pos_moon[2]);
        lambda_moon = atan2(pos_moon[1], pos_moon[0]);
        phi_moon = atan2(pos_moon[2], sqrt(pos_moon[0]*pos_moon[0]+pos_moon[1]*pos_moon[1]));

        n_max = 2;
        if(m_max_degree_Earth_gravity >= 3)
            n_max = 3;   

        for(int n = 2; n <= n_max; n++)
        {   

            /** getting Love number **/ 
            if(n == 2)
                kn = 0.300;
            else if(n == 3)
                kn = 0.093;    

            for(int m = 0; m < n + 1; m++)
            {

                prod = 1.0;
                for(int l = 1; l < 2 * m + 1; l++)
                    prod *= double(n - m + l);

                // factor = 4.0 * kn * sqrt(double(n + 2)) * pow(prod, -1.5); // Montenbruck version
                factor = kn / prod; // IERS version
                if(m > 0)
                    factor *= 2.0;

                /** adding Sun's contribution **/
                factor_sun = factor * ratio_sun * pow(re / r_sun, n + 1) * pow(-1.0, m) * smartmath::Legendre(n, m, sin(phi_sun));
                C[(m_max_degree_Earth_gravity+1)*m+n] += factor_sun * cos(double(m) * lambda_sun);
                S[(m_max_degree_Earth_gravity+1)*m+n] += factor_sun * sin(double(m) * lambda_sun);

                /** adding Moon's contribution **/
                factor_moon = factor * ratio_moon * pow(re / r_moon, n + 1) * pow(-1.0, m) * smartmath::Legendre(n, m, sin(phi_moon));
                C[(m_max_degree_Earth_gravity+1)*m+n] += factor_moon * cos(double(m) * lambda_moon);
                S[(m_max_degree_Earth_gravity+1)*m+n] += factor_moon * sin(double(m) * lambda_moon);

                if((n == 2) && (m_max_degree_Earth_gravity >= 4)) 
                { // effects of tides n=2 on geopotential term for n=4

                    /** getting Love number **/ 
                    if(m == 0)
                        k2m_plus = -0.00087;
                    else if(m == 1)
                        k2m_plus = -0.00079;
                    else if(m == 2)
                        k2m_plus = -0.00057;

                    factor = k2m_plus / prod;
                    if(m > 0)
                        factor *= 2.0;

                    factor_sun = factor * ratio_sun * pow(re / r_sun, n + 1) * pow(-1.0, m) * smartmath::Legendre(n, m, sin(phi_sun));
                    C[(m_max_degree_Earth_gravity+1)*m+4] += factor_sun * cos(double(m) * lambda_sun);
                    S[(m_max_degree_Earth_gravity+1)*m+4] += factor_sun * sin(double(m) * lambda_sun);

                    factor_moon = factor * ratio_moon * pow(re / r_moon, n + 1) * pow(-1.0, m) * smartmath::Legendre(n, m, sin(phi_moon));
                    C[(m_max_degree_Earth_gravity+1)*m+4] += factor_moon * cos(double(m) * lambda_moon);
                    S[(m_max_degree_Earth_gravity+1)*m+4] += factor_moon * sin(double(m) * lambda_moon);
                }

            }

            // C[(m_max_degree_Earth_gravity+1)*0+2] -= -4.108e-9 * sqrt(5.0); // removing time independent term counted twice for zero tide gravitational model, not needed for EGM96
        }
    }

    return 0;
}


template<class T> 
int base_dearth<T>::Earth_gravity(const double &jd, const std::vector<T> &x, std::vector<T> &acc) const{

    /** sanity checks **/
    if(x.size()<3)
        smartastro_throw("EARTH_GRAVITY: state must have at least 3 components");
    if(acc.size()!=3)
        smartastro_throw("EARTH_GRAVITY: acceleration must be a 3 dimensional vector");

    /** declarations **/
    std::vector<T> aux, pos;
    std::vector<double> C = m_C, S = m_S;
    std::vector<double> M(9);
    double mu = constants::mu_earth*pow(m_T_scale,2)/pow(this->m_L_scale,3); // Earth gravitational constant
    double Re=constants::R_earth/this->m_L_scale;  

    /** computing rotation matrix **/
    astrocore::conversion_frames::inertial_to_bf(jd,M);

    /** computing position in Earth-fixed frame **/
    pos.push_back(M[0]*x[0]+M[3]*x[1]+M[6]*x[2]);
    pos.push_back(M[1]*x[0]+M[4]*x[1]+M[7]*x[2]);
    pos.push_back(M[2]*x[0]+M[5]*x[1]+M[8]*x[2]);

    /** computing acceleration in true of date **/
    if(m_flags.size() >= 5) // solid tides
    {
        if(m_flags[4] == 1) 
            tidal_corrections(jd,C,S);
    }
    harmonics(pos,m_max_degree_Earth_gravity,Re,C,S,acc); 
    acc[0]= acc[0]*mu/(Re*Re);
    acc[1]= acc[1]*mu/(Re*Re);
    acc[2]= acc[2]*mu/(Re*Re);

    /** computing acceleration in inertial frame **/
    aux=acc;  
    acc[0]=M[0]*aux[0]+M[1]*aux[1]+M[2]*aux[2];
    acc[1]=M[3]*aux[0]+M[4]*aux[1]+M[5]*aux[2];
    acc[2]=M[6]*aux[0]+M[7]*aux[1]+M[8]*aux[2];

    return 0;
}


template<class T>
int base_dearth<T>::geodetic(const std::vector<T> &pos, T &Z, T &phi) const{

    /** sanity checks **/
    if(pos.size()!=3)
        smartastro_throw("GEODETIC: position must be a 3 dimensional vector");

    /** computing altitude iteratively */
    double Re=constants::R_earth/this->m_L_scale;
    double e2=1.0-(1.0-constants::f_earth)*(1.0-constants::f_earth);
    T dz=e2*pos[2];
    T sinphi=(pos[2]+dz)/sqrt(pos[0]*pos[0]+pos[1]*pos[1]+(pos[2]+dz)*(pos[2]+dz));   
    T N=Re/sqrt(1.0-e2*sinphi*sinphi);
    T N0=N, dz0=dz, sinphi0=sinphi;
    for(int i=0; i<2; i++)
    {
        dz0=dz;N0=N;sinphi0=sinphi; // saving previous step
        sinphi=(pos[2]+dz0)/sqrt(pos[0]*pos[0]+pos[1]*pos[1]+(pos[2]+dz0)*(pos[2]+dz0));
        N=Re/sqrt(1.0-e2*sinphi0*sinphi0);       
        dz=N0*e2*sinphi0;
    }
    T h=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+(pos[2]+dz)*(pos[2]+dz))-N; // scaled altitude
    Z=h*this->m_L_scale*1.0e-3; // altitude in km
    phi=atan2(pos[2]+dz,sqrt(pos[0]*pos[0]+pos[1]*pos[1])); // geodetic latitude in rad    

    return 0;
}


template<class T>
int base_dearth<T>::density(const double &jd, const T &Z, const T &phi, const T &alpha, T &rho) const{

    T zero=0.0*Z;

    double pi,psi,tsa;
    std::vector<double> M(9), pos_sun(3), inter(3);
    pi=constants::pi; 

    /** computing standard temperature at given position and time */
    /** mean solar activity */
    T F_107=zero; // mean solar flux
    if(m_drag.size()>=3)
        F_107+=m_drag[2];
    else
        F_107+=155.0;
    T F = F_107;
    if(m_F10dot7[1].size() > 0)
        var_flux(jd, F_107, F);
    F_107 = F;
    if((-F_107 + 70.0 > 0.0) || (F_107 - 240.0 > 0.0))
        smartastro_throw("DENSITY: mean F10.7 should always be between 70 and 240");
    T Tinf=379.0+3.24*F_107; // standard exospheric temperature

    /** diurnal variation */
    double eps=constants::obliquity_ecliptic_2000;
    double r_sun, lambda_sun;
    Sun(jd,r_sun,lambda_sun);
    r_sun *= 1.0e3 / m_L_scale;
    inter[0]=r_sun*cos(lambda_sun);
    inter[1]=r_sun*sin(lambda_sun)*cos(eps);
    inter[2]=r_sun*sin(lambda_sun)*sin(eps);
    pos_sun[0]=M[0]*inter[0]+M[3]*inter[1]+M[6]*inter[2];
    pos_sun[1]=M[1]*inter[0]+M[4]*inter[1]+M[7]*inter[2];
    pos_sun[2]=M[2]*inter[0]+M[5]*inter[1]+M[8]*inter[2];
    double alpha_sun=atan2(pos_sun[1],pos_sun[0]);
    double delta_sun=atan2(pos_sun[2],sqrt(pos_sun[0]*pos_sun[0]+pos_sun[1]*pos_sun[1]));
    T HA=(alpha-alpha_sun)*180.0/pi;
    T tau=(HA-37.0+6.0*sin((HA+43.0)*(pi/180.0)))*pi/180.0;
    T cos_2theta=cos(phi+delta_sun);
    T cos_2eta=cos(phi-delta_sun);
    Tinf*=1.0+0.3*(pow(0.5*(1.0-cos_2theta),1.1)+(pow(0.5*(1.0+cos_2eta),1.1)-pow(0.5*(1.0-cos_2theta),1.1))*pow(cos(tau/2.0),3));

    /** geomagnetic contribution */
    T Kp=zero; // geomagnetic index
    if(m_drag.size()>=4)
        Kp+=m_drag[3];
    else
        Kp+=4.0; // default value
    T f=0.5*(1.0+tanh(0.04*(Z-350.0)));
    T expKp=exp(Kp);
    Tinf+=f*(28.0*Kp+0.03*expKp)+(1.0-f)*(14.0*Kp+0.02*expKp);

    /** computing density from Jacchia 71 (bi-polynomial representation) */
    std::vector<T> var;
    var.push_back(Tinf);
    var.push_back(Z);
    T log_rho=density_log(var);

    /** computing density corrections at given altitude */  
    log_rho+=(1.0-f)*(0.012*Kp+1.2e-5*expKp); // geomagnetic term

    /** semi-annual variation */   
    psi=(jd-constants::mjd_zero_in_jd-36204.0)/365.2422;
    tsa=psi+0.09544*(pow(0.5+0.5*sin(2.0*pi*psi+6.035),1.65)-0.5);
    log_rho+=(0.06328+5.876e-7*pow(Z,2.331))*exp(-0.002868*Z)
     *(0.02835+(0.3817+0.17829*sin(2.0*pi*tsa+4.137))*sin(4.0*pi*tsa+4.259));

    /** seasonal-latitudinal variation (meaningful for altitude < 120km) NOT ALWAYS WORKING WITH INTRUSIVE */
    // T Delta_Z=Z-90.0;
    // log_rho+=0.014*Delta_Z*exp(-0.0013*Delta_Z*Delta_Z)*sin(2.0*pi*psi+1.72)*signed_square(sin(phi));

    /** Helium density variation ADDING THIS TERM CAN CREATE NEGATIVE VALUES FOR DENSITY*/
    T d_rho=zero;
    T log_n_He=He_number_density_log(var);
    double sign_delta=1.0;   
    if(delta_sun<0.0)
        sign_delta=-1.0;
    T d_log_n_He=0.65*sign_delta*(delta_sun/eps)*(pow(sin(pi/4.0-sign_delta*phi/2.0),3)-0.35355);
    d_rho=exp(log_n_He*log(10.0))*(4.0026/(6.02214129e23))*(exp(d_log_n_He*log(10.0))-1.0);   
    //std::cout << d_log_n_He << std::endl;

    rho=exp(log_rho*log(10.0))+0.0*d_rho;// final atmospheric dentisity

    return 0;
}


template<class T>
int base_dearth<T>::Sun(const double &jd, double &r, double &lambda){

    double pi,angle,M,TT;
    pi=constants::pi;
    TT=(jd-2451545.0)/36525.0;
    angle=282.9400*pi/180.0;
    M=(357.5256+35999.049*TT)*pi/180.0;
    r=(149.619-2.499*cos(M)-0.021*cos(2.0*M))*1.0e6;
    lambda=angle+M+(6892.0*sin(M)+72.0*sin(2.0*M))*pi/(180.0*3600.0);

    return 0;
}


template<class T>
int base_dearth<T>::moon_pos(const double &jd, std::vector<double> &pos_moon){

    /** sanity checks **/
    if(pos_moon.size()!=3)
        smartastro_throw("MOON_ECL: position vector must be of dimension 3");

    /* constants */
    double pi=constants::pi;

    /** precomputations **/
    double TT=(jd-2451545.0)/36525.0;
    double L0=(218.31617+481267.88088*TT-1.3972*TT)*pi/180.0;
    double l=(134.96292+477198.86753*TT)*pi/180.0;
    double lp=(357.52543+35999.04944*TT)*pi/180.0;
    double F=(93.27283+483202.01873*TT)*pi/180.0;
    double D=(297.85207+445267.11135*TT)*pi/180.0;
    double r_moon=(385000.0-20905.0*cos(l)-3699.0*cos(2.0*D-l)-2956.0*cos(2.0*D)-570.0*cos(2.0*l)+246.0*cos(2.0*l-2.0*D)-205.0*cos(lp-2.0*D)-171.0*cos(l+2.0*D)-152.0*cos(l+lp-2.0*D));
    double lambda_moon=L0+(22640.0*sin(l)+769.0*sin(2.0*l)-4586.0*sin(l-2.0*D)+2370.0*sin(2.0*D)-668.0*sin(lp)-412.0*sin(2.0*F)-212.0*sin(2.0*l-2.0*D)-206.0*sin(l+lp-2.0*D)+192.0*sin(l+2.0*D)-165.0*sin(lp-2.0*D)+148.0*sin(l-lp)-125.0*sin(D)-110.0*sin(l+lp)-55.0*sin(2.0*F-2.0*D))*pi/(180.0*3600.0);
    double beta=(18520.0*sin(F+lambda_moon-L0+412.0*sin(2.0*F)*pi/(180.0*3600.0)+541.0*sin(lp)*pi/(180.0*3600.0))-526.0*sin(F-2.0*D)+44.0*sin(l+F-2.0*D)-31.0*sin(-l+F-2.0*D)-25.0*sin(-2.0*l+F)-23.0*sin(lp+F-2.0*D)+21.0*sin(-l+F)+11.0*sin(-lp+F-2.0*D))*pi/(180.0*3600.0);

    /** computing Moon's position with x-y plane being the eccliptic **/
    pos_moon[0]=r_moon*cos(lambda_moon)*cos(beta);
    pos_moon[1]=r_moon*sin(lambda_moon)*cos(beta);
    pos_moon[2]=r_moon*sin(beta);

    return 0;
}


template<class T>
int base_dearth<T>::Moon(const double &jd, std::vector<double> &pos_moon) const{

    /** sanity checks **/
    if(pos_moon.size()!=3)
        smartastro_throw("MOON: position vector must be of dimension 3");

    std::vector<double> pos(3);
    double eps=constants::obliquity_ecliptic_2000;

    /** computing Moon's position from ephemerides in ecliptic plane **/
    moon_pos(jd, pos);

    /** computing Moon's scaled position in reference inertial frame **/
    pos_moon[0]=pos[0]*1.0e3 / this->m_L_scale;
    pos_moon[1]=(cos(-eps)*pos[1]+sin(-eps)*pos[2])*1.0e3 / this->m_L_scale;
    pos_moon[2]=(-sin(-eps)*pos[1]+cos(-eps)*pos[2])*1.0e3 / this->m_L_scale;

    return 0;
}


template<class T>
int base_dearth<T>::lunisolar(const double &jd, const std::vector<T> &x, std::vector<T> &acc) const{

    /** sanity checks **/
    if(x.size()<3)
        smartastro_throw("LUNISOLAR: state must have at least 3 components");
    if(acc.size()!=3)
        smartastro_throw("LUNISOLAR: acceleration must be a 3 dimensional vector");

    /** declarations **/
    std::vector<double> pos_sun(3), pos_moon(3), pos(3);
    double eps,r_sun,mu_sun,lambda_sun,r_moon,mu_moon;

    /** constants **/
    mu_sun=constants::mu_sun*pow(m_T_scale,2)/pow(this->m_L_scale,3);
    mu_moon=constants::mu_moon*pow(m_T_scale,2)/pow(this->m_L_scale,3);
    eps=constants::obliquity_ecliptic_2000;

    /** Sun's ephemerides **/
    Sun(jd,r_sun,lambda_sun);
    r_sun *= 1.0e3 / m_L_scale;

    /** computing Sun's position in inertial frame **/
    pos_sun[0]=r_sun*cos(lambda_sun);
    pos_sun[1]=r_sun*sin(lambda_sun)*cos(eps);
    pos_sun[2]=r_sun*sin(lambda_sun)*sin(eps);

    /** computing relative distance **/
    T r_rel=sqrt((pos_sun[0]-x[0])*(pos_sun[0]-x[0])+(pos_sun[1]-x[1])*(pos_sun[1]-x[1])+(pos_sun[2]-x[2])*(pos_sun[2]-x[2]));

    /** computing solar perturbation in inertial frame **/
    for(int i=0; i<3; i++)
        acc[i]=mu_sun*((pos_sun[i]-x[i])/pow(r_rel,3)-pos_sun[i]/pow(r_sun,3));

    /** computing Moon's position from ephemerides **/
    Moon(jd, pos_moon);

    /** computing Moon's distance to Earth **/
    r_moon = sqrt(pos_moon[0]*pos_moon[0]+pos_moon[1]*pos_moon[1]+pos_moon[2]*pos_moon[2]);

    /** computing relative distance **/
    r_rel=sqrt((pos_moon[0]-x[0])*(pos_moon[0]-x[0])+(pos_moon[1]-x[1])*(pos_moon[1]-x[1])+(pos_moon[2]-x[2])*(pos_moon[2]-x[2]));

    /** adding lunar perturbation in inertial frame **/
    for(int i=0; i<3; i++)
        acc[i]+=mu_moon*((pos_moon[i]-x[i])/pow(r_rel,3)-pos_moon[i]/pow(r_moon,3));

    return 0;
}


template<class T>
T base_dearth<T>::shadow(const double &jd, const std::vector<T> &x) const{

    /** sanity checks **/
    if(x.size()<3)
        smartastro_throw("SHADOW: state must have at least 3 components");

    T zero = 0.0 * x[0];
    T output = zero + 1.0; // default value is 1

    double eps=constants::obliquity_ecliptic_2000; 
    double r_sun, lambda_sun; 
    double Re=constants::R_earth/this->m_L_scale;
  
    /** calling Sun's ephemerides **/
    Sun(jd,r_sun,lambda_sun);
    r_sun *= 1.0e3 / m_L_scale;

    /** computing Sun's position in inertial frame **/
    std::vector<double> dir_sun(3);
    dir_sun[0]=cos(lambda_sun);
    dir_sun[1]=sin(lambda_sun)*cos(eps);
    dir_sun[2]=sin(lambda_sun)*sin(eps);

    /** precomputations **/
    T s0=-(x[0] * dir_sun[0] + x[1] * dir_sun[1] + x[2] * dir_sun[2]);    
    T r_sq = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

    /** computing eclipse conditions **/    
    if(m_flags[2] == 1)
    {
        /* conical shadow */
        double Rs = constants::R_sun / this->m_L_scale;      
        double sinf1 = (Rs+Re) / r_sun;
        double sinf2 = (Rs-Re) / r_sun;
        T l1 = (s0 * sinf1 + Re) / sqrt(1.0 - sinf1 * sinf1);
        T l2 = (s0 * sinf2 - Re) / sqrt(1.0 - sinf2 * sinf2);
        T l = sqrt(r_sq - s0 * s0);
        if(s0 > 0.0)
        { 
            T l2_abs = sqrt(l2 * l2);
            if(l2_abs - l > 0.0) // umbra    
                output = 0.0;           
            else if(l1 - l > 0.0) // penumbra
            {
                // output = 0.5 * (1.0 + tanh((l - 0.5 * (l1 + l2_abs)) / (2.0e-3 * l2_abs))); // sigmoid (old in-house version)
                std::vector<T> pos_rel(3, zero);
                for(int i = 0; i < 3; i++)
                    pos_rel[i] = r_sun * dir_sun[i] - x[i];
                T dist_rel = sqrt(pos_rel[0] * pos_rel[0] + pos_rel[1] * pos_rel[1] + pos_rel[2] * pos_rel[2]);
                T a = asin(Rs / dist_rel);
                T a2 = a * a;
                T r = sqrt(r_sq);
                T b = asin(Re / r);
                T c = acos(-(x[0] * pos_rel[0] + x[1] * pos_rel[1] + x[2] * pos_rel[2]) / (r * dist_rel));
                T x = (c * c + a2 - b * b) / (2.0 * c);
                T y = sqrt(a2 - x * x);
                output = 1.0 - (acos(x / a) + (b * b * acos((c - x) / b) - c * y) / a2) / constants::pi;
            }
        }
    }
    else if(m_flags[2] == 2)
    {
    /* cylindrical shadow */    
        output = 1.0; 
        if((s0 > 0.0)&&(Re * Re - r_sq + s0 * s0 > 0.0))
           output = 0.0;
    }
    else
        smartastro_throw("SHADOW: input for shadow model incorrect");

    return output;
}


template<class T>
int base_dearth<T>::Gaussian_coeff(const double &jd, const int &n_max, std::vector<double> &g, std::vector<double> &h) const{

    if(g.size()!=(n_max+1)*(n_max+1))
        smartastro_throw("GAUSSIAN_COEFF: number of g coefficients must match asked expansion order");    
    if(h.size()!=(n_max+1)*(n_max+1))
        smartastro_throw("GAUSSIAN_COEFF: number of h coefficients must match asked expansion order"); 

    std::vector<double> date(6,0.0);
    date[0]=2015.0;date[1]=1.0;date[2]=1.0;
    double jd2015;
    astrocore::conversion_time::date2jd(date,jd2015); 

    double prod1, g_year, h_year;
    for(int m=0; m<n_max+1; m++)
    {
        for(int n=0; n<n_max+1; n++)
        {
          prod1=1.0;
          if(n>=m)
          {
            for(int k=n-m+1; k<n+m+1; k++)
              prod1*=double(k);      
          }            
          g_year=constants::g_Earth_2015[(constants::n_Earth_magnetic+1)*m+n]+((jd-jd2015)/365.25)*constants::g_dot_Earth_2015[(constants::n_Earth_magnetic+1)*m+n];
          h_year=constants::h_Earth_2015[(constants::n_Earth_magnetic+1)*m+n]+((jd-jd2015)/365.25)*constants::h_dot_Earth_2015[(constants::n_Earth_magnetic+1)*m+n];           
          if(m==0)
          {
            g[(n_max+1)*m+n]=1.0e-9*pow(m_T_scale,2)*g_year*sqrt(1.0/prod1);
            h[(n_max+1)*m+n]=1.0e-9*pow(m_T_scale,2)*h_year*sqrt(1.0/prod1);
          }
          else
          {
            g[(n_max+1)*m+n]=1.0e-9*pow(m_T_scale,2)*g_year*sqrt(2.0/prod1);
            h[(n_max+1)*m+n]=1.0e-9*pow(m_T_scale,2)*h_year*sqrt(2.0/prod1);
          }
          //std::cout << n << " "  << m << " "  << g_year << " " << h_year << std::endl;
        }
    }  

    return 0;
}


template<class T>
int base_dearth<T>::IRI90(const double &jd, const T &Z, const T &phi, const T &t, T &Te) const{

    double pi=constants::pi;
    T theta=pi/2.0-phi; // co magnetic latitude
    T phi_deg=phi*180.0/pi;

    /* Computing the reference values for Booker's profile */
    double DTE[7] = {0.0, 5.0, 5.0, 10.0, 20.0, 20.0, 0.0};

    T zero=0.0*Z;
    std::vector<T> ST;
    for(int i = 0; i < 7; i++)    
      ST.push_back(zero);

    std::vector<T> AHH=ST, ATE=ST;
    T AHH1d=AHH[1], AHH1n=AHH[1];
    T ATE1d=ATE[1], ATE1n=ATE[1];
    T ATE4d=ATE[4], ATE4n=ATE[4];

    AHH[0]=120.0;
    AHH1d=210.0 + 60.0 * exp(- pow((phi_deg/22.4), 2)); // day-time value
    AHH1n=150.0; // night-time value
    AHH[1]=AHH1n + (AHH1d-AHH1n)/(1.0+1.0/exp((t*12.0/pi-6.0)/1.0)) + (AHH1n-AHH1d)/(1.0+1.0/exp((t*12.0/pi-18.0)/1.0));
    AHH[2]= 300.0;
    AHH[3]=400.0;
    AHH[4]=600.0;
    AHH[5]=1400.0;
    AHH[6]=3000.0;   

    ATE[0] = 365.0;
    ATE1d = 1500.0 + 800.0 * exp(- pow((phi_deg/33.0), 2)); // day-time value
    ATE1n = 555.0; // night-time value
    ATE[1]=ATE1n + (ATE1d-ATE1n)/(1.0+1.0/exp((t*12.0/pi-6.0)/1.0)) + (ATE1n-ATE1d)/(1.0+1.0/exp((t*12.0/pi-18.0)/1.0));
    T inter0=exp(phi_deg/11.35);
    T aux0=inter0/pow(1.0 + inter0,2);    
    ATE4d = 2900.0 - 5600.0 * aux0; // 3pm value
    ATE4n = 839.0 + (2000.0-839.0)/(1.0+1.0/exp((sqrt(phi_deg*phi_deg)-45.0)/5.0)); // 3am value    
    ATE[4]=ATE4n + (ATE4d-ATE4n)/(1.0+1.0/exp((t*12.0/pi-6.0)/1.0)) + (ATE4n-ATE4d)/(1.0+1.0/exp((t*12.0/pi-18.0)/1.0));    

    ATE[2] = exp((0.3100e1-0.3215e-2*cos(theta)-sin(theta)*(0.2791e-1*sin(t)-0.1853e0*cos(t)))*log(10.0)); // first order value at equinox
    ATE[3] = exp((0.3136e1+0.6498e-2*cos(theta)-sin(theta)*(0.6555e-1*sin(t)-0.2115e0*cos(t)))*log(10.0)); // first order value at equinox
    ATE[5] = exp((0.3372e1+0.1006e-1*cos(theta)-sin(theta)*(-0.4713e-1*sin(t)-0.2836e0*cos(t)))*log(10.0)); // first order value at equinox 
    ATE[6] = exp((0.3574e1+0.0*cos(theta)-sin(theta)*(-0.5321e-1*sin(t)-0.1768e0*cos(t)))*log(10.0)); // first order value at equinox

    // ATE[2] = exp((0.3137e1+0.6796e-2*cos(theta)-sin(theta)*(0.5284e-1*sin(t)-0.2575e0*cos(t)))*log(10.0)); // first order value at june solstice
    // ATE[3] = exp((0.3144e1+0.8571e-2*cos(theta)-sin(theta)*(0.5232e-1*sin(t)-0.2019e0*cos(t)))*log(10.0)); // first order value at june solstice
    // ATE[5] = exp((0.3367e1+0.1038e-1*cos(theta)-sin(theta)*(-0.5799e-1*sin(t)-0.3110e0*cos(t)))*log(10.0)); // first order value at june solstice
    // ATE[6] = exp((0.3574e1-0.5639e-2*cos(theta)-sin(theta)*(-0.5966e-1*sin(t)-0.1783e0*cos(t)))*log(10.0)); // first order value at june solstice

    ST[0] = 0.0;
    for(int i = 1; i < 7; i++)
      ST[i] = (ATE[i] - ATE[i-1]) / (AHH[i] - AHH[i-1]);

    /* Computing the electron temperature Te */
    T base = (Z - AHH[0]) * ST[1];
    T sum = zero;
    for(int j = 1; j < 5+1; j++)
    {
      T aux1=log(1.0+exp((Z- AHH[j])/ DTE[j]));
      T aux2=log(1.0+exp((AHH[0]- AHH[j])/ DTE[j]));      
      sum += (ST[j+1] - ST[j]) * DTE[j] * (aux1-aux2);
    }    
    Te=ATE[0]+base+sum; // electron temperature

    return 0;
}

template<class T>
T base_dearth<T>::surface_potential_LEO(const T &Te, const double &shape_factor, const T &Vsc, const T &cos_vel) const{

    if((shape_factor!=0.0)&&(shape_factor!=0.5)&&(shape_factor!=1.0))
        smartastro_throw("SURFACE_POTENTIAL_LEO: shape factor can only be 0 (plate), 1/2 (cylinder) or 1 (sphere)");

    double me=constants::me, qe=constants::qe, kb=constants::kb, pi=constants::pi;

    T A_tot_surf=0.0*Te; 
    T A_proj_vel=0.0*Te; 
    if(shape_factor==0.0)
    { // flat plate
        A_tot_surf=2.0*m_drag[0]; // edge effects neglected
        A_proj_vel=A_tot_surf*cos_vel/2.0;
    }
    if(shape_factor==0.5) // cylinder
    {
        A_tot_surf=2.0*pi*sqrt(m_drag[0]/pi)*m_drag[1]; // edge effects neglected
        A_proj_vel=A_tot_surf*sqrt(1.0-cos_vel*cos_vel)/pi;
    }
    if(shape_factor==1.0)
    { // sphere
        A_tot_surf=4.0*m_drag[0];
        A_proj_vel=A_tot_surf/4.0;
    }

    return (kb/qe)*Te*log(sqrt(2.0*pi*me/(kb*Te))*Vsc*(A_proj_vel/A_tot_surf));    
    // T F_phi=0.0*Te, DF_phi=F_phi;    
    // for(int i=0; i<100; i++){
    //     F_phi=sqrt(kb*Te/(2.0*pi*me))*exp(qe*phi/(kb*Te))-Vsc*A_proj_vel/A_tot_surf;
    //     DF_phi=sqrt(kb*Te/(2.0*pi*me))*exp(qe*phi/(kb*Te))*qe/(kb*Te);
    //     phi-=F_phi/DF_phi;
    // }
}

template<class T>
int base_dearth<T>::var_flux(const double &jd, const T &F0, T &F) const{
           
    if(m_F10dot7[1].size() + 2 != m_F10dot7[0].size())
        smartastro_throw("VAR_FLUX: number of interpolated middle values plus 2 must equal the number of interpolation points");          

    unsigned int n = m_F10dot7[0].size(); // number of interpolation points
    double Delta_tau = m_F10dot7[0][n - 1] - m_F10dot7[0][0]; // period
    double tau = m_F10dot7[0][0] + double(floor((jd - m_F10dot7[0][0]) / Delta_tau)) * Delta_tau; // phasis

    std::vector<T> Fs(n, F0);
    std::vector<double> dates(n, tau);
    unsigned int k = 1; // counter 
    for(unsigned int j = 1; j < n - 1; j++)
    {
        dates[j] += m_F10dot7[0][j] - m_F10dot7[0][0];
        Fs[j] = m_F10dot7[1][j - 1];
        if((k == j) && (dates[k] < jd))
            k++;
    }
    dates[n - 1] += Delta_tau;

    F = Fs[k - 1] + ((jd - dates[k - 1]) / (dates[k] - dates[k - 1])) * (Fs[k] - Fs[k - 1]); // linear interpolation

    return 0;
}


template<class T>
double base_dearth<T>::density_log(std::vector<float> X){
  double Tinf=X[0], Z=X[1];
  double Tpower, Zpower;

  int range_T=0, range_Z=0;
      if(Z>1000.0)
        range_Z=3;
      else if(Z>500.0)
        range_Z=2;
      else if(Z>180.0)
        range_Z=1; 

    if(Tinf>850.0)
        range_T=1;

      if(Z>2500.0)
      {
        //std::cout << "altitude is too high: " << Z << "km, there should be no drag" << std::endl;
      }
      if(Z<90.0)
      {
        //std::cout << "altitude is out-of-range (too small): " << Z << "km at time " << setprecision(16) << d << " (modified Julian date)" << std::endl;
        //smartastro_throw("propagation stopped");
      }
      if((Tinf>1900.0)||(Tinf<500.0))
        std::cout << "temperature is out-of-range: " << Tinf << "K" << std::endl;   

    double log_rho=0.0;
    Z*=1.0e-3; // scaled Z
    Tinf*=1.0e-3; // scaled Tinf
    for(int i=0; i<=5; i++)
    {
      if(i==0)
        Zpower=1.0;
      else if(i==1)
        Zpower=Z;
      else 
        Zpower=pow(Z,i);

      for(int j=0; j<=4;j++)
      {
          if(j==0)
            Tpower=1.0;
          else if(j==1)
            Tpower=Tinf;
          else 
            Tpower=pow(Tinf,j);
        log_rho+=constants::coeff_density_JacchiaGill[range_Z][range_T][6*j+i]*Zpower*Tpower;
        }
    }

  return log_rho;
}
template<class T>
double base_dearth<T>::density_log(std::vector<double> X){
  double Tinf=X[0], Z=X[1];
  double Tpower, Zpower;

  int range_T=0, range_Z=0;
      if(Z>1000.0)
        range_Z=3;
      else if(Z>500.0)
        range_Z=2;
      else if(Z>180.0)
        range_Z=1; 

    if(Tinf>850.0)
        range_T=1;

      if(Z>2500.0)
      {
        //std::cout << "altitude is too high: " << Z << "km, there should be no drag" << std::endl;
      }
      if(Z<90.0)
      {
        //std::cout << "altitude is out-of-range (too small): " << Z << "km at time " << setprecision(16) << d << " (modified Julian date)" << std::endl;
        //smartastro_throw("propagation stopped");
      }
      if((Tinf>1900.0)||(Tinf<500.0))
        std::cout << "temperature is out-of-range: " << Tinf << "K" << std::endl;    

    double log_rho=0.0;
    Z*=1.0e-3; // scaled Z
    Tinf*=1.0e-3; // scaled Tinf
    for(int i=0; i<=5; i++)
    {
      if(i==0)
        Zpower=1.0;
      else if(i==1)
        Zpower=Z;
      else 
        Zpower=pow(Z,i);
      for(int j=0; j<=4;j++)
      {
          if(j==0)
            Tpower=1.0;
          else if(j==1)
            Tpower=Tinf;
          else 
            Tpower=pow(Tinf,j);
        log_rho+=constants::coeff_density_JacchiaGill[range_Z][range_T][6*j+i]*Zpower*Tpower;
       }
    }

  return log_rho;
}
template<class T>
double base_dearth<T>::density_log(std::vector<long double> X){
  double Tinf=X[0], Z=X[1];
  double Tpower, Zpower;

  int range_T=0, range_Z=0;
      if(Z>1000.0)
        range_Z=3;
      else if(Z>500.0)
        range_Z=2;
      else if(Z>180.0)
        range_Z=1;

      if(Tinf>850.0)
        range_T=1;

      if(Z>2500.0)
      {
        //std::cout << "altitude is too high: " << Z << "km, there should be no drag" << std::endl;
      }
      if(Z<90.0)
      {
        //std::cout << "altitude is out-of-range (too small): " << Z << "km at time " << setprecision(16) << d << " (modified Julian date)" << std::endl;
        //smartastro_throw("propagation stopped");
      }
      if((Tinf>1900.0)||(Tinf<500.0))
        std::cout << "temperature is out-of-range: " << Tinf << "K" << std::endl;   

    double log_rho=0.0;
    Z*=1.0e-3; // scaled Z
    Tinf*=1.0e-3; // scaled Tinf
    for(int i=0; i<=5; i++)
    {
      if(i==0)
        Zpower=1.0;
      else if(i==1)
        Zpower=Z;
      else 
        Zpower=pow(Z,i);
      for(int j=0; j<=4;j++)
      {
          if(j==0)
            Tpower=1.0;
          else if(j==1)
            Tpower=Tinf;
          else 
            Tpower=pow(Tinf,j);
       log_rho+=constants::coeff_density_JacchiaGill[range_Z][range_T][6*j+i]*Zpower*Tpower;
       }
    }

  return log_rho;
}
#ifdef ENABLE_SMARTUQ
template<class T>
smartuq::polynomial::taylor_polynomial base_dearth<T>::density_log(const std::vector<smartuq::polynomial::taylor_polynomial> &X){
    smartuq::polynomial::taylor_polynomial Tinf=X[0], Z=X[1];
    smartuq::polynomial::taylor_polynomial Tpower=Tinf, Zpower=Z;

    int range_T=0, range_Z=0;
    if(Z.get_coeffs()[0]>1000.0)
        range_Z=3;
    else if(Z.get_coeffs()[0]>500.0)
        range_Z=2;
    else if(Z.get_coeffs()[0]>180.0)
        range_Z=1;

    if(Tinf.get_coeffs()[0]>850.0)
        range_T=1;

    if(Z.get_coeffs()[0]>2500.0)
        std::cout << "nominal altitude is too high: " << Z.get_coeffs()[0] << "km, there should be no drag" << std::endl;
    if(Z.get_coeffs()[0]<90.0)
    {
        std::cout << "nominal altitude is out-of-range (too small): " << Z.get_coeffs()[0] << std::endl;
        smartastro_throw("propagation stopped"); 
    }
    if((Tinf.get_coeffs()[0]>1900.0)||(Tinf.get_coeffs()[0]<500.0))
        std::cout << "nominal temperature is out-of-range: " << Tinf.get_coeffs()[0] << "K" << std::endl;     

    smartuq::polynomial::taylor_polynomial log_rho=0.0*Tinf;
    Z*=1.0e-3;
    Tinf*=1.0e-3;
    for(int i=0; i<=5; i++)
    {
        if(i<=0)
            Zpower=1.0;
        else if(i<=1)
            Zpower=Z;
        else 
            Zpower=pow(Z,i);

        for(int j=0; j<=4;j++)
        {
            if(j<=0)
                Tpower=1.0;
            else if(j<=1)
                Tpower=Tinf;
            else 
                Tpower=pow(Tinf,j);
            log_rho+=constants::coeff_density_JacchiaGill[range_Z][range_T][6*j+i]*Zpower*Tpower;
        }
    }

    return log_rho;

}
template<class T>
smartuq::polynomial::chebyshev_polynomial base_dearth<T>::density_log(const std::vector<smartuq::polynomial::chebyshev_polynomial> &X){
    std::vector<double> range=X[1].get_range();
    double Z_min=range[0], Z_max=range[1];
    std::vector<double> range2=X[0].get_range();
    double T_min=range2[0], T_max=range2[1];

    int range_Z=-1;
    if((T_min>500.0)&&(T_max<1900.0))
    {
        if((Z_min>90.0)&&(Z_max<180.0))
            range_Z=0;       
        else if((Z_min>180.0)&&(Z_max<500.0))
            range_Z=1;
        else if((Z_min>500.0)&&(Z_max<1000.0))
            range_Z=2;
        else if((Z_min>1000.0)&&(Z_max<2500.0))
            range_Z=3;
        else 
        {
        //std::cout << "in-between ranges: " << "altitudes spread from " << Z_min*1.0e-3 << "km to "<< Z_max*1.0e-3 << "km" << std::endl;
        }
    }

    // if((Z_max>2500.0e3)||(Z_min<90.0e3))
    // {
    //   std::cout << "altitude is out-of-range: " << "spreading from " << Z_min*1.0e-3 << "km to "<< Z_max*1.0e-3 << "km" << std::endl;
    // }
    // if((T_max>1900.0)||(T_min<500.0))
    // {
    //   std::cout << "temperature is out-of-range: " << "spreading from " << T_min << "K to " << T_max << "K" << std::endl;
    // }      

    if(range_Z>=0)
    {
        smartuq::polynomial::chebyshev_polynomial Tinf=X[0]*1.0e-3, Z=X[1]*1.0e-3;;
        smartuq::polynomial::chebyshev_polynomial Tpower=Tinf, Zpower=Z;

        smartuq::polynomial::chebyshev_polynomial log_rho=0.0*Tinf;
        smartuq::polynomial::chebyshev_polynomial log_rho_cold=log_rho, log_rho_hot=log_rho;
        for(int i=0; i<=5; i++)
        {
            if(i<=0)
                Zpower=1.0;
            else if(i<=1)
                Zpower=Z;
            else 
                Zpower=pow(Z,i);

            for(int j=0; j<=4;j++)
            {
                if(j<=0)
                    Tpower=1.0;
                else if(j<=1)
                    Tpower=Tinf;
                else 
                    Tpower=pow(Tinf,j);
                log_rho_cold+=constants::coeff_density_JacchiaGill[range_Z][0][6*j+i]*Zpower*Tpower;
                log_rho_hot+=constants::coeff_density_JacchiaGill[range_Z][1][6*j+i]*Zpower*Tpower;
            }
        }
        smartuq::polynomial::chebyshev_polynomial g_T=0.5*(1.0+tanh((Tinf*1.0e3-850.0)/1.0)); 
        log_rho=g_T*log_rho_hot+(1.0-g_T)*log_rho_cold;
        return log_rho;
    }
    else
    {
        std::vector<std::vector<double> > ranges;
        std::vector<double> limits_T, limits_Z;
        limits_T.push_back(500.0);limits_T.push_back(1900.0);
        limits_Z.push_back(90.0);limits_Z.push_back(2500.0);
        ranges.push_back(limits_T);ranges.push_back(limits_Z);
        return smartuq::polynomial::chebyshev_polynomial::approximation(static_cast <double (*)(std::vector<double>)> (&density_log_sigmoid_Tinf), X, ranges);
    }
}

#endif

template<class T>
double base_dearth<T>::density_log_sigmoid_Tinf(std::vector<double> X){
  double Tinf=X[0], Z=X[1];
  double Tpower, Zpower;

  int range_Z=0;
    if(Z>1000.0)
        range_Z=3;
    else if(Z>500.0)
        range_Z=2;
    else if(Z>180.0)
        range_Z=1;

      if((Z>2500.0)||(Z<90.0))
      {
        std::cout << "altitude is out-of-range (sigmoid version): " << Z << "km" << std::endl;
        smartastro_throw("propagation stopped"); 
      }
      if((Tinf>1900.0)||(Tinf<500.0))
      {
        std::cout << "temperature is out-of-range (sigmoid version): " << Tinf << "K" << std::endl;
        smartastro_throw("propagation stopped"); 
      }      

    double log_rho_cold=0.0, log_rho_hot=0.0;
    Z*=1.0e-3; // scaled Z
    Tinf*=1.0e-3; // scaled Tinf
    for(int i=0; i<=5; i++)
    {
      if(i<=0)
        Zpower=1.0;
      else if(i<=1)
        Zpower=Z;
      else 
        Zpower=pow(Z,i);

      for(int j=0; j<=4;j++)
      {
          if(j<=0)
            Tpower=1.0;
          else if(j<=1)
            Tpower=Tinf;
          else 
            Tpower=pow(Tinf,j);
      log_rho_cold+=constants::coeff_density_JacchiaGill[range_Z][0][6*j+i]*Zpower*Tpower;
      log_rho_hot+=constants::coeff_density_JacchiaGill[range_Z][1][6*j+i]*Zpower*Tpower;
       }
    }

    double g_T=0.5*(1.0+tanh((Tinf*1.0e3-850.0)/1.0));
    double log_rho=g_T*log_rho_hot+(1.0-g_T)*log_rho_cold;

  return log_rho;
}


template<class T>
float base_dearth<T>::He_number_density_log(std::vector<float> X){
  double Tinf=X[0], Z=X[1];
  double Tpower, Zpower;

  int range_Z=0;
      if(Z>1000.0)
        range_Z=2;
      else if(Z>500.0)
        range_Z=1;

    float log_rho=0.0;
    for(int i=0; i<=5; i++)
    {
      if(i==0)
        Zpower=1.0;
      else if(i==1)
        Zpower=Z;
      else 
        Zpower=pow(Z,i);

      for(int j=0; j<=4;j++)
      {
          if(j==0)
            Tpower=1.0;
          else if(j==1)
            Tpower=Tinf;
          else 
            Tpower=pow(Tinf,j); 
      log_rho+=constants::coeff_number_density_Helium[range_Z][6*j+i]*Zpower*Tpower;
       }
    }

  return log_rho;
}
template<class T>
double base_dearth<T>::He_number_density_log(std::vector<double> X){
  double Tinf=X[0], Z=X[1];
  double Tpower, Zpower;

  int range_Z=0;
      if(Z>1000.0)
        range_Z=2;
      else if(Z>500.0)
        range_Z=1;

    double log_rho=0.0;
    for(int i=0; i<=5; i++)
    {
      if(i==0)
        Zpower=1.0;
      else if(i==1)
        Zpower=Z;
      else 
        Zpower=pow(Z,i);
      for(int j=0; j<=4;j++)
      {
          if(j==0)
            Tpower=1.0;
          else if(j==1)
            Tpower=Tinf;
          else 
            Tpower=pow(Tinf,j);
      log_rho+=constants::coeff_number_density_Helium[range_Z][6*j+i]*Zpower*Tpower;
       }
    }

  return log_rho;
}
template<class T>
long double base_dearth<T>::He_number_density_log(std::vector<long double> X){
  double Tinf=X[0], Z=X[1];
  double Tpower, Zpower;

  int range_Z=0;
      if(Z>1000.0)
        range_Z=2;
      else if(Z>500.0)
        range_Z=1;

    long double log_rho=0.0;
    for(int i=0; i<=5; i++)
    {
      if(i==0)
        Zpower=1.0;
      else if(i==1)
        Zpower=Z;
      else 
        Zpower=pow(Z,i);

      for(int j=0; j<=4;j++)
      {
          if(j==0)
            Tpower=1.0;
          else if(j==1)
            Tpower=Tinf;
          else 
            Tpower=pow(Tinf,j); 
      log_rho+=constants::coeff_number_density_Helium[range_Z][6*j+i]*Zpower*Tpower;
       }
    }

  return log_rho;
}
#ifdef ENABLE_SMARTUQ
template<class T>
smartuq::polynomial::taylor_polynomial base_dearth<T>::He_number_density_log(const std::vector<smartuq::polynomial::taylor_polynomial> &X){
    smartuq::polynomial::taylor_polynomial Tinf=X[0], Z=X[1];
    smartuq::polynomial::taylor_polynomial Tpower=Tinf, Zpower=Z;

    int range_Z=0;
    if(Z.get_coeffs()[0]>1000.0)
        range_Z=2;
    else if(Z.get_coeffs()[0]>500.0)
        range_Z=1;

    smartuq::polynomial::taylor_polynomial log_n_He=0.0*Tinf;
    for(int i=0; i<=5; i++)
    {
        if(i<=0)
            Zpower=1.0;
        else if(i<=1)
            Zpower=Z;
        else 
            Zpower=pow(Z,i);
        for(int j=0; j<=4;j++)
        {
            if(j<=0)
                Tpower=1.0;
            else if(j<=1)
                Tpower=Tinf;
            else 
                Tpower=pow(Tinf,j);
            log_n_He+=constants::coeff_number_density_Helium[range_Z][6*j+i]*Zpower*Tpower;
        }
    }

  return log_n_He;

}

template<class T>
smartuq::polynomial::chebyshev_polynomial base_dearth<T>::He_number_density_log(const std::vector<smartuq::polynomial::chebyshev_polynomial> &X){
  std::vector<double> range=X[1].get_range();
  double Z_min=range[0], Z_max=range[1];
  std::vector<double> range2=X[0].get_range();
  double T_min=range2[0], T_max=range2[1];

  int range_Z=-1;
  if((T_min>500.0)&&(T_max<1900.0))
  {
    if((Z_min>90.0)&&(Z_max<500.0))
        range_Z=0;
    else if((Z_min>500.0)&&(Z_max<1000.0))
        range_Z=1;
    else if((Z_min>1000.0)&&(Z_max<2500.0))
        range_Z=2;   
    else 
    {
  //std::cout << "in between ranges for Helium: " << "altitudes spread from " << Z_min*1.0e-3 << "km to "<< Z_max*1.0e-3 << "km" << std::endl;
    }
  }

  if(range_Z>=0)
  {
    smartuq::polynomial::chebyshev_polynomial Tinf=X[0], Z=X[1];
    smartuq::polynomial::chebyshev_polynomial Tpower=Tinf, Zpower=Z;

    smartuq::polynomial::chebyshev_polynomial log_n_He=0.0*Tinf;
    for(int i=0; i<=5; i++)
    {
        if(i<=0)
            Zpower=1.0;
        else if(i<=1)
            Zpower=Z;
        else 
            Zpower=pow(Z,i);
        for(int j=0; j<=4;j++)
        {
            if(j<=0)
                Tpower=1.0;
            else if(j<=1)
                Tpower=Tinf;
            else 
                Tpower=pow(Tinf,j);
        log_n_He+=constants::coeff_number_density_Helium[range_Z][6*j+i]*Zpower*Tpower;
        }
    } 
    return log_n_He;
    }
  else
  {
      std::vector<std::vector<double> > ranges;
      std::vector<double> limits_T, limits_Z;
      limits_T.push_back(500.0);limits_T.push_back(1900.0);
      limits_Z.push_back(90.0);limits_Z.push_back(2500.0);
      ranges.push_back(limits_T);ranges.push_back(limits_Z);
      return smartuq::polynomial::chebyshev_polynomial::approximation(static_cast <double (*)(std::vector<double>)> (&He_number_density_log), X, ranges);
  }
}
#endif


template<class T>
float base_dearth<T>::signed_square(float X){
    float ss=X*X;
    if(X<0.0)
        ss=-ss;
    return ss;
}
template<class T>
double base_dearth<T>::signed_square(double X){
    double ss=X*X;
    if(X<0.0)
        ss=-ss;
    return ss;
}
template<class T>
long double base_dearth<T>::signed_square(long double X){
    long double ss=X*X;
    if(X<0.0)
        ss=-ss;
    return ss;
}
#ifdef ENABLE_SMARTUQ
    template<class T>
    smartuq::polynomial::taylor_polynomial base_dearth<T>::signed_square(const smartuq::polynomial::taylor_polynomial &X){
        smartuq::polynomial::taylor_polynomial ss=X*X;
        if(X.get_coeffs()[0]<0.0)
          ss=-ss;
        return ss;
    }
    template<class T>
    smartuq::polynomial::chebyshev_polynomial base_dearth<T>::signed_square(const smartuq::polynomial::chebyshev_polynomial &X){
        smartuq::polynomial::chebyshev_polynomial ss=X;
        std::vector<double> range2=X.get_range();
        double X_min=range2[0], X_max=range2[1];
        if(X_min>=0.0)
          {
            ss=X*X;  
            return ss;
          }
        else if(X_max<=0.0)
          {
            ss=-X*X;
            return ss;
          } 
        else 
        { 
        //ss=tanh(X/0.1)*X*X;return ss;
        return smartuq::polynomial::chebyshev_polynomial::approximation(static_cast <double (*)(double)> (&signed_square), X);
        }
    }
#endif


template class base_dearth<double>;
template class base_dearth<float>;
template class base_dearth<long double>;
#ifdef ENABLE_SMARTUQ
template class base_dearth<smartuq::polynomial::chebyshev_polynomial>;
template class base_dearth<smartuq::polynomial::taylor_polynomial>;
#endif