/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astro-Core/conversion_frames.h"

using namespace smartastro;
using namespace astrocore;

	bool conversion_frames::inertial_to_tod(const double &JD, std::vector<double> &M){
   
        /** sanity checks **/
        if(M.size()!=9)
            smartastro_throw("INERTIAL_TO_TOD: rotation matrix must have 9 components");
   
        /** declarations **/
        std::vector<double> P(9), N(9);
        double T,TT,TTT,z,xi,theta,dpsi,e,ep,de,pi;

        pi=constants::pi;
        T=(JD-2451545.0)/36525.0;
        TT = T*T;
        TTT = TT*T;

        /** computing precession parameters **/
        xi=(2306.2181*T+0.30188*TT+0.017998*TTT)*(1.0/3600.0)*(pi/180.0);
        theta=(2004.3109*T-0.42665*TT-0.041833*TTT)*(1.0/3600.0)*(pi/180.0);
        z=xi+(0.79280*TT+0.000205*TTT)*(1.0/3600.0)*(pi/180.0);

        /** computing precession matrix **/
        P[0]=-sin(z)*sin(xi)+cos(z)*cos(theta)*cos(xi); 
        P[1]=cos(z)*sin(xi)+sin(z)*cos(theta)*cos(xi); 
        P[2]=sin(theta)*cos(xi); 
        P[3]=-sin(z)*cos(xi)-cos(z)*cos(theta)*sin(xi); 
        P[4]=cos(z)*cos(xi)-sin(z)*cos(theta)*sin(xi);  
        P[5]=-sin(theta)*sin(xi);  
        P[6]=-cos(z)*sin(theta); 
        P[7]=-sin(z)*sin(theta); 
        P[8]=cos(theta); 

        /** computing nutation parameters **/
        e=(23.43929111+(-46.8150*T-0.00059*TT+0.001813*TTT)/3600.0)*(pi/180.0);
        //double Omegaa=(125.0+2.0/60.0+40.280/3600.0-(1934.0+8.0/60.0+10.539/3600.0)*T+7.455*T*T/3600.0+0.008*T*T*T/3600.0)*(pi/180.0);
        //de=(9.203*cos(Omegaa))*(1.0/3600.0)*(pi/180.0);
        //ep=e+de;
        //dpsi=(-17.200*sin(Omegaa))*(1.0/3600.0)*(pi/180.0);

        double l=134.0+57.0/60.0+46.733/3600.0+T*(477198.0+52.0/60.0+2.633/3600.0)+TT*(31.310/3600.0)+TTT*(0.064/3600.0);
        double lp=357.0+31.0/60.0+39.804/3600.0+T*(35999.0+3.0/60.0+1.224/3600.0)+TT*(-0.577/3600.0)+TTT*(-0.012/3600.0);
        double F=93.0+16.0/60.0+18.877/3600.0+T*(483202.0+1.0/60.0+3.137/3600.0)+TT*(-13.257/3600.0)+TTT*(0.011/3600.0);
        double D=297.0+51.0/60.0+1.307/3600.0+T*(445267.0+6.0/60.0+41.328/3600.0)+TT*(-6.891/3600.0)+TTT*(0.019/3600.0);
        double Omega=125.0+2.0/60.0+40.280/3600.0-(1934.0+8.0/60.0+10.539/3600.0)*T+7.455*TT/3600.0+TTT*(0.008/3600.0);

        dpsi=0.0;
        de=0.0;
        double phi_i;
        for(int i=0; i<constants::n_nutation; i++)
        { 
            phi_i=(constants::p_l[i]*l+constants::p_lp[i]*lp+constants::p_F[i]*F+constants::p_D[i]*D+constants::p_Omega[i]*Omega)*pi/180.0;
            dpsi+=(constants::d_psi0[i]+T*constants::d_psi1[i])*sin(phi_i);
            de+=(constants::d_eps0[i]+T*constants::d_eps1[i])*cos(phi_i);
        }
        dpsi*=(0.0001/3600.0)*pi/180.0;
        de*=(0.0001/3600.0)*pi/180.0;
        ep=e+de;

        //std::cout << dpsi << "," << dpsip << " " << de << "," << dep << std::endl;

        /** computing nutation matrix **/
        N[0]=cos(dpsi); 
        N[1]=cos(ep)*sin(dpsi); 
        N[2]=sin(ep)*sin(dpsi); 
        N[3]=-cos(e)*sin(dpsi); 
        N[4]=cos(e)*cos(ep)*cos(dpsi)+sin(e)*sin(ep);  
        N[5]=cos(e)*sin(ep)*cos(dpsi)-sin(e)*cos(ep);   
        N[6]=-sin(e)*sin(dpsi); 
        N[7]=sin(e)*cos(ep)*cos(dpsi)-cos(e)*sin(ep);  
        N[8]=sin(e)*sin(ep)*cos(dpsi)+cos(e)*cos(ep);  
 
        /** computing nutation+precession rotation matrix **/
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
                M[i*3+j]=P[i*3]*N[j]+P[i*3+1]*N[j+3]+P[i*3+2]*N[j+6];
        }
        //M[0]=1.0;M[1]=0.0;M[2]=0.0;M[3]=0.0;M[4]=1.0;M[5]=0.0;M[6]=0.0;M[7]=0.0;M[8]=1.0;       
        //std::cout << "M " << M[0] << "," << M[1] << "," << M[2] << "," << M[3] << "," << M[4] << "," << M[5] << "," << M[6] << "," << M[7] << "," << M[8] << std::endl;

        return 0;
    }

        bool conversion_frames::tod_to_bf(const double &JD, std::vector<double> &M){

        /** sanity checks **/
        if(M.size()!=9)
            smartastro_throw("TOD_TO_BF: rotation matrix must have 9 components");

        /** declarations **/
        double T,TT,TTT,T0,e,Omega,dpsi,GAST,GMST,pi,UT1,JD0,fJD;
        std::vector<double> Pi(9), Th(9);

        pi=constants::pi;
        T=(JD-2451545.0)/36525.0;
        TT=T*T;
        TTT=TT*T;
        fJD=double(floor(JD));
        if (JD-fJD>=0.5) 
            JD0=fJD+0.5;
        else 
            JD0=fJD-0.5;
        T0=(JD0-2451545.0)/36525.0;
        UT1=24.0*(JD-JD0); // approximation of UT1 with UTC         
        GMST=(24110.54841+8640184.812866*T0+1.002737909350795*UT1*3600.0+0.093104*TT-0.0000062*TTT)*(pi/12.0)/3600.0;
        e=(23.43929111+(-46.8150*T-0.00059*TT+0.001813*TTT)/3600.0)*(pi/180.0);
        Omega=(125.0+2.0/60.0+40.280/3600.0-(1934.0+8.0/60.0+10.539/3600.0)*T+7.455*TT/3600.0+0.008*TTT/3600.0)*(pi/180.0);
        e+=(9.203*cos(Omega))*(1.0/3600.0)*(pi/180.0);
        dpsi=(-17.200*sin(Omega)/3600.0)*(pi/180.0);

        /** computing parameters **/
        GAST=GMST+dpsi*cos(e);
        GAST -= 2.0*pi*floor(GAST/(2.0*pi));

        /** computing Earth's rotation matrix **/
        Th[0]=cos(GAST); 
        Th[1]=-sin(GAST); 
        Th[2]=0.0; 
        Th[3]=sin(GAST); 
        Th[4]=cos(GAST);  
        Th[5]=0.0;  
        Th[6]=0.0; 
        Th[7]=0.0; 
        Th[8]=1.0;

        /** computing polar motion matrix **/
        double xp=(0.043549/3600.0)*pi/180.0, yp=(0.336691/3600.0)*pi/180.0; // mean polar coordinates in 2000 in rad
        Pi[0]=1.0; 
        Pi[1]=0.0; 
        Pi[2]=-xp; 
        Pi[3]=0.0; 
        Pi[4]=1.0;  
        Pi[5]=yp;  
        Pi[6]=xp; 
        Pi[7]=-yp; 
        Pi[8]=1.0;       

        /** computing global rotation matrix **/
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
                M[i*3+j]=Th[i*3]*Pi[j]+Th[i*3+1]*Pi[j+3]+Th[i*3+2]*Pi[j+6];
        }
        
        return 0;
    }


    bool conversion_frames::inertial_to_bf(const double &JD, std::vector<double> &M){

        /** sanity checks **/
        if(M.size()!=9)
            smartastro_throw("INERTIAL_TO_BF: rotation matrix must have 9 components");

        std::vector<double> A(9), B(9);

        inertial_to_tod(JD,A);
        tod_to_bf(JD,B);

        /** computing rotation matrix **/
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
                M[i*3+j]=A[i*3]*B[j]+A[i*3+1]*B[j+3]+A[i*3+2]*B[j+6];
        }
        return 0;
    }
