/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/dearth_6dof.h"
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include "Astrodynamics/data/aero_database.h"
#include "Astrodynamics/data/aero_database_heavy.h"

#include "IRI_wrapper.h"
#include "SCCal.h"

using namespace smartastro;
using namespace smartastro::astrodynamics;
using namespace smartastro::aero_database_heavy;

template<class T>
dearth_6dof<T>::dearth_6dof(const std::vector<T> &I,const double &r_scale,const double &t_scale,const int &n_max, const std::vector<int> &flags, const std::vector<T> &p, const int &object, const std::vector<T> &coord, const std::vector<T> &other, const std::vector< std::vector<double> > &F10dot7) :
base_dearth<T>("Earth-centered dynamics with 6DoF",r_scale,t_scale,n_max,flags,p,other,F10dot7),m_I(I),m_object(object),m_CoS(coord)
{
  /** sanity checks **/
  if((m_I[0]<m_I[1])||(m_I[0]<m_I[2]))
      smartastro_throw("DEARTH_6DOF: Major principal axis must be the first axis");
  if((m_object!=0)&&(m_object!=1)&&(m_object!=2))
      smartastro_throw("DEARTH_6DOF: There are only three possibilites for the aerodynamic effects: plate (0), GOCE-like cylinder (1) or sphere (2)");  
  if(m_CoS.size()!=3)
      smartastro_throw("DEARTH_6DOF: there must be 3 coordinates for center of symetric geometry");     

}


template<class T>
dearth_6dof<T>::~dearth_6dof(){

}

template<class T>
int dearth_6dof<T>::evaluate(const double &t, const std::vector<T> &x,std::vector<T> &dx) const{

  std::vector<T> acc;
  T zero=0.0*x[0];
  acc.push_back(zero);
  acc.push_back(zero);
  acc.push_back(zero);
  std::vector<T> acc1=acc, acc2=acc, acc3=acc, acc4=acc, acc5=acc;
  std::vector<T> tor=acc, tor2=acc, tor4=acc, tor5=acc;

  double JD=t*m_T_scale/(3600.0*24.0); // Julian date

  /** computing accelerations and torques **/
  Earth_gravity(JD,x,acc1); // Earth gravity

  if(m_flags[0]==1)
   aero(JD,x,acc2,tor2); // atmospheric drag

  if(m_flags[1]==1)
   lunisolar(JD,x,acc3);

  if(m_flags[2]!=0)
   SRP(JD,x,acc4,tor4); 

  if(m_flags[3]==1)
   magnetism(JD,x,acc5,tor5);

  /** summing up contributions **/
  for(int i=0; i<3; i++)
  {
      acc[i] = acc1[i]+acc2[i]+acc3[i]+acc4[i]+acc5[i];
      tor[i] = tor2[i]+tor4[i]+tor5[i];
  }

  /** constructing state derivative **/
  /** position **/
  for(int i=0; i<3; i++)
      dx[i] = x[3+i];
  /** velocity **/
  for(int i=0; i<3; i++)
     dx[3+i] = acc[i];
  /** mass **/
  dx[6]=0.0;
  /** quaternions **/
  dx[7]=0.5*(x[13]*x[8]-x[12]*x[9]+x[11]*x[10]);
  dx[8]=0.5*(-x[13]*x[7]+x[11]*x[9]+x[12]*x[10]);
  dx[9]=0.5*(x[12]*x[7]-x[11]*x[8]+x[13]*x[10]);
  dx[10]=0.5*(-x[11]*x[7]-x[12]*x[8]-x[13]*x[9]);

  /** computations for torque **/
  std::vector<T> q, CM, dir;
  for(int i=0; i<4; i++)
    q.push_back(x[7+i]);
  for(int i=0; i<9; i++)
    CM.push_back(zero);

  orientation(q,CM);
  T r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); // distance to Earth center
  /** computing spacecraft direction in body-fixed frame **/
  dir.push_back((CM[0]*x[0]+CM[3]*x[1]+CM[6]*x[2])/r);
  dir.push_back((CM[1]*x[0]+CM[4]*x[1]+CM[7]*x[2])/r);
  dir.push_back((CM[2]*x[0]+CM[5]*x[1]+CM[8]*x[2])/r);
  /** computing 1st order gravity-gradient torque per inertial moment **/
  double mu=constants::mu_earth*pow(m_T_scale,2)/pow(m_L_scale,3); // Earth gravitational constant
  T factor=3.0*mu/pow(r,3); 
  /** rotation rates **/
  dx[11]=(x[12]*x[13]-factor*dir[1]*dir[2])*(m_I[1]-m_I[2])/m_I[0]+tor[0]*(pow(m_L_scale,2)/m_I[0]);
  dx[12]=(x[11]*x[13]-factor*dir[2]*dir[0])*(m_I[2]-m_I[0])/m_I[1]+tor[1]*(pow(m_L_scale,2)/m_I[1]);
  dx[13]=(x[11]*x[12]-factor*dir[0]*dir[1])*(m_I[0]-m_I[1])/m_I[2]+tor[2]*(pow(m_L_scale,2)/m_I[2]);

  if(std::isnan(dx[13])==1)
     smartastro_throw("EVALUATE: Nan value in the dynamics");

  return 0;
}

template<class T>
int dearth_6dof<T>::orientation(const std::vector<T> &q,std::vector<T> &C) const{

  /** sanity checks **/
  if(q.size()!=4)
      smartastro_throw("ORIENTATION: there must be 4 quaternions");
  if(C.size()!=9)
     smartastro_throw("ORIENTATION: orientation matrix should be 3x3");

  T q00=q[0]*q[0];
  T q11=q[1]*q[1];
  T q22=q[2]*q[2];
  T q33=q[3]*q[3]; 
  T q01=q[0]*q[1];
  T q02=q[0]*q[2];
  T q03=q[0]*q[3];
  T q12=q[1]*q[2];
  T q13=q[1]*q[3];
  T q23=q[2]*q[3];

  C[0]=q00-q11-q22+q33;
  C[1]=2.0*(q01-q23);
  C[2]=2.0*(q02+q13);
  C[3]=2.0*(q01+q23);
  C[4]=-q00+q11-q22+q33;
  C[5]=2.0*(q12-q03);
  C[6]=2.0*(q02-q13);
  C[7]=2.0*(q12+q03);
  C[8]=-q00-q11+q22+q33;

  return 0;
}


template<class T>
int dearth_6dof<T>::aero(const double &JD, const std::vector<T> &x, std::vector<T> &acc, std::vector<T> &tor) const{

  /** sanity checks **/
  if(x.size()!=14)
     smartastro_throw("AERO: state must be a 14 dimensional vector");
  if(acc.size()!=3)
     smartastro_throw("AERO: acceleration must be a 3 dimensional vector");
  if(tor.size()!=3)
     smartastro_throw("AERO: torque must be a 3 dimensional vector");

  double omega=constants::omega_earth*m_T_scale; // Earth spinning rate
  double pi=constants::pi;
  std::vector<double> M(9);
  std::vector<T> pos, aux, q, CM;
  T zero=0.0*x[0];

  astrocore::conversion_frames::inertial_to_tod(JD,M); // rotation matrix

  /** computing position in true of date */
  pos.push_back(M[0]*x[0]+M[3]*x[1]+M[6]*x[2]);
  pos.push_back(M[1]*x[0]+M[4]*x[1]+M[7]*x[2]);
  pos.push_back(M[2]*x[0]+M[5]*x[1]+M[8]*x[2]);

  /** computing altitude iteratively */
  T phi=zero, Z=zero;
  geodetic(pos,Z,phi);
  if(Z>2500.0)
  {
    acc[0]=0.0;acc[1]=0.0;acc[2]=0.0;tor[0]=0.0;tor[1]=0.0;tor[2]=0.0;
    return 0;
  }

  T rho=0.0*x[0];
  density(JD,Z,phi,atan2(pos[1],pos[0]),rho);
  if(rho<0.0)
     smartastro_throw("AERO: negative atmospheric density");

  std::vector<T> v_rot;
  v_rot.push_back(-omega*x[1]);
  v_rot.push_back(+omega*x[0]);
  v_rot.push_back(0.0*x[0]);

  /* computing air velocity in inertial frame */
  std::vector<T> v_atmo=v_rot;
  v_atmo[0]=M[0]*v_rot[0]+M[1]*v_rot[1]+M[2]*v_rot[2];
  v_atmo[1]=M[3]*v_rot[0]+M[4]*v_rot[1]+M[5]*v_rot[2];
  v_atmo[2]=M[6]*v_rot[0]+M[7]*v_rot[1]+M[8]*v_rot[2];

  /* computing relative velocity of atmosphere wrt object */
  std::vector<T> vrel=v_atmo;
  vrel[0]=-x[3]+v_atmo[0];
  vrel[1]=-x[4]+v_atmo[1];
  vrel[2]=-x[5]+v_atmo[2];

  T v = sqrt(vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2]); // norm of relative velocity

  for(int i=0; i<4; i++)
    q.push_back(x[7+i]);
  for(int i=0; i<9; i++)
    CM.push_back(zero);

  orientation(q,CM);

  /** computing relative velocity in body-fixed frame **/
  aux.push_back(CM[0]*vrel[0]+CM[3]*vrel[1]+CM[6]*vrel[2]);
  aux.push_back(CM[1]*vrel[0]+CM[4]*vrel[1]+CM[7]*vrel[2]);
  aux.push_back(CM[2]*vrel[0]+CM[5]*vrel[1]+CM[8]*vrel[2]);

  if(m_object==0)
  {// square plate case

    /** computing aerodynamic angles **/
    T ss=atan2(-aux[2],aux[0]); // side slip angle
    if(ss>pi)
       ss-=2.0*pi;
    T aa=0.0*x[0];
    if(aux[1]/v>0.9999)
       aa=pi/2.0;
    else if(aux[1]/v<-0.9999)
       aa=-pi/2.0;
    else
       aa=asin(aux[1]/v);

    if(std::isnan(aa)==1)
       smartastro_throw("AERO: getting NaN as the arcsin for aerodynamic angle");

    int i=0, j=0;
    while(aa*180.0/pi>AA_square_plate[j+1])
       j+=1;
    while(ss*180.0/pi>SS_square_plate[i+1])
       i+=1;      
    T Fx=0.0*x[0], Fy=Fx, Fz=Fx, Mx=Fx, My=Fx, Mz=Fz; // normalized forces and torques components 

    /*  bilinear interpolation */
    T dA=aa*180.0/pi-AA_square_plate[j];
    T dS=ss*180.0/pi-SS_square_plate[i];
    double DeltaA=AA_square_plate[j+1]-AA_square_plate[j];
    double DeltaS=SS_square_plate[i+1]-SS_square_plate[i];

    T DeltaFxA=coeff_Fx_square_plate[i][j+1]-coeff_Fx_square_plate[i][j];
    T DeltaFxS=coeff_Fx_square_plate[i+1][j]-coeff_Fx_square_plate[i][j];
    T DeltaFxAS=coeff_Fx_square_plate[i][j]+coeff_Fx_square_plate[i+1][j+1]-coeff_Fx_square_plate[i+1][j]-coeff_Fx_square_plate[i][j+1];    
    Fx = coeff_Fx_square_plate[i][j] + DeltaFxA*dA/DeltaA + DeltaFxS*dS/DeltaS + DeltaFxAS*(dA/DeltaA)*(dS/DeltaS); 

    T DeltaFyA=coeff_Fy_square_plate[i][j+1]-coeff_Fy_square_plate[i][j];
    T DeltaFyS=coeff_Fy_square_plate[i+1][j]-coeff_Fy_square_plate[i][j];
    T DeltaFyAS=coeff_Fy_square_plate[i][j]+coeff_Fy_square_plate[i+1][j+1]-coeff_Fy_square_plate[i+1][j]-coeff_Fy_square_plate[i][j+1];    
    Fy = coeff_Fy_square_plate[i][j] + DeltaFyA*dA/DeltaA + DeltaFyS*dS/DeltaS + DeltaFyAS*(dA/DeltaA)*(dS/DeltaS); 

    T DeltaFzA=coeff_Fz_square_plate[i][j+1]-coeff_Fz_square_plate[i][j];
    T DeltaFzS=coeff_Fz_square_plate[i+1][j]-coeff_Fz_square_plate[i][j];
    T DeltaFzAS=coeff_Fz_square_plate[i][j]+coeff_Fz_square_plate[i+1][j+1]-coeff_Fz_square_plate[i+1][j]-coeff_Fz_square_plate[i][j+1];    
    Fz = coeff_Fz_square_plate[i][j] + DeltaFzA*dA/DeltaA + DeltaFzS*dS/DeltaS + DeltaFzAS*(dA/DeltaA)*(dS/DeltaS); 

    T DeltaMxA=coeff_Mx_square_plate[i][j+1]-coeff_Mx_square_plate[i][j];
    T DeltaMxS=coeff_Mx_square_plate[i+1][j]-coeff_Mx_square_plate[i][j];
    T DeltaMxAS=coeff_Mx_square_plate[i][j]+coeff_Mx_square_plate[i+1][j+1]-coeff_Mx_square_plate[i+1][j]-coeff_Mx_square_plate[i][j+1];    
    Mx = coeff_Mx_square_plate[i][j] + DeltaMxA*dA/DeltaA + DeltaMxS*dS/DeltaS + DeltaMxAS*(dA/DeltaA)*(dS/DeltaS); 

    T DeltaMyA=coeff_My_square_plate[i][j+1]-coeff_My_square_plate[i][j];
    T DeltaMyS=coeff_My_square_plate[i+1][j]-coeff_My_square_plate[i][j];
    T DeltaMyAS=coeff_My_square_plate[i][j]+coeff_My_square_plate[i+1][j+1]-coeff_My_square_plate[i+1][j]-coeff_My_square_plate[i][j+1];    
    My = coeff_My_square_plate[i][j] + DeltaMyA*dA/DeltaA + DeltaMyS*dS/DeltaS + DeltaMyAS*(dA/DeltaA)*(dS/DeltaS); 

    T DeltaMzA=coeff_Mz_square_plate[i][j+1]-coeff_Mz_square_plate[i][j];
    T DeltaMzS=coeff_Mz_square_plate[i+1][j]-coeff_Mz_square_plate[i][j];
    T DeltaMzAS=coeff_Mz_square_plate[i][j]+coeff_Mz_square_plate[i+1][j+1]-coeff_Mz_square_plate[i+1][j]-coeff_Mz_square_plate[i][j+1];    
    Mz = coeff_Mz_square_plate[i][j] + DeltaMzA*dA/DeltaA + DeltaMzS*dS/DeltaS + DeltaMzAS*(dA/DeltaA)*(dS/DeltaS); 

    /*  nearest neighbourgh interpolation  */
    // T inter1=sqrt(pow(aa*180.0/pi-AA_square_plate[j],2)+pow(ss*180.0/pi-SS_square_plate[i],2));
    // T inter2=sqrt(pow(aa*180.0/pi-AA_square_plate[j+1],2)+pow(ss*180.0/pi-SS_square_plate[i],2));
    // T inter3=sqrt(pow(aa*180.0/pi-AA_square_plate[j+1],2)+pow(ss*180.0/pi-SS_square_plate[i+1],2));
    // T inter4=sqrt(pow(aa*180.0/pi-AA_square_plate[j],2)+pow(ss*180.0/pi-SS_square_plate[i+1],2));
    // if((inter2<inter1)&&(inter2<inter3)&&(inter2<inter4))
    //   j+=1;
    // if((inter3<inter1)&&(inter3<inter2)&&(inter3<inter4))
    // {
    //    i+=1;j+=1;     
    // }
    // if((inter4<inter1)&&(inter4<inter2)&&(inter4<inter3))
    //   i+=1;
    //Fx=coeff_Fx_square_plate[i][j];
    //Fy=coeff_Fy_square_plate[i][j];
    //Fz=coeff_Fz_square_plate[i][j];
    //Mx=coeff_Mx_square_plate[i][j];
    //My=coeff_My_square_plate[i][j];
    //Mz=coeff_Mz_square_plate[i][j];

    /** computing acceleration and torque in body-fixed frame **/   
    T factor1=0.5*rho*m_drag[0]*m_L_scale*pow(v,2)/x[6];
    double Ax=factor1*Fx;
    double Ay=factor1*Fy;
    double Az=factor1*Fz;
    T factor2=0.5*rho*m_drag[0]*m_drag[1]*pow(v,2);
    tor[0]=factor2*Mx+(m_CoS[1]*Az-m_CoS[2]*Ay)*x[6]/m_L_scale;
    tor[1]=factor2*My+(m_CoS[2]*Ax-m_CoS[0]*Az)*x[6]/m_L_scale;
    tor[2]=factor2*Mz+(m_CoS[0]*Ay-m_CoS[1]*Ax)*x[6]/m_L_scale;

    /** computing acceleration in inertial frame **/
    acc[0]=CM[0]*Ax+CM[1]*Ay+CM[2]*Az;
    acc[1]=CM[3]*Ax+CM[4]*Ay+CM[5]*Az;
    acc[2]=CM[6]*Ax+CM[7]*Ay+CM[8]*Az;
  }

  if(m_object==1)
  {// cylinder case

    T ss=atan2(-aux[2],aux[0]);
    if(ss>pi)
      ss-=2.0*pi;

    T aa=0.0*x[0];
    if(aux[1]/v>0.9999)
      aa=pi/2.0;
    else if(aux[1]/v<-0.9999)
      aa=-pi/2.0;
    else
       aa=asin(aux[1]/v);

    if(std::isnan(aa)==1)
      smartastro_throw("AERO: getting NaN as the arcsin for aerodynamic angle");

    T alpha=-ss, beta=-aa;

    std::vector<double> Wind2Body(9);
    Wind2Body[0]=cos(alpha)*cos(beta);
    Wind2Body[1]=-sin(beta);
    Wind2Body[2]=sin(alpha)*cos(beta);
    Wind2Body[3]=cos(alpha)*sin(beta);
    Wind2Body[4]=cos(beta);
    Wind2Body[5]=sin(alpha)*sin(beta);
    Wind2Body[6]=-sin(alpha);
    Wind2Body[7]=0.0;
    Wind2Body[8]=cos(alpha);

    int i=0, j=0, k=0;
    while(alpha*180.0/pi>alpha_GOCE_cylinder[j+1])
       j+=1;
    while(beta*180.0/pi>beta_GOCE_cylinder[i+1])
       i+=1;  
    while((Z>alt_GOCE_cylinder[k+1])&&(k+1<nb_altitudes_GOCE_cylinder-1))
       k+=1;   

    if((i>92)||(j>180)||(k>nb_altitudes_GOCE_cylinder-1))
      smartastro_throw("AERO: problem with the inputs for the database");

    /** computing acceleration and torque wrt center of symmetry in wind frame **/ 
    T factor1=0.5*rho*m_drag[0]*m_L_scale*pow(v,2)/x[6];
    double A1=factor1*coeff_CD_GOCE_cylinder[k][i][j];
    double A2=factor1*coeff_CL_GOCE_cylinder[k][i][j];
    double A3=factor1*coeff_CS_GOCE_cylinder[k][i][j];
    T factor2=0.5*rho*m_drag[0]*m_drag[1]*pow(v,2);    
    double M1=factor2*coeff_M1_GOCE_cylinder[k][i][j];
    double M2=factor2*coeff_M2_GOCE_cylinder[k][i][j];
    double M3=factor2*coeff_M3_GOCE_cylinder[k][i][j];

    /** computing acceleration and torque wrt center of gravity in body-fixed frame **/   
    double Ax=A1*Wind2Body[0]+A2*Wind2Body[3]+A3*Wind2Body[6];
    double Ay=A1*Wind2Body[1]+A2*Wind2Body[4]+A3*Wind2Body[7];
    double Az=A1*Wind2Body[2]+A2*Wind2Body[5]+A3*Wind2Body[8];
    tor[0]=M1*Wind2Body[0]+M2*Wind2Body[3]+M3*Wind2Body[6]+(m_CoS[1]*Az-m_CoS[2]*Ay)*x[6]/m_L_scale;
    tor[1]=M1*Wind2Body[1]+M2*Wind2Body[4]+M3*Wind2Body[7]+(m_CoS[2]*Ax-m_CoS[0]*Az)*x[6]/m_L_scale;
    tor[2]=M1*Wind2Body[2]+M2*Wind2Body[5]+M3*Wind2Body[8]+(m_CoS[0]*Ay-m_CoS[1]*Ax)*x[6]/m_L_scale;

    /** computing acceleration in inertial frame **/
    acc[0]=CM[0]*Ax+CM[1]*Ay+CM[2]*Az;
    acc[1]=CM[3]*Ax+CM[4]*Ay+CM[5]*Az;
    acc[2]=CM[6]*Ax+CM[7]*Ay+CM[8]*Az;
  }

  if(m_object==2)
  {// sphere case

    /** computing acceleration in inertial frame **/
    double Cd=2.0;
    T factor1=0.5*Cd*rho*m_drag[0]*m_L_scale*pow(v,2)/x[6];
    for(int i=0; i<3; i++)
      acc[i]=vrel[i]*factor1;

    /** computing acceleration in body fixed frame **/    
    std::vector<T> acc2=acc;
    acc2[0]=CM[0]*acc[0]+CM[3]*acc[1]+CM[6]*acc[2];
    acc2[1]=CM[1]*acc[0]+CM[4]*acc[1]+CM[7]*acc[2];
    acc2[2]=CM[2]*acc[0]+CM[5]*acc[1]+CM[8]*acc[2];

    /** computing torque in body-fixed frame **/   
    tor[0]=(m_CoS[1]*acc2[2]-m_CoS[2]*acc2[1])*x[6]/m_L_scale;
    tor[1]=(m_CoS[2]*acc2[0]-m_CoS[0]*acc2[2])*x[6]/m_L_scale;
    tor[2]=(m_CoS[0]*acc2[1]-m_CoS[1]*acc2[0])*x[6]/m_L_scale;
  }

  return 0;
}


template<class T>
int dearth_6dof<T>::SRP(const double &jd, const std::vector<T> &x, std::vector<T> &acc, std::vector<T> &tor) const{

  /** sanity checks **/
  if(x.size()!=14)
      smartastro_throw("SRP: state must be a 14 dimensional vector");
  if(acc.size()!=3)
      smartastro_throw("SRP: acceleration must be a 3 dimensional vector");
  if(tor.size()!=3)
      smartastro_throw("SRP: torque must be a 3 dimensional vector");

  /** declarations **/
  std::vector<double> dir_sun(3);
  double r_sun,eps,lambda_sun;
  eps=constants::obliquity_ecliptic_2000;
  double AU=constants::au*1.0e3/m_L_scale; // scaled astronautical unit
  T C_dr=0.0*x[0];

  /** Sun's ephemerides **/
  Sun(jd,r_sun,lambda_sun);
  r_sun *= 1.0e3 / m_L_scale;

  /** computing Sun's position in inertial frame **/
  dir_sun[0]=cos(lambda_sun);
  dir_sun[1]=sin(lambda_sun)*cos(eps);
  dir_sun[2]=sin(lambda_sun)*sin(eps);

  T PA=constants::SRP*pow(m_T_scale,2)/m_L_scale; 
  T epsilon=m_other[0]; // specular reflectivity
  if(m_other.size()>1)
  {
    C_dr=m_other[1]; // coefficient of diffusive reflection
    if((C_dr>1.0)||(C_dr<0.0))
      smartastro_throw("SRP: diffusive reflection must be between 0.0 and 1.0");
  }

  if((epsilon>1.0)||(epsilon<0.0))
    smartastro_throw("SRP: specular reflectivity must be between 0.0 and 1.0");

  std::vector<T> q, CM, aux;
  T zero=0.0*x[0];
  for(int i=0; i<4; i++)
    q.push_back(x[7+i]);
  for(int i=0; i<9; i++)
    CM.push_back(zero);

  orientation(q,CM);  

  if(m_object==0){ // plate

    /** computing normal to the plate in inertial frame (principal axis is I[0]) **/
    aux.push_back(CM[0]);
    aux.push_back(CM[3]);
    aux.push_back(CM[6]);

    T costheta=aux[0]*dir_sun[0]+aux[1]*dir_sun[1]+aux[2]*dir_sun[2];
    if(costheta<0.0)
    {
      costheta=-costheta;
      aux[0]=-aux[0];aux[1]=-aux[1];aux[2]=-aux[2];
    }

    T factor = -PA*m_drag[0]*costheta*pow(AU/r_sun,2)/x[6];

    /** computing eclipse conditions **/
    factor*=shadow(jd,x);

    /** computing acceleration in inertial frame **/
    for(int i=0; i<3; i++)
      acc[i]=factor*((1.0-epsilon)*dir_sun[i]+2.0*(epsilon*costheta+C_dr/3.0)*aux[i]);
  }  

  if(m_object==1)
  { // cylinder
    T diameter = 2.0*sqrt(m_drag[0]/constants::pi);
    T factor = -PA*pow(AU/r_sun,2)/x[6];

    /** computing eclipse conditions **/
    factor*=shadow(jd,x);

    int index_minor_axis=1;
    if(m_I[2]<m_I[index_minor_axis])
      index_minor_axis=2;
    std::vector<double> dir_cylinder=dir_sun;
    dir_cylinder[0]=CM[index_minor_axis];
    dir_cylinder[1]=CM[index_minor_axis+3];
    dir_cylinder[2]=CM[index_minor_axis+6];
    T phi_sun=acos(dir_sun[0]*dir_cylinder[0]+dir_sun[1]*dir_cylinder[1]+dir_sun[2]*dir_cylinder[2]);

    for(int i=0; i<3; i++)
      acc[i]=factor*(((sin(phi_sun)*(1.0+epsilon/3.0)+constants::pi*C_dr/6.0)*diameter*m_drag[1]+cos(phi_sun)*(1.0-epsilon)*m_drag[0])*dir_sun[i]+(-(4.0*sin(phi_sun)*epsilon/3.0+constants::pi*C_dr/6.0)*diameter*m_drag[1]+2.0*(epsilon*cos(phi_sun)+C_dr/3.0)*m_drag[0])*cos(phi_sun)*dir_cylinder[i]);

  }

  if(m_object==2)
  { // sphere
    T factor = -(1.0+4.0*C_dr/9.0)*PA*m_drag[0]*pow(AU/r_sun,2)/x[6];

    /** computing eclipse conditions **/
    factor*=shadow(jd,x);

    /** computing acceleration in inertial frame  **/
    for(int i=0; i<3; i++)
      acc[i]=factor*dir_sun[i];
  }

  /** computing acceleration in body-fixed frame  **/ 
  std::vector<T> acc2=acc;
  acc2[0]=CM[0]*acc[0]+CM[3]*acc[1]+CM[6]*acc[2];
  acc2[1]=CM[1]*acc[0]+CM[4]*acc[1]+CM[7]*acc[2];
  acc2[2]=CM[2]*acc[0]+CM[5]*acc[1]+CM[8]*acc[2];

  /** computing torque  **/ 
  tor[0]=(m_CoS[1]*acc2[2]-m_CoS[2]*acc2[1])*x[6]/m_L_scale;
  tor[1]=(m_CoS[2]*acc2[0]-m_CoS[0]*acc2[2])*x[6]/m_L_scale;
  tor[2]=(m_CoS[0]*acc2[1]-m_CoS[1]*acc2[0])*x[6]/m_L_scale;

  //std::cout << x[6]*sqrt(acc[0]*acc[0]+acc[1]*acc[1]+acc[2]*acc[2])*m_L_scale/pow(m_T_scale,2) << std::endl;
  return 0;
}


template<class T>
int dearth_6dof<T>::magnetism(const double &jd, const std::vector<T> &x, std::vector<T> &acc, std::vector<T> &tor) const{

  /** sanity checks **/
  if(x.size()!=14)
      smartastro_throw("MAGNETISM: state must be a 14 dimensional vector");
  if(acc.size()!=3)
      smartastro_throw("MAGNETISM: acceleration must be a 3 dimensional vector");
  if(tor.size()!=3)
      smartastro_throw("MAGNETISM: torque must be a 3 dimensional vector");     

  std::vector<double> M(9);
  std::vector<T> posvel, aux=acc;
  double omega=constants::omega_earth*m_T_scale;
  double Re=constants::R_earth/m_L_scale;
   
  /*  position and velocity in Earth fixed frame  */
  smartastro::astrocore::conversion_frames::inertial_to_bf(jd,M); // rotation matrix
  posvel.push_back(M[0]*x[0]+M[3]*x[1]+M[6]*x[2]);
  posvel.push_back(M[1]*x[0]+M[4]*x[1]+M[7]*x[2]);
  posvel.push_back(M[2]*x[0]+M[5]*x[1]+M[8]*x[2]);
  posvel.push_back(M[0]*x[3]+M[3]*x[4]+M[6]*x[5]+omega*posvel[1]);
  posvel.push_back(M[1]*x[3]+M[4]*x[4]+M[7]*x[5]-omega*posvel[0]);
  posvel.push_back(M[2]*x[3]+M[5]*x[4]+M[8]*x[5]);

  T Q=0.0*x[0];
  if(m_other.size()>=3)
    Q=m_other[2]/m_T_scale;
  else
  {
    T C=0.0*x[0], V=0.0*x[0];
    charging(jd, x, C, V); // computing electrostatic capacitance and surface potential in SI units
    Q=C*V/m_T_scale; // scaled charge 
  }

  /*  Computing the Lorentz acceleration  */    
  std::vector<T> pos;
  pos.push_back(posvel[0]);
  pos.push_back(posvel[1]);
  pos.push_back(posvel[2]);

  // T inv_r=1.0/sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
  // std::vector<T> dir;
  // dir.push_back(pos[0]*inv_r);
  // dir.push_back(pos[1]*inv_r);
  // dir.push_back(pos[2]*inv_r);
  // T factor=pow(Re*inv_r,3)*Q/x[6]; 
  // double g11=(constants::g_Earth_2015[(constants::n_Earth_magnetic+1)*1+1]+(date[0]-2015.0)*constants::g_dot_Earth_2015[(constants::n_Earth_magnetic+1)*1+1])*1.0e-9*pow(m_T_scale,2);
  // double h11=(constants::h_Earth_2015[(constants::n_Earth_magnetic+1)*1+1]+(date[0]-2015.0)*constants::h_dot_Earth_2015[(constants::n_Earth_magnetic+1)*1+1])*1.0e-9*pow(m_T_scale,2);
  // double g10=(constants::g_Earth_2015[(constants::n_Earth_magnetic+1)*0+1]+(date[0]-2015.0)*constants::g_dot_Earth_2015[(constants::n_Earth_magnetic+1)*0+1])*1.0e-9*pow(m_T_scale,2);
  // T inter=g11*dir[0]+h11*dir[1]+g10*dir[2];
  // std::vector<T> Bqm2;
  // Bqm2.push_back(factor*(3.0*dir[0]*inter-g11));
  // Bqm2.push_back(factor*(3.0*dir[1]*inter-h11));
  // Bqm2.push_back(factor*(3.0*dir[2]*inter-g10));

  std::vector<double> g((m_max_degree_Earth_magnetic+1)*(m_max_degree_Earth_magnetic+1),0.0),h((m_max_degree_Earth_magnetic+1)*(m_max_degree_Earth_magnetic+1),0.0);
  Gaussian_coeff(jd,m_max_degree_Earth_magnetic,g,h);

  std::vector<T> Bqm=pos;
  harmonics(pos,m_max_degree_Earth_magnetic,Re,g,h,Bqm);
  Bqm[0]*=-Q/x[6];
  Bqm[1]*=-Q/x[6];
  Bqm[2]*=-Q/x[6];
  //std::cout << 100.0*(1.0-sqrt(Bqm2[0]*Bqm2[0]+Bqm2[1]*Bqm2[1]+Bqm2[2]*Bqm2[2])/sqrt(Bqm[0]*Bqm[0]+Bqm[1]*Bqm[1]+Bqm[2]*Bqm[2])) << "%" << std::endl;
  aux[0]=posvel[4]*Bqm[2]-posvel[5]*Bqm[1];
  aux[1]=posvel[5]*Bqm[0]-posvel[3]*Bqm[2];
  aux[2]=posvel[3]*Bqm[1]-posvel[4]*Bqm[0];
  /*  Computing the Lorentz acceleration in inertial frame */ 
  acc[0]=M[0]*aux[0]+M[1]*aux[1]+M[2]*aux[2];
  acc[1]=M[3]*aux[0]+M[4]*aux[1]+M[5]*aux[2];
  acc[2]=M[6]*aux[0]+M[7]*aux[1]+M[8]*aux[2];

  /*  Computing the torque */ 
  std::vector<T> q, CM;
  T zero=0.0*x[0];
  for(int i=0; i<4; i++)
    q.push_back(x[7+i]);
  for(int i=0; i<9; i++)
    CM.push_back(zero);

  orientation(q,CM);
  // computing acceleration in body-fixed frame 
  std::vector<T> acc2=acc;
  acc2[0]=CM[0]*acc[0]+CM[3]*acc[1]+CM[6]*acc[2];
  acc2[1]=CM[1]*acc[0]+CM[4]*acc[1]+CM[7]*acc[2];
  acc2[2]=CM[2]*acc[0]+CM[5]*acc[1]+CM[8]*acc[2];
  tor[0]=(m_CoS[1]*acc2[2]-m_CoS[2]*acc2[1])*x[6]/m_L_scale;
  tor[1]=(m_CoS[2]*acc2[0]-m_CoS[0]*acc2[2])*x[6]/m_L_scale;
  tor[2]=(m_CoS[0]*acc2[1]-m_CoS[1]*acc2[0])*x[6]/m_L_scale;

  return 0;
}


template<class T>
int dearth_6dof<T>::charging(const double &jd, const std::vector<T> &x, T &C, T &V) const {

  /** sanity checks **/
  if(x.size()!=14)
      smartastro_throw("CHARGING: state must be a 14 dimensional vector");

  double qe=constants::qe;
  double kb=constants::kb;
  double eps0=constants::eps0;
  double omega=constants::omega_earth*m_T_scale;
  double pi=constants::pi;
  std::vector<T> pos, vel;
  std::vector<double> M(9,0.0);

  smartastro::astrocore::conversion_frames::inertial_to_bf(jd,M); // rotation matrix
  pos.push_back(M[0]*x[0]+M[3]*x[1]+M[6]*x[2]);
  pos.push_back(M[1]*x[0]+M[4]*x[1]+M[7]*x[2]);
  pos.push_back(M[2]*x[0]+M[5]*x[1]+M[8]*x[2]);
  vel.push_back(M[0]*x[3]+M[3]*x[4]+M[6]*x[5]+omega*pos[1]);
  vel.push_back(M[1]*x[3]+M[4]*x[4]+M[7]*x[5]-omega*pos[0]);
  vel.push_back(M[2]*x[3]+M[5]*x[4]+M[8]*x[5]);

  T zero=0.0*pos[0];
  T Z=zero, phi=zero;
  geodetic(pos,Z,phi);

  /** computing local time */
  double eps=constants::obliquity_ecliptic_2000;
  double r_sun, lambda_sun;
  Sun(jd,r_sun,lambda_sun);
  r_sun *= 1.0e3 / m_L_scale;
  std::vector<double> pos_sun(3);
  pos_sun[0]=r_sun*cos(lambda_sun);
  pos_sun[1]=r_sun*sin(lambda_sun)*cos(eps);
  pos_sun[2]=r_sun*sin(lambda_sun)*sin(eps);
  double alpha_sun=atan2(pos_sun[1],pos_sun[0]);
  T alpha=atan2(pos[1],pos[0]); 
  T local_time=alpha-alpha_sun+pi;

  T Te=0.0*x[0];
  IRI90(jd,Z,phi,local_time,Te); // computing electron temperature
  T Ti=Te; // approximated ion temperature
  double ne=1.0e11, ni=ne;  

  /*  Computing the capacitance C */     
  T L=sqrt(eps0*kb/(ne/Te +ni/Ti))/qe; // approximated Debye length
  if(m_object==0)
      C=0.5*eps0*m_drag[0]/L; // electrostatic capacitance for a plate
  if(m_object==1)
      C=2.0*pi*eps0*sqrt(m_drag[0]/pi)*m_drag[1]/L; // electrostatic capacitance for a cylinder
  if(m_object==2)
      C=4.0*eps0*m_drag[0]*(1.0/L+1.0/sqrt(m_drag[0]/pi)); // electrostatic capacitance for a sphere

  /*  Computing the angles relative to SC velocity and solar flux */
  T cos_vel=zero;
  T cos_Sun=zero;
  if(m_object==0)
  {
    std::vector<T> q, CM;
    for(int i=0; i<4; i++)
      q.push_back(x[7+i]);
    for(int i=0; i<9; i++)
      CM.push_back(zero);
 
    orientation(q,CM);

    cos_Sun=(CM[0]*pos_sun[0]+CM[3]*pos_sun[1]+CM[6]*pos_sun[2])/r_sun;
    if(cos_Sun<0.0)
      cos_Sun=-cos_Sun;

    cos_vel=(CM[0]*vel[0]+CM[3]*vel[1]+CM[6]*vel[2])/sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
    if(cos_vel<0.0)
      cos_vel=-cos_vel;

  }
  if(m_object==1)
  {
    std::vector<T> q, CM;
    for(int i=0; i<4; i++)
      q.push_back(x[7+i]);
    for(int i=0; i<9; i++)
      CM.push_back(zero);

    orientation(q,CM);  

    int index_minor_axis=1;
    if(m_I[2]<m_I[index_minor_axis])
      index_minor_axis=2;          
    cos_Sun = (CM[index_minor_axis]*pos_sun[0]+CM[index_minor_axis+3]*pos_sun[1]+CM[index_minor_axis+6]*pos_sun[2])/r_sun;
    cos_vel = (CM[index_minor_axis]*vel[0]+CM[index_minor_axis+3]*vel[1]+CM[index_minor_axis+6]*vel[2])/sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);        
  }  

  /*  Computing the surface potential V */ 
  T Vsc=(m_L_scale/m_T_scale)*sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]); // SC speed relative to ambient plasma (in standard units)
  double shape_factor;
  if(m_object==0)
     shape_factor=0.0;
  if(m_object==1)
     shape_factor=0.5;
  if(m_object==2)
     shape_factor=1.0;   
  V=surface_potential_LEO(Te,shape_factor,Vsc,cos_vel);

  #ifndef ENABLE_SMARTUQ // WARNING: because of this implementation trick, the fortran code cannot be called if SMART-UQ is ON
    /*  comment if not using IRI12 model in fortran */
    std::vector<double> date(6);
    astrocore::conversion_time::jd2date(jd, date);
    std::map<std::string, bool> iri_options;
    iri_options["geomagnetic"] = false;
    iri_options["universal_time"] = true;
    iri_options["auroral_boundary"] = false;
    iri_options["foE_storm"] = false;
    iri_options["hmF2_foF2_storm"] = false;
    iri_options["topside_foF2_storm"] = false;
    IRI::IRIresult res = IRI::getIRIProfileAt(Z, phi*180.0/pi, atan2(pos[1],pos[0])*180.0/pi, int(date[0]), int(100.0*date[1]+date[2]), float(date[3]+date[4]/60.0+date[5]/3600.0), iri_options);
    Te=double(res.Te)*11605.0; Ti=double(res.Ti)*11605.0; 
    ne=double(res.Ne); ni=ne;   
    SCCal::Spacecraft<double> spacecraft;
    spacecraft.velocity = (m_L_scale/m_T_scale)*sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
    spacecraft.potential = 0.0;    
     if(m_object==0)
     {    
      spacecraft.cos_angle_Sun = cos_Sun;
      spacecraft.cos_angle_vel = cos_vel;
      spacecraft.surface_area = 2.0*m_drag[0]; // edge effects neglected
      spacecraft.shape_factor = 0.0; // shape factor for plate
     }
     if(m_object==1)
     {
        spacecraft.cos_angle_Sun = cos_Sun;
        spacecraft.cos_angle_vel = cos_vel;        
        spacecraft.surface_area = 2.0*pi*sqrt(m_drag[0]/pi)*m_drag[1]; // edge effects neglected
        spacecraft.shape_factor = 0.5; // shape factor for cylinder 
     }   
     if(m_object==2)
     {
      spacecraft.surface_area = 4.0*m_drag[0];
      spacecraft.shape_factor = 1.0; // shape factor for sphere 
     }        
    double ShadowFunction = 0.0;
    if(m_flags[2]!=0)
      ShadowFunction=shadow(jd,x); 
    bool Auroral=false;      
    V = SCCal::solve<double>(res, spacecraft, ShadowFunction, Auroral);
    //  Recomputing the capacitance with IRI12
    L=sqrt(eps0*kb/(ne/Te +ni/Ti))/qe; // approximated Debye length
    if(Z>2000.0)
      L = 10.0;    
    if(m_object==0)
      C=0.5*eps0*m_drag[0]/L; // electrostatic capacitance for a plate
    if(m_object==1)
      C=2.0*pi*eps0*sqrt(m_drag[0]/pi)*m_drag[1]/L; // electrostatic capacitance for a cylinder
    if(m_object==2)
      C=4.0*eps0*m_drag[0]*(1.0/L+1.0/sqrt(m_drag[0]/pi)); // electrostatic capacitance for a sphere
  #endif

  return 0;
}


template class dearth_6dof<double>;
template class dearth_6dof<float>;
template class dearth_6dof<long double>;
#ifdef ENABLE_SMARTUQ
//template class dearth_6dof<smartuq::polynomial::chebyshev_polynomial>;
//template class dearth_6dof<smartuq::polynomial::taylor_polynomial>;
#endif
