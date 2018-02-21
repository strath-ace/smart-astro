/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/dearth_3dof.h"
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

#include "IRI_wrapper.h"
#include "SCCal.h"

using namespace smartastro;
using namespace smartastro::astrodynamics;

template<class T>
dearth_3dof<T>::dearth_3dof(const double &r_scale, const double &t_scale, const int &n_gravity, const std::vector<int> &flags, const std::vector<T> &p1, const std::vector<T> &p2, const std::vector< std::vector<double> > &F10dot7) :
    base_dearth<T>("Earth-centered dynamics with 3DoF",r_scale,t_scale,n_gravity,flags,p1,p2,F10dot7)
{   
    if(m_other.size() < 2)
      smartastro_throw("DEARTH_3DOF: there must be at least 2 'other' parameters");
}


template<class T>
dearth_3dof<T>::~dearth_3dof(){

}


template<class T>
int dearth_3dof<T>::evaluate(const double &t, const std::vector<T> &x,std::vector<T> &dx) const{

  /** sanity checks **/
  if(x.size()!=7)
    smartastro_throw("EVALUATE: state must be a 7 dimensional vector");
  if(dx.size()!=7)
    smartastro_throw("EVALUATE: state derivative must be a 7 dimensional vector");

  std::vector<T> acc1, acc2, acc3, acc4, acc5;
  T zero=0.0*x[0];
  acc1.push_back(zero);
  acc1.push_back(zero);
  acc1.push_back(zero);
  acc2=acc1;acc3=acc1;acc4=acc1;acc5=acc1;

  double jd=t*m_T_scale/(3600.0*24.0); // Julian date

  /** computing accelerations **/
  Earth_gravity(jd,x,acc1);

  if(m_flags[0]==1)
    aero(jd,x,acc2);

  if(m_flags[1]==1)
    lunisolar(jd,x,acc3);
  
  if(m_flags[2]!=0)
    SRP(jd,x,acc4);
    
  if(m_flags[3]==1)
    magnetism(jd,x,acc5);

  /** computing state derivative **/
  /** position **/
  for(int i=0; i<3; i++)
      dx[i] = x[3+i];
  /** velocity **/
  for(int i=0; i<3; i++)
     dx[3+i] = acc1[i]+acc2[i]+acc3[i]+acc4[i]+acc5[i];
  /** mass **/
  dx[6]=0.0;

  return 0;
}


template<class T>
int dearth_3dof<T>::aero(const double &JD, const std::vector<T> &x, std::vector<T> &acc) const{

  /** sanity checks **/
  if(x.size()!=7)
      smartastro_throw("AERO: state must be a 7 dimensional vector");
  if(acc.size()!=3)
      smartastro_throw("AERO: acceleration must be a 3 dimensional vector");

  /** declarations **/
  std::vector<double> M(9);
  std::vector<T> pos, vrel, v_atmo, v_rot; 
  double omega=constants::omega_earth*m_T_scale; // Earth spinning rate

  astrocore::conversion_frames::inertial_to_tod(JD,M); // rotation matrix

  /** computing position in true of date */
  pos.push_back(M[0]*x[0]+M[3]*x[1]+M[6]*x[2]);
  pos.push_back(M[1]*x[0]+M[4]*x[1]+M[7]*x[2]);
  pos.push_back(M[2]*x[0]+M[5]*x[1]+M[8]*x[2]);

  T Z=0.0*x[0], phi=0.0*x[0];
  geodetic(pos,Z,phi); // computing altitude and geodetic latitude
  if(Z>2500.0)
  {
    acc[0]=0.0;acc[1]=0.0;acc[2]=0.0;
    return 0;  
  }

  T rho=0.0*x[0];
  density(JD,Z,phi,atan2(pos[1],pos[0]),rho);

  /* computing air velocity in true of date */
  v_rot=pos;
  v_rot[0]=-omega*pos[1];
  v_rot[1]=+omega*pos[0];
  v_rot[2]=0.0;

  /* computing air velocity in inertial frame */
  v_atmo=v_rot;
  v_atmo[0]=M[0]*v_rot[0]+M[1]*v_rot[1]+M[2]*v_rot[2];
  v_atmo[1]=M[3]*v_rot[0]+M[4]*v_rot[1]+M[5]*v_rot[2];
  v_atmo[2]=M[6]*v_rot[0]+M[7]*v_rot[1]+M[8]*v_rot[2];

  /* computing relative velocity wrt atmosphere */
  vrel=v_atmo;
  vrel[0]=x[3]-v_atmo[0];
  vrel[1]=x[4]-v_atmo[1];
  vrel[2]=x[5]-v_atmo[2];

  T v = sqrt(vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2]); // norm of relative velocity

  T Cd=0.0*x[0]; // drag coefficient
  if(m_drag.size()>=2)
     Cd+=m_drag[1];
  else
     Cd+=2.0; // default value

  /** computing acceleration in inertial frame **/
  T factor=-0.5*Cd*rho*m_drag[0]*m_L_scale*v/x[6];
  for(int i=0; i<3; i++)
      acc[i]=vrel[i]*factor;

  return 0;
}


template<class T>
int dearth_3dof<T>::SRP(const double &jd, const std::vector<T> &x, std::vector<T> &acc) const{

  /** sanity checks **/
  if(x.size()!=7)
      smartastro_throw("SRP: state must be a 7 dimensional vector");
  if(acc.size()!=3)
      smartastro_throw("SRP: acceleration must be a 3 dimensional vector");

  /** declarations **/
  std::vector<double> dir_sun(3);
  double r_sun,eps,lambda_sun;
  eps=constants::obliquity_ecliptic_2000;
  T PA=constants::SRP*m_other[0]*pow(m_T_scale,2)/m_L_scale; // quantity with dimension of a force
  double AU=constants::au*1.0e3/m_L_scale; // scaled astronautical unit

  /** Sun's ephemerides **/
  Sun(jd,r_sun,lambda_sun);
  r_sun *= 1.0e3 / m_L_scale;

  /** computing Sun's position in inertial frame **/
  dir_sun[0]=cos(lambda_sun);
  dir_sun[1]=sin(lambda_sun)*cos(eps);
  dir_sun[2]=sin(lambda_sun)*sin(eps);

  T factor = -m_other[1]*PA*pow(AU/r_sun,2)/x[6];

  /** computing eclipse conditions **/
  factor*=shadow(jd,x);

  /** computing acceleration **/
  for(int i=0; i<3; i++)
      acc[i]=factor*dir_sun[i];
  
  return 0;
}


template<class T>
int dearth_3dof<T>::magnetism(const double &jd, const std::vector<T> &x, std::vector<T> &acc) const{

  /** sanity checks **/
  if(x.size()!=7)
      smartastro_throw("MAGNETISM: state must be a 7 dimensional vector");
  if(acc.size()!=3)
      smartastro_throw("MAGNETISM: acceleration must be a 3 dimensional vector");

  std::vector<double> M(9);
  std::vector<T> posvel, aux=acc;
  //double B0=constants::dipole_Earth*pow(m_T_scale,2)/pow(m_L_scale,3);
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
  if(m_other.size()>2)
    Q=m_other[2]/m_T_scale; // scaled charge
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

  std::vector<double> g((m_max_degree_Earth_magnetic+1)*(m_max_degree_Earth_magnetic+1),0.0),h((m_max_degree_Earth_magnetic+1)*(m_max_degree_Earth_magnetic+1),0.0);
  Gaussian_coeff(jd,m_max_degree_Earth_magnetic,g,h);

  std::vector<T> Bqm=pos;
  harmonics(pos,m_max_degree_Earth_magnetic,Re,g,h,Bqm);
  Bqm[0]*=-Q/x[6];
  Bqm[1]*=-Q/x[6];
  Bqm[2]*=-Q/x[6];
  
  aux[0]=posvel[4]*Bqm[2]-posvel[5]*Bqm[1];
  aux[1]=posvel[5]*Bqm[0]-posvel[3]*Bqm[2];
  aux[2]=posvel[3]*Bqm[1]-posvel[4]*Bqm[0];

  /** computing acceleration in inertial frame **/
  acc[0]=M[0]*aux[0]+M[1]*aux[1]+M[2]*aux[2];
  acc[1]=M[3]*aux[0]+M[4]*aux[1]+M[5]*aux[2];
  acc[2]=M[6]*aux[0]+M[7]*aux[1]+M[8]*aux[2];

  return 0;
}

template<class T>
int dearth_3dof<T>::charging(const double &jd, const std::vector<T> &x, T &C, T &V) const {

  /** sanity checks **/
  if(x.size()!=7)
      smartastro_throw("CHARGING: state must be a 7 dimensional vector");

  double qe=constants::qe;
  double kb=constants::kb;
  double eps0=constants::eps0;
  double omega=constants::omega_earth*m_T_scale;
  double pi=constants::pi;
  std::vector<T> pos, vel;
  std::vector<double> M(9,0.0);

  /*  position and velocity in Earth fixed frame  */
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

  /*  Computing the capacitance */   
  double ne=1.0e11, ni=ne;   
  T L=sqrt(eps0*kb/(ne/Te +ni/Ti))/qe; // approximated Debye length
  C=eps0*4.0*m_drag[0]*(1.0/L+1.0/sqrt(m_drag[0]/constants::pi)); // capacitance for spherical SC
  T Vsc=(m_L_scale/m_T_scale)*sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]); // relative velocity w.r.t. plasma
  V=surface_potential_LEO(Te,1.0,Vsc,zero); // anlytical model using IRI90

  #ifndef ENABLE_SMARTUQ // WARNING: because of this implementation trick, the fortran code cannot be called if SMART-UQ is ON
    /*  comment if not using numerical solver with IRI12 model in fortran */
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
    spacecraft.surface_area = 4.0*m_drag[0];
    spacecraft.shape_factor = 1.0;
    spacecraft.potential = 0.0;
    spacecraft.velocity = Vsc;
    double ShadowFunction = 0.0;
    if(m_flags[2]!=0)
        ShadowFunction=shadow(jd,x); 
    bool Auroral=false;      
    V = SCCal::solve<double>(res, spacecraft, ShadowFunction, Auroral);
    L=sqrt(eps0*kb/(ne/Te +ni/Ti))/qe; // recomputing Debye length with IRI12
    if(Z>2000.0)
      L = 10.0;
    C=eps0*4.0*m_drag[0]*(1.0/L+1.0/sqrt(m_drag[0]/constants::pi)); // recomputing capacitance with IRI12  
  #endif

  return 0;
}


template class dearth_3dof<double>;
template class dearth_3dof<float>;
template class dearth_3dof<long double>;
#ifdef ENABLE_SMARTUQ
template class dearth_3dof<smartuq::polynomial::chebyshev_polynomial>;
template class dearth_3dof<smartuq::polynomial::taylor_polynomial>;
#endif