/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "Astrodynamics/dearth_orb.h"
#include "Astro-Core/conversion_coordinates.h"
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
dearth_orb<T>::dearth_orb(const double &L_scale, const double &t_scale, const int &n_gravity, const T &mass, const std::vector<int> &flags, const std::vector<T> &p1, const std::vector<T> &p2, const std::vector< std::vector<double> > &F10dot7) :
    base_dearth<T>("Earth-centered dynamics using orbital elements", L_scale, t_scale, n_gravity, flags, p1, p2, F10dot7), m_mass(mass){

  m_mu = smartastro::constants::mu_earth * pow(t_scale, 2) / pow(L_scale, 3);
  m_dyn = new smartastro::astrodynamics::dearth_3dof<double>(m_L_scale, m_T_scale, m_max_degree_Earth_gravity, m_flags, m_drag, m_other, F10dot7);

};


template<class T>
dearth_orb<T>::~dearth_orb(){

}


template<class T>
int dearth_orb<T>::evaluate(const double &t, const std::vector<T> &x, std::vector<T> &dx) const{

  /** sanity checks **/
  if(x.size()!=6)
    smartastro_exception("EVALUATE: state must be a 6 dimensional vector");
  if(dx.size()!=6)
    smartastro_exception("EVALUATE: state derivative must be a 6 dimensional vector");

  /** computing position-velocity-mass vector **/
  double pi = constants::pi;
  std::vector<T> dir_rad(3), dir_vel(3), dir_h(3), dir_ortho(3);
  std::vector<T> kep(6), car(6), x_copy = x;
  x_copy[5] -= 2.0 * pi * double(floor(x_copy[5] / (2.0 * pi)));
  astrocore::conversion_coordinates::modeq2car(x_copy, m_mu, car);
  std::vector<T> posvelmass = car;
  posvelmass.push_back(m_mass);

  /** getting perturbed acceleration from Cartesian model **/
  std::vector<T> dstate(7), acc(3);
  m_dyn->evaluate(t, posvelmass, dstate);
  for(int i = 0; i < 3; i++)
    acc[i] = dstate[i + 3];

  /** computing state derivative from Gauss equations **/
  std::vector<T> RTH(3);
  astrocore::conversion_coordinates::car2rth(acc, car, RTH);
  RTH[0] +=  m_mu / (car[0] * car[0] + car[1] * car[1] + car[2] * car[2]);
  gauss_modeq(x, m_mu, RTH, dx);

  return 0;
}

template<class T>
int dearth_orb<T>::gauss_modeq(const std::vector<T> &x, const double &mu, const std::vector<T> &perturb, std::vector<T> &dx){

  /** sanity checks **/
  if(x.size()!=6)
    smartastro_exception("GAUSS_MODEQ: state must be a 6 dimensional vector");
  if(dx.size()!=6)
    smartastro_exception("GAUSS_MODEQ: state derivative must be a 6 dimensional vector");
  if(perturb.size()!=3)
    smartastro_exception("GAUSS_MODEQ: acceleration from perturbations must be a 3 dimensional vector");  

  T S = perturb[0], C = perturb[1], N =  perturb[2]; // getting components of non-Keplerian acceleration

  T L = x[5] - 2.0 * constants::pi * double(floor(x[5] / (2.0 * constants::pi)));  
  T cL = cos(L);
  T sL = sin(L);
  T w = 1.0 + x[1] * cL + x[2] * sL;
  T s2 = 1.0 + x[3] * x[3] + x[4] * x[4];
  T factor = sqrt(x[0] / mu);
  dx[0] = factor * 2.0 * x[0] * C / w;
  T inter1 = w + 1.0;
  T inter2 = (x[3] * sL - x[4] * cL);
  dx[1] = factor * (S * sL + ((inter1 * cL + x[1]) * C - inter2 * x[2] * N) / w);
  dx[2] = factor * (- S * cL + ((inter1 * sL + x[2]) * C + inter2 * x[1] * N) / w);
  T aux = factor * 0.5 * s2 * N / w;
  dx[3] = aux * cL;
  dx[4] = aux * sL;
  dx[5] = factor * (mu * (w / x[0]) * (w / x[0]) + inter2 * N / w);

  return 0;
}

template<class T>
int dearth_orb<T>::charging(const double &jd, const std::vector<T> &x, T &C, T &V) const{

  /** sanity checks **/
  if(x.size()!=6)
      smartastro_throw("CHARGING: state must be a 6 dimensional vector");

  std::vector<double> posvelmass(6), kep(6);
  smartastro::astrocore::conversion_coordinates::modeq2car(x, m_mu, posvelmass); 
  posvelmass.push_back(m_mass);
  m_dyn->charging(jd, posvelmass, C, V);

  return 0;
}


template class dearth_orb<double>;
