/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/

#include "Astrodynamics/cr3bp_3dof.h"
#include <vector>
#include <typeinfo>

using namespace smartastro;
using namespace smartastro::astrodynamics;

template<class T>
cr3bp_3dof_test<T>::cr3bp_3dof_test(const double &mu1) :
    base_astrodynamics<T>("Circular three body problem in rotating frame"), m_mu1(mu1)
{

}


template<class T>
cr3bp_3dof_test<T>::~cr3bp_3dof_test(){

}


template<class T>
int cr3bp_3dof_test<T>::evaluate(const double &t, const std::vector<T> &x,std::vector<T> &dx) const{

   /** sanity checks **/
   if(x.size()!=7)
       smartastro_throw("state must be a 7 dimensional vector");
   if(dx.size()!=7)
       smartastro_throw("state must be a 7 dimensional vector");

   std::vector<T> acc_grav;
   T zero=0.0*x[0];
   double Omega=2.0*smartastro::constants::pi;

   acc_grav.push_back(zero);
   acc_grav.push_back(zero);
   acc_grav.push_back(zero);
   gravity(x,acc_grav);

   /** constructing state derivative **/
   /** position **/
   for(int i=0; i<3; i++){
       dx[i] = x[3+i];
   }
   /** velocity **/
   dx[3]=Omega*Omega*x[0]+2.0*Omega*x[4];
   dx[4]=Omega*Omega*x[1]-2.0*Omega*x[3];  
   dx[5]=0.0*x[0];
   for(int i=0; i<3; i++){
      dx[3+i] += acc_grav[i];
   }
   /** mass **/
   dx[6]=0.0*x[0];

   return 0;
}

template<class T>
int cr3bp_3dof_test<T>::gravity(const std::vector<T> &x, std::vector<T> &acc) const{

   /** sanity checks **/
   if(x.size()!=7)
       smartastro_throw("state must be a 7 dimensional vector");
   if(acc.size()!=3)
       smartastro_throw("acceleration must be a 3 dimensional vector");

     double m1=m_mu1;
     double x1=-1.0+m1;
     double m2=1.0-m1;
     double x2=m1;
     double G=4.0*smartastro::constants::pi*smartastro::constants::pi;

     T r1=sqrt((x[0]-x1)*(x[0]-x1)+x[1]*x[1]+x[2]*x[2]);
     T r2=sqrt((x[0]-x2)*(x[0]-x2)+x[1]*x[1]+x[2]*x[2]);
     T inv_r1_cube=1.0/pow(r1,3);
     T inv_r2_cube=1.0/pow(r2,3);

     acc[0]=-G*(m1*(x[0]-x1)*inv_r1_cube+m2*(x[0]-x2)*inv_r2_cube);
     acc[1]=-x[1]*G*(m1*inv_r1_cube+m2*inv_r2_cube);
     acc[2]=-x[2]*G*(m1*inv_r1_cube+m2*inv_r2_cube);

   return 0;
}


template class cr3bp_3dof_test<double>;
template class cr3bp_3dof_test<float>;
template class cr3bp_3dof_test<long double>;
#ifdef ENABLE_SMARTUQ
//template class cr3bp_3dof_test<smartuq::polynomial::chebyshev_polynomial>;
//template class cr3bp_3dof_test<smartuq::polynomial::taylor_polynomial>;
#endif

