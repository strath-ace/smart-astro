/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

/* This example shows how to use the not-so-well-tested semianalytical dynamics for Earth orbits */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "../include/smartastro.h"

using namespace std;

double sign(double a){

    if(a > 1.0e-6)
        return 1.0;
    else if(a < -1.0e-6)
        return -1.0;
    else
        return 0.0;
}


double atan3(double a, double b){

double epsilon = 0.0000000001;
double y, c;
double pidiv2 = 0.5 * smartastro::constants::pi;

if (fabs(a) < epsilon)
{
   y = (1.0 - sign(b)) * pidiv2;
   return y;
}
else
   c = (2.0 - sign(a)) * pidiv2;

if (fabs(b) < epsilon){
   y = c;
   return y;
}
else
   y = c + sign(a) * sign(b) * (fabs(atan(a / b)) - pidiv2);

return y;
}

bool mean2oscul2(const std::vector<double> &mean, const double &mu, const double & Re, std::vector<double> &oscul){

    double E0, M0;
    double si, ci, ti, si2, ci2, f, sf, cf, s2f, u, e, e2, esf, am;
    double d1, d2, d3, d4, d42, d5, d6, d7, d8, p, r, rdot;
    double twou, twow, s2u, c2u, sf2w, s3f2w, cf2w, c3f2w, q1, di, dp;
    double dummy1, dn, dr, drdot, du, pnw, ainw, annw, rnw, rdotnw, unw, aa, bb, enw2, enw, xfnw, anw, wnw;
    double pi = smartastro::constants::pi;
    double J2 = 0.00108262617385222;

    E0 = 2.0 * atan(sqrt((1.0 - mean[1]) / (1.0 + mean[1])) * tan(mean[5] / 2.0));
    if(E0 < 0.0)
        E0 += 2.0 * pi;
    M0 = E0 - mean[1] * sin(E0);
    am = M0;

    e = mean[1];

    f = mean[5];

    si = sin(mean[2]);
    ci = cos(mean[2]);
    ti = si / ci;
    si2 = si * si;
    ci2 = ci * ci;
     
    sf = sin(f);
    cf = cos(f);
    s2f = sin(2.0 * f);
    u = f + mean[4];

    e2 = e * e;
    esf = e * sf;

    d1 = 1.0 - e2;
    d2 = sqrt(d1);
    d3 = e * cf;
    d4 = 1.0 + d3;
    d42 = d4 * d4;
    d5 = 1.0 + d2;
    d6 = (3.0 * ci2 - 1.0) / d5;
    p = mean[0] * d1;
    d7 = sqrt(mu / p);
    r = p / d4;

    rdot = d7 * esf;

    twou = 2.0 * u;
    twow = 2.0 * mean[4];
    s2u = sin(twou);
    c2u = cos(twou);
    sf2w = sin(f + twow);
    d8 = 3.0 * f + twow;
    s3f2w = sin(d8);
    cf2w = cos(f + twow);
    c3f2w = cos(d8);

    q1 = J2 * (Re / p) * (Re / p);
    di = 0.75 * q1 * si * ci * (c2u + e * cf2w + e / 3.0 * c3f2w);
    dp = 2.0 * p * ti * di;
    dummy1 = f - am + esf - 0.5 * s2u - 0.5 * e * sf2w - e * s3f2w / 6.0;
    dn = -1.5 * q1 * ci * dummy1;

    dr = -0.25 * p * q1 * ((3.0 * ci2 - 1.0)
        * (2.0 * d2 / d4 + d3 / d5 + 1.0) - si2 * c2u);

    drdot = 0.25 * d7 * q1 * (d6 * esf * (d2 * d5 + d42)
        - 2.0 * si2 * d42 * s2u);

    du = -0.125 * q1 * (6.0 * (1.0 - 5.0 * ci2) 
         * (f - am) + 4.0 * esf * ((1.0 - 6.0 * ci2) - d6) 
         - d6 * e2 * s2f + 2.0 * (5.0 * ci2 - 2.0) * e * sf2w 
         + (7.0 * ci2 - 1.0) * s2u + 2.0 * ci2 * e * s3f2w);
    
    pnw = p + dp;

    ainw = mean[2] + di;
    annw = mean[3] + dn;

    rnw = r + dr;

    rdotnw = rdot + drdot;

    unw = u + du;

    aa = pnw / rnw - 1.0;
    bb = sqrt(pnw / mu) * rdotnw;

    enw2 = aa * aa + bb * bb;
    enw = sqrt(enw2);

    xfnw = atan3(bb, aa); 

    anw = pnw / (1.0 - enw2);
    wnw = unw - xfnw;

    oscul[0] = anw;
    oscul[1] = enw;
    oscul[2] = ainw;
    oscul[3] = annw - 2.0 * pi * (double)floor(annw / (2.0 * pi));
    oscul[4] = wnw - 2.0 * pi * (double)floor(wnw / (2.0 * pi));
    oscul[5] = xfnw;

    return 0;
}


int main(){

  std::ofstream semianalytical, symplectic;
  clock_t begin, end;

  double pi=smartastro::constants::pi;

  int steps1, steps2;
  int n_max = 2;

  /** scaling parameters */
  double R_ref=smartastro::constants::R_earth; // distance in meters
  double T_ref=smartastro::constants::T_Earth; // time in seconds

  std::vector<bool> flag_hamilton(1, true), flags_semi(4, false);
  flags_semi[1] = flag_hamilton[0];

  std::vector<double> x0(6), x00(6), x0bis(6), xfbis(6), xf(6), y0(6), yf(6), yfbis(6), y0bis(6), kep(6), car(6), modeq(6);// z0(5), zf(5);
  x0[0] = 7000.0e3 / R_ref;
  x0[1] = 1.0e-2;
  x0[2] = 80.0 * pi / 180.0;
  x0[3] = 10.0 * pi / 180.0;
  x0[4] = 60.0 * pi / 180.0;
  x0[5] = 200.0 * pi / 180.0;

  // x0[0] = 29600.0e3 / R_ref;
  // x0[1] = 0.0;
  // x0[2] = 56.0 * pi / 180.0;
  // x0[3] = 0.0 * pi / 180.0;
  // x0[4] = 0.0 * pi / 180.0;
  // x0[5] = 0.0 * pi / 180.0;

  smartastro::astrocore::conversion_coordinates::kep2car(x0, 1.0, y0);
  kep = x0;
  x0bis = x0;
  if(n_max > 1)
    smartastro::astrocore::conversion_coordinates::oscul2mean(x0, 1.0, 1.0, kep);

  // smartastro::astrocore::conversion_coordinates::mean2oscul(kep, 1.0, 1.0, x0);
  // std::cout << x0[0] * R_ref << " " << x0[1] << " " << x0[2] * 180.0 /pi << " " << x0[3] * 180.0 /pi  << " " << x0[4] * 180.0 /pi << " " << x0[5] * 180.0 /pi << std::endl;

  smartastro::astrocore::conversion_coordinates::kep2modeq(kep, x0bis);
  // std::cout << "Initial mean modified equinoctial: " << x0bis[0] * R_ref << " " << x0bis[1] << " " << x0bis[2] << " " << x0bis[3]  << " " << x0bis[4] << " " << x0bis[5] * 180.0 / pi << std::endl;

  smartastro::astrocore::conversion_coordinates::kep2modeq(x0, xf);
  x0 = xf;

  double mjd = 51624.5;
  double t = (mjd + smartastro::constants::mjd_zero_in_jd) * 3600.0 * 24.0 / T_ref;

  // smartastro::astrodynamics::dearth_mixedvar<double> *dyn = new smartastro::astrodynamics::dearth_mixedvar<double>(n_max, flag_hamilton);
  // smartmath::integrator::forest_mixedvar<double> prop_symplectic(dyn);
  smartastro::astrodynamics::dearth_hamiltonian<double> *dyn = new smartastro::astrodynamics::dearth_hamiltonian<double>(n_max, flag_hamilton);
  smartmath::integrator::yoshida6<double> prop_symplectic(dyn);  
  // smartmath::integrator::ABM<double> prop_symplectic(dyn);    

  smartastro::astrodynamics::dearth_semianalytical *dyn2 = new smartastro::astrodynamics::dearth_semianalytical(flags_semi, n_max);
  smartmath::integrator::rk4<double> prop_semianalytical(dyn2);

  semianalytical.open ("semianalytical.txt");
  symplectic.open ("symplectic.txt");

  double time_elapsed1 = 0.0, time_elapsed2 = 0.0;

  double dt = 24.0 * 3600.0 / T_ref;
  steps1 = 1 + floor(dt * T_ref / 30.0);
  steps2 = 1 + floor(dt * T_ref / (12.0 * 3600.0)); 

  int propagations = 1e2;
  for(int k = 0; k < propagations; k++){
    begin=clock();
    prop_symplectic.integrate(t, t + dt, steps1, y0, yf);
    end=clock();
    time_elapsed1 += (double (end-begin))/CLOCKS_PER_SEC;

    smartastro::astrocore::conversion_coordinates::car2kep(yf, 1.0, kep);
    smartastro::astrocore::conversion_coordinates::kep2modeq(kep, modeq);

    // modeq = kep;
    // smartastro::astrocore::conversion_coordinates::oscul2mean(modeq, 1.0, 1.0, kep);

    // symplectic << setprecision(16) << (t + dt) * T_ref / (3600.0 * 24.0) << " " << modeq[0] * R_ref << " " << modeq[1] << " " << modeq[2] << " " << modeq[3]  << " " << modeq[4] << " " << modeq[5] * 180.0 / pi << " ";
    symplectic << setprecision(16) << (t + dt) * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi << " " << kep[5] * 180.0 /pi << " ";

    begin=clock();
    prop_semianalytical.integrate(t, t + dt, steps2, x0bis, xfbis);
    end=clock();
    time_elapsed2 += (double (end-begin)) / CLOCKS_PER_SEC;
    x0bis = xfbis;

    smartastro::astrocore::conversion_coordinates::modeq2kep(x0bis, kep);

    // xfbis = kep;
    // xfbis[5] = modeq[5];
    // smartastro::astrocore::conversion_coordinates::mean2oscul(xfbis, 1.0, 1.0, kep);

    // semianalytical << setprecision(16) << (t + dt) * T_ref / (3600.0 * 24.0) << " " << x0bis[0] * R_ref << " " << x0bis[1] << " " << x0bis[2] << " " << x0bis[3]  << " " << x0bis[4] << " " << (x0bis[5] - 2.0 * pi * double(floor(x0bis[5] / (2.0 * pi)))) * 180.0 / pi << " ";
    semianalytical << setprecision(16) << (t + dt) * T_ref / (3600.0 * 24.0) << " " << kep[0] * R_ref << " " << kep[1] << " " << kep[2] * 180.0 /pi << " " << kep[3] * 180.0 /pi  << " " << kep[4] * 180.0 /pi << " " << kep[5] * 180.0 /pi << " ";

    symplectic << "\n";
    semianalytical << "\n";

    t += dt;
    y0 = yf;

    // std::cout << k + 1 << std::endl;
  }

  cout << "Time to obtain final state with symplectic : " << time_elapsed1 << " VS semianalytical propagator: " << time_elapsed2 << endl;

  symplectic.close();
  semianalytical.close();

  delete dyn;
  delete dyn2;

  return 0;

}

