/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: Romain Serra (GitHub: Serrof)-----------------------
-------- e-mail: serra.romain@gmail.com ------------------------------
*/

#include "../../include/Ephemerides/analytical_planets.h"

using namespace smartastro;
using namespace ephemerides;


bool analytical_planets::get_orbel(const double &mjd, const int &index_planet, std::vector<double> &kep){

    /** sanity checks **/
    if((index_planet<1)||(index_planet>8))
        smartastro_throw("ANALYTICAL_PLANETS: index of planet must be between 1 and 8");
    if(kep.size()!=6)
        smartastro_throw("ANALYTICAL_PLANETS: vector of orbital elements must have 6 components");    

    double XM;

    // Convert input date to MJD2000
    double mjd2000;
    astrocore::conversion_time::mjd2mjd2000(mjd,mjd2000);

    //  T = JULIAN CENTURIES SINCE 31/12/1899 at 12:00
    double T   = (mjd2000 + 36525.0)/36525.0;
    double TT  = T*T;
    double TTT = T*TT;

    // CLASSICAL PLANETARY ELEMENTS ESTIMATION IN MEAN ECLIPTIC OF DATE
    switch(index_planet){
    
    //  MERCURY
    case 1:
        kep[0] = 0.38709860;
        kep[1] = 0.205614210 + 0.000020460*T - 0.000000030*TT;
        kep[2] = 7.002880555555555560 + 1.86083333333333333e-3*T - 1.83333333333333333e-5*TT;
        kep[3] = 4.71459444444444444e+1 + 1.185208333333333330*T + 1.73888888888888889e-4*TT;
        kep[4] = 2.87537527777777778e+1 + 3.70280555555555556e-1*T +1.20833333333333333e-4*TT;
        XM   = 1.49472515288888889e+5 + 6.38888888888888889e-6*T;
        kep[5] = 1.02279380555555556e2 + XM*T;
    break;
    //  VENUS
    case 2:
        kep[0] = 0.72333160;
        kep[1] = 0.006820690 - 0.000047740*T + 0.0000000910*TT;
        kep[2] = 3.393630555555555560 + 1.00583333333333333e-3*T - 9.72222222222222222e-7*TT;
        kep[3] = 7.57796472222222222e+1 + 8.9985e-1*T + 4.1e-4*TT;
        kep[4] = 5.43841861111111111e+1 + 5.08186111111111111e-1*T -1.38638888888888889e-3*TT;
        XM   = 5.8517803875e+4 + 1.28605555555555556e-3*T;
        kep[5] = 2.12603219444444444e2 + XM*T;
    break;    
    //  EARTH
    case 3:
        kep[0] = 1.000000230;
        kep[1] = 0.016751040 - 0.000041800*T - 0.0000001260*TT;
        kep[2] = 0.00;
        kep[3] = 0.00;
        kep[4] = 1.01220833333333333e+2 + 1.7191750*T + 4.52777777777777778e-4*TT + 3.33333333333333333e-6*TTT;
        XM   = 3.599904975e+4 - 1.50277777777777778e-4*T - 3.33333333333333333e-6*TT;
        kep[5] = 3.58475844444444444e2 + XM*T;
    break;   
    //  MARS
    case 4:
        kep[0] = 1.5236883990;
        kep[1] = 0.093312900 + 0.0000920640*T - 0.0000000770*TT;
        kep[2] = 1.850333333333333330 - 6.75e-4*T + 1.26111111111111111e-5*TT;
        kep[3] = 4.87864416666666667e+1 + 7.70991666666666667e-1*T - 1.38888888888888889e-6*TT - 5.33333333333333333e-6*TTT;
        kep[4] = 2.85431761111111111e+2 + 1.069766666666666670*T +  1.3125e-4*TT + 4.13888888888888889e-6*TTT;
        XM   = 1.91398585e+4 + 1.80805555555555556e-4*T + 1.19444444444444444e-6*TT;
        kep[5] = 3.19529425e2 + XM*T;
    break;    
    //  JUPITER
    case 5:
        kep[0] = 5.2025610;
        kep[1] = 0.048334750 + 0.000164180*T  - 0.00000046760*TT -0.00000000170*TTT;
        kep[2] = 1.308736111111111110 - 5.69611111111111111e-3*T +  3.88888888888888889e-6*TT;
        kep[3] = 9.94433861111111111e+1 + 1.010530*T + 3.52222222222222222e-4*TT - 8.51111111111111111e-6*TTT;
        kep[4] = 2.73277541666666667e+2 + 5.99431666666666667e-1*T + 7.0405e-4*TT + 5.07777777777777778e-6*TTT;
        XM   = 3.03469202388888889e+3 - 7.21588888888888889e-4*T + 1.78444444444444444e-6*TT;
        kep[5] = 2.25328327777777778e2 + XM*T;
    break;    
    //  SATURN
    case 6:
        kep[0] = 9.5547470;
        kep[1] = 0.055892320 - 0.00034550*T - 0.0000007280*TT + 0.000000000740*TTT;
        kep[2] = 2.492519444444444440 - 3.91888888888888889e-3*T - 1.54888888888888889e-5*TT + 4.44444444444444444e-8*TTT;
        kep[3] = 1.12790388888888889e+2 + 8.73195138888888889e-1*T -1.52180555555555556e-4*TT - 5.30555555555555556e-6*TTT;
        kep[4] = 3.38307772222222222e+2 + 1.085220694444444440*T + 9.78541666666666667e-4*TT + 9.91666666666666667e-6*TTT;
        XM   = 1.22155146777777778e+3 - 5.01819444444444444e-4*T - 5.19444444444444444e-6*TT;
        kep[5] = 1.75466216666666667e2 + XM*T;
    break;    
    //  URANUS
    case 7:
        kep[0] = 19.218140;
        kep[1] = 0.04634440 - 0.000026580*T + 0.0000000770*TT;
        kep[2] = 7.72463888888888889e-1 + 6.25277777777777778e-4*T + 3.95e-5*TT;
        kep[3] = 7.34770972222222222e+1 + 4.98667777777777778e-1*T + 1.31166666666666667e-3*TT;
        kep[4] = 9.80715527777777778e+1 + 9.85765e-1*T - 1.07447222222222222e-3*TT - 6.05555555555555556e-7*TTT;
        XM   = 4.28379113055555556e+2 + 7.88444444444444444e-5*T + 1.11111111111111111e-9*TT;
        kep[5] = 7.26488194444444444e1 + XM*T;
    break;    
    //  NEPTUNE
    case 8:
        kep[0] = 30.109570;
        kep[1] = 0.008997040 + 0.0000063300*T - 0.0000000020*TT;
        kep[2] = 1.779241666666666670 - 9.54361111111111111e-3*T - 9.11111111111111111e-6*TT;
        kep[3] = 1.30681358333333333e+2 + 1.0989350*T + 2.49866666666666667e-4*TT - 4.71777777777777778e-6*TTT;
        kep[4] = 2.76045966666666667e+2 + 3.25639444444444444e-1*T + 1.4095e-4*TT + 4.11333333333333333e-6*TTT;
        XM   = 2.18461339722222222e+2 - 7.03333333333333333e-5*T;
        kep[5] = 3.77306694444444444e1 + XM*T;
    break;
    }    
  

    // CONVERSION OF AU INTO KM, DEG INTO RAD AND MEAN TO TRUE ANOMALY
    kep[0]*=constants::au;    
    for(int i=2;i<6;i++){
        kep[i]*=constants::deg2rad;
    }
    double mean_anomaly = kep[5] - 2.0*constants::pi*floor(kep[5]/(2.0*constants::pi)), true_anomaly;
    astrocore::conversion_time::mean2true_anomaly(mean_anomaly,kep[1],true_anomaly);
    kep[5] = true_anomaly;

    return 0;
}
