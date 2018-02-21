/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/

#include "Astro-Core/conversion_time.h"

using namespace smartastro;
using namespace astrocore;

bool conversion_time::hms2fracday(const std::vector<double> &hms, double &fracday){
    /*Sanity check*/
    if(hms.size()!=3)
        smartastro_throw("HMS2FRACDAY: time vector must be three dimensional (h,m,s)");
    if(hms[0]<0 || hms[0] >=24)
        smartastro_throw("HMS2FRACDAY: hours in 24-hours format");
    if(hms[1]<0 || hms[1] >=60)
        smartastro_throw("HMS2FRACDAY: minutes between 0 and 60");
    if(hms[2]<0 || hms[2] >=60)
        smartastro_throw("HMS2FRACDAY: seconds between 0 and 60");

    fracday = (hms[0] + (hms[1] + hms[2]/60)/60) / 24;

    if (fracday <0)
        return 1;

    return 0;
}


bool conversion_time::fracday2hms(const double &fracday, std::vector<double> &hms){

    /*Sanity check*/
    if(fracday<0 || fracday>1)
        smartastro_throw("FRACDAY2HMS: fraction of day must be a positive number between 0 and 1");
    if(hms.size()!=3)
        smartastro_throw("FRACDAY2HMS: time vector must be three dimensional (h,m,s)");

    double temp;

    temp = fracday*24;

    /* Hours */
    if(temp>=0)
        hms[0] = floor(temp);
    else
        hms[0] = floor(temp)+1;

    /* Minutes */
    temp -= hms[0];
    if(temp>=0)
        hms[1] = floor(temp*60);
    else
        hms[1] = floor(temp*60)+1;

    /* Seconds */
    hms[2] = (temp-hms[1]/60)*3600;

    if(hms[0]<0 || hms[0] >=24)
        return 1;
    if(hms[1]<0 || hms[1] >=60)
        return 1;
    if(hms[2]<0 || hms[2] >=60)
        return 1;

    return 0;
}


bool conversion_time::jd2mjd(const double &jd, double &mjd){

    /*Sanity check*/
    if(jd<0)
        smartastro_throw("JD2MJD: julian date must be positive");

    mjd = jd - constants::mjd_zero_in_jd;

    return 0;
}


bool conversion_time::mjd2jd(const double &mjd, double &jd){

    jd = mjd + constants::mjd_zero_in_jd;

    if(jd<0)
        return 1;

    return 0;

}


bool conversion_time::jd2mjd2000(const double &jd, double &mjd2000){

    /*Sanity check*/
    if(jd<0)
        smartastro_throw("JD2MJD2000: julian date must be positive");

    /*                      mjd2000
     *    |-----------------------------------------|
     *             mjd
     *    |-------------------|                       */
    mjd2000 = jd - constants::mjd_zero_in_jd - constants::mjd2000_zero_in_mjd;

    return 0;
}


bool conversion_time::mjd20002jd(const double &mjd2000, double &jd){

    /*                           jd
     *    |----------------------------------------------|
     *                  mjd
     *    |-----------------------------|                  */
    jd = mjd2000 + constants::mjd2000_zero_in_mjd + constants::mjd_zero_in_jd;

    if(jd<0)
        return 1;

    return 0;
}


bool conversion_time::mjd20002mjd(const double &mjd2000, double &mjd){

    mjd = mjd2000 + constants::mjd2000_zero_in_mjd;

    return 0;
}


bool conversion_time::mjd2mjd2000(const double &mjd, double &mjd2000){

    mjd2000 = mjd - constants::mjd2000_zero_in_mjd;

    return 0;
}


bool conversion_time::date2jd(const std::vector<double> &date, double &jd){

    /*Sanity check*/
    if(date.size()!=6)
        smartastro_throw("DATE2JD: date must be a 6-dimensional vector (Y,M,D,hrs,mn,sec)");
    if(date[0]<=0)
        smartastro_throw("DATE2JD: year must be a positive value");
    if(date[1]<=0 || date[1]>12)
        smartastro_throw("DATE2JD: month must be a positive value between 1 and 12");
    if(date[2]<=0 || date[2]>31)
        smartastro_throw("DATE2JD: day must be a positive value between 1 and 31");
    if(date[3]<0 || date[3] >=24)
        smartastro_throw("DATE2JD: hours in 24-hours format");
    if(date[4]<0 || date[4] >=60)
        smartastro_throw("DATE2JD: minutes between 0 and 60");
    if(date[5]<0 || date[5] >=60)
        smartastro_throw("DATE2JD: seconds between 0 and 60");

    double Y,M,D;
    std::vector<double> time(date.begin()+3,date.end());
    double frac;

    /* Manage the input */
    Y   = date[0];
    M   = date[1];
    D   = date[2];

    /* Formula converting Gregorian date into JD */

    hms2fracday(time,frac);

    jd = 367*Y - floor(7*(Y+floor((M+9)/12))/4)
         - floor(3*floor((Y+(M-9)/7)/100+1)/4)
         + floor(275*M/9)
         + D + 1721028.5 + frac;

    /* Check if the input date was valid */
    if (jd<0)
        return 1;

    return 0;

}


bool conversion_time::jd2date(const double &jd, std::vector<double> &date){

    if(jd<0)
        smartastro_throw("JD2DATE: julian date must be positive");
    if(date.size()!=6)
        smartastro_throw("JD2DATE: date must be a 6-dimensional vector (Y,M,D,hrs,mn,sec)");

    double j, g, dg, c, dc, b, db, a, da, y, m, d;

    /* Calculations */
    j = floor(jd+0.5) + 32044;
    g = floor(j/146097);
    dg = fmod(j,146097);
    c = floor((floor(dg/36524)+1) * 3/4);
    dc = dg - c*36524;
    b = floor(dc/1461);
    db = fmod(dc,1461);
    a = floor((floor(db/365)+1) * 3/4);
    da = db - a*365;
    y = g*400 + c*100 + b*4 + a;
    m = floor((da*5 + 308)/153) - 2;
    d = da - floor((m+4)*153/5) + 122;

    /* Year, Month and Day */
    date[0] = y-4800 + floor((m+2)/12);
    date[1] = fmod((m+2),12) + 1;
    date[2] = floor(d+1);

    /* Hour, Minute and Second */
    std::vector<double> time(3);
    fracday2hms(fmod(jd+0.5,floor(jd+0.5)),time);

    for(size_t i=0; i<3;i++)
        date[3+i]=time[i];

    if(date[0]<=0)
        return 1;
    if(date[1]<=0 || date[1]>12)
        return 1;
    if(date[2]<=0 || date[2]>31)
        return 1;
    if(date[3]<0 || date[3] >=24)
        return 1;
    if(date[4]<0 || date[4] >=60)
        return 1;
    if(date[5]<0 || date[5] >=60)
        return 1;

    return 0;
}


bool conversion_time::date2mjd(const std::vector<double> &date, double &mjd){

    double jd;
    bool error;

    error = conversion_time::date2jd(date, jd);

    if(error==0)
        error = jd2mjd(jd,mjd);

    return error;
}


bool conversion_time::mjd2date(const double &mjd, std::vector<double> &date){

    double jd;
    conversion_time::mjd2jd(mjd,jd);

    return conversion_time::jd2date(jd,date);
}


bool conversion_time::date2mjd2000(const std::vector<double> &date, double &mjd2000){

    double jd;
    bool error;

    /* Converte date to JD */
    error = conversion_time::date2jd(date,jd);
    if(error == 0)
        mjd2000 = jd - constants::mjd_zero_in_jd - constants::mjd2000_zero_in_mjd;

    return error;
}


bool conversion_time::mjd20002date(const double &mjd2000, std::vector<double> &date){

    double jd;
    conversion_time::mjd20002jd(mjd2000,jd);

    return conversion_time::jd2date(jd,date);
}

bool conversion_time::mean2eccentric_anomaly(const double &M, const double &eccentricity, double &E){

    if((eccentricity<0.0)||(eccentricity>=1.0))
        smartastro_throw("MEAN2ECCENTRIC_ANOMALY: eccentricity must be between 0 and 1");

    double pi = constants::pi;
    int n_rev = int(floor(M/(2.0*pi)));
    double M0 = M - 2.0*pi*double(n_rev); // mean anomaly modulo 2PI

    E = M0; // initial guess
    double F_E, DF_E; 

    /* Newton-Raphson */
    for(int k=0; k<10; k++)
    {
        F_E = E - eccentricity*sin(E) - M0;
        DF_E = 1.0 - eccentricity*cos(E);
        E -= F_E/DF_E;      
    }

    if((E<0.0)||(E>=2.0*pi))
        return false;

    /* adding the revolutions */
    E += 2.0*pi*double(n_rev);

    return true;
}

bool conversion_time::eccentric2mean_anomaly(const double &E, const double &eccentricity, double &M){

    if((eccentricity<0.0)||(eccentricity>=1.0))
        smartastro_throw("ECCENTRIC2MEAN_ANOMALY: eccentricity must be between 0 and 1");

    M = E - eccentricity*sin(E); // Kepler equation

    return true;
}

bool conversion_time::true2eccentric_anomaly(const double &theta, const double &eccentricity, double &E){

    if((eccentricity<0.0)||(eccentricity>=1.0))
        smartastro_throw("TRUE2ECCENTRIC_ANOMALY: eccentricity must be between 0 and 1");

    double pi = constants::pi;
    int n_rev = floor(theta/(2.0*pi));
    double theta0 = theta - 2.0*pi*double(n_rev); // true anomaly modulo 2PI
    if(theta0 > pi)
        theta0 -= 2.0*pi;

    E = 2.0 * atan( sqrt((1.0 - eccentricity)/(1.0 + eccentricity)) * tan(theta0/2.0));

    if(E<0.0)
        E += 2.0*pi;

    E += 2.0*pi*double(n_rev); //  adding the revolutions

    return true;
}

bool conversion_time::eccentric2true_anomaly(const double &E, const double &eccentricity, double &theta){

    if((eccentricity<0.0)||(eccentricity>=1.0))
        smartastro_throw("ECCENTRIC2TRUE_ANOMALY: eccentricity must be between 0 and 1");

    double pi = constants::pi;
    int n_rev = floor(E/(2.0*pi));
    double E0 = E - 2.0*pi*double(n_rev); // eccentric anomaly modulo 2PI
    if(E0 > pi)
        E0 -= 2.0*pi;

    theta = 2.0 * atan( sqrt((1.0 + eccentricity)/(1.0 - eccentricity)) * tan(E0/2.0));

    if(theta<0.0)
        theta += 2.0*pi;

    theta += 2.0*pi*double(n_rev); // adding the revolutions

    return true;
}

bool conversion_time::true2mean_anomaly(const double &theta, const double &eccentricity, double &M){

    if((eccentricity<0.0)||(eccentricity>=1.0))
        smartastro_throw("TRUE2MEAN_ANOMALY: eccentricity must be between 0 and 1");

    double E;
    true2eccentric_anomaly(theta,eccentricity,E);
    eccentric2mean_anomaly(E,eccentricity,M);

    return true;
}


bool conversion_time::mean2true_anomaly(const double &M, const double &eccentricity, double &theta){

    if((eccentricity<0.0)||(eccentricity>=1.0))
        smartastro_throw("MEAN2TRUE_ANOMALY: eccentricity must be between 0 and 1");

    double E;
    mean2eccentric_anomaly(M,eccentricity,E);
    eccentric2true_anomaly(E,eccentricity,theta);

    return true;
}
