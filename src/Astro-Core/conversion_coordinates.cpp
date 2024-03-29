/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2018 University of Strathclyde and Authors-------
--------- Author: SMART team----------------------
-------- e-mail: smart@strath.ac.uk ------------------------------
*/

#include "Astro-Core/conversion_coordinates.h"

using namespace smartastro;
using namespace astrocore;

bool conversion_coordinates::car2kep(const std::vector<double> &car, const double &mu, std::vector<double> &kep){

    /*Sanity check*/
    if(kep.size()!=6){
        smartastro_throw("CAR2KEP: keplerian vector must be 6-dimensional");
    }
    if(car.size()!=6){
        smartastro_throw("CAR2KEP: cartesian vector must be 6-dimensional");
    }

    /* Declaration */
    Eigen::Vector3d r, v;/* Cartesian position and velocity */
    double nr;                  /* Norm of r */
    Eigen::Vector3d h;          /* cross(r,v): Angular momentum vector */
    double nh;                  /* Norm of h */
    double plocal;              /* local value of p */
    Eigen::Vector3d e;          /* Eccentricity vector */
    double ne;                  /* Eccentricity */
    Eigen::Vector3d e_conv;     /* Arbitrary eccentricity vector */
    double a;                   /* Semi-major axis */
    double inc;                 /* inclination */
    Eigen::Vector3d n;          /* mean motion */
    double nn;                  /* norm of mean motion */
    double Om;                  /* Argument of the ascending node */
    double om;                  /* Argument of the pericentre */
    double th;                  /* True anomaly */

    r << car[0],car[1],car[2];
    v << car[3],car[4],car[5];

    /* Norm of r */
    nr = r.norm();

    /* cross(r,v): Angular momentum vector */
    h=r.cross(v);

    /* Norm of h */
    nh = h.norm();

    plocal = pow(nh,2.0) / mu;

    /* Eccentricity vector */
    e[0] = (v[1]*h[2]-v[2]*h[1]) / mu - r[0]/nr;
    e[1] = (v[2]*h[0]-v[0]*h[2]) / mu - r[1]/nr;
    e[2] = (v[0]*h[1]-v[1]*h[0]) / mu - r[2]/nr;

    /* Eccentricity */
    ne = e.norm();

    /* Semi-major axis */
    if ((ne<(1+constants::tol_ecc)) & (ne>(1-constants::tol_ecc)))        /* Parabola */
    {
        a = HUGE_VAL;
    }
    else
        a = plocal / (1-pow(ne,2.0));

    /* Inclination */
    inc = acos(h[2]/nh);

    /* Line of nodes vector */
    if (inc != 0.0 && inc != acos(-1)) /* Orbit is out of xy plane */
    {
        /* n=cross([0 0 1],h); n=n/norm(n); */
        nn = h.head(2).norm();   /* 2 is correct even if h has 3 elements */
        n[0] = -h[1] / nn;
        n[1] =  h[0] / nn;
        n[2] = 0.0;
    }
    else /* Orbit is in xy plane: n is not defined */
    {
        /* Arbitrary choice */
        n[0] = 1.0;
        n[1] = 0.0;
        n[2] = 0.0;
    }

    /* Argument of the ascending node */
    Om = acos(n[0]);
    if (n[1]<0)
        Om = constants::pi2 - Om;

    /* Argument of the pericentre and True anomaly */
    if (ne>=constants::tol_ecc) /* Non circular orbit */
    {
        /* Argument of the pericentre */
        om = acos(n.dot(e)/ne);
        if(inc != 0.0 && inc != acos(-1)){ /* Orbit is out of xy plane */
            if(e[2]<0)
                om = constants::pi2 - om;
        }
        else /* Orbit is in xy plane */
            if((h[2] > 0 && e[1] < 0) || (h[2] < 0 && e[1] > 0))
                om = constants::pi2 - om;

        /* True anomaly */
        th = acos(std::min(std::max(e.dot(r)/ne/nr,-1.0),1.0));
        if(r.dot(v)<0)
            th = constants::pi2 - th;
    }
    else /* Circular orbit: e is not defined */
    {
        /* Arbitrary eccentricity vector */
        e_conv = n;
        om = 0.0;

        /* True anomaly */
        th = acos(e_conv.dot(r)/nr);
        if(r.dot(v)<0)
            th = constants::pi2 - th;
    }

    /* Assigne the elements of the output kep */
    kep[0] = a;
    kep[1] = ne;
    kep[2] = inc;
    kep[3] = Om;
    kep[4] = om;
    kep[5] = th;

    return 0;
}

bool conversion_coordinates::car2bpl(const std::vector<double> &x_car, const std::vector<double> &u_car, const std::vector<double> &vp_car, std::vector<double> &x_bpl){

    /*Sanity check*/
    if(x_car.size()!=3){
        smartastro_throw("CAR2BPL: cartesian position vector must be 3-dimensional");
    }
    if(u_car.size()!=3){
        smartastro_throw("CAR2BPL: cartesian velocity vector must be 3-dimensional");
    }
    if(vp_car.size()!=3){
        smartastro_throw("CAR2BPL: planet velocity vector must be 3-dimensional");
    }
    if(x_bpl.size()!=3){
        smartastro_throw("CAR2BPL: position vector in the BPL frame must be 3-dimensional");
    }

    Eigen::Vector3d nn;      /* normalised vector of U_cart */
    Eigen::Vector3d ee;
    Eigen::Vector3d cc;
    Eigen::Matrix3d mat;     /* Transformation matrix */
    Eigen::Vector3d temp;    /* temporal vector */
    Eigen::Vector3d x_bpl_tmp;
    Eigen::Vector3d vp_car_tmp(vp_car.data());
    Eigen::Vector3d u_car_tmp(u_car.data());
    Eigen::Vector3d x_car_tmp(x_car.data());

    /* nn */
    nn = u_car_tmp/u_car_tmp.norm();

    /* ee */
    temp = vp_car_tmp.cross(nn);
    ee = temp/temp.norm();

    /* cc */
    cc = ee.cross(nn);

    /* Transformation matrix */
    mat.row(0) = ee;
    mat.row(1) = nn;
    mat.row(2) = cc;

    /* Apply transformation */
    x_bpl_tmp = mat*x_car_tmp;

    for(size_t i=0; i<3; i++)
        x_bpl[i] = x_bpl_tmp[i];

    return 0;

}

bool conversion_coordinates::rth2bpl(const std::vector<double> &x_rth,
                                     const std::vector<double> &s_car,
                                     const std::vector<double> &vp_car,
                                     std::vector<double> &x_bpl)
{
    bool rc;

    //Sanity checks
    if (x_rth.size() != 3)
    {
        smartastro_throw("RTH2BPL: x_car must be a 3-dimensional vector");
    }
    if (s_car.size()!=6)
    {
        smartastro_throw("RTH2BPL: s_car must be a 6-dimensional vector");
    }
    if (vp_car.size()!=3)
    {
        smartastro_throw("RTH2BPL: vp_car must be a 3-dimensional vector");
    }
    if (x_bpl.size()!=3)
    {
        smartastro_throw("RTH2BPL: x_bpl must be a 3-dimensional vector");
    }

    //Transform x rth to cartesian coordinates
    std::vector<double> x_car(3), v_car(3), u_car(3);

    rth2car(x_rth,s_car,x_car);

    v_car[0] = s_car[3];
    v_car[1] = s_car[4];
    v_car[2] = s_car[5];

    for (size_t i = 0; i< 3; i++)
    {
        u_car[i] = v_car[i] - vp_car[i];
    }

    rc = car2bpl(x_car,u_car,vp_car,x_bpl);

    return rc;
}

bool conversion_coordinates::car2rth(const std::vector<double> &x_car, const std::vector<double> &s_car, std::vector<double> &x_rth){

    /*Sanity check*/
    if(x_car.size()!=3){
        smartastro_throw("CAR2RTH: input position vector must be 3-dimensional");
    }
    if(s_car.size()!=6){
        smartastro_throw("CAR2RTH: input position and velocity vector must be 6-dimensional");
    }
    if(x_rth.size()!=3){
        smartastro_throw("CAR2RTH: output position vector must be 3-dimensional");
    }

    Eigen::Vector3d r,v;     /* Position and velocity vector */
    Eigen::Vector3d h;
    Eigen::Vector3d rn,hn,tn;   /* normalised vectors */
    Eigen::Matrix3d mat;            /* transformation matrix */
    double rnorm = 0.0, hnorm = 0.0;
    Eigen::Vector3d x_rth_tmp, x_car_tmp(x_car.data());

    /* Copy the position and velocity vector */
    r[0] = s_car[0];
    r[1] = s_car[1];
    r[2] = s_car[2];
    v[0] = s_car[3];
    v[1] = s_car[4];
    v[2] = s_car[5];

    /* Compute rn */
    rnorm = r.norm();
    if(rnorm==0)
        return 1;
    rn = r/rnorm;

    /* Compute h and hn */
    h = r.cross(v);
    hnorm = h.norm();
    if(hnorm==0)
        return 1;
    hn = h/hnorm;

    /* Compute tn */
    tn = hn.cross(rn);

    /* Transformation matrix */
    mat.row(0) = rn;
    mat.row(1) = tn;
    mat.row(2) = hn;

    /* Apply transformation */
    x_rth_tmp = mat*x_car_tmp;

    for(size_t i=0; i<3; i++)
        x_rth[i] = x_rth_tmp[i];

    return 0;
}

bool conversion_coordinates::car2tnh(const std::vector<double> &x_car, const std::vector<double> &s_car, std::vector<double> &x_tnh){

    /*Sanity check*/
    if(x_car.size()!=3){
        smartastro_throw("CAR2TNH: input position vector must be 3-dimensional");
    }
    if(s_car.size()!=6){
        smartastro_throw("CAR2TNH: input position and velocity vector must be 6-dimensional");
    }
    if(x_tnh.size()!=3){
        smartastro_throw("CAR2TNH: output position vector must be 3-dimensional");
    }

    Eigen::Vector3d r,v;           /* Position and velocity vector */
    Eigen::Vector3d h;
    Eigen::Vector3d tn,nn,hn;   /* normalised vectors */
    Eigen::Matrix3d mat;            /* transformation matrix */
    Eigen::Vector3d x_tnh_tmp, x_car_tmp(x_car.data());

    /* Copy the position and velocity vector */
    for(size_t i=0;i<s_car.size(); i++){
        if(i<3){
            r[i]=s_car[i];
        }
        else{
            v[i-3]=s_car[i];
        }
    }

    /* Compute tn */
    double vnorm = v.norm();
    if(vnorm==0)
        return 1;
    tn = v/vnorm;

    /* Compute h and hn */
    h = r.cross(v);
    double hnorm = h.norm();
    if(hnorm==0)
        return 1;
    hn = h/hnorm;

    /* Compute nn */
    nn = hn.cross(tn);

    /* Transformation matrix */
    mat.row(0) = tn;
    mat.row(1) = nn;
    mat.row(2) = hn;

    /* Apply transformation */
    x_tnh_tmp = mat*x_car_tmp;

    for(size_t i=0; i<3; i++)
        x_tnh[i] = x_tnh_tmp[i];

    return 0;
}

bool conversion_coordinates::car2radec(const std::vector<double> &x_car, const std::vector<double> &s_car, std::vector<double> &x_radec){

    /*Sanity check*/
    if(x_car.size()!=3){
        smartastro_throw("CAR2RADEC: input position vector must be 3-dimensional");
    }
    if(s_car.size()!=6){
        smartastro_throw("CAR2RADEC: input position and velocity vector must be 6-dimensional");
    }
    if(x_radec.size()!=3){
        smartastro_throw("CAR2RADEC: output position vector must be 3-dimensional");
    }

    std::vector<double> x_tnh(3);
    int error;
    error = conversion_coordinates::car2tnh(x_car, s_car, x_tnh);

    if (!error)
    {
        error = conversion_coordinates::tnh2radec(x_tnh, x_radec);
    }

    return error;
}

bool conversion_coordinates::kep2car(const std::vector<double> &kep, const double &mu, std::vector<double> &car){

    /*Sanity check*/
    if(kep.size()!=6){
        smartastro_throw("KEP2CAR: keplerian vector must be 6-dimensional");
    }
    if(kep[0]<=0.0){
        smartastro_throw("KEP2CAR: semimajor axis must be positive");
    }
    if((kep[1]>1.0)||(kep[1]<0.0)){
        smartastro_throw("KEP2CAR: eccentricity must be between 0 and 1");
    }
    if((kep[2]>constants::pi)||(kep[2]<0.0)){
        smartastro_throw("KEP2CAR: inclination must be between 0 and PI");
    }
    if( (kep[3] >= constants::pi2) || (kep[3] < 0.0) ){
        smartastro_throw("KEP2CAR: right ascension of the ascending node must be between 0 and 2*PI");
    }
    if((kep[4]>=constants::pi2)||(kep[4]<0.0)){
        smartastro_throw("KEP2CAR: argument of perigee must be between 0 and 2*PI");
    }
    if((kep[5]>=constants::pi2)||(kep[5]<0.0)){
        smartastro_throw("KEP2CAR: true anomaly must be between 0 and 2*PI");
    }
    if(car.size()!=6){
        smartastro_throw("KEP2CAR: cartesian vector must be 6-dimensional");
    }

    /* Declaration */
    double a   = kep[0];    /* In case of parabola, a corresponds to rp */
    double e   = kep[1];
    double i   = kep[2];
    double Om  = kep[3];
    double om  = kep[4];
    double tho = kep[5];

    std::vector<std::vector<double> > rotmat(3,std::vector<double>(3));

    double p, r, wom_dot, r_dot;
    double xp, yp, vxp, vyp;

    /* Rotation matrix */
    rotmat[0][0] = cos(om)*cos(Om)-sin(om)*cos(i)*sin(Om);
    rotmat[1][0] = cos(om)*sin(Om)+sin(om)*cos(i)*cos(Om);
    rotmat[2][0] = sin(om)*sin(i);

    rotmat[0][1] = -sin(om)*cos(Om)-cos(om)*cos(i)*sin(Om);
    rotmat[1][1] = -sin(om)*sin(Om)+cos(om)*cos(i)*cos(Om);
    rotmat[2][1] = cos(om)*sin(i);

    rotmat[0][2] = sin(i)*sin(Om);
    rotmat[1][2] = -sin(i)*cos(Om);
    rotmat[2][2] = cos(i);

    /* In plane Parameters */
    if ((e<(1+constants::tol_ecc)) & (e>(1-constants::tol_ecc))) /* Parabola */
    {
        p = 2.0*a;     /* In the case of a parabola, kep(1) should be rp */
    }
    else
    {
        p = a*(1.0-pow(e,2.0));
    }

    r       = p/(1.0+e*cos(tho));
    xp      = r*cos(tho);
    yp      = r*sin(tho);
    wom_dot = sqrt(mu*p)/pow(r,2.0);
    r_dot   = sqrt(mu/p)*e*sin(tho);
    vxp     = r_dot*cos(tho)-r*sin(tho)*wom_dot;
    vyp     = r_dot*sin(tho)+r*cos(tho)*wom_dot;

    /* 3D cartesian vector */
    car[0] = rotmat[0][0]*xp + rotmat[0][1]*yp;
    car[1] = rotmat[1][0]*xp + rotmat[1][1]*yp;
    car[2] = rotmat[2][0]*xp + rotmat[2][1]*yp;

    car[3] = rotmat[0][0]*vxp + rotmat[0][1]*vyp;
    car[4] = rotmat[1][0]*vxp + rotmat[1][1]*vyp;
    car[5] = rotmat[2][0]*vxp + rotmat[2][1]*vyp;

    return 0;

}

bool conversion_coordinates::rth2car(const std::vector<double> &x_rth, const std::vector<double> &s_car, std::vector<double> &x_car){

    /*Sanity check*/
    if(x_rth.size()!=3){
        smartastro_throw("RTH2CAR: input position vector must be 3-dimensional");
    }
    if(s_car.size()!=6){
        smartastro_throw("RTH2CAR: position and velocity vector must be 6-dimensional");
    }
    if(x_car.size()!=3){
        smartastro_throw("RTH2CAR: output position vector must be 3-dimensional");
    }

    Eigen::Vector3d r,v;           /* Position and velocity vector */
    Eigen::Vector3d h;
    Eigen::Vector3d rn,hn,tn;   /* normalised vectors */
    Eigen::Matrix3d mat;            /* transformation matrix */
    Eigen::Vector3d x_car_tmp, x_rth_tmp(x_rth.data());

    /* Copy the position and velocity vector */
    for(size_t i=0; i<s_car.size(); i++){
        if(i<3){
            r[i] = s_car[i];
        }
        else{
            v[i-3] = s_car[i];
        }
    }

    /* Compute rn */
    double rnorm = r.norm();
    if(rnorm==0)
        return 1;
    rn = r/rnorm;

    /* Compute h and hn */
    h = r.cross(v);
    double hnorm = h.norm();
    if(hnorm==0)
        return 1;
    hn = h/hnorm;

    /* Compute tn */
    tn = hn.cross(rn);

    /* Transformation matrix */
    mat.row(0) = rn;
    mat.row(1) = tn;
    mat.row(2) = hn;

    /* Apply transformation */
    x_car_tmp = x_rth_tmp.transpose()*mat;

    for(size_t i=0; i<3; i++)
        x_car[i] = x_car_tmp[i];

    return 0;
}

bool conversion_coordinates::rth2tnh(const std::vector<double> &x_rth, const double &a, const double &e, const double &f, const double &mu, std::vector<double> &x_tnh){

    /*Sanity check*/
    if(x_rth.size()!=3){
        smartastro_throw("RTH2TNH: input position vector must be 3-dimensional");
    }
    if(a<=0){
        smartastro_throw("RTH2TNH: semimajor axis must be positive");
    }
    if(e>1 || e<0){
        smartastro_throw("RTH2TNH: eccentricity must be between 0 and 1");
    }
    if(f>=constants::pi2 || f<0){
        smartastro_throw("RTH2TNH: true anomaly must be between 0 and 2*PI");
    }
    if(mu<=0){
        smartastro_throw("RTH2TNH: planetary constant must be positive");
    }
    if(x_tnh.size()!=3){
        smartastro_throw("RTH2TNH: output position vector must be 3-dimensional");
    }

    double p, n, h, r, v;
    double sinb, cosb;
    Eigen::Matrix3d rot;
    Eigen::Vector3d x_rth_tmp(x_rth.data()), x_tnh_tmp;

    if ((e>(1-constants::tol_ecc)) & (e<(1+constants::tol_ecc)))
        return 1;

    /* Parameter */
    p = a*(1-pow(e,2.0));
    n = sqrt(mu/pow(a,3.0));
    h = n*pow(a,2.0);
    h = h*sqrt(1-pow(e,2.0));
    r = p/(1+e*cos(f));
    v = sqrt(2*mu/r - mu/a);

    /* Sinus and cosinus */
    sinb = h*e/(p*v)*sin(f);
    cosb = h/(p*v)*(1+e*cos(f));

    /* Rotation matrix  specified by column*/
    rot(0,0) = sinb;
    rot(0,1) = -cosb;
    rot(0,2) = 0;

    rot(1,0) = cosb;
    rot(1,1) = sinb;
    rot(1,2) = 0;

    rot(2,0) = 0;
    rot(2,1) = 0;
    rot(2,2) = 1;

    /* Apply the rotation matrix */
    x_tnh_tmp = x_rth_tmp.transpose()*rot;

    for(size_t i=0; i<3; i++)
        x_tnh[i] = x_tnh_tmp[i];

    return 0;
}

bool conversion_coordinates::tnh2car(const std::vector<double> &x_tnh, const std::vector<double> &s_car, std::vector<double> &x_car){

    /*Sanity check*/
    if(x_tnh.size()!=3){
        smartastro_throw("TNH2CAR: input position vector must be 3-dimensional");
    }
    if(s_car.size()!=6){
        smartastro_throw("TNH2CAR: position and velocity vector must be 6-dimensional");
    }
    if(x_car.size()!=3){
        smartastro_throw("TNH2CAR: output position vector must be 3-dimensional");
    }

    Eigen::Vector3d r,v;           /* Position and velocity vector */
    Eigen::Vector3d h;
    Eigen::Vector3d tn,nn,hn;   /* normalised vectors */
    Eigen::Matrix3d mat;            /* transformation matrix */
    Eigen::Vector3d x_car_tmp, x_tnh_tmp(x_tnh.data());

    /* Copy the position and velocity vector */

    for(size_t i=0; i<s_car.size(); i++){
        if(i<3){
            r[i] = s_car[i];
        }
        else{
            v[i-3] = s_car[i];
        }
    }

    /* Compute tn */
    double vnorm = v.norm();
    if(vnorm == 0)
        return 1;
    tn = v/vnorm;

    /* Compute h and hn */
    h = r.cross(v);
    double hnorm = h.norm();
    if(hnorm == 0)
        return 1;
    hn = h/hnorm;

    /* Compute nn */
    nn = hn.cross(tn);

    /* Transformation matrix */
    mat.row(0) = tn;
    mat.row(1) = nn;
    mat.row(2) = hn;

    /* Apply transformation */
    x_car_tmp = x_tnh_tmp.transpose()*mat;

    for(size_t i=0; i<3; i++)
        x_car[i] = x_car_tmp[i];

    return 0;
}

bool conversion_coordinates::tnh2radec(const std::vector<double> &x_tnh, std::vector<double> &x_radec){

    /*Sanity check*/
    if(x_tnh.size()!=3){
        smartastro_throw("TNH2RADEC: input position vector must be 3-dimensional");
    }
    if(x_radec.size()!=3){
        smartastro_throw("TNH2RADEC: output position vector must be 3-dimensional");
    }

    Eigen::Vector3d x_tnh_tmp(x_tnh.data());
    x_radec[0] = x_tnh_tmp.norm();

    if (x_radec[0] > 0)
    {
        x_radec[1] = atan2(x_tnh[1],x_tnh[0]);
        x_radec[2] = asin(x_tnh[2]/x_radec[0]);
        return 0;
    }
    else
    {
        x_radec[1] = 0;
        x_radec[2] = 0;
        return 1;
    }

    return 0;
}

bool conversion_coordinates::tnh2rth(const std::vector<double> &x_tnh, const double &a, const double &e, const double &f, const double &mu, std::vector<double> &x_rth){

    /*Sanity check*/
    if(x_tnh.size()!=3){
        smartastro_throw("TNH2RTH: input position vector must be 3-dimensional");
    }
    if(a<=0){
        smartastro_throw("TNH2RTH: semimajor axis must be positive");
    }
    if(e>1 || e<0){
        smartastro_throw("TNH2RTH: eccentricity must be between 0 and 1");
    }
    if(f>=constants::pi2 || f<0){
        smartastro_throw("TNH2RTH: true anomaly must be between 0 and 2*PI");
    }
    if(mu<=0){
        smartastro_throw("TNH2RTH: planetary constant must be positive");
    }
    if(x_rth.size()!=3){
        smartastro_throw("TNH2RTH: output position vector must be 3-dimensional");
    }

    double p, n, h, r, v;
    double sinb, cosb;
    Eigen::Matrix3d rot;
    Eigen::Vector3d x_tnh_tmp(x_tnh.data()), x_rth_tmp;

    if ((e>(1-constants::tol_ecc)) & (e<(1+constants::tol_ecc)))
        return 1;

    /* Parameter */
    p = a*(1-pow(e,2.0));
    n = sqrt(mu/pow(a,3.0));
    h = n*pow(a,2.0);
    h = h*sqrt(1-pow(e,2.0));
    r = p/(1+e*cos(f));
    v = sqrt(2*mu/r - mu/a);

    /* Sinus and cosinus */
    sinb = h*e/(p*v)*sin(f);
    cosb = h/(p*v)*(1+e*cos(f));

    /* Rotation matrix  specified by column*/
    rot(0,0) = sinb;
    rot(0,1) = -cosb;
    rot(0,2) = 0;

    rot(1,0) = cosb;
    rot(1,1) = sinb;
    rot(1,2) = 0;

    rot(2,0) = 0;
    rot(2,1) = 0;
    rot(2,2) = 1;

    /* Apply the rotation matrix */
    x_rth_tmp = rot*x_tnh_tmp;

    for(size_t i=0; i<3; i++)
        x_rth[i] = x_rth_tmp[i];

    return 0;
}

bool conversion_coordinates::radec2tnh(const std::vector<double> &x_radec, std::vector<double> &x_tnh){

    /*Sanity check*/
    if(x_radec.size()!=3){
        smartastro_throw("RADEC2TNH: input position vector must be 3-dimensional");
    }
    if(x_tnh.size()!=3){
        smartastro_throw("RADEC2TNH: output position vector must be 3-dimensional");
    }

    x_tnh[0] = x_radec[0] * cos(x_radec[2])*cos(x_radec[1]);
    x_tnh[1] = x_radec[0] * cos(x_radec[2])*sin(x_radec[1]);
    x_tnh[2] = x_radec[0] * sin(x_radec[2]);

    return 0;
}

bool conversion_coordinates::radec2car(const std::vector<double> &x_radec, const std::vector<double> &s_car, std::vector<double> &x_car){

    /*Sanity check*/
    if(x_radec.size()!=3){
        smartastro_throw("RADEC2CAR: input position vector must be 3-dimensional");
    }
    if(s_car.size()!=6){
        smartastro_throw("RADEC2CAR: position and velocity vector must be 6-dimensional");
    }
    if(x_car.size()!=3){
        smartastro_throw("RADEC2CAR: output position vector must be 3-dimensional");
    }

    std::vector<double> x_tnh(3);
    int error;

    conversion_coordinates::radec2tnh(x_radec, x_tnh);
    error = conversion_coordinates::tnh2car(x_tnh, s_car, x_car);

    return error;
}

bool conversion_coordinates::car2spher(const std::vector<double> &car, std::vector<double> &spher){
    /*Sanity check*/
    if(car.size()!=6){
        smartastro_throw("CAR2SPHER: input Cartesian coordinates must be 6-dimensional");
    }
    if(spher.size()!=6){
        smartastro_throw("CAR2SPHER: output spherical coordinates must be 6-dimensional");
    }    

    /* Position */
    spher[0] = sqrt(car[0]*car[0]+car[1]*car[1]+car[2]*car[2]);
    spher[1] = atan2(car[1],car[0]);
    spher[2] = asin(car[2]/spher[0]);

    /* Velocity */
    double phi = spher[2];
    double theta = spher[1];
    std::vector<double> rotMat(9);
    rotMat[0]=cos(phi)*cos(theta);
    rotMat[1]=-sin(theta);
    rotMat[2]=-sin(phi)*cos(theta);
    rotMat[3]=cos(phi)*sin(theta);
    rotMat[4]=cos(theta);
    rotMat[5]=-sin(phi)*sin(theta);
    rotMat[6]=sin(phi);
    rotMat[7]=0.0;
    rotMat[8]=cos(phi);

    spher[3]=rotMat[0]*car[3]+rotMat[3]*car[4]+rotMat[6]*car[5];
    spher[4]=rotMat[1]*car[3]+rotMat[4]*car[4]+rotMat[7]*car[5];
    spher[5]=rotMat[2]*car[3]+rotMat[5]*car[4]+rotMat[8]*car[5];

    return 0;
}

bool conversion_coordinates::kep2spher(const std::vector<double> &kep, const double &mu, std::vector<double> &spher){
    /*Sanity check*/
    if(kep.size()!=6){
        smartastro_throw("KEP2SPHER: input Cartesian coordinates must be 6-dimensional");
    }
    if(spher.size()!=6){
        smartastro_throw("KEP2SPHER: output spherical coordinates must be 6-dimensional");
    }    

    std::vector<double> car = kep;
    kep2car(kep,mu,car); // conversion to Cartesian coordinates
    car2spher(car,spher); // conversion to spherical coordinates

    return 0;
}

bool conversion_coordinates::modeq2car(const std::vector<double> &modeq, const double &mu, std::vector<double> &car){
    /*Sanity check*/
    if(modeq.size()!=6){
        smartastro_throw("MODEQ2KEP: input orbital elements must be 6-dimensional");
    }
    if(car.size()!=6){
        smartastro_throw("MODEQ2KEP: output Keplerian coordinates must be 6-dimensional");
    }    

    double pi = constants::pi;
    double L = modeq[5];
    L -= 2.0 * pi * double(floor(L / (2.0 * pi)));
    double cL = cos(L), sL = sin(L);
    double r = modeq[0] / (1.0 + modeq[2] * sL + modeq[1] * cL);
    double denom = 1 + modeq[3] * modeq[3] + modeq[4] * modeq[4];
    std::vector<double> R(9);
    R[0] = 1.0 - modeq[4] * modeq[4] + modeq[3] * modeq[3];
    R[1] = 2.0 * modeq[3] * modeq[4];
    R[2] = -2.0 * modeq[4];
    R[3] = R[1];
    R[4] = 1.0 + modeq[4] * modeq[4] - modeq[3] * modeq[3];
    R[5] = 2.0 * modeq[3];
    R[6] = -R[2];
    R[7] = -R[5];
    R[8] = 1.0 - modeq[4] * modeq[4] - modeq[3] * modeq[3];
    for(int n = 0; n < 9; n++)
        R[n] /= denom;
    car[0] = (R[0] * cL + R[3] * sL) * r;
    car[1] = (R[1] * cL + R[4] * sL) * r;
    car[2] = (R[2] * cL + R[5] * sL) * r;
    car[3] = (-R[0] * (modeq[2] + sL) + R[3] * (modeq[1] + cL)) * sqrt(mu / modeq[0]);
    car[4] = (-R[1] * (modeq[2] + sL) + R[4] * (modeq[1] + cL)) * sqrt(mu / modeq[0]);
    car[5] = (-R[2] * (modeq[2] + sL) + R[5] * (modeq[1] + cL)) * sqrt(mu / modeq[0]);
    return 0;
}

bool conversion_coordinates::modeq2kep(const std::vector<double> &modeq, std::vector<double> &kep){
    /*Sanity check*/
    if(modeq.size()!=6){
        smartastro_throw("MODEQ2KEP: input orbital elements must be 6-dimensional");
    }
    if(kep.size()!=6){
        smartastro_throw("MODEQ2KEP: output Keplerian coordinates must be 6-dimensional");
    }    

    double pi = constants::pi;
    kep[0] = modeq[0] / (1.0 - (modeq[1] * modeq[1] + modeq[2] * modeq[2]));
    kep[1] = sqrt(modeq[1] * modeq[1] + modeq[2] * modeq[2]);
    kep[2] = 2.0 * atan(sqrt(modeq[3] * modeq[3] + modeq[4] * modeq[4]));      
    kep[3] = atan2(modeq[4], modeq[3]);
    if(kep[3] < 0.0)
        kep[3] += 2.0 * pi;      
    kep[4] = atan2(modeq[2], modeq[1]) - kep[3];
    kep[4] -= 2.0 * pi * double(floor(kep[4] / (2.0 * pi)));
    kep[5] = modeq[5] - atan2(modeq[2], modeq[1]);
    kep[5] -= 2.0 * pi * double(floor(kep[5] / (2.0 * pi)));

    return 0;
}

bool conversion_coordinates::kep2modeq(const std::vector<double> &kep, std::vector<double> &modeq){
    /*Sanity check*/
    if(modeq.size()!=6){
        smartastro_throw("KEP2MODEQ: output orbital elements must be 6-dimensional");
    }
    if(kep.size()!=6){
        smartastro_throw("KEP2MODEQ: input Keplerian coordinates must be 6-dimensional");
    }    

    modeq[0] = kep[0] * (1.0 - kep[1] * kep[1]);
    modeq[1] = kep[1] * cos(kep[3] + kep[4]);
    modeq[2] = kep[1] * sin(kep[3] + kep[4]);
    modeq[3] = tan(kep[2] / 2.0) * cos(kep[3]);
    modeq[4] = tan(kep[2] / 2.0) * sin(kep[3]);
    modeq[5] = kep[5] + kep[4] + kep[3];
    modeq[5] -= 2.0 * constants::pi * double(floor(modeq[5] / (2.0 * constants::pi)));

    return 0;
}


bool conversion_coordinates::kep2delaunay(const std::vector<double> &kep, const double &mu, std::vector<double> &delaunay){
    /*Sanity check*/
    if(delaunay.size()!=6){
        smartastro_throw("KEP2DELAUNAY: output orbital elements must be 6-dimensional");
    }
    if(kep.size()!=6){
        smartastro_throw("KEP2DELAUNAY: input Keplerian coordinates must be 6-dimensional");
    }    

    double M;
    astrocore::conversion_time::true2mean_anomaly(kep[5], kep[1], M);
    delaunay[0] = M; // l
    delaunay[1] = kep[4]; // g
    delaunay[2] = kep[3]; // h
    delaunay[3] = sqrt(mu * kep[0]); // L
    delaunay[4] = delaunay[3] * sqrt(1.0 - kep[1] * kep[1]); // G
    delaunay[5] = delaunay[4] * cos(kep[2]); // H

    return 0;
}

bool conversion_coordinates::delaunay2kep(const std::vector<double> &delaunay, const double &mu, std::vector<double> &kep){
    /*Sanity check*/
    if(delaunay.size()!=6){
        smartastro_throw("DELAUNAY2KEP: input orbital elements must be 6-dimensional");
    }
    if(kep.size()!=6){
        smartastro_throw("DELAUNAY2KEP: output Keplerian coordinates must be 6-dimensional");
    }  

    kep[0] = delaunay[3] * delaunay[3] / mu;
    kep[1] = sqrt(1.0 - (delaunay[4] * delaunay[4]) / (delaunay[3] * delaunay[3]));
    kep[2] = acos(delaunay[5] / delaunay[4]);
    kep[3] = delaunay[2];
    kep[4] = delaunay[1];
    double theta;
    astrocore::conversion_time::mean2true_anomaly(delaunay[0], kep[1], theta);
    kep[5] =  theta;

    return 0;
}

/* This is an adaptation of the code from David Eagle */
bool conversion_coordinates::oscul2mean(const std::vector<double> &oscul, const double &mu, const double & Re, std::vector<double> &mean){
    /*Sanity check*/
    if(oscul.size()!=6)
        smartastro_throw("OSCUL2MEAN: input osculating Keplerian coordinates must be 6-dimensional");
    if(mean.size()!=6)
        smartastro_throw("OSCUL2MEAN: output mean Keplerian coordinates must be 6-dimensional");

    double E0, M0, M;
    // double mu = 1.0;
    // double Re = 1.0;
    double inter, theta;
    double aa, bb, aos, eos, ios, apos, ranos, maos;
    double lamos, zos, etaos, sl, cl, s2l, c2l, s3l, c3l, s4l, c4l, s2i, ci;
    double am, im, ranm, mam, lamm, zm, etam;
    double asp, isp, ransp, lamsp, zsp, etasp, em, apm;
    double pm, tam, um, hm, esp, eqoc, apsp, masp;
    double E_inter;
    double pi = constants::pi;
    double J2 = 0.00108262617385222;

    aos = oscul[0];
    eos = oscul[1];
    ios = oscul[2];
    apos = oscul[4];
    ranos = oscul[3];

    E0 = 2.0 * atan(sqrt((1.0 - oscul[1]) / (1.0 + oscul[1])) * tan(oscul[5] / 2.0));
    if(E0 < 0.0)
        E0 += 2.0 * pi;
    M0 = E0 - oscul[1] * sin(E0);
    maos = M0;

    bb = 0.5 * sin(oscul[2]) * sin(oscul[2]);
    aa = 1.0 / 3.0 - bb;

    if (oscul[1] < 0.01){

        lamos = maos + apos;
        lamos -= 2.0* pi * (double)floor(lamos / (2.0 * pi));

        zos = eos * cos(apos);
        etaos = eos * sin(apos);
        sl = sin(lamos);
        cl = cos(lamos);
        s2l = sin(2.0 * lamos);
        c2l = cos(2.0 * lamos);
        s3l = sin(3.0 * lamos);
        c3l = cos(3.0 * lamos);
        s4l = sin(4.0 * lamos);
        c4l = cos(4.0 * lamos);
        s2i = sin(2.0 * ios);
        ci = cos(ios);

        am = aos;
        im = ios;
        ranm = ranos;
        mam = maos;
        lamm = lamos;
        zm = zos;
        etam = etaos;

        for(int n = 0; n < 10; n++){
                 asp = 3.0 * J2 * Re * Re / am * (bb * c2l + (1.0 - 3.5 * bb) 
                    * zm * cl + (1.0 - 2.5 * bb) * etam * sl + 3.5 * bb 
                    * (zm * c3l + etam * s3l));

                am = aos-asp;

                isp = 3.0 * J2 / 8.0 * Re * Re / (am * am) * s2i * (c2l - zm * cl 
                    + etam * sl + 7.0/3.0 * (zm * c3l + etam * s3l));

                im = ios - isp;

                ci = cos(im);

                s2i = sin(2.0 * im);

                bb = 0.5 * sin(im) * sin(im);

                ransp = 1.5 * J2 * Re * Re / (am * am) * ci 
                    * (0.5 * s2l - 3.5 * zm * sl + 2.5 * etam * cl 
                    + 7.0/6.0 * (zm * s3l - etam * c3l));

                ranm = ranos - ransp;

                lamsp = 1.5 * J2 * Re * Re / (am * am) 
                    * (-0.5 * (1.0 - 5.0 * bb) * s2l + (7.0 - 77.0/4.0 * bb) 
                    * zm * sl - (6.0 - 55.0/4.0 * bb) * etam * cl - (7.0/6.0 - 77.0/12.0 * bb) 
                    * (zm * s3l - etam * c3l));

                lamm = lamos - lamsp;

                sl = sin(lamm);
                cl = cos(lamm);
                s2l = sin(2.0*lamm);
                c2l = cos(2.0*lamm);
                s3l = sin(3.0*lamm);
                c3l = cos(3.0*lamm);
                s4l = sin(4.0*lamm);
                c4l = cos(4.0*lamm);

                zsp = 1.5 * J2 * Re * Re / (am * am) * ((1.0 - 2.5 * bb) * cl 
                    + 7.0/6.0 * bb * c3l + (1.5 - 5.0 * bb) * zm * c2l + (2.0 - 3.0 * bb) 
                    * etam * s2l + 17/4 * bb * (zm * c4l + etam * s4l));

                zm = zos - zsp;

                etasp = 1.5 * J2 * Re * Re / (am * am) * ((1.0 - 3.5 * bb) * sl
                    + 7.0/6.0 * bb * s3l + (1.0 - 6.0 * bb) * zm * s2l - (1.5 - 4.0 * bb) 
                    * etam * c2l + 17.0/4.0 * bb * (zm * s4l - etam * c4l));

                etam = etaos - etasp;   
        }

        em = sqrt(etam * etam + zm * zm);

        apm = 0.0;

        if (em > 1.0e-8)
            apm = atan2(etam, zm);

        mam = lamm - apm; 

    }
    else{

    pm = aos*(1.0-eos*eos);
    am = aos;
    em = eos;
    im = ios;
    apm = apos;
    ranm = ranos;
    mam = maos;

    tam = oscul[5];

    um = apm+tam;
    um -= 2.0 * pi * (double)floor(um / (2.0 * pi));

    hm = pm/(1.0+em*cos(tam));

        for(int n = 0; n < 10; n++){

            inter = sqrt(1.0/((1.0-em*em)*(1.0-em*em)*(1.0-em*em)));

            asp = 3.0*J2*Re*Re/am*((am/hm)*(am/hm)*(am/hm)*(aa+bb*cos(2.0*um)) 
                -aa*inter);

            am = aos-asp;

            isp = 3.0/8.0*J2*(Re/pm)*(Re/pm)*sin(2.0*im)*(cos(2.0*um) 
                +em*cos(tam+2.0*apm)+1.0/3.0*em*cos(3.0*tam+2.0*apm));

            im = ios-isp;

            aa = 1.0/3.0-0.5*sin(im)*sin(im);

            bb = 0.5*sin(im)*sin(im);

            inter = sqrt(1.0/((1.0-em*em)*(1.0-em*em)*(1.0-em*em)));

            esp = 1.5*J2*Re*Re/(am*am)*((1.0-em*em)/em*((am/hm)*(am/hm)*(am/hm)
                *(aa+bb*cos(2.0*um))-aa*inter)-bb/(em*(1.0-em*em))
                *(cos(2.0*um)+em*cos(tam+2.0*apm)+em*cos(3.0*tam+2.0*apm)/3.0));

            em = eos-esp;

            pm = am*(1.0-em*em);

            hm = pm/(1.0+em*cos(tam));

            tam -= 2.0 * pi * (double)floor(tam / (2.0 * pi));

            mam -= 2.0 * pi * (double)floor(mam / (2.0 * pi));

            if ((fabs(tam-pi) <= 1.0e-10) || (fabs(mam-pi) <= 1.0e-10) || (fabs(tam) <= 1.0e-10) || (fabs(mam) <= 1.0e-10))
                eqoc = 0;
            else
                eqoc = tam-mam;

            ransp = -1.5*J2*(Re/pm)*(Re/pm)*cos(im)*(eqoc+em*sin(tam) 
                -0.5*sin(2.0*um)-0.5*em*sin(tam+2.0*apm)-1.0/6.0*em*sin(3.0*tam+2.0*apm));

            ranm = ranos-ransp;

            apsp = 1.5*J2*(Re/pm)*(Re/pm)*((2.0-5.0*bb)*(eqoc+em*sin(tam)) 
                +(1.0-3.0*bb)*((1.0-0.25*em*em)*sin(tam)/em+0.5*sin(2.0*tam) 
                +em*sin(3.0*tam)/12.0)-(0.5*bb+(0.5-15.0/8.0*bb)*em*em)/em*sin(tam+2.0*apm) 
                +em/8.0*bb*sin(tam-2.0*apm)-0.5*(1-5*bb)*sin(2.0*um)+(7.0/6.0*bb 
                -1.0/6.0*em*em*(1.0-19.0/4.0*bb))/em*sin(3.0*tam+2.0*apm) 
                +0.75*bb*sin(4.0*tam+2.0*apm)+em/8*bb*sin(5.0*tam+2.0*apm));

            apm = apos-apsp;

            masp = 1.5*J2*(Re/pm)*(Re/pm)*sqrt(1.0-em*em)/em*(-(1.0-3.0*bb)*((1.0-em*em/4.0) 
                *sin(tam)+em/2.0*sin(2.0*tam)+em*em/12.0*sin(3.0*tam)) 
                +bb*(0.5*(1.0+1.25*em*em)*sin(tam+2*apm)-em*em/8.0*sin(tam-2.0*apm) 
                -7.0/6.0*(1.0-em*em/28.0)*sin(3*tam+2*apm)-0.75*em*sin(4.0*tam+2.0*apm) 
                -em*em/8.0*sin(5.0*tam+2.0*apm)));

            mam = maos-masp;

            astrocore::conversion_time::mean2eccentric_anomaly(mam, em, E_inter);
            tam = 2.0 * atan(sqrt((1.0 + em) / (1.0 - em)) * tan(E_inter / 2.0));
            um = apm + tam;
            um -= 2.0* pi * (double)floor(um / (2.0 * pi));            

        }   

    }

    M = mam;
    M -= 2.0 * pi * (double)floor(M / (2.0 * pi));

    astrocore::conversion_time::mean2eccentric_anomaly(M, em, E_inter);
    theta = 2.0 * atan(sqrt((1.0 + em) / (1.0 - em)) * tan(E_inter / 2.0));
    if(theta < 0.0)
        theta += 2.0 * pi;

    mean[0] = am;
    mean[1] = em;
    mean[2] = im;
    mean[3] = ranm - 2.0* pi * (double)floor(ranm / (2.0 * pi));
    mean[4] = apm - 2.0* pi * (double)floor(apm / (2.0 * pi));
    mean[5] = theta;

    return 0;
}

bool conversion_coordinates::mean2oscul(const std::vector<double> &mean, const double &mu, const double & Re, std::vector<double> &oscul){
    /*Sanity check*/
    if(oscul.size()!=6)
        smartastro_throw("MEAN2OSCUL: output osculating Keplerian coordinates must be 6-dimensional");
    if(mean.size()!=6)
        smartastro_throw("MEAN2OSCUL: input mean Keplerian coordinates must be 6-dimensional");

    double E0, M0;
    double si, ci, ti, si2, ci2, f, sf, cf, s2f, u, e, e2, esf, am;
    double d1, d2, d3, d4, d42, d5, d6, d7, d8, p, r, rdot;
    double twou, twow, s2u, c2u, sf2w, s3f2w, cf2w, c3f2w, q1, di, dp;
    double dummy1, dn, dr, drdot, du, pnw, ainw, annw, rnw, rdotnw, unw, aa, bb, enw2, enw, xfnw, anw, wnw;
    double pi = constants::pi;
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

    xfnw = atan2(bb, aa); //atan3???
    if(xfnw < 0.0)
        xfnw += 2.0 * pi;

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


bool conversion_coordinates::euler_axis_angle(const std::vector<double> &v, const std::vector<double> &n, const double &theta, std::vector<double> &v1)
{

    /*Sanity check*/
    if(v.size()!=3)
        smartastro_throw("euler_axis_angle: vector to be rotated must be 3-dimensional");
    if(n.size()!=3)
        smartastro_throw("euler_axis_angle: axis vector must be 3-dimensional");

    std::vector<double> nn=n;       /* normalised vector of n */
    double norm = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    for(size_t i=0; i<3; i++)
        nn[i]/=norm;
    
    std::vector<double> rot(9);
    /* Define the rotation matrix by column */
    rot[0] = cos(theta)+(1.0-cos(theta))*nn[0]*nn[0];
    rot[1] = (1.0-cos(theta))*nn[0]*nn[1]-sin(theta)*nn[2];
    rot[2] = (1.0-cos(theta))*nn[0]*nn[2]+sin(theta)*nn[1];
    
    rot[3] = (1.0-cos(theta))*nn[0]*nn[1]+sin(theta)*nn[2];
    rot[4] = cos(theta)+(1.0-cos(theta))*nn[1]*nn[1];
    rot[5] = (1.0-cos(theta))*nn[1]*nn[2]-sin(theta)*nn[0];
    
    rot[6] = (1.0-cos(theta))*nn[0]*nn[2]-sin(theta)*nn[1];
    rot[7] = (1.0-cos(theta))*nn[1]*nn[2]+sin(theta)*nn[0];
    rot[8] = cos(theta)+(1.0-cos(theta))*nn[2]*nn[2];

    v1=v;
    v1[0]=rot[0]*v[0]+rot[3]*v[1]+rot[6]*v[2];
    v1[1]=rot[1]*v[0]+rot[4]*v[1]+rot[7]*v[2];
    v1[2]=rot[2]*v[0]+rot[5]*v[1]+rot[8]*v[2];
      
    return 0;
}




/**
 * @brief Function converting equinoctial elements into Keplerian coordinates
 *
 * Function converting modified equinoctial elements into Keplerian coordinates (scaling must be consistent)
 *
 * @param[in] vector of equinoctial elements [a,p1,p2,q1,q2,L]
 * @param[out] vector of Keplerian coordinates [a,e,i,RAAN,w,true anomaly]
 *
 * @return Error code
 *
 * @author Cristian Greco 2018
 */
bool conversion_coordinates::eq2kep(const std::vector<double> &eq, std::vector<double> &kep) {

    /*Sanity check*/
    if(eq.size()!=6){
        smartastro_throw("Input orbital elements must be 6-dimensional");
    }
    if(kep.size()!=6){
        smartastro_throw("Output Keplerian coordinates must be 6-dimensional");
    }

    double pi = constants::pi;

    double a  = eq[0];
    double P1 = eq[1];
    double P2 = eq[2];
    double Q1 = eq[3];
    double Q2 = eq[4];
    double L  = eq[5];

    double e,i,RAAN,w,teta;

    e = std::sqrt(P1*P1 + P2*P2);
    i = 2.0 * std::atan(sqrt((Q1*Q1 + Q2*Q2)));
    while (i<0.0)
    {
        i += 2.0*pi;
    }
    i = std::fmod(i,2.0*pi);
    if (std::abs(i) <= 1.0e-6)
    {
        RAAN = 0.0;
    }
    else
    {
        double sin_RAAN = Q1 / sqrt(Q1*Q1 + Q2*Q2);
        double cos_RAAN = Q2 / sqrt(Q1*Q1 + Q2*Q2);
        RAAN  = std::atan2(sin_RAAN, cos_RAAN);
        while (RAAN<0.0)
        {
            RAAN += 2.0*pi;
        }
        RAAN = std::fmod(RAAN,2.0*pi);
    }

    double sin_zeta = P1 / std::sqrt(P1*P1 + P2*P2);
    double cos_zeta = P2 / std::sqrt(P1*P1 + P2*P2);
    double zeta = std::fmod(std::atan2(sin_zeta, cos_zeta),2.0*pi);
    if (std::isnan(zeta))
    {
        zeta = 0.0;
    }

    w = zeta-RAAN;
    while (w<0.0)
    {
        w += 2.0*pi;
    }
    w = std::fmod(w,2.0*pi);

    teta = L-RAAN-w;
    while (teta<0.0)
    {
        teta += 2.0*pi;
    }
    teta = std::fmod(teta,2.0*pi);

    kep[0] = a;
    kep[1] = e;
    kep[2] = i;
    kep[3] = RAAN;
    kep[4] = w;
    kep[5] = teta;

    return true;

}


/**
 * @brief Function converting Keplerian coordinates into equinoctial elements
 *
 * Function converting Keplerian coordinates into equinoctial elements (scaling must be consistent)
 *
 * @param[in] vector of Keplerian coordinates [a,e,i,RAAN,w,true anomaly]
 * @param[out] vector of equinoctial elements [a,p1,p2,q1,q2,L]
 *
 * @return Error code
 *
 * @author Cristian Greco 2018
 */
bool conversion_coordinates::kep2eq(const std::vector<double> &kep, std::vector<double> &eq) {

    /*Sanity check*/
    if(eq.size()!=6){
        smartastro_throw("Input orbital elements must be 6-dimensional");
    }
    if(kep.size()!=6){
        smartastro_throw("Output Keplerian coordinates must be 6-dimensional");
    }

    double pi = constants::pi;

    // List Keplerian elements
    double a    = kep[0];
    double e    = kep[1];
    double i    = kep[2];
    double RAAN = kep[3];
    double w    = kep[4];
    double teta = kep[5];

    // Equinoctial elements
    double P1 = e*std::sin(w+RAAN);
    double P2 = e*std::cos(w+RAAN);
    double Q1 = std::tan(i/2.0)*std::sin(RAAN);
    double Q2 = std::tan(i/2.0)*std::cos(RAAN);
    double L  = teta+w+RAAN;
    while (L<0.0)
    {
        L += 2.0*pi;
    }
    L = std::fmod(L,2.0*pi);

    // Set output vector
    eq[0] = a;
    eq[1] = P1;
    eq[2] = P2;
    eq[3] = Q1;
    eq[4] = Q2;
    eq[5] = L;

    return true;

}



/**
 * @brief Function converting equinoctial elements into Cartesian coordinates
 *
 * Function converting equinoctial elements into Cartesian coordinates (scaling must be consistent)
 *
 * @param[in] vector of equinoctial elements         [a,p1,p2,q1,q2,L]
 * @param[in] gravitational constant of central body [mu]
 * @param[out] vector of Cartesian coordinates       [a,e,i,RAAN,w,true anomaly]
 *
 * @return Error code
 *
 * @author Cristian Greco 2018
 */
bool conversion_coordinates::eq2car(const std::vector<double> &eq, const double &mu, std::vector<double> &car)
{
    // Middle passage through keplerian
    std::vector<double> kep(6);
    eq2kep(eq,kep);

    // Convert kepler to cartesian
    kep2car(kep,mu,car);

    return true;
}

/**
 * @brief Function converting Cartesian coordinates into equinoctial elements
 *
 * Function converting Cartesian coordinates into equinoctial elements (scaling must be consistent)
 *
 * @param[out] vector of Cartesian coordinates       [a,e,i,RAAN,w,true anomaly]
 * @param[in] gravitational constant of central body [mu]
 * @param[in] vector of equinoctial elements         [a,p1,p2,q1,q2,L]
 *
 * @return Error code
 *
 * @author Cristian Greco 2018
 */
bool conversion_coordinates::car2eq(const std::vector<double> &car, const double &mu, std::vector<double> &eq)
{
    // Middle passage through Keplerian
    std::vector<double> kep(6);
    car2kep(car,mu,kep);

    // Convert kepler to equinoctial
    kep2eq(kep,eq);

    return true;
}
