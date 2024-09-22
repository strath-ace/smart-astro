"""
    celestialbodiesephemeris_kep(mjd2000,ibody)
    Analytical ephemerides for celestial bodies.
    Planetay orbital elements are restituited in a Sun-centred ecliptic 
    system. These ephemerides were succesfully compared with JPL/NAIF/SPICE
    ephemerides using de405.bps. Lunar orbital elements are restituited
    in an Earth-centered equatorial reference frame.

    INPUT :
    mjd2000[1]  Time, modified Julian day since 01/01/2000, 12:00 noon
                (MJD2000 = MJD-51544.5)
    ibody[1]    Integer number identifying the celestial body (< 11)
                    1:   Mercury
                    2:   Venus
                    3:   Earth
                    4:   Mars
                    5:   Jupiter
                    6:   Saturn
                    7:   Uranus
                    8:   Neptune
                    9:   Pluto
                    10:  Sun
                    11:  Moon

    OUTPUT:
    kep[6]    	Mean Keplerian elements of date
                    kep = [a e i Om om nu] [km, rad] (nu=true anomaly)

    MATLAB VERSION (uplanet.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AUTHOR:
    P. Dysli, 1977 

    PREVIOUS VERSION:
    P. Dysli, 1977 
        - Header and function name in accordance with guidlines.

    CHANGELOG:
    28/12/06, Camilla Colombo: tidied up
    10/01/2007, REVISION, Matteo Ceriotti 
    03/05/2008, Camilla Colombo: Case 11 deleted.
    11/09/2008, Matteo Ceriotti, Camilla Colombo:
        - All ephemerides shifted 0.5 days back in time. Now mjd2000 used
            in this function is referred to 01/01/2000 12:00. In the old
            version it was referred to 02/01/2000 00:00.
        - Corrected ephemeris of Pluto.
    04/10/2010, Camilla Colombo: Header and function name in accordance
        with guidlines.
    
    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language
    
    
"""
function celestialbodiesephemeris_kep(mjd2000,ibody)

    if ibody>11
        Base.error("no celestial body in list")
    elseif ibody==11
        car = eph_moon(mjd2000);
        kep = car2kep(car,398600.4418); 
    else
        DEG2RAD =pi/180;
        KM = 149597870.66;
        
        #T is julian centuries since 31/12/1899 at 12:00
        T = (mjd2000 + 36525)/36525.00;
        TT  = T .* T;
        TTT = T .* TT;
    
        kep = zeros(6);
        if ibody == 1 #Mercury

            kep[1] = 0.38709860;
            kep[2] = 0.205614210 + 0.000020460*T - 0.000000030*TT;
            kep[3] = 7.002880555555555560 + 1.86083333333333333e-3*T - 1.83333333333333333e-5*TT;
            kep[5] = 4.71459444444444444e+1 + 1.185208333333333330*T + 1.73888888888888889e-4*TT;
            kep[4] = 2.87537527777777778e+1 + 3.70280555555555556e-1*T +1.20833333333333333e-4*TT;
            XM   = 1.49472515288888889e+5 + 6.38888888888888889e-6*T;
            kep[6] = 1.02279380555555556e2 + XM.*T;

        elseif ibody == 2 #Venus

            kep[1] = 0.72333160;
            kep[2] = 0.006820690 - 0.000047740*T + 0.0000000910*TT;
            kep[3] = 3.393630555555555560 + 1.00583333333333333e-3*T - 9.72222222222222222e-7*TT;
            kep[5] = 7.57796472222222222e+1 + 8.9985e-1*T + 4.1e-4*TT;
            kep[4] = 5.43841861111111111e+1 + 5.08186111111111111e-1*T -1.38638888888888889e-3*TT;
            XM   = 5.8517803875e+4 + 1.28605555555555556e-3*T;
            kep[6] = 2.12603219444444444e2 + XM.*T;

        elseif ibody == 3 # Earth
    
            kep[1] = 1.000000230;
            kep[2] = 0.016751040 - 0.000041800*T - 0.0000001260*TT;
            kep[3] = 0;
            kep[5] = 0;
            kep[4] = 1.01220833333333333e+2 + 1.7191750*T + 4.52777777777777778e-4*TT + 3.33333333333333333e-6*TTT;
            XM       = 3.599904975e+4 - 1.50277777777777778e-4*T - 3.33333333333333333e-6*TT;
            kep[6] = 3.58475844444444444e2 + XM .* T;
    
        elseif ibody == 4 #Mars

            kep[1] = 1.5236883990;
            kep[2] = 0.093312900 + 0.0000920640*T - 0.0000000770*TT;
            kep[3] = 1.850333333333333330 - 6.75e-4*T + 1.26111111111111111e-5*TT;
            kep[5] = 4.87864416666666667e+1 + 7.70991666666666667e-1*T - 1.38888888888888889e-6*TT - 5.33333333333333333e-6*TTT;
            kep[4] = 2.85431761111111111e+2 + 1.069766666666666670*T +  1.3125e-4*TT + 4.13888888888888889e-6*TTT;
            XM   = 1.91398585e+4 + 1.80805555555555556e-4*T + 1.19444444444444444e-6*TT;
            kep[6] = 3.19529425e2 + XM.*T;

        elseif ibody == 5 #Jupiter

            kep[1] = 5.2025610;
            kep[2] = 0.048334750 + 0.000164180*T  - 0.00000046760*TT -0.00000000170*TTT;
            kep[3] = 1.308736111111111110 - 5.69611111111111111e-3*T +  3.88888888888888889e-6*TT;
            kep[5] = 9.94433861111111111e+1 + 1.010530*T + 3.52222222222222222e-4*TT - 8.51111111111111111e-6*TTT;
            kep[4] = 2.73277541666666667e+2 + 5.99431666666666667e-1*T + 7.0405e-4*TT + 5.07777777777777778e-6*TTT;
            XM   = 3.03469202388888889e+3 - 7.21588888888888889e-4*T + 1.78444444444444444e-6*TT;
            kep[6] = 2.25328327777777778e2 + XM.*T;

        elseif ibody == 6 #Saturn

            kep[1] = 9.5547470;
            kep[2] = 0.055892320 - 0.00034550*T - 0.0000007280*TT + 0.000000000740*TTT;
            kep[3] = 2.492519444444444440 - 3.91888888888888889e-3*T - 1.54888888888888889e-5*TT + 4.44444444444444444e-8*TTT;
            kep[5] = 1.12790388888888889e+2 + 8.73195138888888889e-1*T -1.52180555555555556e-4*TT - 5.30555555555555556e-6*TTT;
            kep[4] = 3.38307772222222222e+2 + 1.085220694444444440*T + 9.78541666666666667e-4*TT + 9.91666666666666667e-6*TTT;
            XM   = 1.22155146777777778e+3 - 5.01819444444444444e-4*T - 5.19444444444444444e-6*TT;
            kep[6] = 1.75466216666666667e2 + XM.*T;

        elseif ibody == 7 #Uranus

            kep[1] = 19.218140;
            kep[2] = 0.04634440 - 0.000026580*T + 0.0000000770*TT;
            kep[3] = 7.72463888888888889e-1 + 6.25277777777777778e-4*T + 3.95e-5*TT;
            kep[5] = 7.34770972222222222e+1 + 4.98667777777777778e-1*T + 1.31166666666666667e-3*TT;
            kep[4] = 9.80715527777777778e+1 + 9.85765e-1*T - 1.07447222222222222e-3*TT - 6.05555555555555556e-7*TTT;
            XM   = 4.28379113055555556e+2 + 7.88444444444444444e-5*T + 1.11111111111111111e-9*TT;
            kep[6] = 7.26488194444444444e1 + XM.*T;

        elseif ibody == 8 #Neptune

            kep[1] = 30.109570;
            kep[2] = 0.008997040 + 0.0000063300*T - 0.0000000020*TT;
            kep[3] = 1.779241666666666670 - 9.54361111111111111e-3*T - 9.11111111111111111e-6*TT;
            kep[5] = 1.30681358333333333e+2 + 1.0989350*T + 2.49866666666666667e-4*TT - 4.71777777777777778e-6*TTT;
            kep[4] = 2.76045966666666667e+2 + 3.25639444444444444e-1*T + 1.4095e-4*TT + 4.11333333333333333e-6*TTT;
            XM   = 2.18461339722222222e+2 - 7.03333333333333333e-5*T;
            kep[6] = 3.77306694444444444e1 + XM.*T;

        elseif ibody == 9

            kep[1] = 39.481686778174627;
            kep[2] = 2.4467e-001;
            kep[3] = 17.150918639446061;
            kep[5] = 110.27718682882954;
            kep[4] = 113.77222937912757;
            XM   = 4.5982945101558835e-008;
            kep[6] = 1.5021e+001 + XM.*mjd2000*86400;

        elseif ibody == 10

            kep = [0,0,0,0,0,0]

        end

        #now, convert these elements from AU -> km and deg -> rad.
        
        kep[1]   = kep[1]*KM;       # a [km]
        kep[3:6] = kep[3:6]*DEG2RAD;    # Transform from deg to rad
        kep[6]   = mod(kep[6],2*pi);
        phi = kep[6]
        for i in 1:5
            g       = kep[6]-(phi-kep[2] * sin(phi)); 
            g_primo = (-1+kep[2] * cos(phi));
            phi     = phi-g./g_primo;   # Computes the eccentric anomaly kep    
        end
        nu=2*atan(sqrt((1+kep[2])./(1-kep[2])) .* tan(phi/2));
    
        kep[6]= nu;
    end
  
    return kep
     
end

"""
    celestialbodiesephemeris_position(ibody,mjd2000)
    Planetay position vectors are restituited in a Sun-centred ecliptic 
    system. Lunar position vector is restituited
    in an Earth-centered equatorial reference frame.
    
    INPUT :
    ibody[1]    Integer number identifying the celestial body (< 11)
                    1:   Mercury
                    2:   Venus
                    3:   Earth
                    4:   Mars
                    5:   Jupiter
                    6:   Saturn
                    7:   Uranus
                    8:   Neptune
                    9:   Pluto
                    10:  Sun
                    11:  Moon
    mjd2000[1]  Time, modified Julian day since 01/01/2000, 12:00 noon
    (MJD2000 = MJD-51544.5)

    OUTPUT:
    rV[3]    	Position vector [km]

    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ORIGINAL VERSION:
       Massimilaino Vasile, 2002, MATLAB, EphSSfun.m
    
    AUTHOR:  
       Matteo Ceriotti, 10/01/2007, MATLAB, EphSS_car.m
    
    PREVIOUS VERSION:
       Matteo Ceriotti, 10/01/2007, MATLAB, EphSS.m
           - Header and function name in accordance with guidlines.
    
    CHANGELOG:
       12/02/2007, Matteo Ceriotti
       03/05/2008, Nicolas Croisard: ephMoon instead of uplanet for the Moon
       03/05/2008, REVISION, Camilla Colombo 
       04/10/2010, Camilla Colombo: Header and function name in accordance
           with guidlines (also changed name of called functions: ephMoon,
           ephNEO, kep2car, astroConstants)


    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language
    10/2023 : Irene Cavallari   : only position vectors as an output (not velocity vector)

"""
function celestialbodiesephemeris_position(ibody,mjd2000)

    if ibody<11 #uplanet needed
        kep = celestialbodiesephemeris_kep(mjd2000,ibody)
    elseif ibody == 11 #ephMoon needed
        rV = lunar_position_equatorialrefframe(mjd2000);
    else #NeoEphemeris needed
        Base.error("no celestial body in the list")
    end

    if ibody != 11 #planet or asteroid, sun centered
        car = kep2car(kep,0.19891000000000E+31*6.67259e-20);
        rV  = car[4:6]
    end

    return rV;
end

"""
    lunar_position_equatorialrefframe(mjd2000)
    Cartesian position of the Moon.
    It gives the position of the Moon at a given epoch in
    Geocentric Equatorial Reference Frame (IAU-76/FK5 J2000, mean equator,
    mean equinox frame) This frame {x,y,z} is characterised by:
        x-axis: on the equatorial plane, along the direction of the gamma
            point
        z-axis: direction of the north pole
        y-axis: on the equatorial plane, completes the reference frame

    INPUT:
    MJD2000[1]	Epoch in Modified Julian Date 2000 (MJD2000 since 12:00
                noon 01/01/2000)

    OUTPUT:
    xP           Position vector of the Moon in cartesian coordinates,
                expressed in the Geocentric Equatorial Reference Frame [km].
 
    
    REFERENCES:
    Algorithm taken from "Fundamentals of Astrodynamics and Applications"
    (3rd edition), D. A. Vallado, p.290 (algorithm 31).

    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AUTHOR:
    Daniel Novak, 04/03/2008, MATLAB, ephMoon.m

    PREVIOUS VERSION:
    Daniel Novak, 04/03/2008, MATLAB, eph_moon.m
        - Header and function name in accordance with guidlines.

    CHANGELOG:
    05/03/2008, REVISION, Matteo Ceriotti
    06/05/2008, Nicolas Croisard: Use of COS and SIN instead of COSD and
        SIND for speed (15 times faster)
    04/10/2010, Camilla Colombo: Header and function name in accordance
        with guidlines.
    16/03/2020, Matteo Manzi: this_style function name

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language
    10/2023 : Irene Cavallari   : only position vectors as an output (not velocity vector)
"""
function lunar_position_equatorialrefframe(mjd2000)

    T_TDB = mjd2000/36525;

    angles = [134.9 + 477198.85*T_TDB;     # L_ecl 1 and p 1
          259.2 - 413335.38*T_TDB;     # L_ecl 2 and p 2
          235.7 + 890534.23*T_TDB;     # L_ecl 3 and p 3
          269.9 +  954397.7*T_TDB;     # L_ecl 4 and p 4
          357.5 +  35999.05*T_TDB;     # L_ecl 5
          186.6 + 966404.05*T_TDB;     # L_ecl 6
           93.3 + 483202.03*T_TDB;     # phi_ecl 1
          228.2 + 960400.87*T_TDB;     # phi_ecl 2
          318.3 +   6003.18*T_TDB;     # phi_ecl 3
          217.6 -  407332.2*T_TDB;     # phi_ecl 4
          ];

    angles = angles * pi/180

    s = sin.(angles[1:10])
    c = cos.(angles[1:4])

    L_ecl = ((218.32 + 481267.883*T_TDB) .+ [6.29 -1.27 0.66 0.21 -0.19 -0.11] * s[1:6])[1];

    phi_ecl = ([5.13 0.28 -0.28 -0.17] * s[7:10])[1]

    P = 0.9508 .+ [0.0518 0.0095 0.0078 0.0028] * c;

    eps = 23.439291 - 0.0130042*T_TDB - 1.64e-7*T_TDB^2 + 5.04e-7*T_TDB^3;

    L_ecl   = L_ecl[1]   * pi/180;
    phi_ecl = phi_ecl[1] * pi/180;
    P       = P[1]       * pi/180;
    eps     = eps     * pi/180;

    r = 1/sin(P) * 6378.16

    xP = r * [cos(L_ecl)*cos(phi_ecl), 
    cos(eps)*cos(phi_ecl)*sin(L_ecl) - sin(eps)*sin(phi_ecl), 
    sin(eps)*cos(phi_ecl)*sin(L_ecl) + cos(eps)*sin(phi_ecl)];

    return xP;

end

"""
    eph_moon(mjd2000)
    Ephemerides (cartesian position and velocity) of the Moon.
    It gives the position and the velocity of the Moon at a given epoch in
    Geocentric Equatorial Reference Frame (IAU-76/FK5 J2000, mean equator,
    mean equinox frame) This frame {x,y,z} is characterised by:
        x-axis: on the equatorial plane, along the direction of the gamma
            point
        z-axis: direction of the north pole
        y-axis: on the equatorial plane, completes the reference frame

    INPUT:
    MJD2000[1]	Epoch in Modified Julian Date 2000 (MJD2000 since 12:00
                noon 01/01/2000)

    OUTPUT:
    xP           Position vector of the Moon in cartesian coordinates,
                expressed in the Geocentric Equatorial Reference Frame [km].
    vP           Velocity vector of the Moon in cartesian coordinates,
                expressed in the Geocentric Equatorial Reference Frame
                [km/s]. The velocity is computed by numerical
                differentiation on a 1 second interval.
 
    
    REFERENCES:
    Algorithm taken from "Fundamentals of Astrodynamics and Applications"
    (3rd edition), D. A. Vallado, p.290 (algorithm 31).

    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AUTHOR:
    Daniel Novak, 04/03/2008, MATLAB, ephMoon.m

    PREVIOUS VERSION:
    Daniel Novak, 04/03/2008, MATLAB, eph_moon.m
        - Header and function name in accordance with guidlines.

    CHANGELOG:
    05/03/2008, REVISION, Matteo Ceriotti
    06/05/2008, Nicolas Croisard: Use of COS and SIN instead of COSD and
        SIND for speed (15 times faster)
    04/10/2010, Camilla Colombo: Header and function name in accordance
        with guidlines.
    16/03/2020, Matteo Manzi: this_style function name

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language
"""
function eph_moon(mjd2000)
    #Returns pos and vel of moon in geocentric equatorial reference frame
    
    xP = lunar_position_equatorialrefframe(mjd2000);
    vP = lunar_position_equatorialrefframe(mjd2000+1/86400) - xP; 

    return append!(vP,xP);

end

