#### interfaces
"""
    getWSRP2(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion 
    (expressed in Sadov variables) with respect to Sadov fast angles
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    m...WTS = terms depending on m pre-computed with the functions in termsWtbWtsForSrpDragW2 in this file
    EL = eccentric longitude (E+om+OM, E=eccentric anomaly, om = argumentum of pericentre, OM = right ascension of the ascending node)
    rSun = position vector of the Sun wrt Earth
    satellite = Dict with the characteristics of the satellite ***
    includeEclipsesEffects = boolean : if true the eclipses effects are considered, if false the body is considered in sun-light (even if it is not)
    passageinshadowoccurs = boolean: if true part of the orbit is in shadow; if false the orbit is in sun-light (if includeEclipsesEffects, set passageinshadowoccurs = false)
    ELin  = value of the eccentric longitude at the entrance of the shadow region (if passageinshadowoccurs, set  Elinshadow = NaN)
    ELout = value of the eccentric longitude at the exit from the shadow region (if passageinshadowoccurs, set  Eloutshadow = NaN)

    **
    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    The reference frame is such that m belongs to in[0,1] and k>0 where
    k = C/A*(B-A)/(C-B) 
    m = k*(C-delta)/(delta-A)
    with
    A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
    delta = G^2/2/Phi, with G the angular momentum of the body and Phi the kinetic energy
    Let GV=diag([A,B,C])omV be angular momentum of the body and G = norm(GV)
    (omV = [p,q,r]^T is the angular velocity) 
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular 
    to GV and the xy plane. 
    sadovvar = [Jl,Jg,Jh,psil,psig,psih] or [Jg,Jh,psil,psig,psih]
            with
                Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m)) [kg m^2 /s^2]
                Jg   = G [kg m^2 /s^2]
                Jh   = Gcos(delta), delta = inclination angle between [kg m^2 /s^2]
                       the XY plane and the plane perpendicular to GV.
                psil = pi/2 ellK(m)ellF(lambda,m) [rad], lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
                psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)) [rad], g = angle between N and N''
                psih = angle between X axis and N [rad]


    ***
    satellite=Dict("MomentsOfInertia"=>IV,"numberOfFacets"=>n,"facets"=>facets);
    IV = [A,B,C] : principal moment of inertia
    IM = intrinsic magnetic moment of the satellite [A m^2]
    n  = number of facets in which the satellite surface is divided
    facets = [facet1...facetn]
        where  
            facetk = [facet_coeff,facet_area,facet_Vinfo,facet_nv]
            facet_coeff  = [cai*psrp,cdi*psrsp,csi*psrp], psrp = solar radiation pressure
            facet_area   = area of the facet [m^2]
            facet_Vinfo  = [rhoiv,vv1,vv2,vv3], rhoiv = centre of mass to facet centroid vector [m], vvi = verteces of the facet [m], i=1..3
            facet_nv     = normal unit vector  
            facet_VVT1   = vector containing combinations of cai*psrp with the components of  normal unit vector and components of rhoiv and area (output of srp_average-->srp_initforaverage_T1(cai,surfi,rhoiv,niv))
            facet_VVT2   = vector containing combinations of cdi*psrp with the components of  normal unit vector and components of rhoiv and area (output of srp_average-->srp_initforaverage_T2(cai,surfi,rhoiv,niv))
            facet_VVT3   = vector containing combinations of csi*psrp with the components of  normal unit vector and components of rhoiv and area (output of srp_average-->srp_initforaverage_T3(cai,surfi,rhoiv,niv))
            with 
                cai = (1-rf*sp)
                cdi = (2/3*(1-sp)*rf+2/3*(1-rf))
                csi = 2*rf*sp*psrp
                with rf reflectivity and sp specular coefficient (see Benson&Sheers2021)                       
        
    OUTPUT
    W2 = vector 6x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,psig,psih] -> [skm',Jg',Jh',psil',psig',psih'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWSRP2(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    W2 = getWSRP2procedure(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,1);
    # W2 = getWSRP2procedure_numerical(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,1);
    return W2;
end

"""
    getWSRP2sadovlike(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion 
    (expressed in Sadov-like variables) with respect to Sadov fast angles
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    EL = eccentric longitude (E+om+OM, E=eccentric anomaly, om = argumentum of pericentre, OM = right ascension of the ascending node)
    m...WTS = terms depending on m pre-computed with the functions in termsWtbWtsForSrpDragW2 in this file
    rSun = position vector of the Sun wrt Earth
    satellite = Dict with the characteristics of the satellite ***
    includeEclipsesEffects = boolean : if true the eclipses effects are considered, if false the body is considered in sun-light (even if it is not)
    passageinshadowoccurs = boolean: if true part of the orbit is in shadow; if false the orbit is in sun-light (if includeEclipsesEffects, set passageinshadowoccurs = false)
    ELin  = value of the eccentric longitude at the entrance of the shadow region (if passageinshadowoccurs, set  Elinshadow = NaN)
    ELout = value of the eccentric longitude at the exit from the shadow region (if passageinshadowoccurs, set  Eloutshadow = NaN)


    **
    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    The reference frame is such that m belongs to in[0,1] and k>0 where
    k = C/A*(B-A)/(C-B) 
    m = k*(C-delta)/(delta-A)
    with
    A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
    delta = G^2/2/Phi, with G the angular momentum of the body and Phi the kinetic energy
    Let GV=diag([A,B,C])omV be angular momentum of the body and G = norm(GV)
    (omV = [p,q,r]^T is the angular velocity) 
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular 
    to GV and the xy plane. 
    sadovvar = [Jl,Jg,Jh,psil,psig,psih] or [Jg,Jh,psil,psig,psih]
            with
                Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m)) [kg m^2 /s^2]
                Jg   = G [kg m^2 /s^2]
                Jh   = Gcos(delta), delta = inclination angle between [kg m^2 /s^2]
                       the XY plane and the plane perpendicular to GV.
                psil = pi/2 ellK(m)ellF(lambda,m) [rad], lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
                psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)) [rad], g = angle between N and N''
                psih = angle between X axis and N [rad]


    ***
    satellite=Dict("MomentsOfInertia"=>IV,"numberOfFacets"=>n,"facets"=>facets);
    IV = [A,B,C] : principal moment of inertia
    IM = intrinsic magnetic moment of the satellite [A m^2]
    n  = number of facets in which the satellite surface is divided
    facets = [facet1...facetn]
        where  
            facetk = [facet_coeff,facet_area,facet_Vinfo,facet_nv]
            facet_coeff  = [cai*psrp,cdi*psrsp,csi*psrp], psrp = solar radiation pressure
            facet_area   = area of the facet [m^2]
            facet_Vinfo  = [rhoiv,vv1,vv2,vv3], rhoiv = centre of mass to facet centroid vector [m], vvi = verteces of the facet [m], i=1..3
            facet_nv     = normal unit vector  
            facet_VVT1   = vector containing combinations of cai*psrp with the components of  normal unit vector and components of rhoiv and area (output of srp_average-->srp_initforaverage_T1(cai,surfi,rhoiv,niv))
            facet_VVT2   = vector containing combinations of cdi*psrp with the components of  normal unit vector and components of rhoiv and area (output of srp_average-->srp_initforaverage_T2(cai,surfi,rhoiv,niv))
            facet_VVT3   = vector containing combinations of csi*psrp with the components of  normal unit vector and components of rhoiv and area (output of srp_average-->srp_initforaverage_T3(cai,surfi,rhoiv,niv))
            with 
                cai = (1-rf*sp)
                cdi = (2/3*(1-sp)*rf+2/3*(1-rf))
                csi = 2*rf*sp*psrp
                with rf reflectivity and sp specular coefficient (see Benson&Sheers2021)                       
        
    OUTPUT
    WSRP2 = vector 7x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,J5,J6,J7] -> [skm',Jg',Jh',psil',J5',J6',J7'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWSRP2sadovlike(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    W2 = getWSRP2procedure(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,2);
    # W2 = getWSRP2procedure_numerical(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,2);
    return W2;
end

"""
    getWDRAG2(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion 
    (expressed in Sadov variables) with respect to Sadov fast angles
    -- the perturbation here is the drag
        
    INPUT
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    m...WTS = terms depending on m pre-computed with the functions in termsWtbWtsForSrpDragW2 in this file
    equi = equinoctial elements (true longitude)
    satellite = Dict with the characteristics of the satellite ***
    cpai = components of the polar axis of the central body expressed in the inertial reference frame
    muPlanet = gravitational parameter of the central planet
    atmosphericinfo = Dict with info relative to the atmospheric model ****
    nE = rotational angular rate of the central planet
    rPlanet = mean radius of the central planet
    inertial2equatorialRM = rotation matrix from inertial reference frame to equatorial reference frame
    mjd2000 = epoch in modified julian day 2000

    **
    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    The reference frame is such that m belongs to in[0,1] and k>0 where
    k = C/A*(B-A)/(C-B) 
    m = k*(C-delta)/(delta-A)
    with
    A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
    delta = G^2/2/Phi, with G the angular momentum of the body and Phi the kinetic energy
    Let GV=diag([A,B,C])omV be angular momentum of the body and G = norm(GV)
    (omV = [p,q,r]^T is the angular velocity) 
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular 
    to GV and the xy plane. 
    sadovvar = [Jl,Jg,Jh,psil,psig,psih] or [Jg,Jh,psil,psig,psih]
            with
                Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m)) [kg m^2 /s^2]
                Jg   = G [kg m^2 /s^2]
                Jh   = Gcos(delta), delta = inclination angle between [kg m^2 /s^2]
                       the XY plane and the plane perpendicular to GV.
                psil = pi/2 ellK(m)ellF(lambda,m) [rad], lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
                psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)) [rad], g = angle between N and N''
                psih = angle between X axis and N [rad]
    

    ***
    satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD,"constantstermsdragtorque"=>constantstermsdragtorque);
    IV = [A,B,C] : principal moment of inertia
    IM = intrinsic magnetic moment of the satellite [A m^2]
    n  = number of facets in which the satellite surface is divided
    facets = [facet1...facetn]
        where  
            facetk = [facet_coeff,facet_area,facet_Vinfo,facet_nv]
            facet_coeff = [cai*psrp,cdi*psrsp,csi*psrp], psrp = solar radiation pressure
            facet_area  = area of the facet [m^2]
            facet_Vinfo  = [rhoiv,vv1,vv2,vv3], rhoiv = centre of mass to facet centroid vector [m], vvi = verteces of the facet [m], i=1..3
            facet_nv    = normal unit vector  
            with 
                cai = (1-rf*sp)
                cdi = (2/3*(1-sp)*rf+2/3*(1-rf))
                csi = 2*rf*sp*psrp
                with rf reflectivity and sp specular coefficient (see Benson&Sheers2021)
    CD: aerodynamic coefficient
    constantstermsdragtorque : array of two vectors collecting constants terms depending on the satellite geometry appearing in the attitude equations

    ****atmosphericinfo
    Dict with keys:
        atmosphericmodel = 1 exponential atmospheric model
        2 nrlmsise00
        other kind of model still not implemented

       
        AtmM =   
            if atmosphericmodel==1
            AtmM = matrix containing the exponential atmospheric model 
            in each row: [min altitude [km], density [kg/m^3], scale height [km]

    
            if atmosphericmodell==2
            AtmM = matrix with solar flux data with format
            In each row [jd-AP1-AP2-AP3-AP4-AP5-AP6-AP7-AP8-AP9-AP_AVG-f107_OBS-f107_OBS_CENTER81]
            with
            AP1 : Planetary Equivalent Amplitude (Ap) for 0000-0300 UT.
            AP2 : Planetary Equivalent Amplitude (Ap) for 0300-0600 UT.
            AP3 : Planetary Equivalent Amplitude (Ap) for 0600-0900 UT.
            AP4 : Planetary Equivalent Amplitude (Ap) for 0900-1200 UT.
            AP5 : Planetary Equivalent Amplitude (Ap) for 1200-1500 UT.
            AP6 : Planetary Equivalent Amplitude (Ap) for 1500-1800 UT.
            AP7 : Planetary Equivalent Amplitude (Ap) for 1800-2100 UT.
            AP8 : Planetary Equivalent Amplitude (Ap) for 2100-0000 UT.
            AP_AVG : Arithmetic average of the 8 Ap indices for the day.
            F10.7_OBS : Observed 10.7-cm Solar Radio Flux (F10.7). Measured at Ottawa at 1700 UT daily from 1947 Feb 14 until 1991 May 31 and measured at Penticton at 2000 UT from 1991 Jun 01 on. Expressed in units of 10-22 W/m2/Hz.
            F10.7_OBS_CENTER81 : Centered 81-day arithmetic average of F10.7 (observed).

        considerthermalCD = boolean (if true in the drag model also friction effects are included)

    OUTPUT
    W2 = vector 6x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,psig,psih] -> [skm',Jg',Jh',psil',psig',psih'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWDRAG2(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    W2 = getWDRAG2procedure(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,1)
    # W2 = getWDRAG2procedure_numerical(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,1,true);
   
    return W2;
end

"""
    getWDRAG2(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    Returns the generating function leading to the transformation 
        of variables allowing to average the attitude equations of motion 
        (expressed in Sadov-like variables) with respect to Sadov fast angles
        -- the perturbation here is the drag
        
    INPUT
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    m...WTS = terms depending on m pre-computed with the functions in termsWtbWtsForSrpDragW2 in this file
    equi = equinoctial elements (true longitude)
    satellite = Dict with the characteristics of the satellite ***
    cpai = components of the polar axis of the central body expressed in the inertial reference frame
    muPlanet = gravitational parameter of the central planet
    atmosphericinfo = Dict with info relative to the atmospheric model ****
    nE = rotational angular rate of the central planet
    rPlanet = mean radius of the central planet
    inertial2equatorialRM = rotation matrix from inertial reference frame to equatorial reference frame
    mjd2000 = epoch in modified julian day 2000

    **
    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    The reference frame is such that m belongs to in[0,1] and k>0 where
    k = C/A*(B-A)/(C-B) 
    m = k*(C-delta)/(delta-A)
    with
    A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
    delta = G^2/2/Phi, with G the angular momentum of the body and Phi the kinetic energy
    Let GV=diag([A,B,C])omV be angular momentum of the body and G = norm(GV)
    (omV = [p,q,r]^T is the angular velocity) 
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular 
    to GV and the xy plane. 
    sadovvar = [Jl,Jg,Jh,psil,psig,psih] or [Jg,Jh,psil,psig,psih]
            with
                Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m)) [kg m^2 /s^2]
                Jg   = G [kg m^2 /s^2]
                Jh   = Gcos(delta), delta = inclination angle between [kg m^2 /s^2]
                       the XY plane and the plane perpendicular to GV.
                psil = pi/2 ellK(m)ellF(lambda,m) [rad], lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
                psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)) [rad], g = angle between N and N''
                psih = angle between X axis and N [rad]
    

    ***
    satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD,"constantstermsdragtorque"=>constantstermsdragtorque);
    IV = [A,B,C] : principal moment of inertia
    IM = intrinsic magnetic moment of the satellite [A m^2]
    n  = number of facets in which the satellite surface is divided
    facets = [facet1...facetn]
        where  
            facetk = [facet_coeff,facet_area,facet_Vinfo,facet_nv]
            facet_coeff = [cai*psrp,cdi*psrsp,csi*psrp], psrp = solar radiation pressure
            facet_area  = area of the facet [m^2]
            facet_Vinfo  = [rhoiv,vv1,vv2,vv3], rhoiv = centre of mass to facet centroid vector [m], vvi = verteces of the facet [m], i=1..3
            facet_nv    = normal unit vector  
            with 
                cai = (1-rf*sp)
                cdi = (2/3*(1-sp)*rf+2/3*(1-rf))
                csi = 2*rf*sp*psrp
                with rf reflectivity and sp specular coefficient (see Benson&Sheers2021)
    CD: aerodynamic coefficient
    constantstermsdragtorque : array of two vectors collecting constants terms depending on the satellite geometry appearing in the attitude equations

    ****atmosphericinfo
    Dict with keys:
        atmosphericmodel = 1 exponential atmospheric model
        2 nrlmsise00
        other kind of model still not implemented

       
        AtmM =   
            if atmosphericmodel==1
            AtmM = matrix containing the exponential atmospheric model 
            in each row: [min altitude [km], density [kg/m^3], scale height [km]

    
            if atmosphericmodell==2
            AtmM = matrix with solar flux data with format
            In each row [jd-AP1-AP2-AP3-AP4-AP5-AP6-AP7-AP8-AP9-AP_AVG-f107_OBS-f107_OBS_CENTER81]
            with
            AP1 : Planetary Equivalent Amplitude (Ap) for 0000-0300 UT.
            AP2 : Planetary Equivalent Amplitude (Ap) for 0300-0600 UT.
            AP3 : Planetary Equivalent Amplitude (Ap) for 0600-0900 UT.
            AP4 : Planetary Equivalent Amplitude (Ap) for 0900-1200 UT.
            AP5 : Planetary Equivalent Amplitude (Ap) for 1200-1500 UT.
            AP6 : Planetary Equivalent Amplitude (Ap) for 1500-1800 UT.
            AP7 : Planetary Equivalent Amplitude (Ap) for 1800-2100 UT.
            AP8 : Planetary Equivalent Amplitude (Ap) for 2100-0000 UT.
            AP_AVG : Arithmetic average of the 8 Ap indices for the day.
            F10.7_OBS : Observed 10.7-cm Solar Radio Flux (F10.7). Measured at Ottawa at 1700 UT daily from 1947 Feb 14 until 1991 May 31 and measured at Penticton at 2000 UT from 1991 Jun 01 on. Expressed in units of 10-22 W/m2/Hz.
            F10.7_OBS_CENTER81 : Centered 81-day arithmetic average of F10.7 (observed).

        considerthermalCD = boolean (if true in the drag model also friction effects are included)

    OUTPUT
    W2 = vector 7x1 of Float64 -
    W2 = vector 7x1 of Float64 -
    generating function of the transformation 
    [skm,Jg,Jh,psil,J5,J6,J7] -> [skm',Jg',Jh',psil',J5',J6',J7'], 
    allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWDRAG2sadovlike(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    W2 = getWDRAG2procedure(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,2);
    # W2 = getWDRAG2procedure_numerical(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,2,true);
    return W2;
end

################## FUNCTIONS FOR SRP and DRAG analytical
function termsWtbWtsForSrpDragW2(skm,k)
    # # coeff matrices
    if skm == 1.0
        # case smk=0.0
        m = 0.0;
        m    = 0.0
        KK   = pi/2;
        EE   = pi/2;
        PP   = pi/2/sqrt(1+k);
        znC1 = 0;
        znC2 = 0;
        znC3 = 0;
        WTb = zeros(2373)
        WTS = zeros(23)
    else
        if k!=0.0
            m     = (1-skm)/skm*k;
            KK    = ellF(pi/2,m);
            EE    = ellE(pi/2,m);
            KP    = ellF(pi/2,1-m);
            PP    = ellP(-k,pi/2,m);
            qn     = exp(-pi*KP/KK);
            sigma = pi*ellF(atan(sqrt(k/m)),1-m)/2/KK;
            M4    = (PP-KK)*sqrt(1+k)/sqrt(skm)/pi;
            WTb = getWTb(k,m,KK,EE,KP,qn,sigma,M4);
            WTS = getWTS(k,m,KK,EE,KP,qn,sigma,M4);
            znC1 = 2*pi*qn/(KK*(-qn^2 + 1));
            znC2 = 2*pi*qn^2/(KK*(-qn^4 + 1));
            znC3 = 2*pi*qn^3/(KK*(-qn^6 + 1));
        else
            m    = 0.0
            KK   = pi/2;
            EE   = pi/2;
            PP   = pi/2;
            WTb = getWTb_k0(skm,1.0-skm);
            WTS = getWTS_k0(skm,1.0-skm);
            znC1 = 0;
            znC2 = 0;
            znC3 = 0;
        end
    end
    return m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS;
end


################## SRP
function getWSRP2procedure(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,typeofvariable)
    
    casesmk0 = false;
    if skm==1.0 
        casesmk0 = true
    end
    # println([KK,EE,PP])

    # vector of constants
    A = satellite["MomentsOfInertia"][1]; C = satellite["MomentsOfInertia"][3];
    numberOfFacets = get(satellite,"numberOfFacets",0);
    facets = get(satellite,"facets",0.0);  
    VVCT = zeros(195,1);
    for kk = 1:numberOfFacets
        VVT1   = facets[kk][5];
        VVT2   = facets[kk][6];
        VVT3   = facets[kk][7];
        VVA    = zeros(195,1);
        VVA[1:30,1]  = VVT1;
        VVA[31:90,1] = VVT2;
        VVA[91:195,1] = VVT3;
        VVCT = VVCT + VVA;
    end
    VVCT = transpose(VVCT);

    # sadov
    if size(sadovvar)[1] == 5
        Jg = sadovvar[1];
        Jh = sadovvar[2];
        psil = sadovvar[3];
        psig = sadovvar[4];
        psih = sadovvar[5];
    else
        Jg = sadovvar[2];
        Jh = sadovvar[3];
        psil = sadovvar[4];
        psig = sadovvar[5];
        psih = sadovvar[6];
    end    
    
    # uVectH --> satellite to Sun unit vector approximated as Earth to Sun unit vector in the ref frame with angular momentum as z axis
    rSunnorm = norm(rSun)
    uVect = rSun/rSunnorm;
    cpsih = cos(psih); 
    spsih = sin(psih);
    cdelta = Jh/Jg;
    sdelta = sqrt(1-cdelta^2);
    Rhdelta = zeros(3,3);
    Rhdelta[1,1] = cpsih;
    Rhdelta[1,2] = spsih;
    Rhdelta[1,3] = 0.0;
    Rhdelta[2,1] = -cdelta*spsih;
    Rhdelta[2,2] = cdelta*cpsih;
    Rhdelta[2,3] = sdelta;
    Rhdelta[3,1] = sdelta*spsih;
    Rhdelta[3,2] = -sdelta*cpsih;
    Rhdelta[3,3] = cdelta;
    uVectH = Rhdelta*uVect;
    
    # vcu (vector depending on psih, delta and the components of vector uVect)
    vcu = zeros(20);
    count = 0;
    for ii=0:3
        for jj=0:3
            for kk=0:3
                if ii+jj+kk>3 
                    break;
                end
                count=count+1;
                vcu[count] = uVectH[1]^ii*uVectH[2]^jj*uVectH[3]^kk;
            end
        end
    end
     
    if !casesmk0
        CMXb11 = getCMXb11(WTb,vcu,1);
        CMYb21 = getCMYb21(WTb,vcu,1);
        CMZb31 = getCMZb31(WTb,vcu,1);
        CMXb12 = getCMXb12(WTb,vcu,1);
        CMYb22 = getCMYb22(WTb,vcu,1);
        CMZb32 = getCMZb32(WTb,vcu,1);
        CMXb13 = getCMXb13(WTb,vcu,1);
        CMYb23 = getCMYb23(WTb,vcu,1);
        CMZb33 = getCMZb33(WTb,vcu,1);
        CMXSx  = getCMXSx(WTb,WTS,vcu,znC1,znC2,znC3,1);
        CMYSy  = getCMYSy(WTb,WTS,vcu,znC1,znC2,znC3,1);
        CMZSz  = getCMZSz(WTb,WTS,vcu,znC1,znC2,znC3,1);
    else
        CMXb11,CMXb12,CMXb13,CMYb21,CMYb22,CMYb23,CMZb31,CMZb32,CMZb33 = CMterms_smk0(vcu,1);
    end

    # VGT and VGGT ---> vectors containing fast variables
    npsig = Jg/A/C*((C-A)PP+A*KK)/KK;
    npsil = -pi/2.0/sqrt(1+k)*sqrt(skm)*Jg/A/C*(C-A)/KK;
    if !casesmk0
        VGT,VGGT = getgeneratinfunctiontermsW2_srp_drag_psil_psig(psil,psig,npsil,npsig);
    else
        VGT,VGGT =  getgeneratinfunctiontermsW2_allpert_psil_psig_smk0(psil,psig,npsil,npsig)
    end

    

    # shadowfun
    if includeEclipsesEffects  && passageinshadowoccurs
        EL1 = ELin;
        EL2 = ELout;
        sf = smoothshadowfun(EL,EL1,EL2);
    else
        sf = 1.0;
    end 

    # terms of the generating function
    W21Block  = -2*skm/Jg*CMXb13 - 2*skm/Jg*(1-m)/(1+k)*CMYb23 + 2*(1-skm)/Jg*CMZb33;
    Tb3Block   = CMXb13+CMYb23+CMZb33;
    Tb2Block   = CMXb12+CMYb22+CMZb32;
    Tb1Block   = CMXb11+CMYb21+CMZb31;

    # gen fun
    dnpsildskm = pi*Jg*(A-C)*(EE-(1-m)*skm*KK)/(4.0*A*C*(1-skm)*sqrt(skm)*sqrt(1+k)-KK^2.0);
    dnpsildJg  = pi*sqrt(skm)*(A-C)/(2.0*A*C*sqrt(1+k)*KK);
    dnpsigdskm = -Jg*(A-C)*((EE-(1-m)*skm*KK)*PP-(1-skm)*KK*EE)/(2*A*C*(1-m)*skm*(1-skm)*KK^2);
    dnpsigdJg  = ((C-A)*PP+A*KK)/(A*C*KK);

    if casesmk0
        dnpsigdskm = -Jg*(A-C)*(-k/sqrt(1+k)/2.0+1.0)/(2*A*C);
    else
        dnpsigdskm = -Jg*(A-C)*((EE-(1-m)*skm*KK)*PP-(1-skm)*KK*EE)/(2*A*C*(1-m)*skm*(1-skm)*KK^2);
    end

    if !casesmk0
        TSBlock    = CMXSx-(1-m)*CMYSy+CMZSz;

        if typeofvariable==1

            W21 = VVCT*W21Block*VGT;
            W22 = VVCT*Tb3Block*VGT;
            W23 = VVCT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = VVCT*(pi/(2*(m-1)*KK*Jg)*TSBlock*VGT + dnpsildskm*W21Block*VGGT + dnpsildJg*Tb3Block*VGGT);
            W25 = VVCT*((PP-(1-skm)*KK)*sqrt(1+k)/(Jg*KK*sqrt(skm)*(1-m))*TSBlock*VGT - cdelta/sdelta/Jg*Tb1Block*VGT + dnpsigdskm*W21Block*VGGT + dnpsigdJg*Tb3Block*VGGT);
            W26 = VVCT*1/sdelta/Jg*Tb1Block*VGT;  

            W2  = append!(W21,W22,W23,W24,W25,W26)*sf;

        else
            W21 = VVCT*W21Block*VGT;
            W22 = VVCT*Tb3Block*VGT;
            W23 = VVCT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = VVCT*(pi/(2*(m-1)*KK*Jg)*TSBlock*VGT + dnpsildskm*W21Block*VGGT + dnpsildJg*Tb3Block*VGGT);
            W25 = VVCT*((PP-(1-skm)*KK)*sqrt(1+k)/(Jg*KK*sqrt(skm)*(1-m))*TSBlock*VGT + psih*sdelta/Jg*Tb2Block*VGT + dnpsigdskm*W21Block*VGGT + dnpsigdJg*Tb3Block*VGGT);
            W26 = VVCT*(-spsih*Tb1Block+sdelta*cpsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
            W27 = VVCT*(cpsih*Tb1Block+sdelta*spsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
        
            W2  = append!(W21,W22,W23,W24,W25,W26,W27)*sf;

        end
    else
        if typeofvariable==1
            W21 = VVCT*W21Block*VGT;
            W22 = VVCT*Tb3Block*VGT;
            W23 = VVCT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = [0.0];
            W25 = VVCT*(- cdelta/sdelta/Jg*Tb1Block*VGT + (dnpsigdskm+dnpsildskm)*W21Block*VGGT + (dnpsigdJg+dnpsildJg)*Tb3Block*VGGT);
            W26 = VVCT*1/sdelta/Jg*Tb1Block*VGT;  
            W2  = append!(W21,W22,W23,W24,W25,W26)*sf;
        else
            W21 = VVCT*W21Block*VGT;
            W22 = VVCT*Tb3Block*VGT;
            W23 = VVCT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = [0.0];
            W25 = VVCT*(psih*sdelta/Jg*Tb2Block*VGT + (dnpsigdskm+dnpsildskm)*W21Block*VGGT + (dnpsigdJg+dnpsildJg)*Tb3Block*VGGT);
            W26 = VVCT*(-spsih*Tb1Block+sdelta*cpsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
            W27 = VVCT*(cpsih*Tb1Block+sdelta*spsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
            W2  = append!(W21,W22,W23,W24,W25,W26,W27)*sf;
        end
    end

    au = astroconstants(0)[2];
    psrpcorr     = (au/rSunnorm)^2.0
    W2 = W2*psrpcorr;
    
    return W2;
end

function getWSRP2procedure_numerical(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,typeofvariable)
    
    casesmk0 = false;
    if skm==1.0 
        casesmk0 = true
    end

    # vector of constants
    A = satellite["MomentsOfInertia"][1]; C = satellite["MomentsOfInertia"][3];
    numberOfFacets = get(satellite,"numberOfFacets",0);
    facets = get(satellite,"facets",0.0);  
    VVT1 = zeros(30);
    VVT2 = zeros(60);
    VVT3 = zeros(105);
    VVT4 = zeros(68);
    for kk = 1:numberOfFacets
        VVT1   = VVT1 + facets[kk][5];
        VVT2   = VVT2 + facets[kk][6];
        VVT3   = VVT3 + facets[kk][7];
    end  
    
    # sadov
    if size(sadovvar)[1] == 5
        Jg = sadovvar[1];
        Jh = sadovvar[2];
        psil = sadovvar[3];
        psig = sadovvar[4];
        psih = sadovvar[5];
    else
        Jg = sadovvar[2];
        Jh = sadovvar[3];
        psil = sadovvar[4];
        psig = sadovvar[5];
        psih = sadovvar[6];
    end    
    cpsih = cos(psih); 
    spsih = sin(psih);
    cdelta = Jh/Jg;
    sdelta = sqrt(1-cdelta^2);
    
    # uVect, uVectH --> satellite to Sun unit vector approximated as Earth to Sun unit vector in the ref frame with angular momentum as z axis
    rSunnorm = norm(rSun);
    uVect = rSun/rSunnorm;
    
    
    # mean motion torque free problem and derivatives
    npsig = Jg/A/C*((C-A)PP+A*KK)/KK;
    npsil = -pi/2.0/sqrt(1+k)*sqrt(skm)*Jg/A/C*(C-A)/KK;
    dnpsildskm = pi*Jg*(A-C)*(EE-(1-m)*skm*KK)/(4.0*A*C*(1-skm)*sqrt(skm)*sqrt(1+k)-KK^2.0);
    dnpsildJg  = pi*sqrt(skm)*(A-C)/(2.0*A*C*sqrt(1+k)*KK);
    if casesmk0
        dnpsigdskm = -Jg*(A-C)*(-k/sqrt(1+k)/2.0+1.0)/(2*A*C);
    else
        dnpsigdskm = -Jg*(A-C)*((EE-(1-m)*skm*KK)*PP-(1-skm)*KK*EE)/(2*A*C*(1-m)*skm*(1-skm)*KK^2);
    end
    dnpsigdJg  = ((C-A)*PP+A*KK)/(A*C*KK);
     
    if !casesmk0
        # constants
        k1r = sqrt(1+k);
        smk = 1.0-skm;
        mM1KK = KK*(m-1.0);

        if m==0.0
            T4 = -pi/2.0/sqrt(1+k);
        elseif m==1.0
            T4 = -atan(sqrt(k))/sqrt(k);
        else
            T4 = (KK*m - PP*(m+k))/k;
        end

        vcu = zeros(20);
        count = 0;
        for ii=0:3
            for jj=0:3
                for kk=0:3
                    if ii+jj+kk>3 
                        break;
                    end
                    count=count+1;
                    vcu[count] = uVect[1]^ii*uVect[2]^jj*uVect[3]^kk;
                end
            end
        end

        # npsil*diff(W_i,psil)+npsig*diff(W_i,psig)=f(psil,psig)-z
        # W_i = W_i_f + W_i_z;  
        # npsil*diff(W_i_z,psil)+npsig*diff(W_i_z,psig)=-z
        # npsil*diff(W_i_f,psil)+npsig*diff(W_i_f,psig)=f(psil,psig)
        
        ######
        # W_i_z
        ZV = getaveragedvectorfield_srp(k,k1r,skm,smk,m,KK,EE,PP,T4,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,vcu,VVT1,VVT2,VVT3,typeofvariable);
        atanpsig = atan(tan(psig/2.0));
        atanpsil = atan(tan(psil/2.0));
        W2_Z = - ZV*(atanpsil/npsil+atanpsig/npsig);
        W2_Z[4] = W2_Z[4] - (dnpsildskm*ZV[1]+dnpsildJg*ZV[2])*((atanpsil^2.0-pi^2.0/12)/(npsil^2) + (atanpsig^2.0-pi^2.0/12)/(npsig^2));
        W2_Z[5] = W2_Z[5] - (dnpsigdskm*ZV[1]+dnpsigdJg*ZV[2])*((atanpsil^2.0-pi^2.0/12)/(npsil^2) + (atanpsig^2.0-pi^2.0/12)/(npsig^2));
    
        f2V,psilV,zetaV,ratio  = getnonaveragedvectorfield_srp(k,k1r,skm,smk,m,KK,EE,PP,T4,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,vcu,VVT1,VVT2,VVT3,typeofvariable,npsil,npsig,1*pi/180.0,1.0*pi/180.0);
        W2_f       = getgenfunterms_psilpsig(npsil,ratio,dnpsildskm,dnpsildJg,dnpsigdskm,dnpsigdJg,psig,psil,psilV,zetaV,f2V)

        # W2        
        W2 = W2_f + W2_Z;

    else
        VVCT = zeros(195,1);
        VVCT[1:30] = transpose(VVT1);
        VVCT[31:90] = transpose(VVT2);
        VVCT[91:end] = transpose(VVT3);
        uVectH = zeros(3);
        uVectH[1] = cpsih*uVect[1] + spsih*uVect[2];
        uVectH[2] = -cdelta*spsih*uVect[1] + cdelta*cpsih*uVect[2] + sdelta*uVect[3]
        uVectH[3] = spsih*sdelta*uVect[1] - cpsih*sdelta*uVect[2] + cdelta*uVect[3]
        vcuH = zeros(20);
        count = 0;
        for ii=0:3
            for jj=0:3
                for kk=0:3
                    if ii+jj+kk>3 
                        break;
                    end
                    count=count+1;
                    vcuH[count] = uVectH[1]^ii*uVectH[2]^jj*uVectH[3]^kk;
                end
            end
        end
        CMXb11,CMXb12,CMXb13,CMYb21,CMYb22,CMYb23,CMZb31,CMZb32,CMZb33 = CMterms_smk0(vcuH,1);
        VGT,VGGT =  getgeneratinfunctiontermsW2_allpert_psil_psig_smk0(psil,psig,npsil,npsig)
        W21Block  = -2*skm/Jg*CMXb13 - 2*skm/Jg*(1-m)/(1+k)*CMYb23 + 2*(1-skm)/Jg*CMZb33;
        Tb3Block   = CMXb13+CMYb23+CMZb33;
        Tb2Block   = CMXb12+CMYb22+CMZb32;
        Tb1Block   = CMXb11+CMYb21+CMZb31;
        if typeofvariable==1
            W21 = VVCT*W21Block*VGT;
            W22 = VVCT*Tb3Block*VGT;
            W23 = VVCT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = [0.0];
            W25 = VVCT*(- cdelta/sdelta/Jg*Tb1Block*VGT + (dnpsigdskm+dnpsildskm)*W21Block*VGGT + (dnpsigdJg+dnpsildJg)*Tb3Block*VGGT);
            W26 = VVCT*1/sdelta/Jg*Tb1Block*VGT;  
            W2  = append!(W21,W22,W23,W24,W25,W26);
        else
            W21 = VVCT*W21Block*VGT;
            W22 = VVCT*Tb3Block*VGT;
            W23 = VVCT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = [0.0];
            W25 = VVCT*(psih*sdelta/Jg*Tb2Block*VGT + (dnpsigdskm+dnpsildskm)*W21Block*VGGT + (dnpsigdJg+dnpsildJg)*Tb3Block*VGGT);
            W26 = VVCT*(-spsih*Tb1Block+sdelta*cpsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
            W27 = VVCT*(cpsih*Tb1Block+sdelta*spsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
            W2  = append!(W21,W22,W23,W24,W25,W26,W27);
        end
    end    

    # shadowfun
    if includeEclipsesEffects  && passageinshadowoccurs
        EL1 = ELin;
        EL2 = ELout;
        sf = smoothshadowfun(EL,EL1,EL2);
    else
        sf = 1.0;
    end 

    au = astroconstants(0)[2];
    psrpcorr     = (au/rSunnorm)^2.0

    W2 = W2*sf*psrpcorr;
    
    return W2;
end

function smoothshadowfun(EL,EL1,EL2)
    sf = 1 + 0.1e1 / pi * sin(-17 * EL2 + 17 * EL) / 17 + 0.1e1 / pi * sin(-16 * EL2 + 16 * EL) / 16 + 0.1e1 / pi * sin(-15 * EL2 + 15 * EL) / 15 + 0.1e1 / pi * sin(-14 * EL2 + 14 * EL) / 14 + 0.1e1 / pi * sin(-13 * EL2 + 13 * EL) / 13 + 0.1e1 / pi * sin(-12 * EL2 + 12 * EL) / 12 + 0.1e1 / pi * sin(-11 * EL2 + 11 * EL) / 11 + 0.1e1 / pi * sin(-10 * EL2 + 10 * EL) / 10 + 0.1e1 / pi * sin(-9 * EL2 + 9 * EL) / 9 + 0.1e1 / pi * sin(-8 * EL2 + 8 * EL) / 8 + 0.1e1 / pi * sin(-7 * EL2 + 7 * EL) / 7 + 0.1e1 / pi * sin(-6 * EL2 + 6 * EL) / 6 + 0.1e1 / pi * sin(-5 * EL2 + 5 * EL) / 5 + 0.1e1 / pi * sin(-4 * EL2 + 4 * EL) / 4 + 0.1e1 / pi * sin(-3 * EL2 + 3 * EL) / 3 + 0.1e1 / pi * sin(-2 * EL2 + 2 * EL) / 2 + 0.1e1 / pi * sin(-EL2 + EL) - 0.1e1 / pi * sin(-20 * EL1 + 20 * EL) / 20 - 0.1e1 / pi * sin(-19 * EL1 + 19 * EL) / 19 - 0.1e1 / pi * sin(-18 * EL1 + 18 * EL) / 18 - 0.1e1 / pi * sin(-17 * EL1 + 17 * EL) / 17 - 0.1e1 / pi * sin(-16 * EL1 + 16 * EL) / 16 - 0.1e1 / pi * sin(-15 * EL1 + 15 * EL) / 15 - 0.1e1 / pi * sin(-14 * EL1 + 14 * EL) / 14 - 0.1e1 / pi * sin(-13 * EL1 + 13 * EL) / 13 - 0.1e1 / pi * sin(-12 * EL1 + 12 * EL) / 12 - 0.1e1 / pi * sin(-11 * EL1 + 11 * EL) / 11 - 0.1e1 / pi * sin(-10 * EL1 + 10 * EL) / 10 - 0.1e1 / pi * sin(-9 * EL1 + 9 * EL) / 9 - 0.1e1 / pi * sin(-8 * EL1 + 8 * EL) / 8 - 0.1e1 / pi * sin(-7 * EL1 + 7 * EL) / 7 - 0.1e1 / pi * sin(-6 * EL1 + 6 * EL) / 6 - 0.1e1 / pi * sin(-5 * EL1 + 5 * EL) / 5 - 0.1e1 / pi * sin(-4 * EL1 + 4 * EL) / 4 - 0.1e1 / pi * sin(-3 * EL1 + 3 * EL) / 3 - 0.1e1 / pi * sin(-2 * EL1 + 2 * EL) / 2 - 0.1e1 / pi * sin(-EL1 + EL) + 0.1e1 / pi * sin(-20 * EL2 + 20 * EL) / 20 + 0.1e1 / pi * sin(-19 * EL2 + 19 * EL) / 19 + 0.1e1 / pi * sin(-18 * EL2 + 18 * EL) / 18 + 0.1e1 / pi * EL1 / 2 - 0.1e1 / pi * EL2 / 2;
    return sf;
end

function getgeneratinfunctiontermsW2_srp_drag_psil_psig(psil,psig,npsil,npsig)
    VG = zeros(158);
    VGG = zeros(158);

    t1 = sin(psig)
    t2 = 0.1e1 / npsig
    VG[1] = t2 * t1
    t3 = sin(psil)
    t4 = 0.1e1 / npsil
    VG[2] = t4 * t3
    t5 = 2 * psig
    t6 = sin(t5)
    VG[3] = t2 * t6 / 2
    t8 = 3 * psig
    t9 = sin(t8)
    VG[4] = t2 * t9 / 3
    t11 = 4 * psig
    t12 = sin(t11)
    VG[5] = t2 * t12 / 4
    t14 = 2 * psil
    t15 = sin(t14)
    VG[6] = t4 * t15 / 2
    t17 = 3 * psil
    t18 = sin(t17)
    VG[7] = t4 * t18 / 3
    t20 = 4 * psil
    t21 = sin(t20)
    VG[8] = t4 * t21 / 4
    t23 = 5 * psil
    t24 = sin(t23)
    VG[9] = t4 * t24 / 5
    t26 = 6 * psil
    t27 = sin(t26)
    VG[10] = t4 * t27 / 6
    t29 = psig - t20
    t30 = sin(t29)
    t31 = 4 * npsil
    t32 = -t31 + npsig
    t33 = 1 / t32
    VG[11] = t33 * t30
    t34 = psig + t20
    t35 = sin(t34)
    t36 = t31 + npsig
    t37 = 1 / t36
    VG[12] = t37 * t35
    t38 = t5 - t20
    t39 = sin(t38)
    t40 = 2 * npsig
    t41 = -t31 + t40
    t42 = 1 / t41
    VG[13] = t42 * t39
    t43 = t5 - t17
    t44 = sin(t43)
    t45 = 3 * npsil
    t46 = -t45 + t40
    t47 = 1 / t46
    VG[14] = t47 * t44
    t48 = -psil + psig
    t49 = 2 * t48
    t50 = sin(t49)
    t51 = -npsil + npsig
    t52 = 2 * t51
    t53 = 1 / t52
    VG[15] = t53 * t50
    t54 = t5 - psil
    t55 = sin(t54)
    t56 = -npsil + t40
    t57 = 1 / t56
    VG[16] = t57 * t55
    t58 = t5 + psil
    t59 = sin(t58)
    t60 = npsil + t40
    t61 = 1 / t60
    VG[17] = t61 * t59
    t62 = psil + psig
    t63 = 2 * t62
    t64 = sin(t63)
    t65 = npsil + npsig
    t66 = 2 * t65
    t67 = 1 / t66
    VG[18] = t67 * t64
    t68 = t5 + t17
    t69 = sin(t68)
    t70 = t45 + t40
    t71 = 1 / t70
    VG[19] = t71 * t69
    t72 = t5 + t20
    t73 = sin(t72)
    t74 = t31 + t40
    t75 = 1 / t74
    VG[20] = t75 * t73
    t76 = t11 - t26
    t77 = sin(t76)
    t78 = 6 * npsil
    t79 = 4 * npsig
    t80 = -t78 + t79
    t81 = 1 / t80
    VG[21] = t81 * t77
    t82 = 4 * t48
    t83 = sin(t82)
    t84 = 4 * t51
    t85 = 1 / t84
    VG[22] = t85 * t83
    t86 = t11 - t14
    t87 = sin(t86)
    t88 = 2 * npsil
    t89 = -t88 + t79
    t90 = 1 / t89
    VG[23] = t90 * t87
    t91 = t11 + t14
    t92 = sin(t91)
    t93 = t88 + t79
    t94 = 1 / t93
    VG[24] = t94 * t92
    t95 = 4 * t62
    t96 = sin(t95)
    t97 = 4 * t65
    t98 = 1 / t97
    VG[25] = t98 * t96
    t99 = t11 + t26
    t100 = sin(t99)
    t101 = t78 + t79
    t102 = 1 / t101
    VG[26] = t102 * t100
    t103 = -t26 + psig
    t104 = sin(t103)
    t105 = -t78 + npsig
    t106 = 1 / t105
    VG[27] = t106 * t104
    t107 = -t26 + t5
    t108 = sin(t107)
    t109 = -t78 + t40
    t110 = 1 / t109
    VG[28] = t110 * t108
    t111 = -t26 + t8
    t112 = sin(t111)
    t113 = 3 * npsig
    t114 = -t78 + t113
    t115 = 1 / t114
    VG[29] = t115 * t112
    t116 = -t23 + psig
    t117 = sin(t116)
    t118 = 5 * npsil
    t119 = -t118 + npsig
    t120 = 1 / t119
    VG[30] = t120 * t117
    t121 = -t23 + t5
    t122 = sin(t121)
    t123 = -t118 + t40
    t124 = 1 / t123
    VG[31] = t124 * t122
    t125 = -t23 + t8
    t126 = sin(t125)
    t127 = -t118 + t113
    t128 = 1 / t127
    VG[32] = t128 * t126
    t129 = -t23 + t11
    t130 = sin(t129)
    t131 = -t118 + t79
    t132 = 1 / t131
    VG[33] = t132 * t130
    t133 = -t20 + t8
    t134 = sin(t133)
    t135 = -t31 + t113
    t136 = 1 / t135
    VG[34] = t136 * t134
    t137 = -t17 + psig
    t138 = sin(t137)
    t139 = -t45 + npsig
    t140 = 1 / t139
    VG[35] = t140 * t138
    t141 = 3 * t48
    t142 = sin(t141)
    t143 = 3 * t51
    t144 = 1 / t143
    VG[36] = t144 * t142
    t145 = -t17 + t11
    t146 = sin(t145)
    t147 = -t45 + t79
    t148 = 1 / t147
    VG[37] = t148 * t146
    t149 = -t14 + psig
    t150 = sin(t149)
    t151 = -t88 + npsig
    t152 = 1 / t151
    VG[38] = t152 * t150
    t153 = -t14 + t8
    t154 = sin(t153)
    t155 = -t88 + t113
    t156 = 1 / t155
    VG[39] = t156 * t154
    t157 = sin(t48)
    t158 = 1 / t51
    VG[40] = t158 * t157
    t159 = -psil + t8
    t160 = sin(t159)
    t161 = -npsil + t113
    t162 = 1 / t161
    VG[41] = t162 * t160
    t163 = -psil + t11
    t164 = sin(t163)
    t165 = -npsil + t79
    t166 = 1 / t165
    VG[42] = t166 * t164
    t167 = sin(t62)
    t168 = 1 / t65
    VG[43] = t168 * t167
    t169 = psil + t8
    t170 = sin(t169)
    t171 = npsil + t113
    t172 = 1 / t171
    VG[44] = t172 * t170
    t173 = psil + t11
    t174 = sin(t173)
    t175 = npsil + t79
    t176 = 1 / t175
    VG[45] = t176 * t174
    t177 = t14 + psig
    t178 = sin(t177)
    t179 = t88 + npsig
    t180 = 1 / t179
    VG[46] = t180 * t178
    t181 = t14 + t8
    t182 = sin(t181)
    t183 = t88 + t113
    t184 = 1 / t183
    VG[47] = t184 * t182
    t185 = t17 + psig
    t186 = sin(t185)
    t187 = t45 + npsig
    t188 = 1 / t187
    VG[48] = t188 * t186
    t189 = 3 * t62
    t190 = sin(t189)
    t191 = 3 * t65
    t192 = 1 / t191
    VG[49] = t192 * t190
    t193 = t17 + t11
    t194 = sin(t193)
    t195 = t45 + t79
    t196 = 1 / t195
    VG[50] = t196 * t194
    t197 = t20 + t8
    t198 = sin(t197)
    t199 = t31 + t113
    t200 = 1 / t199
    VG[51] = t200 * t198
    t201 = t23 + psig
    t202 = sin(t201)
    t203 = t118 + npsig
    t204 = 1 / t203
    VG[52] = t204 * t202
    t205 = t23 + t5
    t206 = sin(t205)
    t207 = t118 + t40
    t208 = 1 / t207
    VG[53] = t208 * t206
    t209 = t23 + t8
    t210 = sin(t209)
    t211 = t118 + t113
    t212 = 1 / t211
    VG[54] = t212 * t210
    t213 = t23 + t11
    t214 = sin(t213)
    t215 = t118 + t79
    t216 = 1 / t215
    VG[55] = t216 * t214
    t217 = t26 + psig
    t218 = sin(t217)
    t219 = t78 + npsig
    t220 = 1 / t219
    VG[56] = t220 * t218
    t221 = t26 + t5
    t222 = sin(t221)
    t223 = t78 + t40
    t224 = 1 / t223
    VG[57] = t224 * t222
    t225 = t26 + t8
    t226 = sin(t225)
    t227 = t78 + t113
    t228 = 1 / t227
    VG[58] = t228 * t226
    t229 = cos(psig)
    VG[59] = -t2 * t229
    t231 = cos(psil)
    VG[60] = -t4 * t231
    t233 = cos(t5)
    VG[61] = -t2 * t233 / 2
    t236 = cos(t8)
    VG[62] = -t2 * t236 / 3
    t239 = cos(t11)
    VG[63] = -t2 * t239 / 4
    t242 = cos(t14)
    VG[64] = -t4 * t242 / 2
    t245 = cos(t17)
    VG[65] = -t4 * t245 / 3
    t248 = cos(t20)
    VG[66] = -t4 * t248 / 4
    t251 = cos(t23)
    VG[67] = -t4 * t251 / 5
    t254 = cos(t26)
    VG[68] = -t4 * t254 / 6
    t257 = cos(t29)
    VG[69] = -t33 * t257
    t259 = cos(t34)
    VG[70] = -t37 * t259
    t261 = cos(t38)
    VG[71] = -t42 * t261
    t263 = cos(t43)
    VG[72] = -t47 * t263
    t265 = cos(t49)
    VG[73] = -t53 * t265
    t267 = cos(t54)
    VG[74] = -t57 * t267
    t269 = cos(t58)
    VG[75] = -t61 * t269
    t271 = cos(t63)
    VG[76] = -t67 * t271
    t273 = cos(t68)
    VG[77] = -t71 * t273
    t275 = cos(t72)
    VG[78] = -t75 * t275
    t277 = cos(t76)
    VG[79] = -t81 * t277
    t279 = cos(t82)
    VG[80] = -t85 * t279
    t281 = cos(t86)
    VG[81] = -t90 * t281
    t283 = cos(t91)
    VG[82] = -t94 * t283
    t285 = cos(t95)
    VG[83] = -t98 * t285
    t287 = cos(t99)
    VG[84] = -t102 * t287
    t289 = cos(t103)
    VG[85] = -t106 * t289
    t291 = cos(t107)
    VG[86] = -t110 * t291
    t293 = cos(t111)
    VG[87] = -t115 * t293
    t295 = cos(t116)
    VG[88] = -t120 * t295
    t297 = cos(t121)
    VG[89] = -t124 * t297
    t299 = cos(t125)
    VG[90] = -t128 * t299
    t301 = cos(t129)
    VG[91] = -t132 * t301
    t303 = cos(t133)
    VG[92] = -t136 * t303
    t305 = cos(t137)
    VG[93] = -t140 * t305
    t307 = cos(t141)
    VG[94] = -t144 * t307
    t309 = cos(t145)
    VG[95] = -t148 * t309
    t311 = cos(t149)
    VG[96] = -t152 * t311
    t313 = cos(t153)
    VG[97] = -t156 * t313
    t315 = cos(t48)
    VG[98] = -t158 * t315
    t317 = cos(t159)
    VG[99] = -t162 * t317
    t319 = cos(t163)
    VG[100] = -t166 * t319
    t321 = cos(t62)
    VG[101] = -t168 * t321
    t323 = cos(t169)
    VG[102] = -t172 * t323
    t325 = cos(t173)
    VG[103] = -t176 * t325
    t327 = cos(t177)
    VG[104] = -t180 * t327
    t329 = cos(t181)
    VG[105] = -t184 * t329
    t331 = cos(t185)
    VG[106] = -t188 * t331
    t333 = cos(t189)
    VG[107] = -t192 * t333
    t335 = cos(t193)
    VG[108] = -t196 * t335
    t337 = cos(t197)
    VG[109] = -t200 * t337
    t339 = cos(t201)
    VG[110] = -t204 * t339
    t341 = cos(t205)
    VG[111] = -t208 * t341
    t343 = cos(t209)
    VG[112] = -t212 * t343
    t345 = cos(t213)
    VG[113] = -t216 * t345
    t347 = cos(t217)
    VG[114] = -t220 * t347
    t349 = cos(t221)
    VG[115] = -t224 * t349
    t351 = cos(t225)
    VG[116] = -t228 * t351
    t353 = 7 * psil
    t354 = sin(t353)
    VG[117] = t4 * t354 / 7
    t356 = 8 * psil
    t357 = sin(t356)
    VG[118] = t4 * t357 / 8
    t359 = t5 - t353
    t360 = sin(t359)
    t361 = 7 * npsil
    t362 = -t361 + t40
    t363 = 1 / t362
    VG[119] = t363 * t360
    t364 = t5 + t353
    t365 = sin(t364)
    t366 = t361 + t40
    t367 = 1 / t366
    VG[120] = t367 * t365
    t368 = t5 + t356
    t369 = sin(t368)
    t370 = 8 * npsil
    t371 = t370 + t40
    t372 = 1 / t371
    VG[121] = t372 * t369
    t373 = -t356 + psig
    t374 = sin(t373)
    t375 = -t370 + npsig
    t376 = 1 / t375
    VG[122] = t376 * t374
    t377 = -t356 + t5
    t378 = sin(t377)
    t379 = -t370 + t40
    t380 = 1 / t379
    VG[123] = t380 * t378
    t381 = -t356 + t8
    t382 = sin(t381)
    t383 = -t370 + t113
    t384 = 1 / t383
    VG[124] = t384 * t382
    t385 = -t353 + psig
    t386 = sin(t385)
    t387 = -t361 + npsig
    t388 = 1 / t387
    VG[125] = t388 * t386
    t389 = -t353 + t8
    t390 = sin(t389)
    t391 = -t361 + t113
    t392 = 1 / t391
    VG[126] = t392 * t390
    t393 = t353 + psig
    t394 = sin(t393)
    t395 = t361 + npsig
    t396 = 1 / t395
    VG[127] = t396 * t394
    t397 = t353 + t8
    t398 = sin(t397)
    t399 = t361 + t113
    t400 = 1 / t399
    VG[128] = t400 * t398
    t401 = t356 + psig
    t402 = sin(t401)
    t403 = t370 + npsig
    t404 = 1 / t403
    VG[129] = t404 * t402
    t405 = t356 + t8
    t406 = sin(t405)
    t407 = t370 + t113
    t408 = 1 / t407
    VG[130] = t408 * t406
    t409 = cos(t353)
    VG[131] = -t4 * t409 / 7
    t412 = cos(t356)
    VG[132] = -t4 * t412 / 8
    t415 = cos(t359)
    VG[133] = -t363 * t415
    t417 = cos(t364)
    VG[134] = -t367 * t417
    t419 = cos(t368)
    VG[135] = -t372 * t419
    t421 = cos(t373)
    VG[136] = -t376 * t421
    t423 = cos(t377)
    VG[137] = -t380 * t423
    t425 = cos(t381)
    VG[138] = -t384 * t425
    t427 = cos(t385)
    VG[139] = -t388 * t427
    t429 = cos(t389)
    VG[140] = -t392 * t429
    t431 = cos(t393)
    VG[141] = -t396 * t431
    t433 = cos(t397)
    VG[142] = -t400 * t433
    t435 = cos(t401)
    VG[143] = -t404 * t435
    t437 = cos(t405)
    VG[144] = -t408 * t437
    t439 = 9 * psil
    t440 = sin(t439)
    VG[145] = t4 * t440 / 9
    t442 = -t439 + psig
    t443 = sin(t442)
    t444 = 9 * npsil
    t445 = -t444 + npsig
    t446 = 1 / t445
    VG[146] = t446 * t443
    t447 = -t439 + t5
    t448 = sin(t447)
    t449 = -t444 + t40
    t450 = 1 / t449
    VG[147] = t450 * t448
    t451 = -t439 + t8
    t452 = sin(t451)
    t453 = -t444 + t113
    t454 = 1 / t453
    VG[148] = t454 * t452
    t455 = t439 + psig
    t456 = sin(t455)
    t457 = t444 + npsig
    t458 = 1 / t457
    VG[149] = t458 * t456
    t459 = t439 + t5
    t460 = sin(t459)
    t461 = t444 + t40
    t462 = 1 / t461
    VG[150] = t462 * t460
    t463 = t439 + t8
    t464 = sin(t463)
    t465 = t444 + t113
    t466 = 1 / t465
    VG[151] = t466 * t464
    t467 = cos(t439)
    VG[152] = -t4 * t467 / 9
    t470 = cos(t442)
    VG[153] = -t446 * t470
    t472 = cos(t447)
    VG[154] = -t450 * t472
    t474 = cos(t451)
    VG[155] = -t454 * t474
    t476 = cos(t455)
    VG[156] = -t458 * t476
    t478 = cos(t459)
    VG[157] = -t462 * t478
    t480 = cos(t463)
    VG[158] = -t466 * t480
    t482 = npsig ^ 2
    t483 = 1 / t482
    VGG[1] = -t229 * t483
    t485 = npsil ^ 2
    t486 = 1 / t485
    VGG[2] = -t231 * t486
    VGG[3] = -t233 * t483 / 4
    VGG[4] = -t236 * t483 / 9
    VGG[5] = -t239 * t483 / 16
    VGG[6] = -t242 * t486 / 4
    VGG[7] = -t245 * t486 / 9
    VGG[8] = -t248 * t486 / 16
    VGG[9] = -t251 * t486 / 25
    VGG[10] = -t254 * t486 / 36
    t504 = t32 ^ 2
    t505 = 1 / t504
    VGG[11] = -t257 * t505
    t507 = t36 ^ 2
    t508 = 1 / t507
    VGG[12] = -t259 * t508
    t510 = t41 ^ 2
    t511 = 1 / t510
    VGG[13] = -t261 * t511
    t513 = t46 ^ 2
    t514 = 1 / t513
    VGG[14] = -t263 * t514
    t516 = t52 ^ 2
    t517 = 1 / t516
    VGG[15] = -t265 * t517
    t519 = t56 ^ 2
    t520 = 1 / t519
    VGG[16] = -t267 * t520
    t522 = t60 ^ 2
    t523 = 1 / t522
    VGG[17] = -t269 * t523
    t525 = t66 ^ 2
    t526 = 1 / t525
    VGG[18] = -t271 * t526
    t528 = t70 ^ 2
    t529 = 1 / t528
    VGG[19] = -t273 * t529
    t531 = t74 ^ 2
    t532 = 1 / t531
    VGG[20] = -t275 * t532
    t534 = t80 ^ 2
    t535 = 1 / t534
    VGG[21] = -t277 * t535
    t537 = t84 ^ 2
    t538 = 1 / t537
    VGG[22] = -t279 * t538
    t540 = t89 ^ 2
    t541 = 1 / t540
    VGG[23] = -t281 * t541
    t543 = t93 ^ 2
    t544 = 1 / t543
    VGG[24] = -t283 * t544
    t546 = t97 ^ 2
    t547 = 1 / t546
    VGG[25] = -t285 * t547
    t549 = t101 ^ 2
    t550 = 1 / t549
    VGG[26] = -t287 * t550
    t552 = t105 ^ 2
    t553 = 1 / t552
    VGG[27] = -t289 * t553
    t555 = t109 ^ 2
    t556 = 1 / t555
    VGG[28] = -t291 * t556
    t558 = t114 ^ 2
    t559 = 1 / t558
    VGG[29] = -t293 * t559
    t561 = t119 ^ 2
    t562 = 1 / t561
    VGG[30] = -t295 * t562
    t564 = t123 ^ 2
    t565 = 1 / t564
    VGG[31] = -t297 * t565
    t567 = t127 ^ 2
    t568 = 1 / t567
    VGG[32] = -t299 * t568
    t570 = t131 ^ 2
    t571 = 1 / t570
    VGG[33] = -t301 * t571
    t573 = t135 ^ 2
    t574 = 1 / t573
    VGG[34] = -t303 * t574
    t576 = t139 ^ 2
    t577 = 1 / t576
    VGG[35] = -t305 * t577
    t579 = t143 ^ 2
    t580 = 1 / t579
    VGG[36] = -t307 * t580
    t582 = t147 ^ 2
    t583 = 1 / t582
    VGG[37] = -t309 * t583
    t585 = t151 ^ 2
    t586 = 1 / t585
    VGG[38] = -t311 * t586
    t588 = t155 ^ 2
    t589 = 1 / t588
    VGG[39] = -t313 * t589
    t591 = t51 ^ 2
    t592 = 1 / t591
    VGG[40] = -t315 * t592
    t594 = t161 ^ 2
    t595 = 1 / t594
    VGG[41] = -t317 * t595
    t597 = t165 ^ 2
    t598 = 1 / t597
    VGG[42] = -t319 * t598
    t600 = t65 ^ 2
    t601 = 1 / t600
    VGG[43] = -t321 * t601
    t603 = t171 ^ 2
    t604 = 1 / t603
    VGG[44] = -t323 * t604
    t606 = t175 ^ 2
    t607 = 1 / t606
    VGG[45] = -t325 * t607
    t609 = t179 ^ 2
    t610 = 1 / t609
    VGG[46] = -t327 * t610
    t612 = t183 ^ 2
    t613 = 1 / t612
    VGG[47] = -t329 * t613
    t615 = t187 ^ 2
    t616 = 1 / t615
    VGG[48] = -t331 * t616
    t618 = t191 ^ 2
    t619 = 1 / t618
    VGG[49] = -t333 * t619
    t621 = t195 ^ 2
    t622 = 1 / t621
    VGG[50] = -t335 * t622
    t624 = t199 ^ 2
    t625 = 1 / t624
    VGG[51] = -t337 * t625
    t627 = t203 ^ 2
    t628 = 1 / t627
    VGG[52] = -t339 * t628
    t630 = t207 ^ 2
    t631 = 1 / t630
    VGG[53] = -t341 * t631
    t633 = t211 ^ 2
    t634 = 1 / t633
    VGG[54] = -t343 * t634
    t636 = t215 ^ 2
    t637 = 1 / t636
    VGG[55] = -t345 * t637
    t639 = t219 ^ 2
    t640 = 1 / t639
    VGG[56] = -t347 * t640
    t642 = t223 ^ 2
    t643 = 1 / t642
    VGG[57] = -t349 * t643
    t645 = t227 ^ 2
    t646 = 1 / t645
    VGG[58] = -t351 * t646
    VGG[59] = -t1 * t483
    VGG[60] = -t3 * t486
    VGG[61] = -t6 * t483 / 4
    VGG[62] = -t9 * t483 / 9
    VGG[63] = -t12 * t483 / 16
    VGG[64] = -t15 * t486 / 4
    VGG[65] = -t18 * t486 / 9
    VGG[66] = -t21 * t486 / 16
    VGG[67] = -t24 * t486 / 25
    VGG[68] = -t27 * t486 / 36
    VGG[69] = -t30 * t505
    VGG[70] = -t35 * t508
    VGG[71] = -t39 * t511
    VGG[72] = -t44 * t514
    VGG[73] = -t50 * t517
    VGG[74] = -t55 * t520
    VGG[75] = -t59 * t523
    VGG[76] = -t64 * t526
    VGG[77] = -t69 * t529
    VGG[78] = -t73 * t532
    VGG[79] = -t77 * t535
    VGG[80] = -t83 * t538
    VGG[81] = -t87 * t541
    VGG[82] = -t92 * t544
    VGG[83] = -t96 * t547
    VGG[84] = -t100 * t550
    VGG[85] = -t104 * t553
    VGG[86] = -t108 * t556
    VGG[87] = -t112 * t559
    VGG[88] = -t117 * t562
    VGG[89] = -t122 * t565
    VGG[90] = -t126 * t568
    VGG[91] = -t130 * t571
    VGG[92] = -t134 * t574
    VGG[93] = -t138 * t577
    VGG[94] = -t142 * t580
    VGG[95] = -t146 * t583
    VGG[96] = -t150 * t586
    VGG[97] = -t154 * t589
    VGG[98] = -t157 * t592
    VGG[99] = -t160 * t595
    VGG[100] = -t164 * t598
    VGG[101] = -t167 * t601
    VGG[102] = -t170 * t604
    VGG[103] = -t174 * t607
    VGG[104] = -t178 * t610
    VGG[105] = -t182 * t613
    VGG[106] = -t186 * t616
    VGG[107] = -t190 * t619
    VGG[108] = -t194 * t622
    VGG[109] = -t198 * t625
    VGG[110] = -t202 * t628
    VGG[111] = -t206 * t631
    VGG[112] = -t210 * t634
    VGG[113] = -t214 * t637
    VGG[114] = -t218 * t640
    VGG[115] = -t222 * t643
    VGG[116] = -t226 * t646
    VGG[117] = -t409 * t486 / 49
    VGG[118] = -t412 * t486 / 64
    t718 = t362 ^ 2
    t719 = 1 / t718
    VGG[119] = -t415 * t719
    t721 = t366 ^ 2
    t722 = 1 / t721
    VGG[120] = -t417 * t722
    t724 = t371 ^ 2
    t725 = 1 / t724
    VGG[121] = -t419 * t725
    t727 = t375 ^ 2
    t728 = 1 / t727
    VGG[122] = -t421 * t728
    t730 = t379 ^ 2
    t731 = 1 / t730
    VGG[123] = -t423 * t731
    t733 = t383 ^ 2
    t734 = 1 / t733
    VGG[124] = -t425 * t734
    t736 = t387 ^ 2
    t737 = 1 / t736
    VGG[125] = -t427 * t737
    t739 = t391 ^ 2
    t740 = 1 / t739
    VGG[126] = -t429 * t740
    t742 = t395 ^ 2
    t743 = 1 / t742
    VGG[127] = -t431 * t743
    t745 = t399 ^ 2
    t746 = 1 / t745
    VGG[128] = -t433 * t746
    t748 = t403 ^ 2
    t749 = 1 / t748
    VGG[129] = -t435 * t749
    t751 = t407 ^ 2
    t752 = 1 / t751
    VGG[130] = -t437 * t752
    VGG[131] = -t354 * t486 / 49
    VGG[132] = -t357 * t486 / 64
    VGG[133] = -t360 * t719
    VGG[134] = -t365 * t722
    VGG[135] = -t369 * t725
    VGG[136] = -t374 * t728
    VGG[137] = -t378 * t731
    VGG[138] = -t382 * t734
    VGG[139] = -t386 * t737
    VGG[140] = -t390 * t740
    VGG[141] = -t394 * t743
    VGG[142] = -t398 * t746
    VGG[143] = -t402 * t749
    VGG[144] = -t406 * t752
    VGG[145] = -t467 * t486 / 81
    t772 = t445 ^ 2
    t773 = 1 / t772
    VGG[146] = -t470 * t773
    t775 = t449 ^ 2
    t776 = 1 / t775
    VGG[147] = -t472 * t776
    t778 = t453 ^ 2
    t779 = 1 / t778
    VGG[148] = -t474 * t779
    t781 = t457 ^ 2
    t782 = 1 / t781
    VGG[149] = -t476 * t782
    t784 = t461 ^ 2
    t785 = 1 / t784
    VGG[150] = -t478 * t785
    t787 = t465 ^ 2
    t788 = 1 / t787
    VGG[151] = -t480 * t788
    VGG[152] = -t440 * t486 / 81
    VGG[153] = -t443 * t773
    VGG[154] = -t448 * t776
    VGG[155] = -t452 * t779
    VGG[156] = -t456 * t782
    VGG[157] = -t460 * t785
    VGG[158] = -t464 * t788
    
    return VG,VGG;
    
end

################# DRAG

function getWDRAG2procedure(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,typeofvariable)
    
    casesmk0 = false;
    if skm==1.0 
        casesmk0 = true
    end

    sma = equi[1]
    P1  = equi[2]
    P2  = equi[3]
    Q1  = equi[4]
    Q2  = equi[5]
    TL  = equi[6]

    # vector of constants
    A = satellite["MomentsOfInertia"][1]; C = satellite["MomentsOfInertia"][3];
    VVT1 = satellite["constantstermsdragtorque"][1]
    VVT2 = satellite["constantstermsdragtorque"][2]
    
    VVdragT = zeros(90)
    VVdragT[1:30]  = VVT1
    VVdragT[31:90] = VVT2    
    VVdragT = transpose(VVdragT)

    # sadov
    if size(sadovvar)[1] == 5
        Jg = sadovvar[1];
        Jh = sadovvar[2];
        psil = sadovvar[3];
        psig = sadovvar[4];
        psih = sadovvar[5];
    else
        Jg = sadovvar[2];
        Jh = sadovvar[3];
        psil = sadovvar[4];
        psig = sadovvar[5];
        psih = sadovvar[6];
    end    
    

    # # uVectH --> satellite to Sun unit vector approximated as Earth to Sun unit vector in the ref frame with angular momentum as z axis
    # Ri2o = equi2rotmatirth(Q1,Q2,sin(TL),cos(TL))
    # Ro2i = transpose(Ri2o)
    # phiL = 1+P1*sin(TL)+P2*cos(TL);
    # eta = sqrt(1-P1^2-P2^2)
    # rc   = sma*eta^2/phiL;
    # VV0  = (sqrt(muPlanet/sma)*1/eta*Ro2i*[P2*sin(TL)-P1*cos(TL), phiL,0] - nE*LinearAlgebra.cross([cpa1,cpa2,cpa3],rc*Ro2i*[1,0,0]))*1000.0;
    # cpsih = cos(psih); 
    # spsih = sin(psih);
    # cdelta = Jh/Jg;
    # sdelta = sqrt(1-cdelta^2);
    # Rhdelta = zeros(3,3);
    # Rhdelta[1,1] = cpsih;
    # Rhdelta[1,2] = spsih;
    # Rhdelta[1,3] = 0.0;
    # Rhdelta[2,1] = -cdelta*spsih;
    # Rhdelta[2,2] = cdelta*cpsih;
    # Rhdelta[2,3] = sdelta;
    # Rhdelta[3,1] = sdelta*spsih;
    # Rhdelta[3,2] = -sdelta*cpsih;
    # Rhdelta[3,3] = cdelta;
    # VV0H = Rhdelta*VV0;
    # V0 = norm(VV0H);
    # voh = zeros(20);
    # count = 0;
    # for ii=0:3
    #     for jj=0:3
    #         for kk=0:3
    #             if ii+jj+kk>3 
    #                 break;
    #             end
    #             count=count+1;
    #             voh[count] = VV0H[1]^ii*VV0H[2]^jj*VV0H[3]^kk;
    #             if ii+jj+kk == 1
    #                 voh[count] = voh[count]*V0;
    #                 println("1 ",count)
    #             end
    #             if ii+jj+kk == 3
    #                 voh[count] = voh[count]/V0;
    #                 println("3 ",count)
    #             end
    #         end
    #     end
    # end
    # # density
    # eta = sqrt(1-P1^2-P2^2)
    # GG  = 1.0+Q1^2+Q2^2;
    # q11 = (1.0-Q1^2+Q2^2)/GG;
    # q12 = (2.0*Q1*Q2)/GG;
    # q13 = (-2.0*Q1)/GG;
    # q21 = (2.0*Q1*Q2)/GG;
    # q22 = (1.0+Q1^2-Q2^2)/GG;
    # q23 = (2.0*Q2)/GG;
    # q31 = (2.0*Q1)/GG;
    # q32 = (-2.0*Q2)/GG;
    # q33 = (1.0-Q1^2-Q2^2)/GG;
    # RQ2i = [q11 q21 q31;
    #     q12 q22 q32;
    #     q13 q23 q33];
    # rc  = sma*(1.0-P1^2-P2^2)/(1.0 + P1 * sin(TL) + P2 * cos(TL));
    # R_flat = rPlanet * planet_flat * (1 -  oblateness* (RQ2i[3,1]*cos(TL)+RQ2i[3,2]*sin(TL))^2.0) +  rPlanet * (1 - planet_flat);
    # density = exponential_atm_model(rc-R_flat,atmosphericinfo["AtmM"]);

    eta = sqrt(1-P1^2-P2^2)
    GG  = 1.0+Q1^2+Q2^2;
    q11 = (1.0-Q1^2+Q2^2)/GG;
    q12 = (2.0*Q1*Q2)/GG;
    q13 = (-2.0*Q1)/GG;
    q21 = (2.0*Q1*Q2)/GG;
    q22 = (1.0+Q1^2-Q2^2)/GG;
    q23 = (2.0*Q2)/GG;
    q31 = (2.0*Q1)/GG;
    q32 = (-2.0*Q2)/GG;
    q33 = (1.0-Q1^2-Q2^2)/GG;
    RQ2i = [q11 q21 q31;
        q12 q22 q32;
        q13 q23 q33];
    MAterms = getMAterms_drag([TL],sma,P1,P2,eta,RQ2i,truelong2meanlong(TL,P1,P2),mjd2000,cpa1,cpa2,cpa3,rPlanet,muPlanet,nE,0.0,planet_flat,oblateness,atmosphericinfo,false,inertial2equatorialRM);
    voh,vO = organizeMAterms_dragforce([TL],MAterms,P1,P2,eta,true,false,false,atmosphericinfo["considerthermalCD"],0);
    voh = voh[1,:];
    
    cpsih = cos(psih); 
    spsih = sin(psih);
    cdelta = Jh/Jg;
    sdelta = sqrt(1-cdelta^2);
    voh[1:20] = vmt2vmtH_order03(voh[1:20],cdelta,sdelta,cpsih,spsih)
    voh[21:40] = vmt2vmtH_order03(voh[21:40],cdelta,sdelta,cpsih,spsih)
   

    if !casesmk0
        CMXb11_P1 = getCMXb11(WTb,voh[1:20],2);
        CMXb11_P2 = getCMXb11(WTb,voh[21:40],2);
        CMYb21_P1 = getCMYb21(WTb,voh[1:20],2);
        CMYb21_P2 = getCMYb21(WTb,voh[21:40],2);
        CMZb31_P1 = getCMZb31(WTb,voh[1:20],2);
        CMZb31_P2 = getCMZb31(WTb,voh[21:40],2);
        CMXb12_P1 = getCMXb12(WTb,voh[1:20],2);
        CMXb12_P2 = getCMXb12(WTb,voh[21:40],2);  
        CMYb22_P1 = getCMYb22(WTb,voh[1:20],2);
        CMYb22_P2 = getCMYb22(WTb,voh[21:40],2);
        CMZb32_P1 = getCMZb32(WTb,voh[1:20],2);
        CMZb32_P2 = getCMZb32(WTb,voh[21:40],2);
        CMXb13_P1 = getCMXb13(WTb,voh[1:20],2);
        CMXb13_P2 = getCMXb13(WTb,voh[21:40],2);
        CMYb23_P1 = getCMYb23(WTb,voh[1:20],2);
        CMYb23_P2 = getCMYb23(WTb,voh[21:40],2);        
        CMZb33_P1 = getCMZb33(WTb,voh[1:20],2);
        CMZb33_P2 = getCMZb33(WTb,voh[21:40],2);
        CMXSx_P1  = getCMXSx(WTb,WTS,voh[1:20],znC1,znC2,znC3,2);
        CMXSx_P2  = getCMXSx(WTb,WTS,voh[21:40],znC1,znC2,znC3,2);
        CMYSy_P1  = getCMYSy(WTb,WTS,voh[1:20],znC1,znC2,znC3,2);
        CMYSy_P2  = getCMYSy(WTb,WTS,voh[21:40],znC1,znC2,znC3,2);
        CMZSz_P1  = getCMZSz(WTb,WTS,voh[1:20],znC1,znC2,znC3,2);
        CMZSz_P2  = getCMZSz(WTb,WTS,voh[21:40],znC1,znC2,znC3,2);

        CMXb11 = zeros(size(CMXb11_P1));
        CMXb11[1:30,:]  = CMXb11_P1[1:30,:];
        CMXb11[61:90,:] = CMXb11_P2[61:90,:];
        CMYb21 = zeros(size(CMYb21_P1));
        CMYb21[1:30,:]  = CMYb21_P1[1:30,:];
        CMYb21[61:90,:] = CMYb21_P2[61:90,:];
        CMZb31 = zeros(size(CMZb31_P1));
        CMZb31[1:30,:]  = CMZb31_P1[1:30,:];
        CMZb31[61:90,:] = CMZb31_P2[61:90,:];
        CMXb12 = zeros(size(CMXb12_P1));
        CMXb12[1:30,:]  = CMXb12_P1[1:30,:];
        CMXb12[61:90,:] = CMXb12_P2[61:90,:];
        CMYb22 = zeros(size(CMYb22_P1));
        CMYb22[1:30,:]  = CMYb22_P1[1:30,:];
        CMYb22[61:90,:] = CMYb22_P2[61:90,:];
        CMZb32 = zeros(size(CMZb32_P1));
        CMZb32[1:30,:]  = CMZb32_P1[1:30,:];
        CMZb32[61:90,:] = CMZb32_P2[61:90,:];
        CMXb13 = zeros(size(CMXb13_P1));
        CMXb13[1:30,:]  = CMXb13_P1[1:30,:];
        CMXb13[61:90,:] = CMXb13_P2[61:90,:];
        CMYb23 = zeros(size(CMYb23_P1));
        CMYb23[1:30,:]  = CMYb23_P1[1:30,:];
        CMYb23[61:90,:] = CMYb23_P2[61:90,:];
        CMZb33 = zeros(size(CMZb33_P1));
        CMZb33[1:30,:]  = CMZb33_P1[1:30,:];
        CMZb33[61:90,:] = CMZb33_P2[61:90,:];
        CMXSx = zeros(size(CMXSx_P1))
        CMXSx[1:30,:]  = CMXSx_P1[1:30,:];
        CMXSx[61:90,:] = CMXSx_P2[61:90,:];
        CMYSy = zeros(size(CMYSy_P1));
        CMYSy[1:30,:]  = CMYSy_P1[1:30,:];
        CMYSy[61:90,:] = CMYSy_P2[61:90,:];
        CMZSz = zeros(size(CMZSz_P1));
        CMZSz[1:30,:]  = CMZSz_P1[1:30,:];
        CMZSz[61:90,:] = CMZSz_P2[61:90,:];

    else
        CMXb11_P1,CMXb12_P1,CMXb13_P1,CMYb21_P1,CMYb22_P1,CMYb23_P1,CMZb31_P1,CMZb32_P1,CMZb33_P1 = CMterms_smk0(voh[1:20],2);
        CMXb11_P2,CMXb12_P2,CMXb13_P2,CMYb21_P2,CMYb22_P2,CMYb23_P2,CMZb31_P2,CMZb32_P2,CMZb33_P2 = CMterms_smk0(voh[21:40],2);

        CMXb11 = zeros(size(CMXb11_P1));
        CMXb11[1:30,:]  = CMXb11_P1[1:30,:];
        CMXb11[61:90,:] = CMXb11_P2[61:90,:];
        CMYb21 = zeros(size(CMYb21_P1));
        CMYb21[1:30,:]  = CMYb21_P1[1:30,:];
        CMYb21[61:90,:] = CMYb21_P2[61:90,:];
        CMZb31 = zeros(size(CMZb31_P1));
        CMZb31[1:30,:]  = CMZb31_P1[1:30,:];
        CMZb31[61:90,:] = CMZb31_P2[61:90,:];
        CMXb12 = zeros(size(CMXb12_P1));
        CMXb12[1:30,:]  = CMXb12_P1[1:30,:];
        CMXb12[61:90,:] = CMXb12_P2[61:90,:];
        CMYb22 = zeros(size(CMYb22_P1));
        CMYb22[1:30,:]  = CMYb22_P1[1:30,:];
        CMYb22[61:90,:] = CMYb22_P2[61:90,:];
        CMZb32 = zeros(size(CMZb32_P1));
        CMZb32[1:30,:]  = CMZb32_P1[1:30,:];
        CMZb32[61:90,:] = CMZb32_P2[61:90,:];
        CMXb13 = zeros(size(CMXb13_P1));
        CMXb13[1:30,:]  = CMXb13_P1[1:30,:];
        CMXb13[61:90,:] = CMXb13_P2[61:90,:];
        CMYb23 = zeros(size(CMYb23_P1));
        CMYb23[1:30,:]  = CMYb23_P1[1:30,:];
        CMYb23[61:90,:] = CMYb23_P2[61:90,:];
        CMZb33 = zeros(size(CMZb33_P1));
        CMZb33[1:30,:]  = CMZb33_P1[1:30,:];
        CMZb33[61:90,:] = CMZb33_P2[61:90,:];

    end

    # VGT and VGGT ---> vectors containing fast variables
    npsig = Jg/A/C*((C-A)PP+A*KK)/KK;
    npsil = -pi/2.0/sqrt(1+k)*sqrt(skm)*Jg/A/C*(C-A)/KK;
    if !casesmk0
        VGT,VGGT = getgeneratinfunctiontermsW2_srp_drag_psil_psig(psil,psig,npsil,npsig);
    else
        VGT,VGGT =  getgeneratinfunctiontermsW2_allpert_psil_psig_smk0(psil,psig,npsil,npsig)
    end

    # terms of the generating function
    W21Block  = -2*skm/Jg*CMXb13 - 2*skm/Jg*(1-m)/(1+k)*CMYb23 + 2*(1-skm)/Jg*CMZb33;
    Tb3Block   = CMXb13+CMYb23+CMZb33;
    Tb2Block   = CMXb12+CMYb22+CMZb32;
    Tb1Block   = CMXb11+CMYb21+CMZb31;
    
    # gen fun
    dnpsildskm = pi*Jg*(A-C)*(EE-(1-m)*skm*KK)/(4.0*A*C*(1-skm)*sqrt(skm)*sqrt(1+k)-KK^2.0);
    dnpsildJg  = pi*sqrt(skm)*(A-C)/(2.0*A*C*sqrt(1+k)*KK);
    dnpsigdJg  = ((C-A)*PP+A*KK)/(A*C*KK);

    if casesmk0
        dnpsigdskm = -Jg*(A-C)*(-k/sqrt(1+k)/2.0+1.0)/(2*A*C);
    else
        dnpsigdskm = -Jg*(A-C)*((EE-(1-m)*skm*KK)*PP-(1-skm)*KK*EE)/(2*A*C*(1-m)*skm*(1-skm)*KK^2);
    end

    if !casesmk0
        TSBlock    = CMXSx-(1-m)*CMYSy+CMZSz;

        if typeofvariable==1
            W21 = VVdragT*W21Block*VGT;
            W22 = VVdragT*Tb3Block*VGT;
            W23 = VVdragT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = VVdragT*(pi/(2*(m-1)*KK*Jg)*TSBlock*VGT + dnpsildskm*W21Block*VGGT + dnpsildJg*Tb3Block*VGGT);
            W25 = VVdragT*((PP-(1-skm)*KK)*sqrt(1+k)/(Jg*KK*sqrt(skm)*(1-m))*TSBlock*VGT - cdelta/sdelta/Jg*Tb1Block*VGT + dnpsigdskm*W21Block*VGGT + dnpsigdJg*Tb3Block*VGGT);
            W26 = VVdragT*1/sdelta/Jg*Tb1Block*VGT;  

            W2  = append!([W21],[W22],[W23],[W24],[W25],[W26]);

        else
            W21 = VVdragT*W21Block*VGT;
            W22 = VVdragT*Tb3Block*VGT;
            W23 = VVdragT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = VVdragT*(pi/(2*(m-1)*KK*Jg)*TSBlock*VGT + dnpsildskm*W21Block*VGGT + dnpsildJg*Tb3Block*VGGT);
            W25 = VVdragT*((PP-(1-skm)*KK)*sqrt(1+k)/(Jg*KK*sqrt(skm)*(1-m))*TSBlock*VGT + psih*sdelta/Jg*Tb2Block*VGT + dnpsigdskm*W21Block*VGGT + dnpsigdJg*Tb3Block*VGGT);
            W26 = VVdragT*(-spsih*Tb1Block+sdelta*cpsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
            W27 = VVdragT*(cpsih*Tb1Block+sdelta*spsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
        
            W2  = append!([W21],[W22],[W23],[W24],[W25],[W26],[W27]);

        end
    else
        if typeofvariable==1

            W21 = VVdragT*W21Block*VGT;
            W22 = VVdragT*Tb3Block*VGT;
            W23 = VVdragT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = 0.0;
            W25 = VVdragT*(- cdelta/sdelta/Jg*Tb1Block*VGT + (dnpsildskm+dnpsigdskm)*W21Block*VGGT + (dnpsildJg+dnpsigdJg)*Tb3Block*VGGT);
            W26 = VVdragT*1/sdelta/Jg*Tb1Block*VGT;  
            
            W2  = [W21,W22,W23,W24,W25,W26];
        else
            W21 = VVdragT*W21Block*VGT;
            W22 = VVdragT*Tb3Block*VGT;
            W23 = VVdragT*(cdelta*Tb3Block + sdelta*Tb2Block)*VGT;
            W24 = [0.0];
            W25 = VVdragT*(psih*sdelta/Jg*Tb2Block*VGT + (dnpsildskm+dnpsigdskm)*W21Block*VGGT + (dnpsildJg+dnpsigdJg)*Tb3Block*VGGT);
            W26 = VVdragT*(-spsih*Tb1Block+sdelta*cpsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
            W27 = VVdragT*(cpsih*Tb1Block+sdelta*spsih*Tb3Block-cdelta*spsih*Tb2Block)*VGT  
           
            W2  = [W21,W22,W23,W24,W25,W26,W27];

        end
    end

    if atmosphericinfo["atmosphericmodel"] == 2 && !casesmk0
        W2_corr =  getWDRAG2procedure_numerical(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,typeofvariable,false);
        W2 = W2 + W2_corr;
    end
    
    return W2;
end

function getWDRAG2procedure_numerical(k,sadovvar,skm,m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi,satellite,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,typeofvariable,isfull)
    
    casesmk0 = false;
    if skm==1.0 
        casesmk0 = true
        Base.error("method still not suitable in the case m=0")
    end

    sma = equi[1]
    P1  = equi[2]
    P2  = equi[3]
    Q1  = equi[4]
    Q2  = equi[5]
    TL  = equi[6]

    # vector of constants
    A = satellite["MomentsOfInertia"][1]; C = satellite["MomentsOfInertia"][3];
    VVT1 = satellite["constantstermsdragtorque"][1]
    VVT2 = satellite["constantstermsdragtorque"][2]
    VVT3 = satellite["constantstermsdragtorque"][3]
    VVT4 = satellite["constantstermsdragtorque"][4]
    
    # sadov
    if size(sadovvar)[1] == 5
        Jg = sadovvar[1];
        Jh = sadovvar[2];
        psil = sadovvar[3];
        psig = sadovvar[4];
        psih = sadovvar[5];
    else
        Jg = sadovvar[2];
        Jh = sadovvar[3];
        psil = sadovvar[4];
        psig = sadovvar[5];
        psih = sadovvar[6];
    end    
    cpsih = cos(psih); 
    spsih = sin(psih);
    cdelta = Jh/Jg;
    sdelta = sqrt(1-cdelta^2);

    GG = 1.0+Q1^2+Q2^2;
    q11 = (1.0-Q1^2+Q2^2)/GG;
    q12 = (2.0*Q1*Q2)/GG;
    q13 = (-2.0*Q1)/GG;
    q21 = (2.0*Q1*Q2)/GG;
    q22 = (1.0+Q1^2-Q2^2)/GG;
    q23 = (2.0*Q2)/GG;
    q31 = (2.0*Q1)/GG;
    q32 = (-2.0*Q2)/GG;
    q33 = (1.0-Q1^2-Q2^2)/GG;

    RQ2i = [q11 q21 q31;
            q12 q22 q32;
            q13 q23 q33];

    eta = sqrt(1.0-P1^2-P2^2)
    
    MAterms = getMAterms_drag([TL],sma,P1,P2,eta,RQ2i,truelong2meanlong(TL,P1,P2),mjd2000,cpa1,cpa2,cpa3,rPlanet,muPlanet,nE,0.0,planet_flat,oblateness,atmosphericinfo,false,inertial2equatorialRM);
    vMAt,vO = organizeMAterms_dragforce([TL],MAterms,P1,P2,eta,true,false,false,atmosphericinfo["considerthermalCD"],0);
    vMAt = vMAt[1,:];

    # mean motion torque free problem and derivatives
    npsig = Jg/A/C*((C-A)PP+A*KK)/KK;
    npsil = -pi/2.0/sqrt(1+k)*sqrt(skm)*Jg/A/C*(C-A)/KK;
    dnpsildskm = pi*Jg*(A-C)*(EE-(1-m)*skm*KK)/(4.0*A*C*(1-skm)*sqrt(skm)*sqrt(1+k)-KK^2.0);
    dnpsildJg  = pi*sqrt(skm)*(A-C)/(2.0*A*C*sqrt(1+k)*KK);
    if casesmk0
        dnpsigdskm = -Jg*(A-C)*(-k/sqrt(1+k)/2.0+1.0)/(2*A*C);
    else
        dnpsigdskm = -Jg*(A-C)*((EE-(1-m)*skm*KK)*PP-(1-skm)*KK*EE)/(2*A*C*(1-m)*skm*(1-skm)*KK^2);
    end
    dnpsigdJg  = ((C-A)*PP+A*KK)/(A*C*KK);
     
    
    # constants
    k1r = sqrt(1+k);
    smk = 1.0-skm;
    mM1KK = KK*(m-1.0);
    if m==1.0
        T4 = -atan(sqrt(k))/sqrt(k);
    else
        if k == 0
            T4 = -pi/2.0;
        else
            T4 = (KK*m - PP*(m+k))/k;
        end
    end
    # npsil*diff(W_i,psil)+npsig*diff(W_i,psig)=f(psil,psig)-z
    # W_i = W_i_f + W_i_z;  
    # npsil*diff(W_i_z,psil)+npsig*diff(W_i_z,psig)=-z
    # npsil*diff(W_i_f,psil)+npsig*diff(W_i_f,psig)=f(psil,psig)
        
    ######
    # W_i_z
    if isfull
        ZV =  getaveragedvectorfield_drag(k,k1r,skm,smk,m,KK,EE,PP,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,psih,vMAt,VVT1,VVT2,VVT3,VVT4,typeofvariable);
    else
        ZV =  getpartialaveragedvectorfield_drag(k,k1r,skm,smk,m,KK,EE,PP,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,psih,vMAt,VVT2,VVT3,VVT4,typeofvariable);
    end
    atanpsig = atan(tan(psig/2.0));
    atanpsil = atan(tan(psil/2.0));
    W2_Z = - ZV*(atanpsil/npsil+atanpsig/npsig);
    W2_Z[4] = W2_Z[4] - (dnpsildskm*ZV[1]+dnpsildJg*ZV[2])*((atanpsil^2.0-pi^2.0/12)/(npsil^2) + (atanpsig^2.0-pi^2.0/12)/(npsig^2));
    W2_Z[5] = W2_Z[5] - (dnpsigdskm*ZV[1]+dnpsigdJg*ZV[2])*((atanpsil^2.0-pi^2.0/12)/(npsil^2) + (atanpsig^2.0-pi^2.0/12)/(npsig^2));
    # println("qui ok")
    if isfull
        f2V,psilV,zetaV,ratio  = getnonaveragedvectorfield_drag(k,k1r,skm,smk,m,KK,EE,PP,T4,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,psih,vMAt,VVT1,VVT2,VVT3,VVT4,typeofvariable,npsil,npsig,2.0*pi/180.0,2.0*pi/180.0);
    else
        f2V,psilV,zetaV,ratio  = getpartialnonaveragedvectorfield_drag(k,k1r,skm,smk,m,KK,EE,PP,T4,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,psih,vMAt,VVT2,VVT3,VVT4,typeofvariable,npsil,npsig,2.0*pi/180.0,2.0*pi/180.0);
    end
    W2_f   = getgenfunterms_psilpsig(npsil,ratio,dnpsildskm,dnpsildJg,dnpsigdskm,dnpsigdJg,psig,psil,psilV,zetaV,f2V)
    
    # W2        
    W2 = W2_f + W2_Z;
    
    return W2;
end

####################################################################################################################################################### numerical

##################### W2 numerical srp
function getaveragedvectorfield_srp(k,k1r,skm,smk,m,KK,EE,PP,T4,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,vcu,VVT1,VVT2,VVT3,typeofvariable)
        
    if m==0.0 
        T1 = 1.0/2.0;
        T2 = -3.0/8.0;
        T3 = 1.0;
        T6 = 0.0;
    else
        T1 = ((m-1.0)*KK+EE)/m/KK;
        T2 = 2.0*(KK-EE)/KK/m^2.0 + (EE-2.0*KK)/m/KK;
        T3 = EE/KK;
        T6 = (KK^2.0*(m-1.0)+EE^2.0)/KK^2.0/k;
    end
    
    cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1 = srpanalyticalaveragedcoeff_T1(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
    cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2 = srpanalyticalaveragedcoeff_T2(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
    cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3 = srpanalyticalaveragedcoeff_T3(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
    
    Mxb11T1,Mxb12T1,Mxb13T1,MxpsilT1,Myb21T1,Myb22T1,Myb23T1,MypsilT1,Mzb31T1,Mzb32T1,Mzb33T1,MzpsilT1 = srp_averagedterms(VVT1,cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1);
    Mxb11T2,Mxb12T2,Mxb13T2,MxpsilT2,Myb21T2,Myb22T2,Myb23T2,MypsilT2,Mzb31T2,Mzb32T2,Mzb33T2,MzpsilT2 = srp_averagedterms(VVT2,cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2);
    Mxb11T3,Mxb12T3,Mxb13T3,MxpsilT3,Myb21T3,Myb22T3,Myb23T3,MypsilT3,Mzb31T3,Mzb32T3,Mzb33T3,MzpsilT3 = srp_averagedterms(VVT3,cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3);
    
    mxb11  = Mxb11T1 + Mxb11T2 + Mxb11T3;
    mxb12  = Mxb12T1 + Mxb12T2 + Mxb12T3;
    mxb13  = Mxb13T1 + Mxb13T2 + Mxb13T3;
    mxpsil = MxpsilT1 + MxpsilT2 + MxpsilT3;
    myb21  = Myb21T1 + Myb21T2 + Myb21T3;
    myb22  = Myb22T1 + Myb22T2 + Myb22T3;
    myb23  = Myb23T1 + Myb23T2 + Myb23T3;
    mypsil = MypsilT1 + MypsilT2 + MypsilT3;
    mzb31  = Mzb31T1 + Mzb31T2 + Mzb31T3;
    mzb32  = Mzb32T1 + Mzb32T2 + Mzb32T3;
    mzb33  = Mzb33T1 + Mzb33T2 + Mzb33T3;
    mzpsil = MzpsilT1 + MzpsilT2 + MzpsilT3;
    
    TBblock1 = mxb11+myb21+mzb31;
    TBblock2 = mxb12+myb22+mzb32;

    Z21 =  -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
    Z22 =  mxb13+myb23+mzb33;
    Z23 =  cdelta*Z22 + sdelta*TBblock2;
    Z24 =  pi/(2*mM1KK*Jg)*(mxpsil+(1-m)*mypsil+mzpsil);

    if typeofvariable==1
        Z26  = TBblock1/Jg/sdelta;
        Z25  = -cdelta*Z26 + 2/pi*sqrt(skm)*T4*k1r*Z24;
        ZV = [Z21,Z22,Z23,Z24,Z25,Z26];
    else
        Z25 = 2/pi*sqrt(skm)*T4*k1r*Z24 + psih*sdelta*TBblock2/Jg;
        Z26 = -spsih*TBblock1 + sdelta*cpsih*Z22 - cdelta*cpsih*TBblock2;
        Z27 =  cpsih*TBblock1 + sdelta*spsih*Z22 - cdelta*spsih*TBblock2;
        ZV = [Z21,Z22,Z23,Z24,Z25,Z26,Z27];
    end
    return ZV;
end

function getnonaveragedvectorfield_srp(k,k1r,skm,smk,m,KK,EE,PP,T4,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,vcu,VVT1,VVT2,VVT3,typeofvariable,npsil,npsig,psilstep,zetastep)
    ratio = npsig/npsil;
    ll = Int(floor(2*pi/psilstep)+1.0)
    if ratio >=0.0
        zetaMin = -ratio*2*pi;
        zetaMax = 2*pi;
    else
        zetaMin = 0.0;
        zetaMax = 2*pi*(1-ratio);
    end
    mm = Int(floor((zetaMax-zetaMin)/zetastep)+1.0);

    psilV = collect(range(0.0,2*pi,ll));
    zetaV =  collect(range(zetaMin,zetaMax,mm));

    nn = 220;
    oo = 12;
    if typeofvariable == 1
        pp = 6;
    else
        pp = 7;
    end
    
    MXT1CM, MYT1CM, MZT1CM = srpcoeffM_T1(vcu,cdelta,sdelta,cpsih,spsih)
    MXT2CM, MYT2CM, MZT2CM = srpcoeffM_T2(vcu,cdelta,sdelta,cpsih,spsih)
    MXT3CM, MYT3CM, MZT3CM = srpcoeffM_T3(vcu,cdelta,sdelta,cpsih,spsih)
    VVT1 = transpose(VVT1);
    VVT2 = transpose(VVT2);
    VVT3 = transpose(VVT3);
    MXCM = VVT1*MXT1CM + VVT2*MXT2CM + VVT3*MXT3CM
    MYCM = VVT1*MYT1CM + VVT2*MYT2CM + VVT3*MYT3CM
    MZCM = VVT1*MZT1CM + VVT2*MZT2CM + VVT3*MZT3CM
    skmR = sqrt(skm);
    smkR = sqrt(smk);
    
    funInit = zeros(ll,mm,pp)
    for ii = 1 : ll
        for jj = 1 : mm
            BV = srptermspsilzeta(psilV[ii],zetaV[jj],ratio,m,KK,PP,EE,k,k1r,skmR,smkR);
            
            TBblock1 = (MXCM*BV[1:nn] + MYCM*BV[4*nn+1:5*nn] + MZCM*BV[8*nn+1:9*nn])[1]
            TBblock2 = (MXCM*BV[nn+1:2*nn] + MYCM*BV[5*nn+1:6*nn] + MZCM*BV[9*nn+1:10*nn])[1]
            mxb13    = (MXCM*BV[2*nn+1:3*nn])[1]
            myb23    = (MYCM*BV[6*nn+1:7*nn])[1]
            mzb33    = (MZCM*BV[10*nn+1:11*nn])[1]
            TSblock  = (MXCM*BV[3*nn+1:4*nn] + (1.0-m)*MYCM*BV[7*nn+1:8*nn] + MZCM*BV[11*nn+1:12*nn])[1]

            funInit[ii,jj,1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
            funInit[ii,jj,2] =  mxb13+myb23+mzb33;
            funInit[ii,jj,3] = cdelta*funInit[ii,jj,2] + sdelta*TBblock2;
            funInit[ii,jj,4] =  pi/(2*mM1KK*Jg)*TSblock;
            
            if typeofvariable == 1
                funInit[ii,jj,6] = TBblock1/Jg/sdelta;
                funInit[ii,jj,5] = -cdelta*funInit[ii,jj,6] + 2/pi*skmR*T4*k1r*funInit[ii,jj,4]; 
            else
                funInit[ii,jj,5] = 2/pi*skmR*T4*k1r*funInit[ii,jj,4]+ psih*sdelta*TBblock2/Jg;
                funInit[ii,jj,6] = -spsih*TBblock1 + sdelta*cpsih*funInit[ii,jj,2] - cdelta*cpsih*TBblock2;
                funInit[ii,jj,7] =  cpsih*TBblock1 + sdelta*spsih*funInit[ii,jj,2] - cdelta*spsih*TBblock2;
            end
        end
    end
            
    return funInit,psilV,zetaV,ratio;
end

function srptermspsilzeta(psil,zeta,ratio,m,KK,PP,EE,k,k1r,skmR,smkR)

    psig = zeta+ratio*psil;
    out = srpBV(psig,psil,m,KK,EE,PP,k,k1r,skmR,smkR)
    return out;

    # u = 2*KK*psil/pi
    # eF  = 2.0*KK*mod(psil,2.0*pi)/pi;
    # lambda = Elliptic.Jacobi.am(eF,m);
    # sn = sin(lambda);
    # cn = cos(lambda);
    # dn = sqrt(1.0-m*sn+2.0);
    # eP  = ellP(-k,lambda,m); 
    # eE  = ellE(lambda,m);
    # zn = eE-eF*EE/KK
    # gA  = zeta + ratio*psil + k1r/skmR*(eF/KK*PP-eP);
    # cgA = cos(gA);
    # sgA = sin(gA);
    # dnkR  = sqrt(1+k*sn^2);

    # bij = zeros(9);
    # bij[1] = (-sgA*cn*dn*skmR-cgA*sn*k1r)/dnkR;
    # bij[2] = ( cgA*cn*dn*skmR-sgA*sn*k1r)/dnkR;
    # bij[3] = smkR*cn
    # bij[4] = ( sgA*sn*dn*skmR*k1r-cgA*cn)/dnkR;
    # bij[5] = (-cgA*sn*dn*skmR*k1r-sgA*cn)/dnkR;
    # bij[6] = -smkR*k1r*sn
    # bij[7] = smkR*sgA*dnkR
    # bij[8] = -smkR*cgA*dnkR
    # bij[9] = skmR*dn

    # BV = zeros(1980);
    # count = 0;
    # for ii =1:7
    #     for jj = ii+1:8
    #         for kk = jj+1:9
    #             count = count + 1;
    #             BV[count] = bij[ii]*bij[jj]*bij[kk];
    #         end
    #     end
    # end
    # for ii =1:9
    #     for kk = 1:9
    #         count = count + 1;
    #         BV[count] = bij[ii]*bij[ii]*bij[kk];
    #     end
    # end
    # for ii =1:9
    #     for jj = ii:9
    #         count = count + 1;
    #         BV[count]  = bij[ii]*bij[jj];
    #     end
    # end
    # for ii =1:9
    #     count = count + 1;
    #     BV[count]  = bij[ii];
    # end
    # count = count + 1;
    # BV[count]  = 1;

    # # BV = vec(BV);
    # if smkR == 0.0
    #     smkR = smkR+1e-15
    # end
    # MV = [bij[1], bij[2], bij[3],(dn*sn-cn*zn)/smkR, bij[4], bij[5], bij[6], (dn*cn+sn*zn)/smkR/k1r, bij[7],bij[8], bij[9], (dn*zn-m*sn*cn)/skmR];
        
    # return BV,MV;
end

##################### W2 numerical drag
function getaveragedvectorfield_drag(k,k1r,skm,smk,m,KK,EE,PP,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,psih,vMAt,VVT1,VVT2,VVT3,VVT4,typeofvariable)
    T1,T2,T3,T4,T6,T7,T8 = getsingularityterms(m,k,smk,skm,KK,EE,PP)

    skmBlock,TBblock1,TBblock2,TBblock3,TSblock = getdragsaaveragedtermsattitudeeq(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,T7,T8,cdelta,sdelta,cpsih,spsih,vMAt,VVT1,VVT2,VVT3,VVT4)

    Z21 =  skmBlock/Jg;
    Z22 =  TBblock3;
    Z23 =  cdelta*Z22 + sdelta*TBblock2;
    Z24 =  pi/(2*mM1KK*Jg)*TSblock;

    if typeofvariable==1
        Z26  = TBblock1/Jg/sdelta;
        Z25  = -cdelta*Z26 + 2/pi*sqrt(skm)*T4*k1r*Z24;
        ZV = [Z21,Z22,Z23,Z24,Z25,Z26];
    else
        Z25 = 2/pi*sqrt(skm)*T4*k1r*Z24 + psih*sdelta*TBblock2/Jg;
        Z26 = -spsih*TBblock1 + sdelta*cpsih*Z22 - cdelta*cpsih*TBblock2;
        Z27 =  cpsih*TBblock1 + sdelta*spsih*Z22 - cdelta*spsih*TBblock2;
        ZV = [Z21,Z22,Z23,Z24,Z25,Z26,Z27];
    end

    return ZV;
end

function getnonaveragedvectorfield_drag(k,k1r,skm,smk,m,KK,EE,PP,T4,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,psih,vMAt,VVT1,VVT2,VVT3,VVT4,typeofvariable,npsil,npsig,psilstep,zetastep)
    
    ratio = npsig/npsil;
    ll1 = Int(floor(2*pi/psilstep)+1.0)
    if ratio >=0.0
        zetaMin = -ratio*2*pi;
        zetaMax = 2*pi;
    else
        zetaMin = 0.0;
        zetaMax = 2*pi*(1-ratio);
    end
    ll2 = Int(floor((zetaMax-zetaMin)/(pi/180.0))+1.0);

    psilV = collect(range(0.0,2*pi,ll1));
    zetaV =  collect(range(zetaMin,zetaMax,ll2));
    
    vMAtH_VO2 = vmt2vmtH_order03(vMAt[1:20],cdelta,sdelta,cpsih,spsih)
    vMAtH_VO3 = vmt2vmtH_order4(vMAt[77:91],cdelta,sdelta,cpsih,spsih)
    vMAtH_VO4 = vmt2vmtH_order5(vMAt[41:61],cdelta,sdelta,cpsih,spsih);
    vMAtH_NO23 = vmt2vmtH_order03(vMAt[21:40],cdelta,sdelta,cpsih,spsih)
    vMAtH_NO4 = vmt2vmtH_order4(vMAt[62:76],cdelta,sdelta,cpsih,spsih)
    
    skmR = sqrt(skm);
    smkR = sqrt(smk);

    if typeofvariable == 1
        funInit = zeros(ll1,ll2,6)
    else
        funInit = zeros(ll1,ll2,7)
    end
    for ii = 1 : ll1
        for jj = 1 : ll2
            b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz = dragtermspsilzeta(psilV[ii],zetaV[jj],ratio,m,KK,PP,EE,k,k1r,skmR,smkR);
                        
            mm = dragtermsV_O2(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_VO2)
            nn = length(VVT1)
            mxb11    = LinearAlgebra.dot(VVT1,mm[1:nn]); 
            mxb12    = LinearAlgebra.dot(VVT1,mm[nn+1:2*nn]); 
            mxb13    = LinearAlgebra.dot(VVT1,mm[2*nn+1:3*nn]); 
            mxpsil   = LinearAlgebra.dot(VVT1,mm[3*nn+1:4*nn]); 
            myb21    = LinearAlgebra.dot(VVT1,mm[4*nn+1:5*nn]); 
            myb22    = LinearAlgebra.dot(VVT1,mm[5*nn+1:6*nn]); 
            myb23    = LinearAlgebra.dot(VVT1,mm[6*nn+1:7*nn]); 
            mypsil   = LinearAlgebra.dot(VVT1,mm[7*nn+1:8*nn]); 
            mzb31    = LinearAlgebra.dot(VVT1,mm[8*nn+1:9*nn]); 
            mzb32    = LinearAlgebra.dot(VVT1,mm[9*nn+1:10*nn]); 
            mzb33    = LinearAlgebra.dot(VVT1,mm[10*nn+1:11*nn]); 
            mzpsil   = LinearAlgebra.dot(VVT1,mm[11*nn+1:12*nn]); 

            mm = dragtermsV_O3(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_VO3)
            nn = length(VVT2)
            mxb11    = mxb11 + LinearAlgebra.dot(VVT2,mm[1:nn]); 
            mxb12    = mxb12 + LinearAlgebra.dot(VVT2,mm[nn+1:2*nn]); 
            mxb13    = mxb13 + LinearAlgebra.dot(VVT2,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT2,mm[3*nn+1:4*nn]); 
            myb21    = myb21 + LinearAlgebra.dot(VVT2,mm[4*nn+1:5*nn]); 
            myb22    = myb22 + LinearAlgebra.dot(VVT2,mm[5*nn+1:6*nn]); 
            myb23    = myb23 + LinearAlgebra.dot(VVT2,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT2,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31 + LinearAlgebra.dot(VVT2,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32 + LinearAlgebra.dot(VVT2,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33 + LinearAlgebra.dot(VVT2,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT2,mm[11*nn+1:12*nn]); 

            mm = dragtermsV_O4(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_VO4)
            nn = length(VVT3)
            mxb11    = mxb11 + LinearAlgebra.dot(VVT3,mm[1:nn]); 
            mxb12    = mxb12 + LinearAlgebra.dot(VVT3,mm[nn+1:2*nn]); 
            mxb13    = mxb13 + LinearAlgebra.dot(VVT3,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT3,mm[3*nn+1:4*nn]); 
            myb21    = myb21 + LinearAlgebra.dot(VVT3,mm[4*nn+1:5*nn]); 
            myb22    = myb22 + LinearAlgebra.dot(VVT3,mm[5*nn+1:6*nn]); 
            myb23    = myb23 + LinearAlgebra.dot(VVT3,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT3,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31 + LinearAlgebra.dot(VVT3,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32 + LinearAlgebra.dot(VVT3,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33 + LinearAlgebra.dot(VVT3,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT3,mm[11*nn+1:12*nn]); 

            mm = dragtermsN_O2(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_NO23)
            nn = length(VVT2)
            mxb11    = mxb11 + LinearAlgebra.dot(VVT2,mm[1:nn]); 
            mxb12    = mxb12 + LinearAlgebra.dot(VVT2,mm[nn+1:2*nn]); 
            mxb13    = mxb13 + LinearAlgebra.dot(VVT2,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT2,mm[3*nn+1:4*nn]); 
            myb21    = myb21 + LinearAlgebra.dot(VVT2,mm[4*nn+1:5*nn]); 
            myb22    = myb22 + LinearAlgebra.dot(VVT2,mm[5*nn+1:6*nn]); 
            myb23    = myb23 + LinearAlgebra.dot(VVT2,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT2,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31 + LinearAlgebra.dot(VVT2,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32 + LinearAlgebra.dot(VVT2,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33 + LinearAlgebra.dot(VVT2,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT2,mm[11*nn+1:12*nn]); 

            mm = dragtermsN_O3(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_NO23)
            nn = length(VVT3)
            mxb11    = mxb11 + LinearAlgebra.dot(VVT3,mm[1:nn]); 
            mxb12    = mxb12 + LinearAlgebra.dot(VVT3,mm[nn+1:2*nn]); 
            mxb13    = mxb13 + LinearAlgebra.dot(VVT3,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT3,mm[3*nn+1:4*nn]); 
            myb21    = myb21 + LinearAlgebra.dot(VVT3,mm[4*nn+1:5*nn]); 
            myb22    = myb22 + LinearAlgebra.dot(VVT3,mm[5*nn+1:6*nn]); 
            myb23    = myb23 + LinearAlgebra.dot(VVT3,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT3,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31 + LinearAlgebra.dot(VVT3,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32 + LinearAlgebra.dot(VVT3,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33 + LinearAlgebra.dot(VVT3,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT3,mm[11*nn+1:12*nn]); 

            mm = dragtermsN_O4(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_NO4)
            nn = length(VVT4)
            mxb11    = mxb11  + LinearAlgebra.dot(VVT4,mm[1:nn]); 
            mxb12    = mxb12  + LinearAlgebra.dot(VVT4,mm[nn+1:2*nn]); 
            mxb13    = mxb13  + LinearAlgebra.dot(VVT4,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT4,mm[3*nn+1:4*nn]); 
            myb21    = myb21  + LinearAlgebra.dot(VVT4,mm[4*nn+1:5*nn]); 
            myb22    = myb22  + LinearAlgebra.dot(VVT4,mm[5*nn+1:6*nn]); 
            myb23    = myb23  + LinearAlgebra.dot(VVT4,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT4,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31  + LinearAlgebra.dot(VVT4,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32  + LinearAlgebra.dot(VVT4,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33  + LinearAlgebra.dot(VVT4,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT4,mm[11*nn+1:12*nn]); 

            TBblock1 = mxb11 + myb21 + mzb31
            TBblock2 = mxb12 + myb22 + mzb32
            TSblock  = mxpsil + (1-m)*mypsil + mzpsil

            funInit[ii,jj,1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
            funInit[ii,jj,2] =  mxb13+myb23+mzb33;
            funInit[ii,jj,3] =  cdelta*funInit[ii,jj,2] + sdelta*TBblock2;
            funInit[ii,jj,4] =  pi/(2*mM1KK*Jg)*TSblock;
            
            if typeofvariable == 1
                funInit[ii,jj,6] = TBblock1/Jg/sdelta;
                funInit[ii,jj,5] = -cdelta*funInit[ii,jj,6] + 2/pi*skmR*T4*k1r*funInit[ii,jj,4]; 
            else
                funInit[ii,jj,5] = 2/pi*skmR*T4*k1r*funInit[ii,jj,4]+ psih*sdelta*TBblock2/Jg;
                funInit[ii,jj,6] = -spsih*TBblock1 + sdelta*cpsih*funInit[ii,jj,2] - cdelta*cpsih*TBblock2;
                funInit[ii,jj,7] =  cpsih*TBblock1 + sdelta*spsih*funInit[ii,jj,2] - cdelta*spsih*TBblock2;
            end
        end
    end
            
    return funInit,psilV,zetaV,ratio;
end

function getpartialaveragedvectorfield_drag(k,k1r,skm,smk,m,KK,EE,PP,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,psih,vMAt,VVT2,VVT3,VVT4,typeofvariable)
    T1,T2,T3,T4,T6,T7,T8 = getsingularityterms(m,k,smk,skm,KK,EE,PP)

    cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2 = dragsemianalyticalaveragedcoeff_fV_O3(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vMAt[71:91]);
    cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3 = dragsemianalyticalaveragedcoeff_T3(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,T7,T8,cdelta,sdelta,cpsih,spsih,vMAt);
    cMxb11T4,cMxb12T4,cMxb13T4,cMxpsilT4,cMyb21T4,cMyb22T4,cMyb23T4,cMypsilT4,cMzb31T4,cMzb32T4,cMzb33T4,cMzpsilT4 = dragsemianalyticalaveragedcoeff_T4(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,T7,T8,cdelta,sdelta,cpsih,spsih,vMAt);
 
    Mxb11T2,Mxb12T2,Mxb13T2,MxpsilT2,Myb21T2,Myb22T2,Myb23T2,MypsilT2,Mzb31T2,Mzb32T2,Mzb33T2,MzpsilT2 = drag_averagedterms(VVT2,cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2);
    Mxb11T3,Mxb12T3,Mxb13T3,MxpsilT3,Myb21T3,Myb22T3,Myb23T3,MypsilT3,Mzb31T3,Mzb32T3,Mzb33T3,MzpsilT3 = drag_averagedterms(VVT3,cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3);
    Mxb11T4,Mxb12T4,Mxb13T4,MxpsilT4,Myb21T4,Myb22T4,Myb23T4,MypsilT4,Mzb31T4,Mzb32T4,Mzb33T4,MzpsilT4 = drag_averagedterms(VVT4,cMxb11T4,cMxb12T4,cMxb13T4,cMxpsilT4,cMyb21T4,cMyb22T4,cMyb23T4,cMypsilT4,cMzb31T4,cMzb32T4,cMzb33T4,cMzpsilT4);

    mxb11  = Mxb11T2 + Mxb11T3 + Mxb11T4; 
    mxb12  = Mxb12T2 + Mxb12T3 + Mxb12T4; 
    mxb13  = Mxb13T2 + Mxb13T3 + Mxb13T4;
    mxpsil = MxpsilT2 + MxpsilT3 + MxpsilT4;

    myb21  = Myb21T2 + Myb21T3 + Myb21T4;
    myb22  = Myb22T2 + Myb22T3 + Myb22T4; 
    myb23  = Myb23T2 + Myb23T3 + Myb23T4;
    mypsil = MypsilT2 + MypsilT3 + MypsilT4; 

    mzb31  = Mzb31T2 + Mzb31T3 + Mzb31T4; 
    mzb32  = Mzb32T2 + Mzb32T3 + Mzb32T4;
    mzb33  = Mzb33T2 + Mzb33T3 + Mzb33T4;
    mzpsil = MzpsilT2 + MzpsilT3 + MzpsilT4; 

    skmBlock = -2*mxb13*skm - 2*myb23*(1-m)/(1+k)*skm + 2*smk*mzb33
    TBblock1 = mxb11+myb21+mzb31;
    TBblock2 = mxb12+myb22+mzb32;
    TBblock3 = mxb13+myb23+mzb33;
    TSblock  = mxpsil+(1-m)*mypsil+mzpsil;

    Z21 =  skmBlock/Jg;
    Z22 =  TBblock3;
    Z23 =  cdelta*Z22 + sdelta*TBblock2;
    Z24 =  pi/(2*mM1KK*Jg)*TSblock;

    if typeofvariable==1
        Z26  = TBblock1/Jg/sdelta;
        Z25  = -cdelta*Z26 + 2/pi*sqrt(skm)*T4*k1r*Z24;
        ZV = [Z21,Z22,Z23,Z24,Z25,Z26];
    else
        Z25 = 2/pi*sqrt(skm)*T4*k1r*Z24 + psih*sdelta*TBblock2/Jg;
        Z26 = -spsih*TBblock1 + sdelta*cpsih*Z22 - cdelta*cpsih*TBblock2;
        Z27 =  cpsih*TBblock1 + sdelta*spsih*Z22 - cdelta*spsih*TBblock2;
        ZV = [Z21,Z22,Z23,Z24,Z25,Z26,Z27];
    end

    return ZV;
end

function getpartialnonaveragedvectorfield_drag(k,k1r,skm,smk,m,KK,EE,PP,T4,mM1KK,Jg,cdelta,sdelta,cpsih,spsih,psih,vMAt,VVT2,VVT3,VVT4,typeofvariable,npsil,npsig,psilstep,zetastep)
    
    ratio = npsig/npsil;
    ll1 = Int(floor(2*pi/psilstep)+1.0)
    if ratio >=0.0
        zetaMin = -ratio*2*pi;
        zetaMax = 2*pi;
    else
        zetaMin = 0.0;
        zetaMax = 2*pi*(1-ratio);
    end
    ll2 = Int(floor((zetaMax-zetaMin)/(pi/180.0))+1.0);

    psilV = collect(range(0.0,2*pi,ll1));
    zetaV =  collect(range(zetaMin,zetaMax,ll2));
    
    vMAtH_VO3   = vmt2vmtH_order4(vMAt[77:91],cdelta,sdelta,cpsih,spsih)
    vMAtH_VO4   = vmt2vmtH_order5(vMAt[41:61],cdelta,sdelta,cpsih,spsih);
    vMAtH_NO23  = vmt2vmtH_order03(vMAt[21:40],cdelta,sdelta,cpsih,spsih)
    vMAtH_NO4   = vmt2vmtH_order4(vMAt[62:76],cdelta,sdelta,cpsih,spsih)
    
    skmR = sqrt(skm);
    smkR = sqrt(smk);

    if typeofvariable == 1
        funInit = zeros(ll1,ll2,6)
    else
        funInit = zeros(ll1,ll2,7)
    end
    for ii = 1 : ll1
        for jj = 1 : ll2

            mxb11    = 0.0
            mxb12    = 0.0
            mxb13    = 0.0
            mxpsil   = 0.0
            myb21    = 0.0
            myb22    = 0.0
            myb23    = 0.0
            mypsil   = 0.0
            mzb31    = 0.0
            mzb32    = 0.0
            mzb33    = 0.0
            mzpsil   = 0.0

            b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz = dragtermspsilzeta(psilV[ii],zetaV[jj],ratio,m,KK,PP,EE,k,k1r,skmR,smkR);

            mm = dragtermsV_O3(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_VO3)
            nn = length(VVT2)
            mxb11    = mxb11 + LinearAlgebra.dot(VVT2,mm[1:nn]); 
            mxb12    = mxb12 + LinearAlgebra.dot(VVT2,mm[nn+1:2*nn]); 
            mxb13    = mxb13 + LinearAlgebra.dot(VVT2,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT2,mm[3*nn+1:4*nn]); 
            myb21    = myb21 + LinearAlgebra.dot(VVT2,mm[4*nn+1:5*nn]); 
            myb22    = myb22 + LinearAlgebra.dot(VVT2,mm[5*nn+1:6*nn]); 
            myb23    = myb23 + LinearAlgebra.dot(VVT2,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT2,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31 + LinearAlgebra.dot(VVT2,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32 + LinearAlgebra.dot(VVT2,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33 + LinearAlgebra.dot(VVT2,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT2,mm[11*nn+1:12*nn]); 

            mm = dragtermsV_O4(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_VO4)
            nn = length(VVT3)
            mxb11    = mxb11 + LinearAlgebra.dot(VVT3,mm[1:nn]); 
            mxb12    = mxb12 + LinearAlgebra.dot(VVT3,mm[nn+1:2*nn]); 
            mxb13    = mxb13 + LinearAlgebra.dot(VVT3,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT3,mm[3*nn+1:4*nn]); 
            myb21    = myb21 + LinearAlgebra.dot(VVT3,mm[4*nn+1:5*nn]); 
            myb22    = myb22 + LinearAlgebra.dot(VVT3,mm[5*nn+1:6*nn]); 
            myb23    = myb23 + LinearAlgebra.dot(VVT3,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT3,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31 + LinearAlgebra.dot(VVT3,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32 + LinearAlgebra.dot(VVT3,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33 + LinearAlgebra.dot(VVT3,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT3,mm[11*nn+1:12*nn]); 

            mm = dragtermsN_O3(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_NO23)
            nn = length(VVT3)
            mxb11    = mxb11 + LinearAlgebra.dot(VVT3,mm[1:nn]); 
            mxb12    = mxb12 + LinearAlgebra.dot(VVT3,mm[nn+1:2*nn]); 
            mxb13    = mxb13 + LinearAlgebra.dot(VVT3,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT3,mm[3*nn+1:4*nn]); 
            myb21    = myb21 + LinearAlgebra.dot(VVT3,mm[4*nn+1:5*nn]); 
            myb22    = myb22 + LinearAlgebra.dot(VVT3,mm[5*nn+1:6*nn]); 
            myb23    = myb23 + LinearAlgebra.dot(VVT3,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT3,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31 + LinearAlgebra.dot(VVT3,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32 + LinearAlgebra.dot(VVT3,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33 + LinearAlgebra.dot(VVT3,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT3,mm[11*nn+1:12*nn]); 

            mm = dragtermsN_O4(b11,b12,b13,b21,b22,b23,b31,b32,b33,Sx,Sy,Sz,vMAtH_NO4)
            nn = length(VVT4)
            mxb11    = mxb11  + LinearAlgebra.dot(VVT4,mm[1:nn]); 
            mxb12    = mxb12  + LinearAlgebra.dot(VVT4,mm[nn+1:2*nn]); 
            mxb13    = mxb13  + LinearAlgebra.dot(VVT4,mm[2*nn+1:3*nn]); 
            mxpsil   = mxpsil + LinearAlgebra.dot(VVT4,mm[3*nn+1:4*nn]); 
            myb21    = myb21  + LinearAlgebra.dot(VVT4,mm[4*nn+1:5*nn]); 
            myb22    = myb22  + LinearAlgebra.dot(VVT4,mm[5*nn+1:6*nn]); 
            myb23    = myb23  + LinearAlgebra.dot(VVT4,mm[6*nn+1:7*nn]); 
            mypsil   = mypsil + LinearAlgebra.dot(VVT4,mm[7*nn+1:8*nn]); 
            mzb31    = mzb31  + LinearAlgebra.dot(VVT4,mm[8*nn+1:9*nn]); 
            mzb32    = mzb32  + LinearAlgebra.dot(VVT4,mm[9*nn+1:10*nn]); 
            mzb33    = mzb33  + LinearAlgebra.dot(VVT4,mm[10*nn+1:11*nn]); 
            mzpsil   = mzpsil + LinearAlgebra.dot(VVT4,mm[11*nn+1:12*nn]); 

            TBblock1 = mxb11 + myb21 + mzb31
            TBblock2 = mxb12 + myb22 + mzb32
            TSblock  = mxpsil + (1-m)*mypsil + mzpsil

            funInit[ii,jj,1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
            funInit[ii,jj,2] =  mxb13+myb23+mzb33;
            funInit[ii,jj,3] =  cdelta*funInit[ii,jj,2] + sdelta*TBblock2;
            funInit[ii,jj,4] =  pi/(2*mM1KK*Jg)*TSblock;
            
            if typeofvariable == 1
                funInit[ii,jj,6] = TBblock1/Jg/sdelta;
                funInit[ii,jj,5] = -cdelta*funInit[ii,jj,6] + 2/pi*skmR*T4*k1r*funInit[ii,jj,4]; 
            else
                funInit[ii,jj,5] = 2/pi*skmR*T4*k1r*funInit[ii,jj,4]+ psih*sdelta*TBblock2/Jg;
                funInit[ii,jj,6] = -spsih*TBblock1 + sdelta*cpsih*funInit[ii,jj,2] - cdelta*cpsih*TBblock2;
                funInit[ii,jj,7] =  cpsih*TBblock1 + sdelta*spsih*funInit[ii,jj,2] - cdelta*spsih*TBblock2;
            end
        end
    end
            
    return funInit,psilV,zetaV,ratio;
end

function dragtermspsilzeta(psil,zeta,ratio,m,KK,PP,EE,k,k1r,skmR,smkR)
    psig = zeta+ratio*psil;
    psilpsigterms  = getPsilPsigterms_drag([psig],[psil],m,KK,EE,PP,k,k1r,skmR,smkR)
    return  psilpsigterms[1,1,1], psilpsigterms[1,1,2], psilpsigterms[1,1,3], psilpsigterms[1,1,4], psilpsigterms[1,1,5], psilpsigterms[1,1,6],psilpsigterms[1,1,7], psilpsigterms[1,1,8], psilpsigterms[1,1,9],psilpsigterms[1,1,10], psilpsigterms[1,1,11], psilpsigterms[1,1,12];
end

##################### W2 numerical computation

function getgenfunterms_psilpsig(npsil,ratio,dnpsildskm,dnpsildJg,dnpsigdskm,dnpsigdJg,psig,psil,psilV,zetaV,funStruct)    
    
    ll = length(psilV);
    mm = length(zetaV);
    
    if size(funStruct)[1]!=ll || size(funStruct)[2]!=mm
        Base.error("wrong size");
    end
    nn = size(funStruct)[3]
    
    if nn!=6 && nn!=7
        Base.error("wrong size");
    end

    genfun = zeros(nn)
    ggenfun = zeros(2)

    for kk = 1 : nn
        fun = funStruct[:,:,kk]

        ## get FUN = int_0^psil fun
        FUN = zeros(ll,mm)
        for ii = 2 : ll
            FUN[ii,:] = FUN[ii-1,:] + (psilV[ii]-psilV[ii-1])*(fun[ii,:]+fun[ii-1,:])/2.0;
        end

        ## correction
        FUNcc = zeros(mm);
        for jj = 1 : mm
            FUNcc[jj] = sum((FUN[2:ll,jj]+FUN[1:ll-1,jj]).*(psilV[2:ll]-psilV[1:ll-1])/2.0);
        end
        FUNcc = -FUNcc/(2*pi);
        
        for ii = 1 : ll
            FUN[ii,:] = FUN[ii,:] + FUNcc;
        end
        
        ## get average of FUN 
        FUNmean = numericalintegral_psilzeta(FUN,psilV,zetaV,ratio)/(4*pi^2.0);
        
        # interpolation
        FUN_psil = zeros(ll)
        for ii = 1 : ll
            spline  = CubicSpline(FUN[ii,:],zetaV);
            FUN_psil[ii]  = spline(psig-ratio*psilV[ii]);
        end
        spline = CubicSpline(FUN_psil,psilV);
        
        # generating function term
        genfun[kk] = (spline(psil)-FUNmean)/npsil;

        if kk==1 || kk ==2
            
            ## get FFUN = int_0^psil FUN
            FFUN = zeros(ll,mm)
            for ii = 2 : ll
                FFUN[ii,:] = FUN[ii-1,:] + (psilV[ii]-psilV[ii-1])*(FUN[ii,:]+FUN[ii-1,:])/2.0;
            end

            ## correction
            FUNcc = zeros(mm);
            for jj = 1 : mm
                FUNcc[jj] = sum((FFUN[2:ll,jj]+FFUN[1:ll-1,jj]).*(psilV[2:ll]-psilV[1:ll-1])/2.0);
            end
            FUNcc = -FUNcc/(2*pi);
            for ii = 1 : ll
                FFUN[ii,:] = FFUN[ii,:] + FUNcc;
            end

            ## get average of FFUN 
            FFUNmean = numericalintegral_psilzeta(FFUN,psilV,zetaV,ratio)/(4*pi^2.0);

            # interpolation
            FFUN_psil = zeros(ll)
            for ii = 1 : ll
                spline  = CubicSpline(FFUN[ii,:],zetaV);
                FFUN_psil[ii]  = spline(psig-ratio*psil);
            end
            spline = CubicSpline(FFUN_psil,psilV);

            # generating function term
            ggenfun[kk] = (spline(psil) - FFUNmean)/(npsil^2.0) - FUNmean/npsil*(atan(tan(psil/2))/npsil+atan(tan(psig/2))/(ratio*npsil));
        end
    end

    W2 = genfun;
    W2[4] = W2[4] + dnpsildskm*ggenfun[1] + dnpsildJg*ggenfun[2];
    W2[5] = W2[5] + dnpsigdskm*ggenfun[1] + dnpsigdJg*ggenfun[2];

    return W2;
end

function numericalintegral_psilzeta(fun_psil_zeta,psilV,zetaV,ratio)
    SS = size(fun_psil_zeta);
    fun_psil = zeros(SS[1])

    for ii = 1 : SS[1]
        zMin = -ratio*psilV[ii]
        zMax = 2.0*pi-ratio*psilV[ii];

        if zMax > zetaV[end]
            if abs(zMax-zetaV[end])<1e-13
                zMax = zetaV[end]
            else
                println( abs(zMax-zetaV[end]))
                Base.error("something wrong in the definition of the integration interval")
            end
        end

        if zMin <zetaV[1]
            if abs(zMin-zetaV[1])<1e-13
                zMin = zetaV[1]
            else
                Base.error("something wrong in the definition of the integration interval")
            end
        end
       
        ffzeta  = fun_psil_zeta[ii,:];
        spline  = CubicSpline(ffzeta,zetaV);
        ffMin = spline(zMin);
        ffMax = spline(zMax);      

        idx = zetaV.>zMin .&& zetaV.<zMax;
        zzV = zetaV[idx];
        ffV = ffzeta[idx];
    
        S2 = size(zzV)[1]
        fun_psil[ii] = (ffMin +  ffV[1])*(zzV[1]-zMin)/2.0 + (ffV[end] +  ffMax)*(zMax-zzV[end])/2.0 + sum((ffV[1:S2-1]+ffV[2:S2]).*(zzV[2:S2]-zzV[1:S2-1])/2.0);

    end

    int = sum((fun_psil[1:SS[1]-1]+fun_psil[2:SS[1]]).*(psilV[2:SS[1]]-psilV[1:SS[1]-1])/2.0);
    return int;
end
