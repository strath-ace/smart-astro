"""
    getDRAG2bis(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion associated to Jg and skm 
    (Sadov variables) with respect to the orbital mean anomaly
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    IV = [A,B,C] principal moments of inertia
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    equi = equinoctial elements (with true longitude)
    satellite = Dict with the characteristics of the satellite ***
    vmt = terms of the equations of motion depending on the orbital mean anomaly averaged
    cpa1,cpa2,cpa3 = planet polar axis direction cosines in the inertial ref frame
    planet_flat = 1 to consider planet oblateness (important with an exponential model of the atmosphere)
    oblateness  = oblateness coefficient
    muPlanet = gravitational parameter of the central body
    atmosphericinfo = Dict with info about atmospheric model ****
    nE   = rotational angular rate of central body
    rPlanet = mean radius of the central body
    inertial2equatorialRM = rotation matrix from inertial reference frame to equatorial reference fra
    mjd2000 = epoch in modiefied julian day 2000

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
                
    ****
    atmosphericinfo
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
    W2bis = vector 6x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,psig,psih] -> [skm',Jg',Jh',psil',psig',psih'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWDRAG2bis(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)

    # constants depending on the satellite
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

    # quantities depending on m and k
    if skm == 1.0
        m = 0.0
    else
        m     = (1-skm)/skm*k;
    end
    KK    = ellF(pi/2,m);
    EE    = ellE(pi/2,m);
    KP    = ellF(pi/2,1-m);
    PP    = ellP(-k,pi/2,m);
    k1r = sqrt(1+k);
    smk = 1-skm;
    T1,T2,T3,T4,T6,T7,T8 = getsingularityterms(m,k,smk,skm,KK,EE,PP)

    # numerically averaged terms -> vmt

    # values of terms depending on mean anomaly before averaging
    sma = equi[1]
    P1  = equi[2]
    P2  = equi[3]
    Q1  = equi[4]
    Q2  = equi[5]
    TL  = equi[6]

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
    eta = sqrt(1.0-P1^2.0-P2^2.0)

    MAterms = getMAterms_drag([TL],sma,P1,P2,eta,RQ2i,truelong2meanlong(TL,P1,P2),mjd2000,cpa1,cpa2,cpa3,rPlanet,muPlanet,nE,0.0,planet_flat,oblateness,atmosphericinfo,false,inertial2equatorialRM);
    vMAt,vO = organizeMAterms_dragforce([TL],MAterms,P1,P2,eta,true,false,false,atmosphericinfo["considerthermalCD"],0);
    vMAt = vMAt[1,:];
    
    # averaged terms
    skmBlockA,TBblock1A,TBblock2A,TBblock3A,TSblockA = getdragsaaveragedtermsattitudeeq(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,T7,T8,cdelta,sdelta,cpsih,spsih,vmt,VVT1,VVT2,VVT3,VVT4)

    # non averaged terms
    skmBlock,TBblock1,TBblock2,TBblock3,TSblock = getdragsaaveragedtermsattitudeeq(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,T7,T8,cdelta,sdelta,cpsih,spsih,vMAt,VVT1,VVT2,VVT3,VVT4)
   
    # averaged term wrt psil and psig but not wrt to M
    Z21 =  skmBlock/Jg;
    Z22 =  TBblock3; 

    # averaged term wrt psil and psig and M
    Z21A =  skmBlockA/Jg;
    Z22A =  TBblock3A;

    # gen fun
    npsig = Jg/A/C*((C-A)PP+A*KK)/KK;
    dnpsildskm = pi*Jg*(A-C)*(EE-(1-m)*skm*KK)/(4.0*A*C*(1-skm)*sqrt(skm)*sqrt(1+k)-KK^2.0);
    dnpsildJg  = pi*sqrt(skm)*(A-C)/(2.0*A*C*sqrt(1+k)*KK);
    dnpsigdJg  = ((C-A)*PP+A*KK)/(A*C*KK);
    if skm == 1.0
        dnpsigdskm = -Jg*(A-C)*(-k/sqrt(1+k)/2.0+1.0)/(2*A*C);
    else
        dnpsigdskm = -Jg*(A-C)*((EE-(1-m)*skm*KK)*PP-(1-skm)*KK*EE)/(2*A*C*(1-m)*skm*(1-skm)*KK^2);
    end

    W21bis = 2*(Z21-Z21A)*atan(tan(psig/2))/npsig;
    W22bis = 2*(Z22-Z22A)*atan(tan(psig/2))/npsig;
    W23bis = 0.0;
    W24bis = 2.0*(dnpsildskm*(Z21-Z21A)+dnpsildJg*(Z22-Z22A))*(atan(tan(psig/2))^2-pi^2/12)/(npsig^2.0);
    W25bis = 2.0*(dnpsigdskm*(Z21-Z21A)+dnpsigdJg*(Z22-Z22A))*(atan(tan(psig/2))^2-pi^2/12)/(npsig^2.0);
    W26bis = 0.0; 
    W2bis  = [W21bis,W22bis,W23bis,W24bis,W25bis,W26bis];

    return W2bis;
end

"""
    getDRAG2bissadovlike(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion associated to Jg and skm 
    (Sadov variables) with respect to the orbital mean anomaly
    -- the perturbation here is the solar radiation pressure

    INPUT
    IV = [A,B,C] principal moments of inertia
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    equi = equinoctial elements (with true longitude)
    satellite = Dict with the characteristics of the satellite ***
    vmt = terms of the equations of motion depending on the orbital mean anomaly averaged
    cpa1,cpa2,cpa3 = planet polar axis direction cosines in the inertial ref frame
    planet_flat = 1 to consider planet oblateness (important with an exponential model of the atmosphere)
    oblateness  = oblateness coefficient
    muPlanet = gravitational parameter of the central body
    atmosphericinfo = Dict with info about atmospheric model ****
    nE   = rotational angular rate of central body
    rPlanet = mean radius of the central body
    inertial2equatorialRM = rotation matrix from inertial reference frame to equatorial reference fra
    mjd2000 = epoch in modiefied julian day 2000

    **
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
                
    ****
    atmosphericinfo
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
    W2bis = vector 7x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,J5,J6,J7] -> [skm',Jg',Jh',psil',J5',J6',J7'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWDRAG2bissadovlike(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    W2bis =  getWDRAG2bis(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000);
    W2bis = append!(W2bis,[0.0]);
    return W2bis;
end


"""
    getWDRAG3(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion associated to Jg and skm 
    (Sadov variables) with respect to the orbital mean anomaly
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    IV = [A,B,C] principal moments of inertia
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    equi = equinoctial elements (with true longitude)
    satellite = Dict with the characteristics of the satellite ***
    vmt = terms of the equations of motion depending on the orbital mean anomaly averaged
    cpa1,cpa2,cpa3 = planet polar axis direction cosines in the inertial ref frame
    planet_flat = 1 to consider planet oblateness (important with an exponential model of the atmosphere)
    oblateness  = oblateness coefficient
    muPlanet = gravitational parameter of the central body
    atmosphericinfo = Dict with info about atmospheric model ****
    nE   = rotational angular rate of central body
    rPlanet = mean radius of the central body
    inertial2equatorialRM = rotation matrix from inertial reference frame to equatorial reference fra
    mjd2000 = epoch in modiefied julian day 2000

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
                
    ****
    atmosphericinfo
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
    W3 = vector 6x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,psig,psih] -> [skm',Jg',Jh',psil',psig',psih'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWDRAG3(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)

    W3 = getWDRAG3procedure(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,1)
   
    return W3;
end

"""
    getWDRAG3sadovlike(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion associated to Jg and skm 
    (Sadov-like variables) with respect to the orbital mean anomaly
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    IV = [A,B,C] principal moments of inertia
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    equi = equinoctial elements (with true longitude)
    satellite = Dict with the characteristics of the satellite ***
    vmt = terms of the equations of motion depending on the orbital mean anomaly averaged
    cpa1,cpa2,cpa3 = planet polar axis direction cosines in the inertial ref frame
    planet_flat = 1 to consider planet oblateness (important with an exponential model of the atmosphere)
    oblateness  = oblateness coefficient
    muPlanet = gravitational parameter of the central body
    atmosphericinfo = Dict with info about atmospheric model ****
    nE   = rotational angular rate of central body
    rPlanet = mean radius of the central body
    inertial2equatorialRM = rotation matrix from inertial reference frame to equatorial reference fra
    mjd2000 = epoch in modiefied julian day 2000

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
                
    ****
    atmosphericinfo
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
    W3 = vector 7x1 of Float64 -
        generating function of the transformation 
            [skm,Jg,Jh,psil,J5,J6,J7] -> [skm',Jg',Jh',psil',J5',J6',J7'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWDRAG3sadovlike(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000)

    W3 = getWDRAG3procedure(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,2)

    return W3;
end

############### function required

function  getWDRAG3procedure(k,sadovvar,skm,equi,satellite,vmt,cpa1,cpa2,cpa3,planet_flat,oblateness,muPlanet,atmosphericinfo,nE,rPlanet,inertial2equatorialRM,mjd2000,typeofvariable)

    # constants depending on the satellite
    # A = satellite["MomentsOfInertia"][1]; C = satellite["MomentsOfInertia"][3];
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

    # quantities depending on m and k
    if skm == 1.0
        m = 0.0
    else
        m     = (1-skm)/skm*k;
    end
    KK    = ellF(pi/2,m);
    EE    = ellE(pi/2,m);
    KP    = ellF(pi/2,1-m);
    PP    = ellP(-k,pi/2,m);
    k1r = sqrt(1+k);
    smk = 1-skm;
    T1,T2,T3,T4,T6,T7,T8 = getsingularityterms(m,k,smk,skm,KK,EE,PP)
    mM1KK = (m-1)*KK; 
    

    # generating function terms
    vmtgenfun = gettermsgeneratingfunction_averageoverM(vmt,equi,cpa1,cpa2,cpa3,atmosphericinfo,nE,muPlanet,rPlanet,planet_flat,oblateness,inertial2equatorialRM,mjd2000);

    skmBlock,TBblock1,TBblock2,TBblock3,TSblock = getdragsaaveragedtermsattitudeeq(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,T7,T8,cdelta,sdelta,cpsih,spsih,vmtgenfun,VVT1,VVT2,VVT3,VVT4)


    if typeofvariable == 1
        W3 = zeros(6);
        W3[3] = cdelta*TBblock3+ sdelta*TBblock2;
        W3[4] = pi/(2*mM1KK*Jg)*TSblock;
        W3[6] = TBblock1/Jg/sdelta;
        W3[5] = -cdelta*W3[6]+ 2/pi*sqrt(skm)*T4*k1r*W3[4];
    else
        W3 = zeros(7);
        W3[3] = cdelta*TBblock3 + sdelta*TBblock2;
        W3[4] =  pi/(2*mM1KK*Jg)*TSblock;
        W3[5] = 2/pi*sqrt(skm)*T4*k1r*W3[4] + psih*sdelta*TBblock2/Jg;
        W3[6] = -spsih*TBblock1 + sdelta*cpsih*TBblock3 - cdelta*cpsih*TBblock2;
        W3[7] =  cpsih*TBblock1 + sdelta*spsih*TBblock3 - cdelta*spsih*TBblock2;
    end

    return W3;

end


function gettermsgeneratingfunction_averageoverM(vmt,equi,cpa1,cpa2,cpa3,atmosphericinfo,nE,muPlanet,rPlanet,planet_flat,oblateness,inertial2equatorialRM,mjd2000)

    sma = equi[1]
    P1  = equi[2]
    P2  = equi[3]
    Q1  = equi[4]
    Q2  = equi[5]
    TLcurrent  = equi[6]
    MLcurrent  = truelong2meanlong(TLcurrent,P1,P2); 

    eta = sqrt(1-P1^2.0-P2^2.0);
    
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

    meanmotion = sqrt(muPlanet/(sma^3.0)); 

    ll = 361;
    TLV = collect(range(TLcurrent,TLcurrent + 2.0*pi,ll));
    MAterms = getMAterms_drag(TLV,sma,P1,P2,eta,RQ2i,truelong2meanlong(TLcurrent,P1,P2),mjd2000,cpa1,cpa2,cpa3,rPlanet,muPlanet,nE,0.0,planet_flat,oblateness,atmosphericinfo,false,inertial2equatorialRM);
    Ztj,useless = organizeMAterms_dragforce(TLV,MAterms,P1,P2,eta,true,false,false,atmosphericinfo["considerthermalCD"],0);
    Ztj     = Ztj[:,1:length(vmt)]

    zM = zeros(length(vmt),ll)
    for ii = 1 : ll
        zM[:,ii] = Ztj[ii,:] - vmt*eta^3/((1.0 .+ P1*sin.(TLV[ii]) + P2*cos.(TLV[ii])).^2.0);
    end

    F0M = zeros(length(vmt),ll); 
    for jj = 2:ll
        F0M[:,jj] = F0M[:,jj-1] + (zM[:,jj] .+ zM[:,jj-1])*(TLV[jj]-TLV[jj-1])/2;
    end


    F0MA = zeros(length(vmt));
    for jj = 2:ll
        F0MA= F0MA + (F0M[:,jj] .+ F0M[:,jj-1])*(TLV[jj]-TLV[jj-1])/2;
    end
    F0MA = F0MA/2/pi;

    W = F0M .- F0MA; 

    idxnu = 1;
    for jj=1:length(TLV)-1
        if TLcurrent>=TLV[jj] && TLcurrent<TLV[jj+1]
            idxnu = jj;
            break; 
        end 
    end

    W1 = W[:,idxnu]; 
    W2 = W[:,idxnu+1];
    W  = (W2-W1)/(TLV[idxnu+1]-TLV[idxnu])*TLcurrent + W2 - (W2-W1)/(TLV[idxnu+1]-TLV[idxnu])*TLV[idxnu+1];


   return W;

end

