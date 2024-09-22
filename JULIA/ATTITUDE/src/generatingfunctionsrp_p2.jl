"""
    getWSRP2bis(IV,k,sadovvar,skm,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion associated to Jg and skm 
    (Sadov variables) with respect to the orbital mean anomaly
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    IV = [A,B,C] principal moments of inertia
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    P1 = ec*sin(om+OM), ec = orbital eccentricity, om = argumentum of pericentre, OM = right ascension of the ascending node  
    P2 = ec*cos(om+OM) 
    EL = eccentric longitude (E+om+OM, E=eccentric anomaly)
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
    W2bis = vector 6x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,psig,psih] -> [skm',Jg',Jh',psil',psig',psih'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWSRP2bis(k,sadovvar,skm,P1,P2,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    
    if includeEclipsesEffects  == false || passageinshadowoccurs == false
        W2bis = zeros(6);
    else
        # constants depending on the satellite
        A = satellite["MomentsOfInertia"][1]; C = satellite["MomentsOfInertia"][3];
        numberOfFacets = get(satellite,"numberOfFacets",0);
        facets = get(satellite,"facets",0.0);  
        
        VVT1T = zeros(30);
        VVT2T = zeros(60);
        VVT3T = zeros(105);
        for kk = 1:numberOfFacets
            VVT1   = facets[kk][5];
            VVT2   = facets[kk][6];
            VVT3   = facets[kk][7];
            VVT1T  = VVT1T + VVT1;
            VVT2T  = VVT2T + VVT2;
            VVT3T  = VVT3T + VVT3;
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

        # uVect --> satellite to Sun unit vector approximated as Earth to Sun unit vector in the inertial ref frame
        earthsundist = norm(rSun);
        uVect = rSun/earthsundist;        
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

        # quantities depending on m and k
        if skm == 1.0
            m = 0.0
        else
            m     = (1-skm)/skm*k;
        end
        KK    = ellF(pi/2,m);
        EE    = ellE(pi/2,m);
        KP    = ellF(pi/2,1-m);
        EP    = ellE(pi/2,1-m);
        PP    = ellP(-k,pi/2,m);
        qn     = exp(-pi*KP/KK);
        k1r = sqrt(1+k);
        smk = 1-skm;
        if m==0.0 
            T1 = 1.0/2.0;
            T2 = -3.0/8.0;
            T3 = 1.0;
            T4 = -pi/2.0/sqrt(1+k);
            T6 = 0.0;
        else
            T1 = ((m-1.0)*KK+EE)/m/KK;
            T2 = 2.0*(KK-EE)/KK/m^2.0 + (EE-2.0*KK)/m/KK;
            T3 = EE/KK;
            if m==1.0
                T4 = -atan(sqrt(k))/sqrt(k);
            else
                T4 = (KK*m - PP*(m+k))/k;
            end
            T6 = (KK^2.0*(m-1.0)+EE^2.0)/KK^2.0/k;
        end

        # averaged terms
        cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1 = srpanalyticalaveragedcoeff_T1(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
        cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2 = srpanalyticalaveragedcoeff_T2(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
        cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3 = srpanalyticalaveragedcoeff_T3(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);

        Mxb11T1,Mxb12T1,Mxb13T1,MxpsilT1,Myb21T1,Myb22T1,Myb23T1,MypsilT1,Mzb31T1,Mzb32T1,Mzb33T1,MzpsilT1 = srp_averagedterms(VVT1T,cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1);
        Mxb11T2,Mxb12T2,Mxb13T2,MxpsilT2,Myb21T2,Myb22T2,Myb23T2,MypsilT2,Mzb31T2,Mzb32T2,Mzb33T2,MzpsilT2 = srp_averagedterms(VVT2T,cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2);
        Mxb11T3,Mxb12T3,Mxb13T3,MxpsilT3,Myb21T3,Myb22T3,Myb23T3,MypsilT3,Mzb31T3,Mzb32T3,Mzb33T3,MzpsilT3 = srp_averagedterms(VVT3T,cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3);
        mxb13  = Mxb13T1 + Mxb13T2 + Mxb13T3;
        myb23  = Myb23T1 + Myb23T2 + Myb23T3;
        mzb33  = Mzb33T1 + Mzb33T2 + Mzb33T3;


        # shadowfun
        EL1 = ELin;
        EL2 = ELout;
        sf  = smoothshadowfun(EL,EL1,EL2)
        msf = meansmoothshadowfun(EL1,EL2,P1,P2)

        # averaged term wrt psil and psig
        Z21 =  -2*skm/Jg*mxb13 - 2*skm/Jg*(1-m)/(1+k)*myb23 + 2*(1-skm)/Jg*mzb33;
        Z22 =  mxb13+myb23+mzb33; 

        # psrp corrections
        au = astroconstants(0)[2];
        psrpcorr     = (au/earthsundist)^2.0

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

        atanpsig = atan(tan(psig/2));
        # if atanpsig<0.0
        #     atanpsig = atanpsig + pi;
        # end
        
        W21bis = 2*Z21*(sf-msf)*atanpsig/npsig;
        W22bis = 2*Z22*(sf-msf)*atanpsig/npsig;
        W23bis = 0.0;
        W24bis = 2.0*(dnpsildskm*Z21+dnpsildJg*Z22)*(sf-msf)*(atanpsig^2-pi^2/12)/(npsig^2.0);
        W25bis = 2.0*(dnpsigdskm*Z21+dnpsigdJg*Z22)*(sf-msf)*(atanpsig^2-pi^2/12)/(npsig^2.0);
        W26bis = 0.0; 
        W2bis  = [W21bis,W22bis,W23bis,W24bis,W25bis,W26bis]*psrpcorr;
    end

    return W2bis;
end

"""
    getWSRP2bis(IV,k,sadovvar,skm,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion associated to Jg and skm 
    (Sadov-like variables) with respect to the orbital mean anomaly
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    IV = [A,B,C] principal moments of inertia
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    P1 = ec*sin(om+OM), ec = orbital eccentricity, om = argumentum of pericentre, OM = right ascension of the ascending node  
    P2 = ec*cos(om+OM) 
    EL = eccentric longitude (E+om+OM, E=eccentric anomaly)
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
    W2bis = vector 6x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,J5,J6,J7] -> [skm',Jg',Jh',psil',J5',J6',J7'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWSRP2bissadovlike(k,sadovvar,skm,P1,P2,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    W2bis = getWSRP2bis(k,sadovvar,skm,P1,P2,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    W2bis = append!(W2bis,[0.0]);
    return W2bis;
end

"""
    getWSRP3(IV,k,sadovvar,skm,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion associated to Jg and skm 
    (Sadov variables) with respect to the orbital mean anomaly
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    IV = [A,B,C] principal moments of inertia
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    P1 = ec*sin(om+OM), ec = orbital eccentricity, om = argumentum of pericentre, OM = right ascension of the ascending node  
    P2 = ec*cos(om+OM) 
    EL = eccentric longitude (E+om+OM, E=eccentric anomaly)
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
    W2bis = vector 6x1 of Float64 -
        generating function of the transformation 
        [skm,Jg,Jh,psil,psig,psih] -> [skm',Jg',Jh',psil',psig',psih'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWSRP3(k,sadovvar,skm,sma,P1,P2,EL,rSun,mu,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    W3 =  getWSRP3procedure(k,sadovvar,skm,sma,P1,P2,EL,rSun,mu,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,1)
    return W3;
end

"""
    getWSRP3sadovlike(IV,k,sadovvar,skm,EL,rSun,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    Returns the generating function leading to the transformation 
    of variables allowing to average the attitude equations of motion associated to Jg and skm 
    (Sadov-like variables) with respect to the orbital mean anomaly
    -- the perturbation here is the solar radiation pressure
    
    INPUT
    IV = [A,B,C] principal moments of inertia
    k  = (B-A)/(C-A)*C/A
    sadovvar,skm = sadov variables and sadov related variable **
    P1 = ec*sin(om+OM), ec = orbital eccentricity, om = argumentum of pericentre, OM = right ascension of the ascending node  
    P2 = ec*cos(om+OM) 
    EL = eccentric longitude (E+om+OM, E=eccentric anomaly)
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
    W2bis = vector 6x1 of Float64 -
        generating function of the transformation 
            [skm,Jg,Jh,psil,J5,J6,J7] -> [skm',Jg',Jh',psil',J5',J6',J7'], 
        allowing to obtain averaged attitude equations of motion with respect to psil and psig
"""
function getWSRP3sadovlike(k,sadovvar,skm,sma,P1,P2,EL,rSun,mu,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout)
    W3 =  getWSRP3procedure(k,sadovvar,skm,sma,P1,P2,EL,rSun,mu,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,2)
    return W3;
end

################## 
function getWSRP3procedure(k,sadovvar,skm,sma,P1,P2,EL,rSun,mu,satellite,includeEclipsesEffects ,passageinshadowoccurs,ELin,ELout,typeofvariables)

    if includeEclipsesEffects  == false || passageinshadowoccurs == false
        W3 = zeros(6);
    else
        # constants depending on the satellite
        A = satellite["MomentsOfInertia"][1]; C = satellite["MomentsOfInertia"][3];
        numberOfFacets = get(satellite,"numberOfFacets",0);
        facets = get(satellite,"facets",0.0);  
        
        VVT1T = zeros(30);
        VVT2T = zeros(60);
        VVT3T = zeros(105);
        for kk = 1:numberOfFacets
            VVT1   = facets[kk][5];
            VVT2   = facets[kk][6];
            VVT3   = facets[kk][7];
            VVT1T  = VVT1T + VVT1;
            VVT2T  = VVT2T + VVT2;
            VVT3T  = VVT3T + VVT3;
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

        # uVect --> satellite to Sun unit vector approximated as Earth to Sun unit vector in the inertial ref frame
        earthsundist = norm(rSun);
        uVect = rSun/earthsundist;        
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

        # quantities depending on m and k
        if skm == 1.0
            m = 0.0
        else
            m     = (1-skm)/skm*k;
        end
        KK    = ellF(pi/2,m);
        EE    = ellE(pi/2,m);
        KP    = ellF(pi/2,1-m);
        EP    = ellE(pi/2,1-m);
        PP    = ellP(-k,pi/2,m);
        qn     = exp(-pi*KP/KK);
        k1r = sqrt(1+k);
        smk = 1-skm;
        if m==0.0 
            T1 = 1.0/2.0;
            T2 = -3.0/8.0;
            T3 = 1.0;
            T4 = -pi/2.0/sqrt(1+k);
            T6 = 0.0;
        else
            T1 = ((m-1.0)*KK+EE)/m/KK;
            T2 = 2.0*(KK-EE)/KK/m^2.0 + (EE-2.0*KK)/m/KK;
            T3 = EE/KK;
            if m==1.0
                T4 = -atan(sqrt(k))/sqrt(k);
            else
                T4 = (KK*m - PP*(m+k))/k;
            end
            T6 = (KK^2.0*(m-1.0)+EE^2.0)/KK^2.0/k;
        end

        # averaged terms
        cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1 = srpanalyticalaveragedcoeff_T1(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
        cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2 = srpanalyticalaveragedcoeff_T2(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
        cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3 = srpanalyticalaveragedcoeff_T3(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);

        Mxb11T1,Mxb12T1,Mxb13T1,MxpsilT1,Myb21T1,Myb22T1,Myb23T1,MypsilT1,Mzb31T1,Mzb32T1,Mzb33T1,MzpsilT1 = srp_averagedterms(VVT1T,cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1);
        Mxb11T2,Mxb12T2,Mxb13T2,MxpsilT2,Myb21T2,Myb22T2,Myb23T2,MypsilT2,Mzb31T2,Mzb32T2,Mzb33T2,MzpsilT2 = srp_averagedterms(VVT2T,cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2);
        Mxb11T3,Mxb12T3,Mxb13T3,MxpsilT3,Myb21T3,Myb22T3,Myb23T3,MypsilT3,Mzb31T3,Mzb32T3,Mzb33T3,MzpsilT3 = srp_averagedterms(VVT3T,cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3);
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
        JhdotMcdeltaJgdotOversdelta = mxb12+myb22+mzb32;
        psihdotsdeltaJg = (mxb11+myb21+mzb31);

        # generating function shadow
        EL1 = ELin;
        EL2 = ELout;
        n = sqrt(mu/sma^3.0);
        phiEA = -P1*cos(EL) + P2*sin(EL);
        gfunsf = genfunsmoothshadowfun(EL,EL1,EL2,phiEA,P1,P2,n); 

        # psrp corrections
        au = astroconstants(0)[2];
        psrpcorr     = (au/earthsundist)^2.0

        if typeofvariables == 1
            Z23 = cdelta*(mxb13+myb23+mzb33)+ sdelta*JhdotMcdeltaJgdotOversdelta;
            Z24 = pi/(2*(m-1)*KK*Jg)*(mxpsil+(1-m)*mypsil+mzpsil);
            Z26 = psihdotsdeltaJg/Jg/sdelta;
            Z25 = -cdelta*Z26+ 2/pi*sqrt(skm)*T4*k1r*Z24;
            W31  = [0.0];
            W32  = [0.0];
            W33  = Z23;
            W34  = Z24;
            W35  = Z25;
            W36  = Z26;
            W3   = append!(W31,W32,W33,W34,W35,W36);
        else
            Jgdot = mxb31myb23+mzb33;
            # terms averaged over psil and psig      
            Z23 = cdelta*Jgdot + sdelta*JhdotMcdeltaJgOversdelta;
            Z24 = pi/(2*mM1KK*Jg)*(mxpsil+(1-m)*mypsil+mzpsil);
            Z25 = 2/pi*sqrt(skm)*T4*k1r*Z24 + psih*sdelta*JhdotMcdeltaJgOversdelta/Jg;
            Z26 = -spsih*psihdotsdeltaJg + sdelta*cpsih*Jgdot - cdelta*cpsih*JhdotMcdeltaJgOversdelta;
            Z27 =  cpsih*psihdotsdeltaJg + sdelta*spsih*Jgdot - cdelta*spsih*JhdotMcdeltaJgOversdelta;
            
            # W3   
            W31  = [0.0];
            W32  = [0.0];
            W33  = Z23;
            W34  = Z24;
            W35  = Z25;
            W36  = Z26;
            W37  = Z27;
            W3   = append!(W31,W32,W33,W34,W35,W36,W37);
        end
        W3   = W3*gfunsf*psrpcorr;

    end
    return W3;
end


function meansmoothshadowfun(EL1,EL2,P1,P2)
    sfmean = P1 * cos(EL1) / pi / 2 - P1 * cos(EL2) / pi / 2 - P2 * sin(EL1) / pi / 2 + P2 * sin(EL2) / pi / 2 + (2 * pi + EL1 - EL2) / pi / 2;
    return sfmean;
end

function genfunsmoothshadowfun(EL,EL1,EL2,phiEA,P1,P2,n)
    t1 = 0.1e1 / pi
    t2 = P2 * t1
    t3 = 0.1e1 / n
    t4 = 3 * EL1
    t5 = 2 * EL
    t6 = -t4 + t5
    t7 = cos(t6)
    t11 = 3 * EL
    t12 = 4 * EL1
    t13 = t11 - t12
    t14 = cos(t13)
    t18 = 4 * EL2
    t19 = t11 - t18
    t20 = cos(t19)
    t24 = 2 * EL2
    t25 = t11 - t24
    t26 = cos(t25)
    t30 = 4 * EL
    t31 = 5 * EL1
    t32 = t30 - t31
    t33 = cos(t32)
    t37 = 5 * EL2
    t38 = t30 - t37
    t39 = cos(t38)
    t43 = 3 * EL2
    t44 = t30 - t43
    t45 = cos(t44)
    t49 = 5 * EL
    t50 = 6 * EL1
    t51 = t49 - t50
    t52 = cos(t51)
    t56 = t49 - t12
    t57 = cos(t56)
    t61 = 6 * EL2
    t62 = t49 - t61
    t63 = cos(t62)
    t67 = t49 - t18
    t68 = cos(t67)
    t72 = 6 * EL
    t73 = 7 * EL1
    t74 = t72 - t73
    t75 = cos(t74)
    t79 = t72 - t31
    t80 = cos(t79)
    t84 = -t7 * t3 * t2 / 12 - t14 * t3 * t2 / 24 + t20 * t3 * t2 / 24 + t26 * t3 * t2 / 12 - t33 * t3 * t2 / 40 + t39 * t3 * t2 / 40 + t45 * t3 * t2 / 24 - t52 * t3 * t2 / 60 - t57 * t3 * t2 / 40 + t63 * t3 * t2 / 60 + t68 * t3 * t2 / 40 - t75 * t3 * t2 / 84 - t80 * t3 * t2 / 60
    t85 = 7 * EL2
    t86 = t72 - t85
    t87 = cos(t86)
    t91 = P1 * t1
    t92 = 15 * EL
    t93 = 16 * EL2
    t94 = t92 - t93
    t95 = sin(t94)
    t99 = 14 * EL2
    t100 = t92 - t99
    t101 = sin(t100)
    t105 = 16 * EL
    t106 = 15 * EL1
    t107 = t105 - t106
    t108 = sin(t107)
    t112 = 16 * EL1
    t113 = t92 - t112
    t114 = sin(t113)
    t118 = 14 * EL1
    t119 = t92 - t118
    t120 = sin(t119)
    t124 = 14 * EL
    t125 = 15 * EL2
    t126 = t124 - t125
    t127 = sin(t126)
    t131 = 13 * EL2
    t132 = t124 - t131
    t133 = sin(t132)
    t137 = t124 - t106
    t138 = sin(t137)
    t142 = 13 * EL1
    t143 = t124 - t142
    t144 = sin(t143)
    t148 = 13 * EL
    t149 = t148 - t99
    t150 = sin(t149)
    t154 = 12 * EL2
    t155 = t148 - t154
    t156 = sin(t155)
    t160 = 12 * EL
    t161 = 11 * EL2
    t162 = t160 - t161
    t163 = sin(t162)
    t167 = t87 * t3 * t2 / 84 - t95 * t3 * t91 / 480 + t101 * t3 * t91 / 420 - t108 * t3 * t91 / 480 + t114 * t3 * t91 / 480 - t120 * t3 * t91 / 420 - t127 * t3 * t91 / 420 + t133 * t3 * t91 / 364 + t138 * t3 * t91 / 420 - t144 * t3 * t91 / 364 - t150 * t3 * t91 / 364 + t156 * t3 * t91 / 312 + t163 * t3 * t91 / 264
    t169 = t148 - t118
    t170 = sin(t169)
    t174 = 12 * EL1
    t175 = t148 - t174
    t176 = sin(t175)
    t180 = t160 - t131
    t181 = sin(t180)
    t185 = t160 - t142
    t186 = sin(t185)
    t190 = 11 * EL1
    t191 = t160 - t190
    t192 = sin(t191)
    t196 = 11 * EL
    t197 = t196 - t154
    t198 = sin(t197)
    t202 = 10 * EL2
    t203 = t196 - t202
    t204 = sin(t203)
    t208 = t196 - t174
    t209 = sin(t208)
    t213 = t5 - t43
    t214 = sin(t213)
    t218 = t5 - EL2
    t219 = sin(t218)
    t223 = 10 * EL1
    t224 = 9 * EL
    t225 = -t223 + t224
    t226 = sin(t225)
    t230 = -t223 + t196
    t231 = sin(t230)
    t235 = sin(t6)
    t239 = t170 * t3 * t91 / 364 - t176 * t3 * t91 / 312 - t181 * t3 * t91 / 312 + t186 * t3 * t91 / 312 - t192 * t3 * t91 / 264 - t198 * t3 * t91 / 264 + t204 * t3 * t91 / 220 + t209 * t3 * t91 / 264 - t214 * t3 * t91 / 12 + t219 * t3 * t91 / 4 + t226 * t3 * t91 / 180 - t231 * t3 * t91 / 220 + t235 * t3 * t91 / 12
    t240 = 21 * EL
    t241 = 20 * EL1
    t242 = t240 - t241
    t243 = sin(t242)
    t247 = 18 * EL1
    t248 = 17 * EL
    t249 = -t247 + t248
    t250 = sin(t249)
    t254 = 19 * EL
    t255 = -t247 + t254
    t256 = sin(t255)
    t260 = 17 * EL1
    t261 = -t260 + t105
    t262 = sin(t261)
    t266 = 18 * EL
    t267 = 19 * EL1
    t268 = t266 - t267
    t269 = sin(t268)
    t273 = 17 * EL2
    t274 = t266 - t273
    t275 = sin(t274)
    t279 = t254 - t241
    t280 = sin(t279)
    t284 = 20 * EL
    t285 = t284 - t267
    t286 = sin(t285)
    t290 = t248 - t112
    t291 = sin(t290)
    t295 = t248 - t93
    t296 = sin(t295)
    t300 = t105 - t125
    t301 = sin(t300)
    t305 = t105 - t273
    t306 = sin(t305)
    t310 = -EL1 + t5
    t311 = cos(t310)
    t315 = 20 * EL2
    t316 = -t315 + t254
    t317 = cos(t316)
    t321 = -t243 * t3 * t91 / 840 + t250 * t3 * t91 / 612 - t256 * t3 * t91 / 684 + t262 * t3 * t91 / 544 + t269 * t3 * t91 / 684 + t275 * t3 * t91 / 612 + t280 * t3 * t91 / 760 - t286 * t3 * t91 / 760 - t291 * t3 * t91 / 544 + t296 * t3 * t91 / 544 + t301 * t3 * t91 / 480 - t306 * t3 * t91 / 544 - t311 * t3 * t2 / 4 + t317 * t3 * t2 / 760
    t324 = -t315 + t240
    t325 = cos(t324)
    t329 = -t4 + t30
    t330 = cos(t329)
    t334 = 2 * EL1
    t335 = -t334 + EL
    t336 = cos(t335)
    t340 = -t334 + t11
    t341 = cos(t340)
    t345 = cos(t225)
    t349 = cos(t230)
    t353 = cos(t255)
    t357 = cos(t261)
    t361 = -t260 + t266
    t362 = cos(t361)
    t366 = cos(t249)
    t370 = cos(t274)
    t374 = cos(t279)
    t378 = cos(t285)
    t382 = t325 * t3 * t2 / 840 - t330 * t3 * t2 / 24 - t336 * t3 * t2 / 4 - t341 * t3 * t2 / 12 - t345 * t3 * t2 / 180 - t349 * t3 * t2 / 220 - t353 * t3 * t2 / 684 - t357 * t3 * t2 / 544 - t362 * t3 * t2 / 612 - t366 * t3 * t2 / 612 + t370 * t3 * t2 / 612 - t374 * t3 * t2 / 760 - t378 * t3 * t2 / 760
    t383 = cos(t242)
    t387 = cos(t305)
    t391 = cos(t300)
    t395 = cos(t290)
    t399 = cos(t295)
    t403 = cos(t268)
    t407 = cos(t126)
    t411 = cos(t132)
    t415 = cos(t113)
    t419 = cos(t119)
    t423 = cos(t94)
    t427 = cos(t100)
    t431 = cos(t107)
    t435 = -t383 * t3 * t2 / 840 + t387 * t3 * t2 / 544 + t391 * t3 * t2 / 480 - t395 * t3 * t2 / 544 + t399 * t3 * t2 / 544 - t403 * t3 * t2 / 684 + t407 * t3 * t2 / 420 + t411 * t3 * t2 / 364 - t415 * t3 * t2 / 480 - t419 * t3 * t2 / 420 + t423 * t3 * t2 / 480 + t427 * t3 * t2 / 420 - t431 * t3 * t2 / 480
    t437 = cos(t149)
    t441 = cos(t155)
    t445 = cos(t137)
    t449 = cos(t143)
    t453 = sin(t13)
    t457 = sin(t19)
    t461 = sin(t25)
    t465 = sin(t32)
    t469 = 10 * EL
    t470 = t469 - t161
    t471 = sin(t470)
    t475 = t3 * t1
    t476 = -EL2 + EL
    t478 = cos(5 * t476)
    t482 = cos(4 * t476)
    t486 = cos(3 * t476)
    t490 = cos(2 * t476)
    t493 = t437 * t3 * t2 / 364 + t441 * t3 * t2 / 312 - t445 * t3 * t2 / 420 - t449 * t3 * t2 / 364 + t453 * t3 * t91 / 24 - t457 * t3 * t91 / 24 + t461 * t3 * t91 / 12 + t465 * t3 * t91 / 40 - t471 * t3 * t91 / 220 - t478 * t475 / 25 - t482 * t475 / 16 - t486 * t475 / 9 - t490 * t475 / 4
    t494 = cos(t476)
    t496 = 8 * EL
    t497 = t496 - t85
    t498 = cos(t497)
    t502 = 8 * EL1
    t503 = t224 - t502
    t504 = cos(t503)
    t508 = t224 - t202
    t509 = cos(t508)
    t513 = 8 * EL2
    t514 = t224 - t513
    t515 = cos(t514)
    t519 = t469 - t190
    t520 = cos(t519)
    t524 = 9 * EL1
    t525 = t469 - t524
    t526 = cos(t525)
    t530 = cos(t470)
    t534 = 9 * EL2
    t535 = t469 - t534
    t536 = cos(t535)
    t540 = cos(t208)
    t544 = cos(t197)
    t548 = cos(t203)
    t552 = cos(t185)
    t556 = cos(t191)
    t560 = -t494 * t475 + t498 * t3 * t2 / 112 - t504 * t3 * t2 / 144 + t509 * t3 * t2 / 180 + t515 * t3 * t2 / 144 - t520 * t3 * t2 / 220 - t526 * t3 * t2 / 180 + t530 * t3 * t2 / 220 + t536 * t3 * t2 / 180 - t540 * t3 * t2 / 264 + t544 * t3 * t2 / 264 + t548 * t3 * t2 / 220 - t552 * t3 * t2 / 312 - t556 * t3 * t2 / 264
    t564 = cos(t180)
    t568 = cos(t162)
    t572 = cos(t169)
    t576 = cos(t175)
    t580 = EL - t24
    t581 = sin(t580)
    t585 = 18 * EL2
    t586 = -t585 + t254
    t587 = cos(t586)
    t591 = 19 * EL2
    t592 = -t591 + t266
    t593 = cos(t592)
    t597 = -t591 + t284
    t598 = cos(t597)
    t602 = -t585 + t248
    t603 = cos(t602)
    t607 = sin(t324)
    t611 = sin(t592)
    t615 = sin(t597)
    t619 = sin(t602)
    t623 = t564 * t3 * t2 / 312 + t568 * t3 * t2 / 264 - t572 * t3 * t2 / 364 - t576 * t3 * t2 / 312 - t581 * t3 * t91 / 4 + t587 * t3 * t2 / 684 + t593 * t3 * t2 / 684 + t598 * t3 * t2 / 760 + t603 * t3 * t2 / 612 + t607 * t3 * t91 / 840 - t611 * t3 * t91 / 684 + t615 * t3 * t91 / 760 - t619 * t3 * t91 / 612
    t624 = sin(t586)
    t628 = sin(t329)
    t632 = sin(t335)
    t636 = sin(t340)
    t640 = sin(t310)
    t644 = sin(t316)
    t648 = sin(t361)
    t652 = sin(t535)
    t657 = 2 * pi + EL1 - EL2
    t662 = sin(t38)
    t666 = sin(t44)
    t670 = sin(t51)
    t674 = sin(t56)
    t678 = t624 * t3 * t91 / 684 - t628 * t3 * t91 / 24 + t632 * t3 * t91 / 4 - t636 * t3 * t91 / 12 - t640 * t3 * t91 / 4 - t644 * t3 * t91 / 760 - t648 * t3 * t91 / 612 + t652 * t3 * t91 / 180 + phiEA * t3 * t1 * t657 / 2 - t662 * t3 * t91 / 40 + t666 * t3 * t91 / 24 + t670 * t3 * t91 / 60 - t674 * t3 * t91 / 40
    t680 = sin(t62)
    t684 = sin(t67)
    t688 = sin(EL1)
    t692 = cos(EL2)
    t696 = P2 ^ 2
    t697 = cos(t334)
    t701 = sin(EL2)
    t705 = cos(t24)
    t709 = cos(EL1)
    t713 = P1 ^ 2
    t720 = sin(t503)
    t724 = sin(t508)
    t728 = sin(t514)
    t732 = -t680 * t3 * t91 / 60 + t684 * t3 * t91 / 40 + t688 * P1 * t475 / 2 - t692 * P2 * t475 / 2 - t697 * t696 * t475 / 8 - t701 * P1 * t475 / 2 + t705 * t696 * t475 / 8 + t709 * P2 * t475 / 2 + t697 * t713 * t475 / 8 - t705 * t713 * t475 / 8 - t720 * t3 * t91 / 144 - t724 * t3 * t91 / 180 + t728 * t3 * t91 / 144
    t733 = sin(t519)
    t737 = sin(t525)
    t741 = t496 - t534
    t742 = sin(t741)
    t746 = sin(t497)
    t750 = 7 * EL
    t751 = t750 - t61
    t752 = sin(t751)
    t756 = t496 - t524
    t757 = sin(t756)
    t761 = t496 - t73
    t762 = sin(t761)
    t766 = t750 - t513
    t767 = sin(t766)
    t771 = sin(t74)
    t775 = sin(t79)
    t779 = sin(t86)
    t783 = t72 - t37
    t784 = sin(t783)
    t788 = t750 - t502
    t789 = sin(t788)
    t793 = t750 - t50
    t794 = sin(t793)
    t798 = t733 * t3 * t91 / 220 - t737 * t3 * t91 / 180 - t742 * t3 * t91 / 144 + t746 * t3 * t91 / 112 + t752 * t3 * t91 / 84 + t757 * t3 * t91 / 144 - t762 * t3 * t91 / 112 - t767 * t3 * t91 / 112 + t771 * t3 * t91 / 84 - t775 * t3 * t91 / 60 - t779 * t3 * t91 / 84 + t784 * t3 * t91 / 60 + t789 * t3 * t91 / 112 - t794 * t3 * t91 / 84
    t801 = cos(t783)
    t805 = cos(t788)
    t809 = cos(t793)
    t813 = cos(t766)
    t817 = cos(t751)
    t821 = cos(t580)
    t825 = cos(t213)
    t829 = cos(t218)
    t833 = cos(t756)
    t837 = cos(t761)
    t841 = cos(t741)
    t845 = -EL1 + EL
    t847 = cos(20 * t845)
    t851 = cos(19 * t845)
    t854 = t801 * t3 * t2 / 60 - t805 * t3 * t2 / 112 - t809 * t3 * t2 / 84 + t813 * t3 * t2 / 112 + t817 * t3 * t2 / 84 + t821 * t3 * t2 / 4 + t825 * t3 * t2 / 12 + t829 * t3 * t2 / 4 - t833 * t3 * t2 / 144 - t837 * t3 * t2 / 112 + t841 * t3 * t2 / 144 + t847 * t475 / 400 + t851 * t475 / 361
    t856 = cos(18 * t845)
    t860 = cos(17 * t845)
    t864 = cos(16 * t845)
    t868 = cos(15 * t845)
    t872 = cos(14 * t845)
    t876 = cos(13 * t845)
    t880 = cos(12 * t845)
    t884 = cos(11 * t845)
    t888 = cos(10 * t845)
    t892 = cos(9 * t845)
    t896 = cos(8 * t845)
    t900 = cos(7 * t845)
    t904 = cos(6 * t845)
    t908 = cos(5 * t845)
    t911 = t856 * t475 / 324 + t860 * t475 / 289 + t864 * t475 / 256 + t868 * t475 / 225 + t872 * t475 / 196 + t876 * t475 / 169 + t880 * t475 / 144 + t884 * t475 / 121 + t888 * t475 / 100 + t892 * t475 / 81 + t896 * t475 / 64 + t900 * t475 / 49 + t904 * t475 / 36 + t908 * t475 / 25
    t914 = cos(4 * t845)
    t918 = cos(3 * t845)
    t922 = cos(2 * t845)
    t925 = cos(t845)
    t928 = cos(20 * t476)
    t932 = cos(19 * t476)
    t936 = cos(18 * t476)
    t940 = cos(17 * t476)
    t944 = cos(16 * t476)
    t948 = cos(15 * t476)
    t952 = cos(14 * t476)
    t956 = cos(13 * t476)
    t960 = cos(12 * t476)
    t963 = t914 * t475 / 16 + t918 * t475 / 9 + t922 * t475 / 4 + t925 * t475 - t928 * t475 / 400 - t932 * t475 / 361 - t936 * t475 / 324 - t940 * t475 / 289 - t944 * t475 / 256 - t948 * t475 / 225 - t952 * t475 / 196 - t956 * t475 / 169 - t960 * t475 / 144
    t965 = cos(11 * t476)
    t969 = cos(10 * t476)
    t973 = cos(9 * t476)
    t977 = cos(8 * t476)
    t981 = cos(7 * t476)
    t985 = cos(6 * t476)
    t992 = P2 * P1
    t993 = sin(t334)
    t997 = sin(t24)
    t1002 = cos(EL)
    t1019 = sin(EL)
    t1023 = -t965 * t475 / 121 - t969 * t475 / 100 - t973 * t475 / 81 - t977 * t475 / 64 - t981 * t475 / 49 - t985 * t475 / 36 - phiEA * t692 * t3 * t91 / 2 - t993 * t992 * t475 / 4 + t997 * t992 * t475 / 4 + t1002 * t475 * t657 * P1 / 2 - phiEA * t688 * t3 * t2 / 2 + phiEA * t701 * t3 * t2 / 2 + phiEA * t709 * t3 * t91 / 2 - t1019 * t475 * t657 * P2 / 2
    genfunsf = t84 + t167 + t239 + t321 + t382 + t435 + t493 + t560 + t623 + t678 + t732 + t798 + t854 + t911 + t963 + t1023
    return genfunsf;

end


