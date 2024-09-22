# EQUATIONS OF MOTION: angular velocity + quaternions
function orbitAttitudeEv_quat(du,u,p,t)

    # if mod(floor(t/86400.0),2) == 0.0 
    #     println(t/86400)
    # end
    # println(t/86400)


    planetsunparams  = p[1];
    inertialrefframeinfo = p[2];
    satellite        = p[3];
    atmosphericinfo  = p[4];
    settings         = p[5];
    mjd2000_0        = p[6];
       
    # equinoctial elements
    mu = planetsunparams["muPlanet"]
    sma = u[1] 
    P1  = u[2]
    P2  = u[3]
    Q1  = u[4]
    Q2  = u[5]
    ML  = u[6] # mean anomaly
    
    eta = sqrt(1 - P1^2 - P2^2)
    bb = sma * eta;
    meanmotion = sqrt(mu / sma^3);
    hh = meanmotion * sma * bb;    

    # true longitude L = Omega + omega + nu (Battin, page. 493)
    TL = meanlong2truelong(ML,P1,P2);
    # if ier != 1
    #     Base.error("something went wrong in the solution of Kepler's equation")
    # end
    sTL = sin(TL); 
    cTL = cos(TL);     
    Phi = 1 + P1 * sTL + P2 * cTL


    # central body - centre of mass distance
    rnorm = sma*eta^2.0/Phi;
        
    # quaternion update
    quat = u[7:10];

    # angular velocity
    om  = u[11:13];
 
    # moments of inertia
    IV = get(satellite,"MomentsOfInertia",[0;0;0]);
    A = IV[1]; B = IV[2]; C = IV[3];

    # perturbations   
    exttorque = zeros(3,1); 
    pertacc   = zeros(3,1);

    # perturbations flags

    # attitude perturbations
    includeGravityTorque  = settings["includeGravityTorque"];
    includeMagneticTorque = settings["includeMagneticTorque"];
    includeSrpTorque      = settings["includeSrpTorque"];
    includeDragTorque     = settings["includeDragTorque"];

    # orbit perturbations
    includeZonalHarmsAcc  = settings["includeZonalHarmsAcc"];   
    includeThirdBodyAcc   = settings["includeThirdBodyAcc"];
    includeSunGravityAcc  = settings["includeSunGravityAcc"];
    includeSrpAcc         = settings["includeSrpAcc"];
    includeDragAcc        = settings["includeDragAcc"];

    # eclipses effects
    includeEclipsesEffectsOnAttitude =  get(settings,"includeEclipsesEffectsOnAttitude",false);
    includeEclipsesEffectsOnOrbit    =  get(settings,"includeEclipsesEffectsOnOrbit",false);
    pertacc = zeros(3)
    exttorque = zeros(3)
    # perturbations
    if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque || includeSrpAcc || includeDragAcc || includeThirdBodyAcc || includeSunGravityAcc
        Ri2b = transpose(quat2birotmat(quat));
        Ri2o = equi2rotmatirth(Q1,Q2,sTL,cTL);
        Ro2i = Transpose(Ri2o);
        rV = rnorm*Ro2i*[1,0.0,0.0]; 
        rUV = Ri2b*rV; rUV=rUV/rnorm;
        mjd2000   = mjd2000_0 + t/24.0/3600.0
        
        # Recef2eci = r_ecef_to_eci(J2000(), PEF(),mjd20002jd(mjd2000))
        # polaraxis = inertialrefframeinfo["equatorial2inertial"]*Recef2eci*[0.0,0.0,1.0];
        polaraxis = inertialrefframeinfo["equatorial2inertial"]*[0.0,0.0,1.0];

        gtorque = [0.0;0.0;0.0];
        if includeGravityTorque
            gtorque = getgravitytorque(mu,rUV,rnorm,A,B,C) 
        end

        mtorque = [0.0;0.0;0.0];
        if includeMagneticTorque
            mtorque = getmagnetictorque(planetsunparams["muM"],Ri2b*polaraxis,rUV,rnorm,satellite["intrinsicMagneticMoment"]);
        end

        dragtorque = [0.0;0.0;0.0];
        dragacc    = [0.0;0.0;0.0];
        if (includeDragTorque ||includeDragAcc) && rnorm-planetsunparams["rPlanet"]<1000.0
            vV =  sqrt(mu/sma/eta^2.0)*Ro2i*[P2*sTL-P1*cTL,P1*sTL+P2*cTL+1.0,0.0];
            atmAngVel = planetsunparams["planetRotation"]*polaraxis;
            if atmosphericinfo["atmosphericmodel"] == 1
                dragtorque,dragforce = getdragtorqueandforce_atmexp(rnorm,atmosphericinfo,planetsunparams["rPlanet"],satellite["CD"],satellite["facets"],satellite["numberOfFacets"],atmAngVel,vV,rV,Ri2b,om);
            else
                dragtorque,dragforce = getdragtorqueandforce_atmnrlmsise00(atmosphericinfo,satellite["facets"],satellite["numberOfFacets"],rV,vV,atmAngVel,Ri2b,transpose(inertialrefframeinfo["equatorial2inertial"]),mjd2000);
            end
            if includeDragAcc
                dragacc = Transpose(Ri2b)*dragforce/1000/satellite["mass"];
             end
             if includeDragTorque==false
                dragtorque = [0.0;0.0;0.0];
             end
        
        end

        srptorque = [0.0;0.0;0.0];
        srpacc    = [0.0;0.0;0.0];
        sungacc = [0.0;0.0;0.0];
        if includeSrpTorque || includeSrpAcc || includeSunGravityAcc
            rSunV = -inertialrefframeinfo["ecliptic2inertial"]*celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000)

            if includeSrpTorque || includeSrpAcc
                rsunnorm = sqrt(rSunV[1]^2.0+rSunV[2]^2.0+rSunV[3]^2.0)
                srptorque,srpforce = getsrptorqueandforce(rSunV,rsunnorm,rV,rnorm,Ri2b,satellite["facets"],satellite["numberOfFacets"],planetsunparams["rPlanet"],planetsunparams["psrp_refdist"],includeEclipsesEffectsOnAttitude,includeEclipsesEffectsOnOrbit);   
                if includeSrpAcc
                    srpacc = Transpose(Ri2b)*srpforce/1000/satellite["mass"];
                end
                if includeSrpTorque==false
                    srptorque = [0.0;0.0;0.0];
                end
            end
            
            if includeSunGravityAcc
                sungacc = getSunGravityPert(rSunV,rV);
            end
        end
   
        tbgacc = [0.0;0.0;0.0]
        if includeThirdBodyAcc 
            for jj = 1 : size(planetsunparams["perturbingbodiesid"])[1]
                idx3b = planetsunparams["perturbingbodiesid"][jj];
                if idx3b == 11
                    if planetsunparams["centralBodyIDX"] == 3
                        rV_3b = inertialrefframeinfo["equatorial2inertial"]*celestialbodiesephemeris_position(idx3b,mjd2000);
                    else
                        continue;
                    end
                else
                    rV_3bWRTS = celestialbodiesephemeris_position(idx3b,mjd2000);
                    rV_CBWRTS = celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000);
                    rV_3b = inertialrefframeinfo["ecliptic2inertial"]*(rV_3bWRTS-rV_CBWRTS);
                end
                tbgacc = tbgacc + get3BodyAcc(rV,planetsunparams["perturbingbodiesgravparam"][jj],rV_3b);
            end
        end

        pertacc = pertacc + Ri2o*(tbgacc+sungacc+dragacc+srpacc);

        exttorque = exttorque + gtorque+mtorque+srptorque+dragtorque;  
        #exttorque[abs.(exttorque).<1e-15].=0.0;     
    end

    if includeZonalHarmsAcc
        J2 = planetsunparams["zonalharmonicscoff"][1]
        J3 = planetsunparams["zonalharmonicscoff"][2]
        J4 = planetsunparams["zonalharmonicscoff"][3]
        J5 = planetsunparams["zonalharmonicscoff"][4]
        potpert = getCentralBodyPotentialAcc(sma,eta,Q1,Q2,cTL,sTL,Phi,mu,planetsunparams["rPlanet"],J2,J3,J4,J5)
        pertacc = pertacc + potpert;
    end  

    # variation in time of equinoctial elements
    aR  = pertacc[1];
    aT  = pertacc[2]; 
    aN  = pertacc[3];   

    du[1] = (2 / eta) * sqrt(sma^3 / mu) * ((P2 * sTL - P1 * cTL) * aR + Phi * aT)
    du[2] = eta * sqrt(sma / mu) * (-cTL * aR + (((P1 + sTL) / Phi) + sTL) * aT - P2 * ((Q1 * cTL - Q2 * sTL) / Phi) * aN)
    du[3] = eta * sqrt(sma / mu) * (sTL * aR + (((P2 + cTL) / Phi) + cTL) * aT + P1 * ((Q1 * cTL - Q2 * sTL) / Phi) * aN)
    du[4] = 0.5 * eta * sqrt(sma / mu) * (1 + Q1^2 + Q2^2) * sTL / Phi * aN
    du[5] = 0.5 * eta * sqrt(sma / mu) * (1 + Q1^2 + Q2^2) * cTL / Phi * aN    
    du[6] = meanmotion - (rnorm / hh) * ((sma / (sma + bb) * Phi * (P1 * sTL + P2 * cTL) + (2 * bb / sma)) * aR +
            sma / (sma + bb) * (1.0 + Phi) * (P1 * cTL - P2 * sTL) * aT +
            (Q1 * cTL - Q2 * sTL) * aN)

    # variation in time of quaternions
    du[7] = -0.5*(om[1]*quat[2]+om[2]*quat[3]+om[3]*quat[4]);
    du[8] =  0.5*(quat[1]*om[1]-om[2]*quat[4]+om[3]*quat[3]);
    du[9] =  0.5*(quat[1]*om[2]-om[3]*quat[2]+om[1]*quat[4]);
    du[10] = 0.5*(quat[1]*om[3]-om[1]*quat[3]+om[2]*quat[2]);

       
    # variation in time of components of angular velocity (Euler's equations)
    du[11] = (exttorque[1] - om[2]*om[3]*(C-B))/A;
    du[12] = (exttorque[2] - om[3]*om[1]*(A-C))/B;
    du[13] = (exttorque[3] - om[1]*om[2]*(B-A))/C;
        
    # println(exttorque,du[11:13])
    # println(t/86400)
    return du;
end

function orbitAttitudeEv_sadov(du,u,p,t)

    # if mod(floor(t/86400.0),2) == 0.0 
    #     println(t/86400)
    # end
    planetsunparams  = p[1];
    inertialrefframeinfo = p[2];
    satellite        = p[3];
    atmosphericinfo  = p[4];
    settings         = p[5];
    mjd2000_0        = p[6];
       
    # equinoctial elements
    mu = planetsunparams["muPlanet"]
    sma = u[1] 
    P1  = u[2]
    P2  = u[3]
    Q1  = u[4]
    Q2  = u[5]
    ML  = u[6] # mean anomaly
    
    eta = sqrt(1 - P1^2 - P2^2)
    bb = sma * eta;
    meanmotion = sqrt(mu / sma^3);
    hh = meanmotion * sma * bb;    

    # true longitude L = Omega + omega + nu (Battin, page. 493)
    TL = meanlong2truelong(ML,P1,P2);
    # if ier != 1
    #     Base.error("something went wrong in the solution of Kepler's equation")
    # end
    sTL = sin(TL); 
    cTL = cos(TL);     
    Phi = 1.0 + P1 * sTL + P2 * cTL


    # central body - centre of mass distance
    rnorm = sma*eta^2.0/Phi;
        
    # sadov c
    IV = get(satellite,"MomentsOfInertia",[0;0;0]);
    k  = get(satellite,"k_constant",0.0);
    k1r = sqrt(1+k);
    A = IV[1]; C = IV[3];
    Jg  = u[8];
    skm = u[7];
    smk = 1.0-skm;
    if skm == 0.0
        m = 0.0;
    else
        m = k*smk/skm;
    end
    if m==1.0
        KK = Inf;
        PP = Inf;
        PPMsmkKK = sqrt(k)*atan(sqrt(k))/(m+k);
    else
        KK  = ellF(pi/2,m); 
        PP  = ellP(-k,pi/2.0,m); 
        PPMsmkKK = (PP-smk*KK);
    end

    # moments of inertia
    IV = get(satellite,"MomentsOfInertia",[0;0;0]);
    A = IV[1]; B = IV[2]; C = IV[3];

    # perturbations   
    exttorque = zeros(3,1); 
    pertacc   = zeros(3,1);

    # perturbations flags

    # attitude perturbations
    includeGravityTorque  = settings["includeGravityTorque"];
    includeMagneticTorque = settings["includeMagneticTorque"];
    includeSrpTorque      = settings["includeSrpTorque"];
    includeDragTorque     = settings["includeDragTorque"];

    # orbit perturbations
    includeZonalHarmsAcc  = settings["includeZonalHarmsAcc"];   
    includeThirdBodyAcc   = settings["includeThirdBodyAcc"];
    includeSunGravityAcc  = settings["includeSunGravityAcc"];
    includeSrpAcc         = settings["includeSrpAcc"];
    includeDragAcc        = settings["includeDragAcc"];

    # eclipses effects
    includeEclipsesEffectsOnAttitude =  get(settings,"includeEclipsesEffectsOnAttitude",false);
    includeEclipsesEffectsOnOrbit    =  get(settings,"includeEclipsesEffectsOnOrbit",false);
    pertacc = zeros(3)
    torquepert = zeros(6)
    # perturbations
    if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque || includeSrpAcc || includeDragAcc || includeThirdBodyAcc || includeSunGravityAcc
        B = IV[2];
        skmR = sqrt(skm);
        smkR = sqrt(smk);

        psil = mod(u[10],2*pi);
        psig = mod(u[11],2*pi);
        eF  = 2.0*KK*psil/pi;
        lambda = Elliptic.Jacobi.am(eF,m);
        sn = cos(lambda);
        cn = sin(lambda)
        dn = sqrt(1.0-m*sn^2.0);
        eP  = ellP(-k,lambda,m) 
        eE  = ellE(lambda,m);
        EE  = ellE(pi/2.0,m);
        zn = eE-eF*EE/KK
        deltag  = k1r/skmR*(eF/KK*PP-eP);
        cgA = cos(deltag+psig);
        sgA = sin(deltag+psig);
        dnkR  = sqrt(1+k*sn^2);

        cdelta = u[9]/Jg;
        sdelta = sqrt(1-cdelta^2);
        cpsih = cos(u[12]);
        spsih = sin(u[12]);
        
        Rpsilpsig = zeros(3,3)
        Rpsilpsig[1,1] = (-sgA*cn*dn*skmR-cgA*sn*k1r)/dnkR;
        Rpsilpsig[1,2] = ( cgA*cn*dn*skmR-sgA*sn*k1r)/dnkR;
        Rpsilpsig[1,3] = smkR*cn
        Rpsilpsig[2,1] = ( sgA*sn*dn*skmR*k1r-cgA*cn)/dnkR;
        Rpsilpsig[2,2] = (-cgA*sn*dn*skmR*k1r-sgA*cn)/dnkR;
        Rpsilpsig[2,3] = -smkR*k1r*sn
        Rpsilpsig[3,1] = smkR*sgA*dnkR
        Rpsilpsig[3,2] = -smkR*cgA*dnkR
        Rpsilpsig[3,3] = skmR*dn   
        Sx = (dn*sn-cn*zn)/smkR;
        Sy = (dn*cn+sn*zn)/smkR/k1r;
        Sz = (dn*zn-m*sn*cn)/skmR;
        
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

        Ri2b = Rpsilpsig*Rhdelta;
        Ri2o = equi2rotmatirth(Q1,Q2,sTL,cTL);
        Ro2i = Transpose(Ri2o);
        rV = rnorm*Ro2i*[1,0.0,0.0]; 
        rUV = Ri2b*rV; rUV=rUV/rnorm;
        mjd2000   = mjd2000_0 + t/24.0/3600.0

        gtorque = [0.0;0.0;0.0];
        if includeGravityTorque
            gtorque = getgravitytorque(mu,rUV,rnorm,A,B,C) 
        end

        mtorque = [0.0;0.0;0.0];
        if includeMagneticTorque
            polaraxis = Ri2b*(inertialrefframeinfo["equatorial2inertial"]*[0.0,0.0,1.0]);
            mtorque = getmagnetictorque(planetsunparams["muM"],polaraxis,rUV,rnorm,satellite["intrinsicMagneticMoment"]);
        end

        dragtorque = [0.0;0.0;0.0];
        dragacc    = [0.0;0.0;0.0];
        if (includeDragTorque ||includeDragAcc) && rnorm-planetsunparams["rPlanet"]<1000.0
            vV =  sqrt(mu/sma/eta^2.0)*Ro2i*[P2*sTL-P1*cTL,P1*sTL+P2*cTL+1.0,0.0];
            atmAngVel = planetsunparams["planetRotation"]*inertialrefframeinfo["equatorial2inertial"]*[0.0,0.0,1.0];
            if atmosphericinfo["atmosphericmodel"] == 1
                dragtorque,dragforce = getdragtorqueandforce_atmexp(rnorm,atmosphericinfo,planetsunparams["rPlanet"],satellite["CD"],satellite["facets"],satellite["numberOfFacets"],atmAngVel,vV,rV,Ri2b,om);
            else
                dragtorque,dragforce = getdragtorqueandforce_atmnrlmsise00(atmosphericinfo,satellite["facets"],satellite["numberOfFacets"],rV,vV,atmAngVel,Ri2b,transpose(inertialrefframeinfo["equatorial2inertial"]),mjd2000);
            end
            if includeDragAcc
                dragacc = Transpose(Ri2b)*dragforce/1000/satellite["mass"];
             end
             if includeDragTorque==false
                dragtorque = [0.0;0.0;0.0];
             end
        
        end

        srptorque = [0.0;0.0;0.0];
        srpacc    = [0.0;0.0;0.0];
        sungacc = [0.0;0.0;0.0];
        if includeSrpTorque || includeSrpAcc || includeSunGravityAcc
            rSunV = -inertialrefframeinfo["ecliptic2inertial"]*celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000)

            if includeSrpTorque || includeSrpAcc
                rsunnorm = sqrt(rSunV[1]^2.0+rSunV[2]^2.0+rSunV[3]^2.0)
                srptorque,srpforce = getsrptorqueandforce(rSunV,rsunnorm,rV,rnorm,Ri2b,satellite["facets"],satellite["numberOfFacets"],planetsunparams["rPlanet"],planetsunparams["psrp_refdist"],includeEclipsesEffectsOnAttitude,includeEclipsesEffectsOnOrbit);   
                if includeSrpAcc
                    srpacc = Transpose(Ri2b)*srpforce/1000/satellite["mass"];
                end
                if includeSrpTorque==false
                    srptorque = [0.0;0.0;0.0];
                end
            end
            
            if includeSunGravityAcc
                sungacc = getSunGravityPert(rSunV,rV);
            end
        end
   
        tbgacc = [0.0;0.0;0.0]
        if includeThirdBodyAcc 
            for jj = 1 : size(planetsunparams["perturbingbodiesid"])[1]
                idx3b = planetsunparams["perturbingbodiesid"][jj];
                if idx3b == 11
                    if planetsunparams["centralBodyIDX"] == 3
                        rV_3b = inertialrefframeinfo["equatorial2inertial"]*celestialbodiesephemeris_position(idx3b,mjd2000);
                    else
                        continue;
                    end
                else
                    rV_3bWRTS = celestialbodiesephemeris_position(idx3b,mjd2000);
                    rV_CBWRTS = celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000);
                    rV_3b = inertialrefframeinfo["ecliptic2inertial"]*(rV_3bWRTS-rV_CBWRTS);
                end
                tbgacc = tbgacc + get3BodyAcc(rV,planetsunparams["perturbingbodiesgravparam"][jj],rV_3b);
            end
        end

        pertacc = pertacc + Ri2o*(tbgacc+sungacc+dragacc+srpacc);

        exttorque = exttorque + gtorque+mtorque+srptorque+dragtorque;  
        exttorque[abs.(exttorque).<1e-15].=0.0;     

        mxb11 = exttorque[1]*Rpsilpsig[1,1];
        myb21 = exttorque[2]*Rpsilpsig[2,1];
        mzb31 = exttorque[3]*Rpsilpsig[3,1];
        mxb12 = exttorque[1]*Rpsilpsig[1,2];
        myb22 = exttorque[2]*Rpsilpsig[2,2];
        mzb32 = exttorque[3]*Rpsilpsig[3,2];
        mxb13 = exttorque[1]*Rpsilpsig[1,3];        
        myb23 = exttorque[2]*Rpsilpsig[2,3];       
        mzb33 = exttorque[3]*Rpsilpsig[3,3];
        mxSx = exttorque[1]*Sx
        mySy = exttorque[2]*Sy
        mzSz = exttorque[3]*Sz
    
        if m==0.0
            T4 = -pi/2/sqrt(1+k);
        elseif m==1.0
            T4 = -atan(sqrt(k))/sqrt(k);
        else
            T4 = (KK*m - PP*(m+k))/k;
        end
        torquepert  = zeros(6);

        torquepert[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
    
        torquepert[2] = mxb13+myb23+mzb33;

        torquepert[3] = cdelta*torquepert[2] + sdelta*(mxb12+myb22+mzb32);

        torquepert[4] = -pi/(2*Jg*KK*(1-m))*(mxSx+(1.0-m)*mySy+mzSz)
       
        torquepert[6] = (mxb11+myb21+mzb31)/Jg/sdelta;
            
        torquepert[5] = -cdelta*torquepert[6] + 2/pi*sqrt(skm)*T4*k1r*torquepert[4] ; 
  

    end

    if includeZonalHarmsAcc
        J2 = planetsunparams["zonalharmonicscoff"][1]
        J3 = planetsunparams["zonalharmonicscoff"][2]
        J4 = planetsunparams["zonalharmonicscoff"][3]
        J5 = planetsunparams["zonalharmonicscoff"][4]
        potpert = getCentralBodyPotentialAcc(sma,eta,Q1,Q2,cTL,sTL,Phi,mu,planetsunparams["rPlanet"],J2,J3,J4,J5)
        pertacc = pertacc + potpert;
    end  

    # variation in time of equinoctial elements
    aR  = pertacc[1];
    aT  = pertacc[2]; 
    aN  = pertacc[3];   

    du[1] = (2 / eta) * sqrt(sma^3 / mu) * ((P2 * sTL - P1 * cTL) * aR + Phi * aT)
    du[2] = eta * sqrt(sma / mu) * (-cTL * aR + (((P1 + sTL) / Phi) + sTL) * aT - P2 * ((Q1 * cTL - Q2 * sTL) / Phi) * aN)
    du[3] = eta * sqrt(sma / mu) * (sTL * aR + (((P2 + cTL) / Phi) + cTL) * aT + P1 * ((Q1 * cTL - Q2 * sTL) / Phi) * aN)
    du[4] = 0.5 * eta * sqrt(sma / mu) * (1 + Q1^2 + Q2^2) * sTL / Phi * aN
    du[5] = 0.5 * eta * sqrt(sma / mu) * (1 + Q1^2 + Q2^2) * cTL / Phi * aN    
    du[6] = meanmotion - (rnorm / hh) * ((sma / (sma + bb) * Phi * (P1 * sTL + P2 * cTL) + (2 * bb / sma)) * aR +
            sma / (sma + bb) * (1.0 + Phi) * (P1 * cTL - P2 * sTL) * aT +
            (Q1 * cTL - Q2 * sTL) * aN)

    # variation in time of sadov        
    du[7] = torquepert[1]
    du[8] = torquepert[2]
    du[9] = torquepert[3]
    du[10] = -pi/k1r*sqrt(skm)*Jg/2/A/C*(C-A)/KK + torquepert[4]
    du[11] = Jg/A/C*(C-A)*PPMsmkKK/KK + Jg/A/C*(C*smk+A*skm) + torquepert[5]
    du[12] = torquepert[6]
        
    # println(exttorque,du[11:13])
    return du;
end

###################################
# torque and orbital perturbations

function getgravitytorque(mu,rUV,rnorm,A,B,C)
    gtorque = [0.0,0.0,0.0];
    gtorque[1] = 3.0*mu/(rnorm^3)*(C-B)*rUV[2]*rUV[3];
    gtorque[2] = 3.0*mu/(rnorm^3)*(A-C)*rUV[3]*rUV[1];
    gtorque[3] = 3.0*mu/(rnorm^3)*(B-A)*rUV[1]*rUV[2];
    return gtorque;
end

function getmagnetictorque(muM,Zeq,rUV,rnorm,IM)
    # println(muM," ",Zeq," ",rUV)
    Hb  = (Zeq - 3.0*(Zeq[1]*rUV[1]+Zeq[2]*rUV[2]+Zeq[3]*rUV[3])*rUV)*muM/rnorm^3;
    mtorque = LinearAlgebra.cross(IM,Hb);
    return mtorque;
end

function getsrptorqueandforce(rVSun,rsunnorm,rV,rnorm,Ri2b,facets,numberOfFacets,rPlanet,refdist,includeEclipsesEffectsOnAttitude,includeEclipsesEffectsOnOrbit)
    srpforce  = [0.0,0.0,0.0];
    srptorque = [0.0,0.0,0.0];
    shadowfunattitude = 1.0;
    shadowfunorbit    = 1.0;

    if includeEclipsesEffectsOnAttitude || includeEclipsesEffectsOnOrbit
        sc = (rVSun[1]*rV[1]+rVSun[2]*rV[2]+rVSun[3]*rV[3])/rsunnorm+sqrt(rnorm^2-rPlanet^2);
        shadowfun = 1/2*(1+tanh(1000*sc));
        if includeEclipsesEffectsOnAttitude 
            shadowfunattitude = shadowfun
        end
        if includeEclipsesEffectsOnOrbit 
            shadowfunorbit = shadowfun
        end
    end    


    uVect = Ri2b*(rVSun-rV);
    ssd = norm(uVect); 
    uVect = uVect/ssd;
    for kk = 1:numberOfFacets 
        surfi = facets[kk][2];
        niv = facets[kk][4];
        rhoiv = facets[kk][3][1];
        niudot = niv[1]*uVect[1]+niv[2]*uVect[2]+niv[3]*uVect[3];
        gifun = max(0,niudot);
        
        cai = facets[kk][1][1]*(refdist/ssd)^2;
        cdi = facets[kk][1][2]*(refdist/ssd)^2;
        csi = facets[kk][1][3]*(refdist/ssd)^2;

        srpforcekk_T1 = - (surfi*gifun*cai*uVect);
        srpforcekk_T2 = - (surfi*gifun*cdi*niv);
        srpforcekk_T3 = - (surfi*gifun*csi*niudot*niv)
        
        srpforce = srpforce + (srpforcekk_T1 + srpforcekk_T2 + srpforcekk_T3) ;

        srptorquekk_T1 = LinearAlgebra.cross(rhoiv,srpforcekk_T1);
        srptorquekk_T2 = LinearAlgebra.cross(rhoiv,srpforcekk_T2);
        srptorquekk_T3 = LinearAlgebra.cross(rhoiv,srpforcekk_T3);
        srptorquekk = (srptorquekk_T1 + srptorquekk_T2 + srptorquekk_T3); 
        srptorque = srptorque + srptorquekk;
    end        
    srptorque = srptorque*shadowfunattitude;
    srpforce  = srpforce*shadowfunorbit;
    
    return srptorque,srpforce;
end

function getdragtorqueandforceVC_atmexp(rnorm,atmosphericinfo,rPlanet,CD,facets,numerOfFacets,nEV,vV,rV,Ri2b,omV)
    dragtorque = [0.0,0.0,0.0];
    dragforce    = [0.0,0.0,0.0];
   
    density =  exponential_atm_model(rnorm-rPlanet,atmosphericinfo["AtmM"])
    
    V0V = Ri2b*(vV-LinearAlgebra.cross(nEV,rV))*1000.0; # km2m
    eV0 = V0V/V0;

    for kk = 1:numerOfFacets
        surfi = facets[kk][2];
        niv = facets[kk][4];
        rhoiv = facets[kk][3][1];
        vv1   = facets[kk][3][2];
        vv2   = facets[kk][3][3];
        vv3   = facets[kk][3][4];

        niveV0 = eV0[1]*niv[1]+eV0[2]*niv[2]+eV0[3]*niv[3];
        g0 = 1/3/pi + niveV0/2 + 4*(niveV0)^2/3/pi;
        g2 = 1/2*(1+16/3/pi*niveV0);
        gV = 2*g0*eV0+g2*niv;
        gVComV = LinearAlgebra.cross(gV,omV);
        
        fi1 = -1/2*CD*density*surfi*V0^2*g0*eV0
        fi2 = -1/2*CD*density*surfi*V0*LinearAlgebra.dot(gVComV,rhoiv)*eV0;
        fi3 = -1/2*CD*density*surfi*V0*g0*LinearAlgebra.cross(omV,rhoiv);
    
        mi1 = LinearAlgebra.cross(rhoiv,fi1); 
        mi2 = 3/4*LinearAlgebra.cross(rhoiv,fi2)-1/24*CD*density*surfi*V0*(vv1*LinearAlgebra.dot(vv1,gVComV)+vv2*LinearAlgebra.dot(vv2,gVComV)+vv2*LinearAlgebra.dot(vv3,gVComV));
        mi3 = 3/4*LinearAlgebra.cross(rhoiv,fi3)-1/24*CD*density*surfi*V0*(LinearAlgebra.cross(vv1,LinearAlgebra.cross(omV,vv1))+LinearAlgebra.cross(vv2,LinearAlgebra.cross(omV,vv2))+LinearAlgebra.cross(vv3,LinearAlgebra.cross(omV,vv3)));
       
        dragtorque = dragtorque + mi1 + mi2 + mi3;
        dragforce    = dragforce + fi1 + fi2 + fi3;
    end      
    return dragtorque,dragforce;
end

function getdragtorqueandforce_atmexp(rnorm,atmosphericinfo,rPlanet,CD,facets,numerOfFacets,nEV,vV,rV,Ri2b,omV)
    dragtorque = [0.0,0.0,0.0];
    dragforce    = [0.0,0.0,0.0];
    density =  exponential_atm_model(rnorm-rPlanet,atmosphericinfo["AtmM"])

    V0V = Ri2b*(vV-LinearAlgebra.cross(nEV,rV))*1000.0; # km2m
    V0 = norm(V0V);
    eV0 = V0V/V0;

    for kk = 1:numerOfFacets
        surfi = facets[kk][2];
        niv = facets[kk][4];
        rhoiv = facets[kk][3][1];
              
        niveV0 = eV0[1]*niv[1]+eV0[2]*niv[2]+eV0[3]*niv[3];
        g0 = 1/3/pi + niveV0/2 + 4*(niveV0)^2/3/pi;
       
        fi1 = -1/2*CD*density*surfi*V0^2*g0*eV0
    
        mi1 = LinearAlgebra.cross(rhoiv,fi1); 
        dragtorque = dragtorque + mi1;
        dragforce    = dragforce + fi1;
    end        

    return dragtorque,dragforce;
end

function getdragtorqueandforce_atmnrlmsise00(atmosphere,facets,numerOfFacets,rV,vV,nEV,Ri2b,Ri2EECI,mjd2000)
    
    dragtorque = zeros(3);
    dragforce  = zeros(3);
    
    rhoinf,Tinf,ainf,ndens,Te = getatmosphereproperties_nrlmsise00(mjd2000,Ri2EECI*rV,atmosphere["AtmM"]);

    #[kg/m^3,K,m/s]
    
    VVinf = (vV-cross(nEV,rV))*1000.0; # km2m
    Vinf  = norm(VVinf);
    eVinf = VVinf/Vinf;
    eVinf = Ri2b*eVinf;
    
    fV = zeros(3)
    fN = zeros(3)
    mV = zeros(3)
    mN = zeros(3)

    Sinf   = Vinf/ainf;
    Tratio = sqrt(300.0/Tinf);
   
    
    bsi0 = besseli(0.0,Sinf^2.0/2.0);
    bsi1 = besseli(1.0,Sinf^2.0/2.0);
    expSM = exp(-Sinf^2.0/2.0)
    expS  = exp(-Sinf^2.0)
    pi32   = pi^(3.0/2.0)
    hyp1  = pFq((1/2,1),(3/2,3/2),-Sinf^2.0)
    hyp2  = pFq((1/2,2),(3/2,5/2),-Sinf^2.0)
    hyp3  = pFq((1/2,3),(3/2,7/2),-Sinf^2.0)
    
    ff1a0 = expSM*((Sinf^2.0+1.0)*bsi0+Sinf^2.0*bsi1)/(2.0*sqrt(pi)*Sinf) 
    ff1a1 = (8.0*Sinf^3.0*hyp2 + 3.0*expS*erfi(Sinf)*sqrt(pi))/(3.0*Sinf^2.0*pi32)
    ff1a2 = expSM*(Sinf^2.0*bsi0+(Sinf^2.0-1.0)*bsi1)/(3.0*sqrt(pi)*Sinf)
    ff1a3 = -8.0*(Sinf^5.0*hyp2 - 16.0*Sinf^5.0*hyp3/15.0 - erfi(Sinf)*sqrt(pi)*(Sinf^2.0+2.0)*expS/8.0 + Sinf/2.0)/(Sinf^4.0*pi32)
    ff1a4 = -expSM*((Sinf^4.0-4.0*Sinf^2.0+12.0)*bsi1+ Sinf^2.0*(Sinf^2.0-3.0)*bsi0)/(15.0*sqrt(pi)*Sinf^3.0)

    ff2a0 = (pi32+4.0*hyp1*Sinf)/(4.0*pi32*Sinf^2.0)
    ff2a1 = (pi*Sinf*expSM*bsi0+pi*Sinf*expSM*bsi1+2.0*sqrt(pi))/(2.0*Sinf^2.0*pi32)
    ff2a2 = (8.0*hyp2-6.0*hyp1)/(3.0*pi32*Sinf)
    ff2a3 = -(pi*expSM*(Sinf^2.0-4.0)*bsi1+Sinf*(pi*Sinf*expSM*bsi0+2.0*sqrt(pi)))/(6.0*Sinf^3.0*pi32)
    ff2a4 = 2.0*(64.0*hyp3-80.0*hyp2+15.0*hyp1)/(15.0*Sinf*pi32)

    A0  = 1.0/3.0/pi + (ff1a0 - ff1a2 + ff1a4);
    A1  = 1/2.0 + (ff1a1 - 3.0*ff1a3);
    A2  = 4.0/pi/3.0 + (2.0*ff1a2 - 8.0*ff1a4);
    A3  = 4.0*ff1a3;
    A4  = 8.0*ff1a4;
 
    B0  = (ff2a0 - ff2a2 + ff2a4);
    B1  = (ff2a1 - 3.0*ff2a3);
    B2  = (2.0*ff2a2 - 8.0*ff2a4);
    B3  = 4.0*ff2a3;
    B4  = 8.0*ff2a4; 


    for kk = 1:numerOfFacets
        surfi = facets[kk][2];
        niv = facets[kk][4];
        rhoiv = facets[kk][3][1];

        ctheta = LinearAlgebra.dot(niv,eVinf);
        CV = A0 + A1*ctheta + A2*ctheta^2 + A3*ctheta^3 + A4*ctheta^4;
        CN = B0 + B1*ctheta + B2*ctheta^2 + B3*ctheta^3 + B4*ctheta^4 + sqrt(pi)/2.0/Sinf*Tratio*CV;
        # CV,CN  = getdragcoefficients(ctheta,Vinf,ainf,Tinf,300.0);

        fiV = -0.5*rhoinf*Vinf^2.0*surfi*CV*eVinf;
        fiN = -0.5*rhoinf*Vinf^2.0*surfi*CN*niv;
        miV = LinearAlgebra.cross(rhoiv,fiV);
        miN = LinearAlgebra.cross(rhoiv,fiN);

        # cnsn = -niv*surfi*(CN+CV*ctheta)
        # ctst = -(eVinf-niv*ctheta)*surfi*CV
        # println([mjd20002jd(mjd2000),cnsn,ctst])
   
        fV = fV + fiV;
        fN = fN + fiN
        mV = mV + miV
        mN = mN + miN
    end
    dragtorque = mV + mN;
    dragforce = fV + fN

    return dragtorque,dragforce;

end

function getdragcoefficients(ctheta,Vinf,ainf,Tinf,TW)    

    Sinf   = Vinf/ainf;
    Tratio = sqrt(TW/Tinf);
    Sn     = Sinf*ctheta;
    erfSn  = SpecialFunctions.erf(Sn);
    delta  = max(ctheta/abs(ctheta),0.0);
    sfun   = (1+tanh(5.0*ctheta))/2.0;
    
    CV  = (1.0+erfSn)*max(ctheta,0.0) + exp(-Sn^2)/sqrt(pi)/Sinf*sfun;
    CN  = (1.0+erfSn)/2/(Sinf^2)*sfun + sqrt(pi)/2.0/Sinf*CV*Tratio;

  
    return CV,CN;


    # Sinf   = Vinf/ainf;
    # Tratio = sqrt(300.0/Tinf);
    # thetaV   = collect(range(0.0,pi/2,180)); 
    # cthetaV  = cos.(thetaV);
    # SnV      = Sinf*cthetaV;
    # ff2   = (1.0 .+  SpecialFunctions.erf.(SnV))/2.0/(Sinf^2.0); 
    # ff2a0 = 1.0/pi*trapezoidalrule_V2(ff2,thetaV);
    # ff2aV = zeros(5);
    # ff2aV[1] = ff2a0
    # for ii = 1 : size(ff2aV)[1]-1
    #     ff2aV[ii+1] = 2.0/pi*trapezoidalrule_V2(ff2.*cos.(thetaV*ii),thetaV);
    # end
    # ff1   = (cthetaV + SpecialFunctions.erf.(SnV).*cthetaV + exp.(-SnV.^2.0)/Sinf/sqrt(pi));
    # ff1a0 = 1.0/pi*trapezoidalrule_V2(ff1,thetaV);
    # ff1aV = zeros(5);
    # ff1aV[1] = ff1a0
    # for ii = 1 : size(ff1aV)[1]-1
    #     ff1aV[ii+1] = 2.0/pi*trapezoidalrule_V2(ff1.*cos.(thetaV*ii),thetaV);
    # end
    # CV = Polynomials.evalpoly(ctheta,Polynomials.ChebyshevT(ff1aV));
    # CN = Polynomials.evalpoly(ctheta,Polynomials.ChebyshevT(ff2aV)) + sqrt(pi)/2.0/Sinf*Tratio*CV;

    ##############################################################################
    # ff2a1 = 2.0/pi*trapezoidalrule_V2(ff2.*cthetaV,thetaV);
    # ff2a2 = 2.0/pi*trapezoidalrule_V2(ff2.*cos.(thetaV*2.0),thetaV);
    # ff2a3 = 2.0/pi*trapezoidalrule_V2(ff2.*cos.(thetaV*3.0),thetaV);
    # ff2a4 = 2.0/pi*trapezoidalrule_V2(ff2.*cos.(thetaV*4.0),thetaV);
    # B0  = (ff2a0 - ff2a2 + ff2a4);
    # B1  = (ff2a1 - 3.0*ff2a3);
    # B2  = (2.0*ff2a2 - 8.0*ff2a4);
    # B3  = 4.0*ff2a3;
    # B4  = 8.0*ff2a4; 

    # ff1   = (cthetaV + SpecialFunctions.erf.(SnV).*cthetaV + exp.(-SnV.^2.0)/Sinf/sqrt(pi));
    # ff1a0 = 1.0/pi*trapezoidalrule_V2(ff1,thetaV);
    # ff1a1 = 2.0/pi*trapezoidalrule_V2(ff1.*cthetaV,thetaV);
    # ff1a2 = 2.0/pi*trapezoidalrule_V2(ff1.*cos.(thetaV*2.0),thetaV);
    # ff1a3 = 2.0/pi*trapezoidalrule_V2(ff1.*cos.(thetaV*3.0),thetaV);
    # ff1a4 = 2.0/pi*trapezoidalrule_V2(ff1.*cos.(thetaV*4.0),thetaV);
    # A0  = (ff1a0 - ff1a2 + ff1a4);
    # A1  = (ff1a1 - 3.0*ff1a3);
    # A2  = (2.0*ff1a2 - 8.0*ff1a4);
    # A3  = 4.0*ff1a3;
    # A4  = 8.0*ff1a4;

    # CV = A0 + A1*ctheta + A2*ctheta^2 + A3*ctheta^3 + A4*ctheta^4;
    # CN = B0 + B1*ctheta + B2*ctheta^2 + B3*ctheta^3 + B4*ctheta^4 + sqrt(pi)/2.0/Sinf*Tratio*CV;
    
    

    # ### prova
    # ll       = 181;
    # thetaV       = collect(range(0.0,pi/2,ll)); 
    # cthetaV  = cos.(thetaV);
    # SnV      = Sinf*cthetaV;
    # ff1      = (SpecialFunctions.erf.(SnV).*cthetaV + exp.(-SnV.^2.0)/sqrt(pi)/Sinf);
    # ff1a0    = 1/pi*trapezoidalrule_V2(ff1,thetaV);
    # ff1a1    = 2/pi*trapezoidalrule_V2(ff1.*cthetaV,thetaV);
    # ff1a2    = 2/pi*trapezoidalrule_V2(ff1.*cos.(2.0*thetaV),thetaV);
    # ff1a3    = 2/pi*trapezoidalrule_V2(ff1.*cos.(3.0*thetaV),thetaV);
    # ff1a4    = 2/pi*trapezoidalrule_V2(ff1.*cos.(4.0*thetaV),thetaV);
    # ff1a5    = 2/pi*trapezoidalrule_V2(ff1.*cos.(5.0*thetaV),thetaV);
    # ff1a6    = 2/pi*trapezoidalrule_V2(ff1.*cos.(6.0*thetaV),thetaV);
    # ff1a7    = 2/pi*trapezoidalrule_V2(ff1.*cos.(7.0*thetaV),thetaV);

    
    # A0  = ff1a0 - ff1a2 + ff1a4 - ff1a6;
    # A1  = ff1a1 - 3.0*ff1a3 + 5.0*ff1a5 - 7.0*ff1a7;
    # A2  = 2.0*ff1a2  - 8.0*ff1a4  + 18.0*ff1a6;
    # A3  = 4.0*ff1a3  - 20.0*ff1a5 + 56.0*ff1a7;
    # A4  = 8.0*ff1a4  - 48.0*ff1a6 
    # A5  = 16.0*ff1a5 - 112.0*ff1a7
    # A6  = 32.0*ff1a6
    # A7  = 64.0*ff1a7

    # CV  = ctheta + A0 + A1*ctheta + A2*ctheta^2 + A3*ctheta^3 + A4*ctheta^4 + A5 * ctheta^5 + A6*ctheta^6 + A7*ctheta^7;
    # CN  = (1+erfSn)/2/(Sinf^2)*delta + sqrt(pi)/2/Sinf*CV*Tratio;

    #     # B0  = ff2a0 - ff2a2 + ff2a4 - ff2a6;
    #     # B1  = ff2a1 - 3.0*ff2a3 + 5.0*ff2a5 - 7.0*ff2a7;
    #     # B2  = 2.0*ff2a2 - 8.0*ff2a4 + 18.0*ff2a6;
    #     # B3  = 4.0*ff2a3 - 20.0*ff2a5 + 56.0*ff2a7;
    #     # B4  = 8.0*ff2a4 - 48.0*ff2a6 
    #     # B5  = 16.0*ff2a5 - 112.0*ff2a7
    #     # B6  = 32.0*ff2a6
    #     # B7  = 64.0*ff2a7

end

function getCentralBodyPotentialAcc(sma,eta,Q1,Q2,cTL,sTL,Phi,mu,rPlanet,J2,J3,J4,J5)
    #Eq of Motion: derivative of the equinoctial elements wrt. to time

    Psi = Q2 * sTL - Q1 * cTL
    G = 1 + Q1^2 + Q2^2

    pert = zeros(3,1);    

    #Accelerations J2
    if J2 != 0.0
        aR_J2 = 3.0  * mu * J2 * rPlanet^2.0 / (2.0 * eta^8.0 * sma^4.0 * G^2.0) * (12.0 * (Q1 * cTL - Q2 * sTL)^2.0 - G^2.0) * Phi^4.0;
        aT_J2 = 12.0 * mu * J2 * rPlanet^2.0 / (eta^8.0 * sma^4.0 * G^2.0) * (Q2 * cTL + Q1 * sTL) * (Q1 * cTL - Q2 * sTL) * Phi^4.0;
        aN_J2 = 6.0  * mu * J2 * rPlanet^2.0 / (eta^8.0 * sma^4.0 * G^2.0) * (Q1 * cTL - Q2 * sTL) * (1.0 - Q1^2.0 - Q2^2.0) * Phi^4.0;
        pert[1] = pert[1] + aR_J2
        pert[2] = pert[2] + aT_J2
        pert[3] = pert[3] + aN_J2
    end

    #==================================#
    #Accelerations J3
    if J3 != 0.0
        aR_J3 = 4.0 * mu * J3 * rPlanet^3.0 * Phi^5.0 / (eta^10.0 * sma^5.0 * G) * (Q2 * sTL - Q1 * cTL) * (20.0/G^2.0 * (Q2 * sTL - Q1 * cTL)^2.0 - 3.0);
        aT_J3 = 3.0 * mu * J3 * rPlanet^3.0 * Phi^5.0 / (eta^10.0 * sma^5.0 * G) * (Q2 * cTL + Q1 * sTL) * (1 - 20.0/G^2 * (Q2 * sTL - Q1 * cTL)^2.0);
        aN_J3 = 1.5 * mu * J3 * rPlanet^3.0 * Phi^5.0 / (eta^10.0 * sma^5.0 * G) * (1.0 - Q1^2.0 - Q2^2.0) * (1.0 - 20.0/G^2.0 * (Q2 * sTL - Q1 * cTL)^2.0);
        pert[1] = pert[1] + aR_J3
        pert[2] = pert[2] + aT_J3
        pert[3] = pert[3] + aN_J3
    end

    #==================================#
    #Accelerations J4
    if J4 != 0.0
        k_J4 = mu*J4*rPlanet^4.0/(eta^12.0*sma^6.0) * Phi^6.0;
        aR_J4 =  k_J4/8 * (2800/G^4*Psi^4 - 600/G^2*Psi^2 + 15);
        aT_J4 = -k_J4/G^2 * (280/G^2*Psi^2 - 30)*(Q2*cTL+Q1*sTL)*Psi;
        aN_J4 = -k_J4/G^2 * (1-Q1^2-Q2^2)*(140/G^2*Psi^2-15)*Psi;
        pert[1] = pert[1] + aR_J4
        pert[2] = pert[2] + aT_J4
        pert[3] = pert[3] + aN_J4
    end

    #==================================#
    #Accelerations J5
    if J5 != 0.0
        k_J5 = mu*J5*rPlanet^5/(G*eta^14*sma^7)*Phi^7;
        aR_J5 =  k_J5/2*Psi*(3024/G^4*Psi^4 - 840/G^2*Psi^2 + 45);
        aT_J5 = -k_J5/4*(Q2*cTL+Q1*sTL)*(5040/G^4*Psi^4 - 840/G^2*Psi^2 + 15);
        aN_J5 = -k_J5/8* (1-Q1^2-Q2^2) * (5040/G^4*Psi^4 - 840/G^2*Psi^2 + 15);
        pert[1] = pert[1] + aR_J5
        pert[2] = pert[2] + aT_J5
        pert[3] = pert[3] + aN_J5
    end

   return pert;    
end

function get3BodyAcc(rV,mu_3b,rV_3b)
    a3bp = mu_3b * ( ( rV_3b - rV) / norm(rV_3b - rV)^3.0 - rV_3b / norm(rV_3b)^3.0);
    return a3bp
end

function getSunGravityPert(rSunV,rV)
    mu_Sun = 0.19891000000000E+31 * 6.67259e-20;
    aSun = mu_Sun * ( (rSunV - rV)/ norm(rSunV - rV)^3.0 - rSunV / norm(rSunV)^3.0);
    return aSun;
end





# function orbitAttitudeEv_quatVC(du,u,p,t)

#     planetsunparams  = p[1];
#     inertialrefframeinfo = p[2];
#     satellite        = p[3];
#     atmosphericinfo  = p[4];
#     settings         = p[5];
#     mjd2000_0        = p[6];
       
#     # equinoctial elements
#     mu = planetsunparams["muPlanet"]
#     sma = u[1] 
#     P1  = u[2]
#     P2  = u[3]
#     Q1  = u[4]
#     Q2  = u[5]
#     ML  = u[6] # mean anomaly
    
#     eta = sqrt(1 - P1^2 - P2^2)
#     bb = sma * eta;
#     meanmotion = sqrt(mu / sma^3);
#     hh = meanmotion * sma * bb;    

#     # true longitude L = Omega + omega + nu (Battin, page. 493)
#     TL = meanlong2truelong(ML,P1,P2);
#     # if ier != 1
#     #     Base.error("something went wrong in the solution of Kepler's equation")
#     # end
#     sTL = sin(TL); 
#     cTL = cos(TL);     
#     Phi = 1 + P1 * sTL + P2 * cTL


#     # central body - centre of mass distance
#     rnorm = sma*eta^2.0/Phi;
        
#     # quaternion update
#     quat = u[7:10];

#     # angular velocity
#     om  = u[11:13];
 
#     # moments of inertia
#     IV = get(satellite,"MomentsOfInertia",[0;0;0]);
#     A = IV[1]; B = IV[2]; C = IV[3];

#     # perturbations   
#     exttorque = zeros(3,1); 
#     pertacc   = zeros(3,1);

#     # perturbations flags

#     # attitude perturbations
#     includeGravityTorque  = settings["includeGravityTorque"];
#     includeMagneticTorque = settings["includeMagneticTorque"];
#     includeSrpTorque      = settings["includeSrpTorque"];
#     includeDragTorque     = settings["includeDragTorque"];

#     # orbit perturbations
#     includeZonalHarmsAcc  = settings["includeZonalHarmsAcc"];   
#     includeThirdBodyAcc   = settings["includeThirdBodyAcc"];
#     includeSunGravityAcc  = settings["includeSunGravityAcc"];
#     includeSrpAcc         = settings["includeSrpAcc"];
#     includeDragAcc        = settings["includeDragAcc"];

#     # eclipses effects
#     includeEclipsesEffectsOnAttitude =  get(settings,"includeEclipsesEffectsOnAttitude",false);
#     includeEclipsesEffectsOnOrbit    =  get(settings,"includeEclipsesEffectsOnOrbit",false);

#     # perturbations
#     if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque || includeSrpAcc || includeDragAcc || includeThirdBodyAcc || includeSunGravityAcc
#         Ri2b = transpose(quat2birotmat(quat));
#         Ri2o = equi2rotmatirth(Q1,Q2,sTL,cTL);
#         Ro2i = Transpose(Ri2o);
#         rV = rnorm*Ro2i*[1,0.0,0.0]; 
#         rUV = Ri2b*rV; rUV=rUV/rnorm;
#         mjd2000   = mjd2000_0 + t/24.0/3600.0

#         gtorque = [0.0;0.0;0.0];
#         if includeGravityTorque
#             gtorque = getgravitytorque(mu,rUV,rnorm,A,B,C) 
#         end

#         mtorque = [0.0;0.0;0.0];
#         if includeMagneticTorque
#             polaraxis = Ri2b*(inertialrefframeinfo["equatorial2inertial"]*[0.0,0.0,1.0]);
#             mtorque = getmagnetictorque(planetsunparams["muM"],polaraxis,rUV,rnorm,satellite["intrinsicMagneticMoment"]);
#         end

#         dragtorque = [0.0;0.0;0.0];
#         dragacc    = [0.0;0.0;0.0];
#         if (includeDragTorque ||includeDragAcc) #&& rnorm-planetsunparams["rPlanet"]<1000.0
#             vV =  sqrt(mu/sma/eta^2.0)*Ro2i*[P2*sTL-P1*cTL,P1*sTL+P2*cTL+1.0,0.0];
#             atmAngVel = planetsunparams["planetRotation"]*inertialrefframeinfo["equatorial2inertial"]*[0.0,0.0,1.0];
#             if atmosphericinfo["atmosphericmodel"] == 1
#                 dragtorque,dragforce = getdragtorqueandforceVC_atmexp(rnorm,atmosphericinfo,planetsunparams["rPlanet"],satellite["CD"],satellite["facets"],satellite["numberOfFacets"],atmAngVel,vV,rV,Ri2b,om);
#             else
#                 dragtorque,dragforce = getdragtorqueandforce_atmnrlmsise00(atmosphericinfo,satellite["facets"],satellite["numberOfFacets"],rV,vV,atmAngVel,Ri2b,transpose(inertialrefframeinfo["equatorial2inertial"]),mjd2000);
#             end
#             if includeDragAcc
#                 dragacc = Transpose(Ri2b)*dragforce/1000/satellite["mass"];
#              end
#              if includeDragTorque==false
#                 dragtorque = [0.0;0.0;0.0];
#              end
        
#         end

#         srptorque = [0.0;0.0;0.0];
#         srpacc    = [0.0;0.0;0.0];
#         sungacc = [0.0;0.0;0.0];
#         if includeSrpTorque || includeSrpAcc || includeSunGravityAcc
#             rSunV = -inertialrefframeinfo["ecliptic2inertial"]*celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000)

#             if includeSrpTorque || includeSrpAcc
#                 rsunnorm = sqrt(rSunV[1]^2.0+rSunV[2]^2.0+rSunV[3]^2.0)
#                 srptorque,srpforce = getsrptorqueandforce(rSunV,rsunnorm,rV,rnorm,Ri2b,satellite["facets"],satellite["numberOfFacets"],planetsunparams["rPlanet"],planetsunparams["psrp_refdist"],includeEclipsesEffectsOnAttitude,includeEclipsesEffectsOnOrbit);   
#                 if includeSrpAcc
#                     srpacc = Transpose(Ri2b)*srpforce/1000/satellite["mass"];
#                 end
#                 if includeSrpTorque==false
#                     srptorque = [0.0;0.0;0.0];
#                 end
#             end
            
#             if includeSunGravityAcc
#                 sungacc = getSunGravityPert(rSunV,rV);
#             end
#         end
   
#         tbgacc = [0.0;0.0;0.0]
#         if includeThirdBodyAcc 
#             for jj = 1 : size(planetsunparams["perturbingbodiesid"])[1]
#                 idx3b = planetsunparams["perturbingbodiesid"][jj];
#                 if idx3b == 11
#                     if planetsunparams["centralBodyIDX"] == 3
#                         rV_3b = inertialrefframeinfo["equatorial2inertial"]*celestialbodiesephemeris_position(idx3b,mjd2000);
#                     else
#                         continue;
#                     end
#                 else
#                     rV_3bWRTS = celestialbodiesephemeris_position(idx3b,mjd2000);
#                     rV_CBWRTS = celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000);
#                     rV_3b = inertialrefframeinfo["ecliptic2inertial"]*(rV_3bWRTS-rV_CBWRTS);
#                 end
#                 tbgacc = tbgacc + get3BodyAcc(rV,planetsunparams["perturbingbodiesgravparam"][jj],rV_3b);
#             end
#         end

#         pertacc = pertacc + Ri2o*(tbgacc+sungacc+dragacc+srpacc);

#         exttorque = exttorque + gtorque+mtorque+srptorque+dragtorque;  
#         exttorque[abs.(exttorque).<1e-15].=0.0;     
#     end

#     if includeZonalHarmsAcc
#         J2 = planetsunparams["zonalharmonicscoff"][1]
#         J3 = planetsunparams["zonalharmonicscoff"][2]
#         J4 = planetsunparams["zonalharmonicscoff"][3]
#         J5 = planetsunparams["zonalharmonicscoff"][4]
#         potpert = getCentralBodyPotentialAcc(sma,eta,Q1,Q2,cTL,sTL,Phi,mu,planetsunparams["rPlanet"],J2,J3,J4,J5)
#         pertacc = pertacc + potpert;
#     end  

#     # variation in time of equinoctial elements
#     aR  = pertacc[1];
#     aT  = pertacc[2]; 
#     aN  = pertacc[3];   

#     du[1] = (2 / eta) * sqrt(sma^3 / mu) * ((P2 * sTL - P1 * cTL) * aR + Phi * aT)
#     du[2] = eta * sqrt(sma / mu) * (-cTL * aR + (((P1 + sTL) / Phi) + sTL) * aT - P2 * ((Q1 * cTL - Q2 * sTL) / Phi) * aN)
#     du[3] = eta * sqrt(sma / mu) * (sTL * aR + (((P2 + cTL) / Phi) + cTL) * aT + P1 * ((Q1 * cTL - Q2 * sTL) / Phi) * aN)
#     du[4] = 0.5 * eta * sqrt(sma / mu) * (1 + Q1^2 + Q2^2) * sTL / Phi * aN
#     du[5] = 0.5 * eta * sqrt(sma / mu) * (1 + Q1^2 + Q2^2) * cTL / Phi * aN    
#     du[6] = meanmotion - (rnorm / hh) * ((sma / (sma + bb) * Phi * (P1 * sTL + P2 * cTL) + (2 * bb / sma)) * aR +
#             sma / (sma + bb) * (1.0 + Phi) * (P1 * cTL - P2 * sTL) * aT +
#             (Q1 * cTL - Q2 * sTL) * aN)

#     # variation in time of quaternions
#     du[7] = -0.5*(om[1]*quat[2]+om[2]*quat[3]+om[3]*quat[4]);
#     du[8] =  0.5*(quat[1]*om[1]-om[2]*quat[4]+om[3]*quat[3]);
#     du[9] =  0.5*(quat[1]*om[2]-om[3]*quat[2]+om[1]*quat[4]);
#     du[10] = 0.5*(quat[1]*om[3]-om[1]*quat[3]+om[2]*quat[2]);

    
        
#     # variation in time of components of angular velocity (Euler's equations)
#     du[11] = (exttorque[1] - om[2]*om[3]*(C-B))/A;
#     du[12] = (exttorque[2] - om[3]*om[1]*(A-C))/B;
#     du[13] = (exttorque[3] - om[1]*om[2]*(B-A))/C;
        
#     # println(exttorque,du[11:13])
#     return du;
# end





