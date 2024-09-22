# new 
##############################################################################################################################################################################
# TRIAXIAL BODIES
##############################################################################################################################################################################
# perturbations

function gravitytorque_girf_triaxial(k,k1r,skm,smk,m,mM1KK,T1,T3,T4,T6,Jg,spsih,cpsih,sdelta,cdelta,sma,eta,q11,q12,q13,q21,q22,q23,A,B,C,mu)
    out = zeros(6,1);
    
    Jhdot_gravityA = (3.0*A*((-3.0*T1 + 1.0)*smk + skm)*mu*(2.0*(cdelta - 1.0)*(cdelta + 1.0)*(q11*q12 + q21*q22)*cpsih^2.0 + ((cdelta - 1.0)*(cdelta + 1.0)*(-q11^2.0 + q12^2.0 - q21^2.0 + q22^2.0)*spsih + cdelta*sdelta*(q11*q13 + q21*q23))*cpsih + cdelta*sdelta*(q12*q13 + q22*q23)*spsih - (cdelta - 1.0)*(cdelta + 1.0)*(q11*q12 + q21*q22)))/(4.0*sma^3.0*eta^3.0);       
    Jhdot_gravityB = -(9.0*mu*B*(-2*(cdelta - 1)*(cdelta + 1)*(q11*q12 + q21*q22)*cpsih^2 + ((cdelta - 1)*(cdelta + 1)*(q11^2 - q12^2 + q21^2 - q22^2)*spsih - cdelta*sdelta*(q11*q13 + q21*q23))*cpsih - cdelta*sdelta*(q12*q13 + q22*q23)*spsih + (cdelta - 1)*(cdelta + 1)*(q11*q12 + q21*q22))*((T1*m - m + 1/3)*skm + smk*(T1 - 2/3)))/(4*sma^3*eta^3);     
    Jhdot_gravityC = -9.0*((T3 - 1/3)*skm - smk/3)*C*mu*(2*(cdelta - 1)*(cdelta + 1)*(q11*q12 + q21*q22)*cpsih^2 + ((cdelta - 1)*(cdelta + 1)*(-q11^2 + q12^2 - q21^2 + q22^2)*spsih + cdelta*sdelta*(q11*q13 + q21*q23))*cpsih + cdelta*sdelta*(q12*q13 + q22*q23)*spsih - (cdelta - 1)*(cdelta + 1)*(q11*q12 + q21*q22))/(4*sma^3*eta^3);
    Jhdot_gravity  = Jhdot_gravityA + Jhdot_gravityB + Jhdot_gravityC;

    psildot_gravityA = -(9*mu*A*(T1^2*m + T6 + m - 1)*sqrt(skm)*((cdelta - 1)*(cdelta + 1)*(q11^2 - q12^2 + q21^2 - q22^2)*cpsih^2 + (2*(cdelta - 1)*(cdelta + 1)*(q11*q12 + q21*q22)*spsih - 2*cdelta*sdelta*(q12*q13 + q22*q23))*cpsih + 2*cdelta*sdelta*(q11*q13 + q21*q23)*spsih + (-q11^2 + q13^2 - q21^2 + q23^2)*cdelta^2 - q23^2/3 + (2*q21^2)/3 + (2*q11^2)/3 - q12^2/3 - q13^2/3 - q22^2/3)*pi)/(16*k1r*Jg*sma^3*eta^3*mM1KK);
    psildot_gravityB =  (9*mu*B*k1r*sqrt(skm)*(T1^2*m + T6 - m + 1)*((cdelta - 1)*(cdelta + 1)*(q11^2 - q12^2 + q21^2 - q22^2)*cpsih^2 + (2*(cdelta - 1)*(cdelta + 1)*(q11*q12 + q21*q22)*spsih - 2*cdelta*sdelta*(q12*q13 + q22*q23))*cpsih + 2*cdelta*sdelta*(q11*q13 + q21*q23)*spsih + (-q11^2 + q13^2 - q21^2 + q23^2)*cdelta^2 - q23^2/3 + (2*q21^2)/3 + (2*q11^2)/3 - q12^2/3 - q13^2/3 - q22^2/3)*pi)/(16*Jg*sma^3*eta^3*mM1KK);
    psildot_gravityC = -(9*C*mu*((cdelta - 1)*(cdelta + 1)*(q11^2 - q12^2 + q21^2 - q22^2)*cpsih^2 + (2*(cdelta - 1)*(cdelta + 1)*(q11*q12 + q21*q22)*spsih - 2*cdelta*sdelta*(q12*q13 + q22*q23))*cpsih + 2*cdelta*sdelta*(q11*q13 + q21*q23)*spsih + (-q11^2 + q13^2 - q21^2 + q23^2)*cdelta^2 - q23^2/3 + (2*q21^2)/3 + (2*q11^2)/3 - q12^2/3 - q13^2/3 - q22^2/3)*((T1^2*m + T6 - m + 1)*k - 2*m + 2)*sqrt(skm)*pi)/(16*k1r*Jg*sma^3*eta^3*mM1KK);
    psildot_gravity  = psildot_gravityA + psildot_gravityB + psildot_gravityC;

    psihdot_gravityA = -(3*A*(cdelta*sdelta*(q11^2 - q12^2 + q21^2 - q22^2)*cpsih^2 + (2*cdelta*sdelta*(q11*q12 + q21*q22)*spsih + (cdelta - sdelta)*(cdelta + sdelta)*(q12*q13 + q22*q23))*cpsih - (cdelta - sdelta)*(cdelta + sdelta)*(q11*q13 + q21*q23)*spsih - cdelta*sdelta*(q11^2 - q13^2 + q21^2 - q23^2))*((-3*T1 + 1)*smk + skm)*mu)/(4*sma^3*eta^3*sdelta*Jg);
    psihdot_gravityB = -(9*mu*B*(cdelta*sdelta*(q11^2 - q12^2 + q21^2 - q22^2)*cpsih^2 + (2*cdelta*sdelta*(q11*q12 + q21*q22)*spsih + (cdelta - sdelta)*(cdelta + sdelta)*(q12*q13 + q22*q23))*cpsih - (cdelta - sdelta)*(cdelta + sdelta)*(q11*q13 + q21*q23)*spsih + cdelta*sdelta*(-q11^2 + q13^2 - q21^2 + q23^2))*((T1*m - m + 1/3)*skm + smk*(T1 - 2/3)))/(4*sma^3*eta^3*sdelta*Jg);
    psihdot_gravityC =  (9*C*(cdelta*sdelta*(q11^2 - q12^2 + q21^2 - q22^2)*cpsih^2 + (2*cdelta*sdelta*(q11*q12 + q21*q22)*spsih + (cdelta - sdelta)*(cdelta + sdelta)*(q12*q13 + q22*q23))*cpsih - (cdelta - sdelta)*(cdelta + sdelta)*(q11*q13 + q21*q23)*spsih - cdelta*sdelta*(q11^2 - q13^2 + q21^2 - q23^2))*mu*((T3 - 1/3)*skm - smk/3))/(4*sma^3*eta^3*sdelta*Jg);
    psihdot_gravity  = psihdot_gravityA + psihdot_gravityB + psihdot_gravityC;

    psigdot_gravityA = -cdelta*psihdot_gravityA + 2.0/pi*sqrt(skm)*T4*k1r*psildot_gravityA;
    psigdot_gravityB = -cdelta*psihdot_gravityB + 2.0/pi*sqrt(skm)*T4*k1r*psildot_gravityB;
    psigdot_gravityC = -cdelta*psihdot_gravityC + 2.0/pi*sqrt(skm)*T4*k1r*psildot_gravityC;
    psigdot_gravity  = psigdot_gravityA + psigdot_gravityB + psigdot_gravityC;
    
    # res
    out[3] = out[3] +  Jhdot_gravity;
    out[4] = out[4] +  psildot_gravity;
    out[5] = out[5] +  psigdot_gravity;
    out[6] = out[6] +  psihdot_gravity;
    out[abs.(out).<1e-15].=0.0
    return out;
end

function magnetictorque_girf_triaxial(k,k1r,skm,m,mM1KK,T1,T4,KK,Jg,spsih,cpsih,sdelta,cdelta,cpa1,cpa2,cpa3,sma,eta,q11,q12,q13,q21,q22,q23,muM,IMz)
    out = zeros(6,1);
    
    Jhdot_magnetic   = -(3.0*pi*IMz*muM*((q11^2*cpa1 + (cpa2*q12 + cpa3*q13)*q11 + cpa1*q21^2 + (cpa2*q22 + cpa3*q23)*q21 - (2*cpa1)/3)*cpsih + spsih*(q12^2*cpa2 + cpa2*q22^2 + cpa1*q21*q22 + q12*q13*cpa3 + q11*q12*cpa1 + cpa3*q22*q23 - 2/3*cpa2))*sdelta*sqrt(skm))/(4*sma^3*eta^3*KK);
    psildot_magnetic = -(3.0*muM*((m + k)*T1 - m + 1)*pi^2*(-((q11*q12 + q21*q22)*cpa1 + (q12^2 + q22^2 - 2/3)*cpa2 + cpa3*(q12*q13 + q22*q23))*sdelta*cpsih + ((q11^2 + q21^2 - 2/3)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*sdelta*spsih + cdelta*((q11*q13 + q21*q23)*cpa1 + (q12*q13 + q22*q23)*cpa2 + cpa3*(-2/3 + q13^2 + q23^2)))*IMz)/(8*KK*k1r*Jg*sma^3*eta^3*mM1KK);
    psihdot_magnetic = -(3.0*pi*muM*(-(q12^2*cpa2 + cpa2*q22^2 + cpa1*q21*q22 + q12*q13*cpa3 + q11*q12*cpa1 + cpa3*q22*q23 - 2/3*cpa2)*cdelta*cpsih + (q11^2*cpa1 + (cpa2*q12 + cpa3*q13)*q11 + cpa1*q21^2 + (cpa2*q22 + cpa3*q23)*q21 - (2*cpa1)/3)*cdelta*spsih - (q11*q13*cpa1 + q12*q13*cpa2 + cpa3*q13^2 + cpa1*q21*q23 + cpa2*q22*q23 + cpa3*q23^2 - 2/3*cpa3)*sdelta)*IMz*sqrt(skm))/(4*sma^3*eta^3*KK*sdelta*Jg);
    psigdot_magnetic = -cdelta*psihdot_magnetic + 2/pi*sqrt(skm)*T4*k1r*psildot_magnetic;

    # checks
    if abs(Jhdot_magnetic)<1e-15
        Jhdot_magnetic   = 0.0;
    end
    if abs(psildot_magnetic)<1e-15
        psildot_magnetic = 0.0;
    end
    if abs(psihdot_magnetic)<1e-15
        psihdot_magnetic = 0.0;
    end
    if abs(psigdot_magnetic)<1e-15
        psigdot_magnetic = 0.0;
    end

    out[3] = out[3] +  Jhdot_magnetic;
    out[4] = out[4] +  psildot_magnetic;
    out[5] = out[5] +  psigdot_magnetic;
    out[6] = out[6] +  psihdot_magnetic;
     
    return out;
end

function srptorque_semianalytical(satellite,sma,P1,P2,eta,q11,q12,q13,q21,q22,q23,rSunV,k,k1r,skm,smk,m,mM1KK,Jg,KK,EE,T1,T2,T3,T4,T6,npsig,cdelta,sdelta,cpsih,spsih,muPlanet,higherordercorrections,includeEclipsesEffects,inshadow,ELin,ELout)

    mxb11 = 0.0; mxb12 = 0.0; mxb13 = 0.0; mxpsil = 0.0;
    myb21 = 0.0; myb22 = 0.0; myb23 = 0.0; mypsil = 0.0;
    mzb31 = 0.0; mzb32 = 0.0; mzb33 = 0.0; mzpsil = 0.0;       
    
    numberOfFacets = get(satellite,"numberOfFacets",0);
    facets = get(satellite,"facets",0.0);  
    
    vcu = getAveragedScSunUnitVector(rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,1,includeEclipsesEffects,inshadow,ELin,ELout);  
    cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1 = srpanalyticalaveragedcoeff_T1(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
    cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2 = srpanalyticalaveragedcoeff_T2(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
    cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3 = srpanalyticalaveragedcoeff_T3(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
    
    # srphoc = zeros(6,1);
    # if higherordercorrections 
    #     rsunnorm = sqrt(rSunV[1]^2.0+rSunV[2]^2.0+rSunV[3]^2.0)
    #     cs1      = rSunV[1]/rsunnorm
    #     cs2      = rSunV[2]/rsunnorm
    #     cs3      = rSunV[3]/rsunnorm
    #     EL1 = 0.0;
    #     EL2 = 0.0;
    #     if includeEclipsesEffects && inshadow
    #         EL1 = ELin
    #         EL2 = ELout 
    #         if ELout<ELin
    #             EL2 = EL2 + 2.0*pi;
    #         end
    #     end
    #     srphoc =  gethigherordertermsgravitysrp(k,skm,m,KK,EE,T4,Jg,cdelta,sdelta,cpsih,spsih,satellite,cs1,cs2,cs3,rsunnorm,sma,P1,P2,eta,q11,q12,q13,q21,q22,q23,EL1,EL2,muPlanet)
    # end     
        
    for kk = 1:numberOfFacets
        
        VVT1 = facets[kk][5];
        VVT2 = facets[kk][6];
        VVT3 = facets[kk][7];
        
        Mxb11T1,Mxb12T1,Mxb13T1,MxpsilT1,Myb21T1,Myb22T1,Myb23T1,MypsilT1,Mzb31T1,Mzb32T1,Mzb33T1,MzpsilT1 = srp_averagedterms(VVT1,cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1);
        Mxb11T2,Mxb12T2,Mxb13T2,MxpsilT2,Myb21T2,Myb22T2,Myb23T2,MypsilT2,Mzb31T2,Mzb32T2,Mzb33T2,MzpsilT2 = srp_averagedterms(VVT2,cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2);
        Mxb11T3,Mxb12T3,Mxb13T3,MxpsilT3,Myb21T3,Myb22T3,Myb23T3,MypsilT3,Mzb31T3,Mzb32T3,Mzb33T3,MzpsilT3 = srp_averagedterms(VVT3,cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3);


        mxb11 = mxb11 + Mxb11T1 + Mxb11T2 + Mxb11T3;
        mxb12 = mxb12 + Mxb12T1 + Mxb12T2 + Mxb12T3;
        mxb13 = mxb13 + Mxb13T1 + Mxb13T2 + Mxb13T3;
        mxpsil = mxpsil + MxpsilT1 + MxpsilT2 + MxpsilT3;

        myb21 = myb21 + Myb21T1 + Myb21T2 + Myb21T3;
        myb22 = myb22 + Myb22T1 + Myb22T2 + Myb22T3;
        myb23 = myb23 + Myb23T1 + Myb23T2 + Myb23T3;
        mypsil = mypsil + MypsilT1 + MypsilT2 + MypsilT3;
        
        mzb31 = mzb31 + Mzb31T1 + Mzb31T2 + Mzb31T3;
        mzb32 = mzb32 + Mzb32T1 + Mzb32T2 + Mzb32T3;
        mzb33 = mzb33 + Mzb33T1 + Mzb33T2 + Mzb33T3;
        mzpsil = mzpsil + MzpsilT1 + MzpsilT2 + MzpsilT3;
        
    end

    au = astroconstants(0)[2];
    earthsundist = norm(rSunV);
    psrpcorr     = (au/earthsundist)^2.0
    
    mxb11 = mxb11*psrpcorr
    mxb12 = mxb12*psrpcorr
    mxb13 = mxb13*psrpcorr
    mxpsil = mxpsil*psrpcorr
    myb21 = myb21*psrpcorr
    myb22 = myb22*psrpcorr
    myb23 = myb23*psrpcorr
    mypsil = mypsil*psrpcorr        
    mzb31 = mzb31*psrpcorr
    mzb32 = mzb32*psrpcorr
    mzb33 = mzb33*psrpcorr
    mzpsil = mzpsil*psrpcorr

    out = zeros(6,1);

    JhdotMcdeltaJgdotOversdelta = mxb12+myb22+mzb32;
    psihdotsdeltaJg = (mxb11+myb21+mzb31);
    
    out[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;

    out[2] = mxb13+myb23+mzb33;     

    out[3] = cdelta*out[2] + sdelta*JhdotMcdeltaJgdotOversdelta;

    out[4] = pi/(2*mM1KK*Jg)*(mxpsil+(1-m)*mypsil+mzpsil);

    out[6] = psihdotsdeltaJg/Jg/sdelta;
    
    out[5] = -cdelta*out[6] + 2/pi*sqrt(skm)*T4*k1r*out[4];

    out[abs.(out).<1e-15] .=0.0;
    
    # srphoc[abs.(srphoc).<1e-15].=0.0;
    
    # out    = srphoc + out

    # if higherordercorrections
    #     n = sqrt(muPlanet/(sma^3.0));
    #     hot = gethigherordertermsgravitysrp(k,k1r,skm,smk,m,KK,mM1KK,EE,T1,T2,T3,T4,T6,Jg,cdelta,sdelta,cpsih,spsih,satellite,rSunV,includeEclipsesEffects,inshadow,ELin,ELout,P1,P2,n,psrpcorr,npsig)
    #     out = out + hot;
    # end
    
    return out;

end

function dragtorque_semianalytical(satellite,vmt,k,k1r,skm,smk,m,mM1KK,Jg,KK,T1,T2,T3,T4,T6,T7,T8,cdelta,sdelta,cpsih,spsih)
    out = zeros(6,1);
  
    VVT1 = satellite["constantstermsdragtorque"][1]
    VVT2 = satellite["constantstermsdragtorque"][2]
    VVT3 = satellite["constantstermsdragtorque"][3]
    VVT4 = satellite["constantstermsdragtorque"][4]

    skmBlock,TBblock1,TBblock2,TBblock3,TSblock = getdragsaaveragedtermsattitudeeq(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,T7,T8,cdelta,sdelta,cpsih,spsih,vmt,VVT1,VVT2,VVT3,VVT4)

    out[1] = skmBlock/Jg;

    out[2] = TBblock3;     

    out[3] = cdelta*out[2] + sdelta*TBblock2;

    out[4] = pi/(2*mM1KK*Jg)*TSblock;

    out[6] = TBblock1/Jg/sdelta;
    
    out[5] = -cdelta*out[6] + 2/pi*sqrt(skm)*T4*k1r*out[4];

    out[abs.(out).<1e-15].= 0.0;
  
    return out;

end

function getsingularityterms(m,k,smk,skm,KK,EE,PP)
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

    if k == 0.0
        T7 = -15/32*(smk+2*skm)*smk^2.0/(skm^3.0);
    else
        T7 = (-10*(m/10 + k + 4/5)*(m - 1)*KK + 5*EE*((2*m^2)/5 + (k + 3/5)*m - 2*k - 8/5))/(2*KK*k^3);
    end
    if m == 0.0
        T8 = (68.0*skm*smk + 75.0*smk^2)/(178.0*skm^2) + 24.0/178.0;
    else
        T8 = ((24.0*k^2 + (88.0*m - 32.0)*k + 157.0*m^2 - 176.0*m + 64.0)*KK - 24.0*EE*(k^2.0 + (3.0*m - 4.0/3.0)*k + (89.0*m^2.0)/24.0 - 6.0*m + 8.0/3.0))/(89.0*m*k^2*KK);
    end

    return T1,T2,T3,T4,T6,T7,T8
end

# ignore
function srptorque_semianalytical_numerical(satellite,vmt,rSunV,k,k1r,skm,smk,m,mM1KK,Jg,KK,T4,cdelta,sdelta,cpsih,spsih)

    mxb11 = 0.0; mxb12 = 0.0; mxb13 = 0.0; mxpsil = 0.0;
    myb21 = 0.0; myb22 = 0.0; myb23 = 0.0; mypsil = 0.0;
    mzb31 = 0.0; mzb32 = 0.0; mzb33 = 0.0; mzpsil = 0.0;       
    
    numberOfFacets = get(satellite,"numberOfFacets",0);
    facets = get(satellite,"facets",0.0);  
    
    cUV = rSunV/norm(rSunV);
    cU1 = cUV[1]; cU2 = cUV[2]; cU3 = cUV[3];
    cUH1 = cpsih*cU1 + spsih*cU2;
    cUH2 = -cdelta*spsih*cU1 + cdelta*cpsih*cU2 + sdelta*cU3
    cUH3 = spsih*sdelta*cU1 - cpsih*sdelta*cU2 + cdelta*cU3
    cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1 = srpnumericalaveragedcoeff_T1(vmt,cUH1,cUH2,cUH3)
    cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2 = srpnumericalaveragedcoeff_T2(vmt,cUH1,cUH2,cUH3)
    cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3 = srpnumericalaveragedcoeff_T3(vmt,cUH1,cUH2,cUH3)
            
    for kk = 1:numberOfFacets
        
        VVT1 = facets[kk][5];
        VVT2 = facets[kk][6];
        VVT3 = facets[kk][7];
        
        Mxb11T1,Mxb12T1,Mxb13T1,MxpsilT1,Myb21T1,Myb22T1,Myb23T1,MypsilT1,Mzb31T1,Mzb32T1,Mzb33T1,MzpsilT1 = srp_averagedterms(VVT1,cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1);
        Mxb11T2,Mxb12T2,Mxb13T2,MxpsilT2,Myb21T2,Myb22T2,Myb23T2,MypsilT2,Mzb31T2,Mzb32T2,Mzb33T2,MzpsilT2 = srp_averagedterms(VVT2,cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2);
        Mxb11T3,Mxb12T3,Mxb13T3,MxpsilT3,Myb21T3,Myb22T3,Myb23T3,MypsilT3,Mzb31T3,Mzb32T3,Mzb33T3,MzpsilT3 = srp_averagedterms(VVT3,cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3);

        mxb11 = mxb11 + Mxb11T1 + Mxb11T2 + Mxb11T3;
        mxb12 = mxb12 + Mxb12T1 + Mxb12T2 + Mxb12T3;
        mxb13 = mxb13 + Mxb13T1 + Mxb13T2 + Mxb13T3;
        mxpsil = mxpsil + MxpsilT1 + MxpsilT2 + MxpsilT3;

        myb21 = myb21 + Myb21T1 + Myb21T2 + Myb21T3;
        myb22 = myb22 + Myb22T1 + Myb22T2 + Myb22T3;
        myb23 = myb23 + Myb23T1 + Myb23T2 + Myb23T3;
        mypsil = mypsil + MypsilT1 + MypsilT2 + MypsilT3;
        
        mzb31 = mzb31 + Mzb31T1 + Mzb31T2 + Mzb31T3;
        mzb32 = mzb32 + Mzb32T1 + Mzb32T2 + Mzb32T3;
        mzb33 = mzb33 + Mzb33T1 + Mzb33T2 + Mzb33T3;
        mzpsil = mzpsil + MzpsilT1 + MzpsilT2 + MzpsilT3;
        
    end

    au = astroconstants(0)[2];
    earthsundist = norm(rSunV);
    psrpcorr     = (au/earthsundist)^2.0
    mxb11 = mxb11*psrpcorr
    mxb12 = mxb12*psrpcorr
    mxb13 = mxb13*psrpcorr
    mxpsil = mxpsil*psrpcorr
    myb21 = myb21*psrpcorr
    myb22 = myb22*psrpcorr
    myb23 = myb23*psrpcorr
    mypsil = mypsil*psrpcorr        
    mzb31 = mzb31*psrpcorr
    mzb32 = mzb32*psrpcorr
    mzb33 = mzb33*psrpcorr
    mzpsil = mzpsil*psrpcorr

    out = zeros(6,1);

    JhdotMcdeltaJgdotOversdelta = mxb12+myb22+mzb32;
    psihdotsdeltaJg = (mxb11+myb21+mzb31);
    
    out[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;

    out[2] = mxb13+myb23+mzb33;     

    out[3] = cdelta*out[2] + sdelta*JhdotMcdeltaJgdotOversdelta;

    out[4] = pi/(2*mM1KK*Jg)*(mxpsil+(1-m)*mypsil+mzpsil);

    out[6] = psihdotsdeltaJg/Jg/sdelta;
    
    out[5] = -cdelta*out[6] + 2/pi*sqrt(skm)*T4*k1r*out[4];

    out[abs.(out).<1e-15] .=0.0;

    
    # out    = srphoc + out
    
    return out;

end

# semi-analyitical model 
function attitudeMeanEv_semianalytical_triaxial(du,u,p,t)
    # println(t/24/3600)
    # inputs
    planetsunparams      = p[1];
    inertialrefframeinfo = p[2];
    satellite            = p[3];
    draginfo             = p[4];
    settings             = p[5];
    mjd2000_0            = p[6];
    
    # unperturbed orbital evolution
    mu  =  planetsunparams["muPlanet"];
    sma = u[7];
    P1  = u[8];
    P2  = u[9];
    Q1 =  u[10];
    Q2 =  u[11];
    du[7]  = 0.0; #sma
    du[8]  = 0.0;  #P1
    du[9]  = 0.0;  #P2
    du[10] = 0.0;  #Q1
    du[11] = 0.0;  #Q2
    du[12] = sqrt(mu/sma^3); 


    # unperturbed attitude evolution
    IV = get(satellite,"MomentsOfInertia",[0.0,0.0,0.0]);
    k  = get(satellite,"k_constant",[0.0]);
    k1r = sqrt(1+k);
    A = IV[1]; C = IV[3]; 
    Jg  = u[2];
    skm = u[1];
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
        if k==0.0
            PPMsmkKK = skm*pi/2.0;
            KK = pi/2.0
            PP = pi/2.0
        else
            KK  = ellF(pi/2,m); #Elliptic.K(m);
            PP  = ellP(-k,pi/2.0,m); #Elliptic.Pi(-k,pi/2.0,m);
            PPMsmkKK = (PP-smk*KK);
        end
    end

    mM1KK = (m-1)*KK;
    if abs(mM1KK)<1e-15
        mM1KK = mM1KK +1e-10;
    end

    npsil = -pi/k1r*sqrt(skm)*Jg/2/A/C*(C-A)/KK;
    npsig = Jg/A/C*(C-A)*PPMsmkKK/KK + Jg/A/C*(C*smk+A*skm)
    
    du[1] = 0.0;
    du[2] = 0.0;
    du[3] = 0.0;
    du[4] = npsil;
    du[5] = npsig;
    du[6] = 0.0;
    # println([2*pi/npsil, 2*pi/npsig])
    # println(t," - ", m, " ", u[1:6])


    # perturbations for torque
    includeGravityTorque   = settings["includeGravityTorque"];
    includeMagneticTorque  = settings["includeMagneticTorque"];
    includeSrpTorque       = settings["includeSrpTorque"];
    includeDragTorque      = settings["includeDragTorque"];
    higherordercorrections = settings["higherordercorrections"];

    # perturbations for orbit
    includeZonalHarmsAcc  = settings["includeZonalHarmsAcc"];
    includeSrpAcc         = settings["includeSrpAcc"];
    includeThirdBodyAcc   = settings["includeThirdBodyAcc"];
    includeSunGravityAcc  = settings["includeSunGravityAcc"];
    includeDragAcc         = settings["includeDragAcc"];

    # constants needed for multiple perturbations
    if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque || includeSrpAcc || includeDragAcc || includeZonalHarmsAcc || includeThirdBodyAcc || includeSunGravityAcc
        B = IV[2];
        cdelta = u[3]/Jg;
        if abs(cdelta)==1 || abs(cdelta)>1
            Base.error("singularity")            
        end
        sdelta = sqrt(1.0-cdelta^2.0);
        cpsih = cos(u[6]);
        spsih = sin(u[6]);
        EE  = Elliptic.E(m);
        T1,T2,T3,T4,T6,T7,T8 = getsingularityterms(m,k,smk,skm,KK,EE,PP)

        eta = sqrt(1.0-P1^2-P2^2);
        GG = 1.0+Q1^2.0+Q2^2.0;
        q11 = (1.0-Q1^2.0+Q2^2.0)/GG;
        q12 = (2.0*Q1*Q2)/GG;
        q13 = (-2.0*Q1)/GG;
        q21 = (2.0*Q1*Q2)/GG;
        q22 = (1.0+Q1^2.0-Q2^2.0)/GG;
        q23 = (2.0*Q2)/GG;
        q31 = (2.0*Q1)/GG;
        q32 = (-2.0*Q2)/GG;
        q33 = (1.0-Q1^2.0-Q2^2.0)/GG;  

        if includeDragTorque || includeDragAcc || includeMagneticTorque
            # Recef2eci = r_ecef_to_eci(J2000(), PEF(),mjd20002jd(mjd2000_0))
            # polaraxis = (inertialrefframeinfo["equatorial2inertial"]*Recef2eci*[0.0,0.0,1.0]);
            polaraxis = [inertialrefframeinfo["cpa1"],inertialrefframeinfo["cpa2"],inertialrefframeinfo["cpa3"]]
        end
        
        if  includeSrpTorque ||    includeSrpAcc || includeThirdBodyAcc || includeSunGravityAcc  || includeDragTorque || includeDragAcc 
            mjd2000 = mjd2000_0 + t/24.0/3600.0;
        end
    end

    # perturbations effects
    if includeGravityTorque 
        gtcontr = gravitytorque_girf_triaxial(k,k1r,skm,smk,m,mM1KK,T1,T3,T4,T6,Jg,spsih,cpsih,sdelta,cdelta,sma,eta,q11,q12,q13,q21,q22,q23,A,B,C,mu);
        du[1:6] = du[1:6] + gtcontr;    
    end

    if includeMagneticTorque
        muM = get(planetsunparams,"muM",0.0);
        IM = get(satellite,"intrinsicMagneticMoment",[0.0;0.0;0.0]);
        IMz = IM[3];
        cpa1 =  polaraxis[1];
        cpa2 =  polaraxis[2];
        cpa3 =  polaraxis[3];
        mtcontr = magnetictorque_girf_triaxial(k,k1r,skm,m,mM1KK,T1,T4,KK,Jg,spsih,cpsih,sdelta,cdelta,cpa1,cpa2,cpa3,sma,eta,q11,q12,q13,q21,q22,q23,muM,IMz);
        du[1:6] = du[1:6] + mtcontr;  
    end

    if higherordercorrections && (includeGravityTorque || includeMagneticTorque)
        # if m!=0.0
        if includeMagneticTorque==false
            muM = 0.0;
            IMz = 0.0;
        end
        hoccontr = gethigherordercorrections(k,m,skm,smk,mM1KK,KK,EE,PP,T2,Jg,cdelta,sdelta,u[6],sma,eta,P1,P2,Q1,Q2,A,B,C,mu,IMz,muM,inertialrefframeinfo["cpa1"],inertialrefframeinfo["cpa2"],inertialrefframeinfo["cpa3"],[includeGravityTorque,includeMagneticTorque],1);
        # hoccontr = gethigherordercorrections(k,m,KK,EE,PP,Jg,cdelta,sdelta,u[6],sma,eta,P1,P2,Q1,Q2,A,B,C,mu,IMz,muM,inertialrefframeinfo["cpa1"],inertialrefframeinfo["cpa2"],inertialrefframeinfo["cpa3"],[includeGravityTorque,includeMagneticTorque]);
        
        du[1:6] = du[1:6]+hoccontr;
        # end
    end

    if includeSrpTorque || includeSrpAcc || includeSunGravityAcc
        # sun ephemeris
        if settings["sunFixed"]
            rSunV = inertialrefframeinfo["ecliptic2inertial"]*settings["sunposition"];
        else
            rSunV = -inertialrefframeinfo["ecliptic2inertial"]*celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000)
        end

        if includeSrpTorque || includeSrpAcc
            # shadow
            includeEclipsesEffectsOnAttitude =  get(settings,"includeEclipsesEffectsOnAttitude",false);
            includeEclipsesEffectsOnOrbit    =  get(settings,"includeEclipsesEffectsOnOrbit",false);

            if includeEclipsesEffectsOnAttitude || includeEclipsesEffectsOnOrbit
                inshadow,TLin,TLout = eclipseinout(rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,planetsunparams["rPlanet"]);
                if inshadow
                    ELin  = truelong2ecclong(TLin,P1,P2);
                    ELout = truelong2ecclong(TLout,P1,P2);
                else
                    ELin = NaN;
                    ELout = NaN;
                end
            else
                inshadow = false;
                ELin = NaN;
                ELout = NaN;
            end
            if includeSrpTorque
                srpTorqueAnalytical = get(settings,"srpTorqueAnalytical",true);
                if srpTorqueAnalytical 
                    srpcontr = srptorque_semianalytical(satellite,sma,P1,P2,eta,q11,q12,q13,q21,q22,q23,rSunV,k,k1r,skm,smk,m,mM1KK,Jg,KK,EE,T1,T2,T3,T4,T6,npsig,cdelta,sdelta,cpsih,spsih,mu,higherordercorrections,includeEclipsesEffectsOnAttitude,inshadow,ELin,ELout);
                else
                    vmt = draginfo["srpmeanv"];
                    srpcontr =  srptorque_semianalytical_numerical(satellite,vmt,rSunV,k,k1r,skm,smk,m,mM1KK,Jg,KK,T4,cdelta,sdelta,cpsih,spsih);
                end
                du[1:6] = du[1:6] + srpcontr;
            end
            if includeSrpAcc
                srpacccontr = srpaveragedgauss_wrapper1_AV(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffectsOnOrbit,inshadow,ELin,ELout);
                du[7:12] = du[7:12] + srpacccontr;
            end
        end

        if  includeSunGravityAcc
            sungcontr = tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rSunV,mu,0.19891000000000E+31 * 6.67259e-20,5); 
            du[7:12] = du[7:12] + sungcontr;
        end

    end

    if includeDragTorque || includeDragAcc

        # # if t>=draginfo["timeupdate"]
        # meanOverMDrag,Idrag = dragaverageOverMUpdate(u[7:12],inertialrefframeinfo["cpa1"],inertialrefframeinfo["cpa2"],inertialrefframeinfo["cpa3"],
        #     planetsunparams["rPlanet"],planetsunparams["muPlanet"],planetsunparams["planetRotation"],planetsunparams["planet_flat"],planetsunparams["oblatenesscoeff"],
        #     planetsunparams["zonalharmonicscoff"][1],includeDragTorque,includeDragAcc,draginfo["atmosphericinfo"],draginfo["correction"],transpose(inertialrefframeinfo["equatorial2inertial"]),mjd2000);
        # # draginfo["averageoverM_attitude"] = meanOverMDrag;
        # draginfo["averageoverM_orbit"] = Idrag;
        # draginfo["averageoverM_coupling"] = eV0mean;
        # draginfo["timeupdate"]   = t + 2.0*pi*sqrt(sma^3.0/mu);
            # println("qui ", t/24/3600)
        # end       

        if includeDragTorque
            dragcontr = dragtorque_semianalytical(satellite,draginfo["averageoverM_attitude"],k,k1r,skm,smk,m,mM1KK,Jg,KK,T1,T2,T3,T4,T6,T7,T8,cdelta,sdelta,cpsih,spsih);
            du[1:6] = du[1:6] + dragcontr;
        end
        if includeDragAcc
            # if draginfo["atmosphericinfo"]["atmosphericmodel"]==1
            #     dragacccontr = dragaveragedgauss_wrapper1_AV(draginfo["averageoverM_orbit"][1:9],sma,P1,P2,Q1,Q2,satellite,k,skm,KK,T1,cdelta,sdelta,cpsih,spsih,
            #     planetsunparams["muPlanet"],draginfo["averageoverM_orbit"][10:20]);
            #     du[7:12] = du[7:12] + dragacccontr;
            # end
            dragacccontr = dragveragedgauss(k,k1r,skm,smk,m,KK,T1,T2,T3,cdelta,sdelta,cpsih,spsih,sma,P1,P2,eta,Q1,Q2,GG,draginfo["averageoverM_orbit"],satellite,mu,draginfo["couplingstrategy"]);
            du[7:12] = du[7:12] + dragacccontr
        end
    end

    if includeZonalHarmsAcc
        J2 = planetsunparams["zonalharmonicscoff"][1]
        J3 = planetsunparams["zonalharmonicscoff"][2]
        J4 = planetsunparams["zonalharmonicscoff"][3]
        J5 = planetsunparams["zonalharmonicscoff"][4]
        rPlanet = planetsunparams["rPlanet"]
        zharmcontr =  zonalharmonicsaveragedgauss(sma,eta,P1,P2,Q1,Q2,mu,rPlanet,J2,J3,J4,J5)
        du[7:12] = du[7:12] + zharmcontr;
    end

    if includeThirdBodyAcc 
        tbgcontr = zeros(6);
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
            tbgcontr = tbgcontr + tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rV_3b,mu,planetsunparams["perturbingbodiesgravparam"][jj],2)
        end
        du[7:12] = du[7:12] + tbgcontr;
    end

    # println(t/24/3600," ",draginfo["averageoverM_attitude"][1])
    # println(t/24/3600)

    return du;
end

# event handling
function conditionpatch(u,t,integrator)
    skm = u[1]
    k = integrator.p[2]["k_constant"];
    if k!=0.0
        fun = skm-1/(1+0.9999/k);
    else
        fun = skm-1e-4;
    end
    return fun;
end

function affectpatch!(integrator)
    terminate!(integrator)
end

###################### still keplerian

####################################################################################################################################################################
# TRIAXIAL BODIES handling singularity sin(delta) = 0
#####################################################################################################################################################################

# perturbations

function gravitytorque_girf_triaxial_sadovlike(k,k1r,skm,smk,m,mM1KK,T1,T3,T4,T6,J2,J3,J6,J7,sma,eta,q11,q12,q13,q21,q22,q23,A,B,C,mu)
    out = zeros(7,1);
    
    if J7==0.0 && J6==0.0
        atanval = 0.0;
    else
        atanval = mod(atan(J7,J6),2.0*pi)
    end

    J3dot_gravityA = -3*((-3*T1 + 1)*smk + skm)*((q11*q12 + q21*q22)*J6^2 + (J7*(-q11^2 + q12^2 - q21^2 + q22^2) - (q11*q13 + q21*q23)*J3)*J6 - J7*(J7*(q11*q12 + q21*q22) + (q12*q13 + q22*q23)*J3))*A*mu/(4*J2^2*sma^3*eta^3);
    J3dot_gravityB = 9*((T1*m - m + 1/3)*skm + smk*(T1 - 2/3))*mu*B*((-q11*q12 - q21*q22)*J6^2 + ((q11*q13 + q21*q23)*J3 + J7*(q11^2 - q12^2 + q21^2 - q22^2))*J6 + J7*(J7*(q11*q12 + q21*q22) + (q12*q13 + q22*q23)*J3))/(4*J2^2*sma^3*eta^3);
    J3dot_gravityC = 9*((q11*q12 + q21*q22)*J6^2 + (J7*(-q11^2 + q12^2 - q21^2 + q22^2) - (q11*q13 + q21*q23)*J3)*J6 - J7*(J7*(q11*q12 + q21*q22) + (q12*q13 + q22*q23)*J3))*C*mu*((T3 - 1/3)*skm - smk/3)/(4*J2^2*sma^3*eta^3);
    J3dot_gravity  = J3dot_gravityA + J3dot_gravityB + J3dot_gravityC;

    J4dot_gravityA = (3*mu*(T1^2*m + T6 + m - 1)*((-3*q12^2 + 3*q13^2 - 3*q22^2 + 3*q23^2)*J6^2 + ((6*q11*q12 + 6*q21*q22)*J7 + 6*(q12*q13 + q22*q23)*J3)*J6 + (-3*q11^2 + 3*q13^2 - 3*q21^2 + 3*q23^2)*J7^2 - 6*J3*(q11*q13 + q21*q23)*J7 + J2^2*(q11^2 + q12^2 - 2*q13^2 + q21^2 + q22^2 - 2*q23^2))*sqrt(skm)*A*pi)/(16*J2^3*k1r*sma^3*eta^3*mM1KK);
    J4dot_gravityB =  -(3*mu*(T1^2*m + T6 - m + 1)*k1r*((-3*q12^2 + 3*q13^2 - 3*q22^2 + 3*q23^2)*J6^2 + ((6*q11*q12 + 6*q21*q22)*J7 + 6*(q12*q13 + q22*q23)*J3)*J6 + (-3*q11^2 + 3*q13^2 - 3*q21^2 + 3*q23^2)*J7^2 - 6*J3*(q11*q13 + q21*q23)*J7 + J2^2*(q11^2 + q12^2 - 2*q13^2 + q21^2 + q22^2 - 2*q23^2))*sqrt(skm)*B*pi)/(16*J2^3*sma^3*eta^3*mM1KK);
    J4dot_gravityC = 3*((T1^2*m + T6 - m + 1)*k - 2*m + 2)*pi*sqrt(skm)*C*((-3*q12^2 + 3*q13^2 - 3*q22^2 + 3*q23^2)*J6^2 + ((6*q11*q12 + 6*q21*q22)*J7 + 6*(q12*q13 + q22*q23)*J3)*J6 + (-3*q11^2 + 3*q13^2 - 3*q21^2 + 3*q23^2)*J7^2 - 6*J3*(q11*q13 + q21*q23)*J7 + J2^2*(q11^2 + q12^2 - 2*q13^2 + q21^2 + q22^2 - 2*q23^2))*mu/(16*J2^3*k1r*sma^3*eta^3*mM1KK);
    J4dot_gravity  = J4dot_gravityA + J4dot_gravityB + J4dot_gravityC;

    J5dot_gravityA = 2.0/pi*sqrt(skm)*T4*k1r*J4dot_gravityA - 3*((-3*T1 + 1)*smk + skm)*((2*q11*q12 + 2*q21*q22)*J6^2 + ((-q11*q13 - q21*q23)*J3 + J7*(-q11^2 + q12^2 - q21^2 + q22^2))*J6 + (q11*q12 + q21*q22)*J3^2 - J3*(q12*q13 + q22*q23)*J7 - J2^2*(q11*q12 + q21*q22))*atanval*A*mu/(4*eta^3*sma^3*J2^3)
    J5dot_gravityB = 2.0/pi*sqrt(skm)*T4*k1r*J4dot_gravityB + 9*((-2*q11*q12 - 2*q21*q22)*J6^2 + ((q11*q13 + q21*q23)*J3 + J7*(q11^2 - q12^2 + q21^2 - q22^2))*J6 + (-q11*q12 - q21*q22)*J3^2 + J3*(q12*q13 + q22*q23)*J7 + J2^2*(q11*q12 + q21*q22))*((T1*m - m + 1/3)*skm + smk*(T1 - 2/3))*mu*B*atanval/(4*eta^3*sma^3*J2^3);
    J5dot_gravityC = 2.0/pi*sqrt(skm)*T4*k1r*J4dot_gravityC + 9*((2*q11*q12 + 2*q21*q22)*J6^2 + ((-q11*q13 - q21*q23)*J3 + J7*(-q11^2 + q12^2 - q21^2 + q22^2))*J6 + (q11*q12 + q21*q22)*J3^2 - J3*(q12*q13 + q22*q23)*J7 - J2^2*(q11*q12 + q21*q22))*atanval*((T3 - 1/3)*skm - smk/3)*C*mu/(4*eta^3*sma^3*J2^3);
    J5dot_gravity  = J5dot_gravityA + J5dot_gravityB + J5dot_gravityC;
    
    J6dot_gravityA = -3*((-3*T1 + 1)*smk + skm)*A*mu*((-2*q11*q13 - 2*q21*q23)*J7^2 + ((q11^2 - q13^2 + q21^2 - q23^2)*J3 + (q12*q13 + q22*q23)*J6)*J7 - J3*J6*(q11*q12 + q21*q22) + (J2 - J6)*(J2 + J6)*(q11*q13 + q21*q23))/(4*J2^2*sma^3*eta^3);
    J6dot_gravityB = -9*((T1*m - m + 1/3)*skm + smk*(T1 - 2/3))*mu*B*((-2*q11*q13 - 2*q21*q23)*J7^2 + ((q11^2 - q13^2 + q21^2 - q23^2)*J3 + (q12*q13 + q22*q23)*J6)*J7 - J3*J6*(q11*q12 + q21*q22) + (J2 - J6)*(J2 + J6)*(q11*q13 + q21*q23))/(4*J2^2*sma^3*eta^3);
    J6dot_gravityC = (9*C*mu*((T3 - 1/3)*skm - smk/3)*((-2*q11*q13 - 2*q21*q23)*J7^2 + ((q11^2 - q13^2 + q21^2 - q23^2)*J3 + (q12*q13 + q22*q23)*J6)*J7 - J3*J6*(q11*q12 + q21*q22) + (J2 - J6)*(J2 + J6)*(q11*q13 + q21*q23)))/(4*J2^2*sma^3*eta^3);
    J6dot_gravity  = J6dot_gravityA + J6dot_gravityB + J6dot_gravityC;
    
    J7dot_gravityA = 3*((-3*T1 + 1)*smk + skm)*A*mu*((-q12*q13 - q22*q23)*J3^2 + ((q12^2 - q13^2 + q22^2 - q23^2)*J6 - J7*(q11*q12 + q21*q22))*J3 + ((q12*q13 + q22*q23)*J6 - J7*(q11*q13 + q21*q23))*J6)/(4*J2^2*sma^3*eta^3);
    J7dot_gravityB = -9*((T1*m - m + 1/3)*skm + smk*(T1 - 2/3))*mu*B*((q12*q13 + q22*q23)*J3^2 + ((-q12^2 + q13^2 - q22^2 + q23^2)*J6 + J7*(q11*q12 + q21*q22))*J3 + J6*((-q12*q13 - q22*q23)*J6 + J7*(q11*q13 + q21*q23)))/(4*J2^2*sma^3*eta^3);
    J7dot_gravityC = -9*((-q12*q13 - q22*q23)*J3^2 + ((q12^2 - q13^2 + q22^2 - q23^2)*J6 - J7*(q11*q12 + q21*q22))*J3 + ((q12*q13 + q22*q23)*J6 - J7*(q11*q13 + q21*q23))*J6)*((T3 - 1/3)*skm - smk/3)*C*mu/(4*J2^2*sma^3*eta^3);
    J7dot_gravity  = J7dot_gravityA + J7dot_gravityB + J7dot_gravityC;

    out[3] = out[3] +  J3dot_gravity;
    out[4] = out[4] +  J4dot_gravity;
    out[5] = out[5] +  J5dot_gravity;
    out[6] = out[6] +  J6dot_gravity;
    out[7] = out[7] +  J7dot_gravity;

    for jj = 1 : 7
        if abs( out[jj])<1e-15
            out[jj] = 0.0;
        end
    end
 

    return out;
end

function magnetictorque_girf_triaxial_sadovlike(k,k1r,skm,m,mM1KK,T1,T4,KK,J2,J3,J6,J7,cpa1,cpa2,cpa3,sma,eta,q11,q12,q13,q21,q22,q23,muM,IMz)
    out = zeros(7,1);

    if J7==0.0 && J6==0.0
        atanval = 0.0;
    else
        atanval = mod(atan(J7,J6),2.0*pi)
    end
    
    J3dot_magnetic = -3*(((-2/3 + q11^2 + q21^2)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*J6 + ((q11*q12 + q21*q22)*cpa1 + (-2/3 + q12^2 + q22^2)*cpa2 + cpa3*(q12*q13 + q22*q23))*J7)*pi*muM*sqrt(skm)*IMz/(4*J2*sma^3*eta^3*KK);
    J4dot_magnetic = -(3*IMz*pi^2*muM*((m + k)*T1 - m + 1)*(((q11*q13 + q21*q23)*J3 + (-q11*q12 - q21*q22)*J6 + J7*(-2/3 + q11^2 + q21^2))*cpa1 + ((q12*q13 + q22*q23)*J3 + (2/3 - q12^2 - q22^2)*J6 + J7*(q11*q12 + q21*q22))*cpa2 + ((q13^2 + q23^2 - 2/3)*J3 + (-q12*q13 - q22*q23)*J6 + J7*(q11*q13 + q21*q23))*cpa3))/(8*J2^2*KK*k1r*sma^3*eta^3*mM1KK);
    J5dot_magnetic = 2/pi*sqrt(skm)*T4*k1r*J4dot_magnetic -(3*muM*(((-2/3 + q11^2 + q21^2)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*J6 + ((q11*q12 + q21*q22)*cpa1 + (-2/3 + q12^2 + q22^2)*cpa2 + cpa3*(q12*q13 + q22*q23))*J7)*pi*atanval*IMz*sqrt(skm))/(4*J2^2*sma^3*eta^3*KK);
    J6dot_magnetic = (3*pi*(((-2/3 + q11^2 + q21^2)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*J3 - J7*((q11*q13 + q21*q23)*cpa1 + (q12*q13 + q22*q23)*cpa2 + cpa3*(q13^2 + q23^2 - 2/3)))*muM*sqrt(skm)*IMz)/(4*J2*sma^3*eta^3*KK);
    J7dot_magnetic = 3*(((q11*q12 + q21*q22)*cpa1 + (-2/3 + q12^2 + q22^2)*cpa2 + cpa3*(q12*q13 + q22*q23))*J3 + J6*((q11*q13 + q21*q23)*cpa1 + (q12*q13 + q22*q23)*cpa2 + cpa3*(q13^2 + q23^2 - 2/3)))*pi*muM*sqrt(skm)*IMz/(4*J2*sma^3*eta^3*KK);

    out[3] = out[3] +  J3dot_magnetic;
    out[4] = out[4] +  J4dot_magnetic;
    out[5] = out[5] +  J5dot_magnetic;
    out[6] = out[6] +  J6dot_magnetic;
    out[7] = out[7] +  J7dot_magnetic;

    # checks
    for jj = 1 : 7
        if abs( out[jj])<1e-15
            out[jj] = 0.0;
        end
    end

    return out;
end

function srptorque_semianalytical_sadovelike(satellite,sma,P1,P2,q11,q12,q13,q21,q22,q23,rSunV,k,k1r,skm,smk,m,mM1KK,KK,T1,T2,T3,T4,T6,J2,J3,J6,J7,includeEclipsesEffects,inshadow,ELin,ELout)
    
    Jg = J2;
    cdelta = J3/Jg;
    if J6==J7 && J6==0.0
        sdelta = 0.0;
        cpsih = 1.0;
        spsih = 0.0;
    else
        sdelta = sqrt(J6^2+J7^2)/Jg;
        cpsih  = J6/sdelta/Jg;
        spsih  = J7/sdelta/Jg;
    end
    psih = mod(atan(spsih,cpsih),2.0*pi)
    
    numberOfFacets = satellite["numberOfFacets"];
    facets = satellite["facets"];  

    vcu = getAveragedScSunUnitVector(rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,0,includeEclipsesEffects,inshadow,ELin,ELout); 
    cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1 = srpanalyticalaveragedcoeff_T1(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
    cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2 = srpanalyticalaveragedcoeff_T2(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);
    cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3 = srpanalyticalaveragedcoeff_T3(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vcu);

    mxb11 = 0.0;     mxb12 = 0.0;     mxb13 = 0.0;     mxpsil = 0.0;
    myb21 = 0.0;     myb22 = 0.0;     myb23 = 0.0;     mypsil = 0.0;
    mzb31 = 0.0;     mzb32 = 0.0;     mzb33 = 0.0;     mzpsil = 0.0;  
    
    for kk = 1:numberOfFacets
        
        VVT1 = facets[kk][5];
        VVT2 = facets[kk][6];
        VVT3 = facets[kk][7];
        
        Mxb11T1,Mxb12T1,Mxb13T1,MxpsilT1,Myb21T1,Myb22T1,Myb23T1,MypsilT1,Mzb31T1,Mzb32T1,Mzb33T1,MzpsilT1 = srp_averagedterms(VVT1,cMxb11T1,cMxb12T1,cMxb13T1,cMxpsilT1,cMyb21T1,cMyb22T1,cMyb23T1,cMypsilT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzpsilT1);
        Mxb11T2,Mxb12T2,Mxb13T2,MxpsilT2,Myb21T2,Myb22T2,Myb23T2,MypsilT2,Mzb31T2,Mzb32T2,Mzb33T2,MzpsilT2 = srp_averagedterms(VVT2,cMxb11T2,cMxb12T2,cMxb13T2,cMxpsilT2,cMyb21T2,cMyb22T2,cMyb23T2,cMypsilT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzpsilT2);
        Mxb11T3,Mxb12T3,Mxb13T3,MxpsilT3,Myb21T3,Myb22T3,Myb23T3,MypsilT3,Mzb31T3,Mzb32T3,Mzb33T3,MzpsilT3 = srp_averagedterms(VVT3,cMxb11T3,cMxb12T3,cMxb13T3,cMxpsilT3,cMyb21T3,cMyb22T3,cMyb23T3,cMypsilT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzpsilT3);

        mxb11 = mxb11 + Mxb11T1 + Mxb11T2 + Mxb11T3;
        mxb12 = mxb12 + Mxb12T1 + Mxb12T2 + Mxb12T3;
        mxb13 = mxb13 + Mxb13T1 + Mxb13T2 + Mxb13T3;
        mxpsil = mxpsil + MxpsilT1 + MxpsilT2 + MxpsilT3;

        myb21 = myb21 + Myb21T1 + Myb21T2 + Myb21T3;
        myb22 = myb22 + Myb22T1 + Myb22T2 + Myb22T3;
        myb23 = myb23 + Myb23T1 + Myb23T2 + Myb23T3;
        mypsil = mypsil + MypsilT1 + MypsilT2 + MypsilT3;
        
        mzb31 = mzb31 + Mzb31T1 + Mzb31T2 + Mzb31T3;
        mzb32 = mzb32 + Mzb32T1 + Mzb32T2 + Mzb32T3;
        mzb33 = mzb33 + Mzb33T1 + Mzb33T2 + Mzb33T3;
        mzpsil = mzpsil + MzpsilT1 + MzpsilT2 + MzpsilT3;
    end
    

    au = astroconstants(0)[2];
    earthsundist = norm(rSunV);
    psrpcorr     = (au/earthsundist)^2.0
    mxb11 = mxb11*psrpcorr
    mxb12 = mxb12*psrpcorr
    mxb13 = mxb13*psrpcorr
    mxpsil = mxpsil*psrpcorr
    myb21 = myb21*psrpcorr
    myb22 = myb22*psrpcorr
    myb23 = myb23*psrpcorr
    mypsil = mypsil*psrpcorr        
    mzb31 = mzb31*psrpcorr
    mzb32 = mzb32*psrpcorr
    mzb33 = mzb33*psrpcorr
    mzpsil = mzpsil*psrpcorr

    out = zeros(7,1);

    JhdotMcdeltaJgOversdelta = mxb12+myb22+mzb32;
    psihdotsdeltaJg = (mxb11+myb21+mzb31);
    
    out[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;

    out[2] = mxb13+myb23+mzb33;     

    out[3] = cdelta*out[2] + sdelta*JhdotMcdeltaJgOversdelta;

    out[4] = pi/(2*mM1KK*Jg)*(mxpsil+(1-m)*mypsil+mzpsil);

    out[5] = 2/pi*sqrt(skm)*T4*k1r*out[4] + psih*sdelta*JhdotMcdeltaJgOversdelta/Jg;

    out[6] = -spsih*psihdotsdeltaJg + sdelta*cpsih*out[2] - cdelta*cpsih*JhdotMcdeltaJgOversdelta;
    
    out[7] =  cpsih*psihdotsdeltaJg + sdelta*spsih*out[2] - cdelta*spsih*JhdotMcdeltaJgOversdelta;

    out[abs.(out).<1e-15] .=0.0;

    return out;

end

function dragtorque_semianalytical_sadovelikeO(satellite,vmt,k,k1r,skm,smk,m,mM1KK,KK,T1,T2,T3,T4,T6,J2,J3,J6,J7)
    
    Jg = J2;
    cdelta = J3/Jg;
    if J6==J7 && J6==0.0
        sdelta = 0.0;
        cpsih = 1.0;
        spsih = 0.0;
    else
        sdelta = sqrt(J6^2+J7^2)/Jg;
        cpsih  = J6/sdelta/Jg;
        spsih  = J7/sdelta/Jg;
    end
    psih = mod(atan(spsih,cspih),2.0*pi)

    mxb11 = 0.0;     mxb12 = 0.0;     mxb13 = 0.0;     mxpsil = 0.0;
    myb21 = 0.0;     myb22 = 0.0;     myb23 = 0.0;     mypsil = 0.0;
    mzb31 = 0.0;     mzb32 = 0.0;     mzb33 = 0.0;     mzpsil = 0.0;  

    numberOfFacets = get(satellite,"numberOfFacets",0);
    facets = get(satellite,"facets",0.0);  
        
    cMxb11,cMxb12,cMxb13,cMxpsil,cMyb21,cMyb22,cMyb23,cMypsil,cMzb31,cMzb32,cMzb33,cMzpsil = dragsemianalyticalaveragedcoeff(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vmt);
    # cMxb11,cMxb12,cMxb13,cMxpsil,cMyb21,cMyb22,cMyb23,cMypsil,cMzb31,cMzb32,cMzb33,cMzpsil = dragsemianalyticalaveragedcoeffV2(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,cdelta,sdelta,cpsih,spsih,vmt);
    for kk = 1:numberOfFacets        
        VVdrag = facets[kk][8];
        Mxb11,Mxb12,Mxb13,Mxpsil,Myb21,Myb22,Myb23,Mypsil,Mzb31,Mzb32,Mzb33,Mzpsil = drag_averagedterms(VVdrag,cMxb11,cMxb12,cMxb13,cMxpsil,cMyb21,cMyb22,cMyb23,cMypsil,cMzb31,cMzb32,cMzb33,cMzpsil);
                        
        mxb11 = mxb11 + Mxb11;
        mxb12 = mxb12 + Mxb12;
        mxb13 = mxb13 + Mxb13;
        mxpsil = mxpsil + Mxpsil;

        myb21 = myb21 + Myb21;
        myb22 = myb22 + Myb22;
        myb23 = myb23 + Myb23;
        mypsil = mypsil + Mypsil;
        
        mzb31 = mzb31 + Mzb31;
        mzb32 = mzb32 + Mzb32;
        mzb33 = mzb33 + Mzb33;
        mzpsil = mzpsil + Mzpsil;
        
    end
    

    out = zeros(7,1);

    JhdotMcdeltaJgOversdelta = mxb12+myb22+mzb32;
    psihdotsdeltaJg = (mxb11+myb21+mzb31);
    
    out[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;

    out[2] = mxb13+myb23+mzb33;     

    out[3] = cdelta*srpout[2] + sdelta*JhdotMcdeltaJgOversdelta;

    out[4] = pi/(2*mM1KK*Jg)*(mxpsil+(1-m)*mypsil+mzpsil);

    out[5] = 2/pi*sqrt(skm)*T4*k1r*out[4] + psih*sdelta*JhdotMcdeltaJgOversdelta/Jg;

    out[6] = -spsih*psihdotsdeltaJg + sdelta*cpsih*out[2] - cdelta*cpsih*JhdotMcdeltaJgOversdelta;
    
    out[7] =  cpsih*psihdotsdeltaJg + sdelta*spsih*out[2] - cdelta*spsih*JhdotMcdeltaJgOversdelta;

    out[abs.(out).<1e-15] .=0.0;

    return out;

end

function dragtorque_semianalytical_sadovelike(satellite,vmt,k,k1r,skm,smk,m,mM1KK,KK,T1,T2,T3,T4,T6,T7,T8,J2,J3,J6,J7)
    
    Jg = J2;
    cdelta = J3/Jg;
    if J6==J7 && J6==0.0
        sdelta = 0.0;
        cpsih = 1.0;
        spsih = 0.0;
    else
        sdelta = sqrt(J6^2+J7^2)/Jg;
        cpsih  = J6/sdelta/Jg;
        spsih  = J7/sdelta/Jg;
    end
    psih = mod(atan(spsih,cpsih),2.0*pi)

    
    VVT1 = satellite["constantstermsdragtorque"][1]
    VVT2 = satellite["constantstermsdragtorque"][2]
    VVT3 = satellite["constantstermsdragtorque"][3]
    VVT4 = satellite["constantstermsdragtorque"][4]

    skmBlock,TBblock1,TBblock2,TBblock3,TSblock = getdragsaaveragedtermsattitudeeq(k,k1r,skm,smk,m,KK,T1,T2,T3,T6,T7,T8,cdelta,sdelta,cpsih,spsih,vmt,VVT1,VVT2,VVT3,VVT4)
    out = zeros(7,1);

    out[1] = skmBlock/Jg;

    out[2] = TBblock3;     

    out[3] = cdelta*out[2] + sdelta*TBblock2;

    out[4] = pi/(2*mM1KK*Jg)*TSblock;

    out[5] = 2/pi*sqrt(skm)*T4*k1r*out[4] + psih*sdelta*TBblock2/Jg;

    out[6] = -spsih*TBblock1 + sdelta*cpsih*out[2] - cdelta*cpsih*TBblock2;
    
    out[7] =  cpsih*TBblock1 + sdelta*spsih*out[2] - cdelta*spsih*TBblock2;

    out[abs.(out).<1e-15] .=0.0;

    return out;

end

# semi-analitycal model with sadovlike variables
function attitudeMeanEv_semianalytical_triaxial_sadovlike(du,u,p,t)
    
    # inputs
    planetsunparams       = p[1];
    inertialrefframeinfo  = p[2];
    satellite             = p[3];
    draginfo              = p[4];
    settings              = p[5];
    mjd2000_0             = p[6];


    # orbit evolution Kepler
    mu  = planetsunparams["muPlanet"];
    sma = u[8];
    P1  = u[9];
    P2  = u[10];
    Q1 =  u[11];
    Q2 =  u[12];
    eta = sqrt(1.0-P1^2-P2^2);
    du[8] = 0.0; #sma
    du[9] = 0.0  #P1
    du[10] = 0.0  #P2
    du[11] = 0.0  #Q1
    du[12] = 0.0  #Q2
    du[13] = sqrt(mu/sma^3);

    # unperturbed attitude evolution
    IV = satellite["MomentsOfInertia"];
    k  = satellite["k_constant"];
    k1r = sqrt(1+k);
    A = IV[1]; C = IV[3]; 
    Jg  = u[2];
    skm = u[1];
    smk = 1-skm;

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
        KK  = ellF(pi/2,m); #Elliptic.K(m);
        PP  = ellP(-k,pi/2.0,m); #Elliptic.Pi(-k,pi/2.0,m);
        PPMsmkKK = (PP-smk*KK);
    end

    mM1KK = (m-1)*KK;
    if abs(mM1KK)<1e-13
        mM1KK = mM1KK +1e-10;
        println("qui mM1KK")
    end

    # println(t," - m: ",m," ",skm)
    
    du[1] = 0.0;
    du[2] = 0.0;
    du[3] = 0.0;
    du[4] = -pi/k1r*sqrt(skm)*Jg/2/A/C*(C-A)/KK;
    du[5] = Jg/A/C*(C-A)*PPMsmkKK/KK + Jg/A/C*(C*smk+A*skm);
    du[6] = 0.0;
    du[7] = 0.0;


    # perturbations for torque
    includeGravityTorque  = get(settings,"includeGravityTorque",false);
    includeMagneticTorque = get(settings,"includeMagneticTorque",false);
    includeSrpTorque      = get(settings,"includeSrpTorque",false);
    includeDragTorque     = get(settings,"includeDragTorque",false);
    higherordercorrections = settings["higherordercorrections"];

    # perturbations for orbit
    includeZonalHarmsAcc  = settings["includeZonalHarmsAcc"];
    includeSrpAcc         = settings["includeSrpAcc"];
    includeThirdBodyAcc   = settings["includeThirdBodyAcc"];
    includeSunGravityAcc  = settings["includeSunGravityAcc"];
    includeDragAcc         = settings["includeDragAcc"];

    if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque ||  includeZonalHarmsAcc || includeThirdBodyAcc || includeSunGravityAcc || includeSrpAcc || includeDragAcc
        B = IV[2];
        EE  = Elliptic.E(m);
        T1,T2,T3,T4,T6,T7,T8 = getsingularityterms(m,k,smk,skm,KK,EE,PP)

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
        
        if  includeSrpTorque ||    includeSrpAcc || includeThirdBodyAcc || includeSunGravityAcc
            mjd2000 = mjd2000_0 + t/24.0/3600.0;
        end

        cdelta = u[3]/Jg;
        if u[6]==u[7] && u[6]==0.0
            sdelta = 0.0;
            cpsih = 1.0;
            spsih = 0.0;
            psih = 0.0;
            # psig = mod(u[5],2.0*pi);
        else
            sdelta = sqrt(u[6]^2+u[7]^2)/Jg;
            cpsih  = u[6]/sdelta/Jg;
            spsih  = u[7]/sdelta/Jg;
            psih = mod(atan(spsih,cpsih),2.0*pi);
            # psig = mod(u[5]- cdelta*psih,2.0*pi);
        end
    end

    if includeGravityTorque 
        gtcontr = gravitytorque_girf_triaxial_sadovlike(k,k1r,skm,smk,m,mM1KK,T1,T3,T4,T6,Jg,u[3],u[6],u[7],sma,eta,q11,q12,q13,q21,q22,q23,A,B,C,mu);
        du[1:7] = du[1:7] + gtcontr;
    end

    if includeMagneticTorque
        muM = get(planetsunparams,"muM",0.0);
        IM = get(satellite,"intrinsicMagneticMoment",[0.0;0.0;0.0]);
        IMz = IM[3];
        cpa1 =  inertialrefframeinfo["cpa1"];
        cpa2 =  inertialrefframeinfo["cpa2"];
        cpa3 =  inertialrefframeinfo["cpa3"];
        mtcontr = magnetictorque_girf_triaxial_sadovlike(k,k1r,skm,m,mM1KK,T1,T4,KK,u[2],u[3],u[6],u[7],cpa1,cpa2,cpa3,sma,eta,q11,q12,q13,q21,q22,q23,muM,IMz);
        du[1:7] = du[1:7] + mtcontr;
    end

    if higherordercorrections && (includeGravityTorque || includeMagneticTorque)
        
        if includeMagneticTorque==false
            muM = 0.0;
            IMz = 0.0;
        end
        hoccontr = gethigherordercorrections(k,m,skm,smk,mM1KK,KK,EE,PP,T2,Jg,cdelta,sdelta,psih,sma,eta,P1,P2,Q1,Q2,A,B,C,mu,IMz,muM,inertialrefframeinfo["cpa1"],inertialrefframeinfo["cpa2"],inertialrefframeinfo["cpa3"],[includeGravityTorque,includeMagneticTorque],2);
        
        du[1:7] = du[1:7]+hoccontr;
   
    end

    if includeSrpTorque || includeSrpAcc || includeSunGravityAcc
        # sun ephemeris
        if settings["sunFixed"]
           
            rSunV = inertialrefframeinfo["ecliptic2inertial"]*settings["sunposition"];
        else
            rSunV = -inertialrefframeinfo["ecliptic2inertial"]*celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000)
        end

        if includeSrpTorque || includeSrpAcc
            # shadow
            includeEclipsesEffectsOnAttitude =  get(settings,"includeEclipsesEffectsOnAttitude",false);
            includeEclipsesEffectsOnOrbit    =  get(settings,"includeEclipsesEffectsOnOrbit",false);
            
            if  includeEclipsesEffectsOnAttitude || includeEclipsesEffectsOnOrbit 
                inshadow,TLin,TLout = eclipseinout(rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,planetsunparams["rPlanet"]);
                
                if inshadow
                    ELin  = truelong2ecclong(TLin,P1,P2);
                    ELout = truelong2ecclong(TLout,P1,P2);
                else
                    ELin = NaN;
                    ELout = NaN;
                end
            else
                inshadow = false;
                ELin = NaN;
                ELout = NaN;
            end
            # srp torque
            if includeSrpTorque
                srpcontr = srptorque_semianalytical_sadovelike(satellite,sma,P1,P2,q11,q12,q13,q21,q22,q23,rSunV,k,k1r,skm,smk,m,mM1KK,KK,T1,T2,T3,T4,T6,Jg,u[3],u[6],u[7], includeEclipsesEffectsOnAttitude,inshadow,ELin,ELout);
                du[1:7] = du[1:7] + srpcontr;
            end
            # srp acc
            if includeSrpAcc
                srpacccontr = srpaveragedgauss_wrapper1_AV(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffectsOnOrbit,inshadow,ELin,ELout);
                du[8:13] = du[8:13] + srpacccontr;
            end
        end

        if  includeSunGravityAcc
            sungcontr = tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rSunV,mu,0.19891000000000E+31 * 6.67259e-20,5); 
            du[8:13] = du[8:13] + sungcontr;
        end

    end

    if includeDragTorque
        
        if includeDragTorque
            dragcontr = dragtorque_semianalytical_sadovelike(satellite,draginfo["averageoverM_attitude"],k,k1r,skm,smk,m,mM1KK,KK,T1,T2,T3,T4,T6,T7,T8,u[2],u[3],u[6],u[7]); 
            du[1:7] = du[1:7] + dragcontr;
        end

        if includeDragAcc
            dragacccontr = dragveragedgauss(k,k1r,skm,smk,m,KK,T1,T2,T3,cdelta,sdelta,cpsih,spsih,sma,P1,P2,eta,Q1,Q2,GG,draginfo["averageoverM_orbit"],satellite,mu,draginfo["couplingstrategy"]);
            du[8:13] = du[8:13] + dragacccontr;
        end
    end

    if includeZonalHarmsAcc
        J2 = planetsunparams["zonalharmonicscoff"][1]
        J3 = planetsunparams["zonalharmonicscoff"][2]
        J4 = planetsunparams["zonalharmonicscoff"][3]
        J5 = planetsunparams["zonalharmonicscoff"][4]
        rPlanet = planetsunparams["rPlanet"]
        zharmcontr =  zonalharmonicsaveragedgauss(sma,eta,P1,P2,Q1,Q2,mu,rPlanet,J2,J3,J4,J5)
        du[8:13] = du[8:13] + zharmcontr;
    end

    if includeThirdBodyAcc 
        tbgcontr = zeros(6);
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
            tbgcontr = tbgcontr + tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rV_3b,mu,planetsunparams["perturbingbodiesgravparam"][jj],6)
        end
        du[8:13] = du[8:13] + tbgcontr;
    end

    return du;
end

# event handling
function conditionjump(u,t,integrator)
    fun = u[7] 
    return fun;
end

function affectjump!(integrator)
    if integrator.u[6]>=0.0 && abs(integrator.t)>1e-4
        terminate!(integrator);
    end
end

function conditionjumppatch(out,u,t,integrator)
    out[1] = u[7]
    skm = u[1]
    k = integrator.p[2]["k_constant"];
    if k!=0.0
        out[2] = skm-k/(k+0.9999);
    else
        out[2] = skm-1e-4;
    end
end

function affectjumppatch!(integrator,idx)
    if idx == 1
        if integrator.u[6]!=0.0 && abs(integrator.t)>1e-4
            terminate!(integrator);
        end
    elseif idx==2
        terminate!(integrator)
    end
end


##############################################################################################################################################################################
# AXISYMMETRIC BODIES (A=B=C)
##############################################################################################################################################################################

# perturbations

function magnetictorque_girf_axisym(csigma,ssigma,GA,cdelta,sdelta,clA,slA,chA,shA,cpa1,cpa2,cpa3,sma,eta,q11,q12,q13,q21,q22,q23,muM,IMx,IMy,IMz)
    out = zeros(6,1);
    
    LAdot_magnetic = -(3*ssigma*(-((q11*q12 + q21*q22)*cpa1 + (q12^2 + q22^2 - 2/3)*cpa2 + cpa3*(q12*q13 + q22*q23))*sdelta*chA + sdelta*((q21^2 + q11^2 - 2/3)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*shA + cdelta*((q11*q13 + q21*q23)*cpa1 + (q12*q13 + q22*q23)*cpa2 + cpa3*(-2/3 + q23^2 + q13^2)))*(IMx*clA - slA*IMy)*muM)/(2*sma^3*eta^3);
    HAdot_magnetic = -3*(((q21^2 + q11^2 - 2/3)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*chA + ((q11*q12 + q21*q22)*cpa1 + (q12^2 + q22^2 - 2/3)*cpa2 + cpa3*(q12*q13 + q22*q23))*shA)*sdelta*muM*(IMx*slA*ssigma + IMy*clA*ssigma + IMz*csigma)/(2*sma^3*eta^3);
    lAdot_magnetic = -3*(IMx*slA*csigma + IMy*clA*csigma - IMz*ssigma)*muM*(-sdelta*((q11*q12 + q21*q22)*cpa1 + (q22^2 + q12^2 - 2/3)*cpa2 + cpa3*(q12*q13 + q22*q23))*chA + sdelta*((q21^2 + q11^2 - 2/3)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*shA + cdelta*((q11*q13 + q21*q23)*cpa1 + (q12*q13 + q22*q23)*cpa2 + cpa3*(-2/3 + q23^2 + q13^2)))/(2*sma^3*eta^3*ssigma*GA);
    hAdot_magnetic = -3*(IMx*slA*ssigma + IMy*clA*ssigma + IMz*csigma)*(-((q11*q12 + q21*q22)*cpa1 + (q22^2 + q12^2 - 2/3)*cpa2 + cpa3*(q12*q13 + q22*q23))*cdelta*chA + ((q21^2 + q11^2 - 2/3)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*cdelta*shA - sdelta*((q11*q13 + q21*q23)*cpa1 + (q12*q13 + q22*q23)*cpa2 + cpa3*(-2/3 + q23^2 + q13^2)))*muM/(2*sma^3*eta^3*sdelta*GA);
    gAdot_magnetic = -csigma*lAdot_magnetic - cdelta*hAdot_magnetic;

    out[1] = out[1] + LAdot_magnetic;
    out[3] = out[3] + HAdot_magnetic;
    out[4] = out[4] + lAdot_magnetic;
    out[5] = out[5] + gAdot_magnetic;
    out[6] = out[6] + hAdot_magnetic;

    # checks
    for jj = 1 : 6
        if abs( out[jj])<1e-15
            out[jj] = 0.0;
        end
    end
   
    return out;
end

function srptorque_semianalytical_axisymm(satellite,sma,P1,P2,q11,q12,q13,q21,q22,q23,rSunV,csigma,ssigma,GA,cdelta,sdelta,clA,slA,chA,shA,includeEclipsesEffects,inshadow,ELin,ELout)
    mxb11 = 0.0;     mxb12 = 0.0;     mxb13 = 0.0;     mxl = 0.0;
    myb21 = 0.0;     myb22 = 0.0;     myb23 = 0.0;     myl = 0.0;
    mzb31 = 0.0;     mzb32 = 0.0;     mzb33 = 0.0;     mzL = 0.0;      
 
    numberOfFacets = get(satellite,"numberOfFacets",0);
    facets = get(satellite,"facets",0.0);  

    vcu = getAveragedScSunUnitVector(rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,0,includeEclipsesEffects,inshadow,ELin,ELout); 
    cMxb11T1,cMxb12T1,cMxb13T1,cMxlT1,cMyb21T1,cMyb22T1,cMyb23T1,cMylT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzT1 = srpanalyticalaveragedcoeff_T1_andoyer_cube(ssigma,csigma,cdelta,sdelta,clA,slA,chA,shA,vcu);
    cMxb11T2,cMxb12T2,cMxb13T2,cMxlT2,cMyb21T2,cMyb22T2,cMyb23T2,cMylT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzT2 = srpanalyticalaveragedcoeff_T2_andoyer_cube(ssigma,csigma,cdelta,sdelta,clA,slA,chA,shA,vcu);
    cMxb11T3,cMxb12T3,cMxb13T3,cMxlT3,cMyb21T3,cMyb22T3,cMyb23T3,cMylT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzT3 = srpanalyticalaveragedcoeff_T3_andoyer_cube(ssigma,csigma,cdelta,sdelta,clA,slA,chA,shA,vcu);


    for kk = 1:numberOfFacets   
        
        VVT1 = facets[kk][5];
        VVT2 = facets[kk][6];
        VVT3 = facets[kk][7];
        
        Mxb11T1,Mxb12T1,Mxb13T1,MxlT1,Myb21T1,Myb22T1,Myb23T1,MylT1,Mzb31T1,Mzb32T1,Mzb33T1,MzT1 = srp_averagedterms(VVT1,cMxb11T1,cMxb12T1,cMxb13T1,cMxlT1,cMyb21T1,cMyb22T1,cMyb23T1,cMylT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzT1);
        Mxb11T2,Mxb12T2,Mxb13T2,MxlT2,Myb21T2,Myb22T2,Myb23T2,MylT2,Mzb31T2,Mzb32T2,Mzb33T2,MzT2 = srp_averagedterms(VVT2,cMxb11T2,cMxb12T2,cMxb13T2,cMxlT2,cMyb21T2,cMyb22T2,cMyb23T2,cMylT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzT2);
        Mxb11T3,Mxb12T3,Mxb13T3,MxlT3,Myb21T3,Myb22T3,Myb23T3,MylT3,Mzb31T3,Mzb32T3,Mzb33T3,MzT3 = srp_averagedterms(VVT3,cMxb11T3,cMxb12T3,cMxb13T3,cMxlT3,cMyb21T3,cMyb22T3,cMyb23T3,cMylT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzT3);


        mxb11 = mxb11  + Mxb11T2 + Mxb11T3 + Mxb11T1
        mxb12 = mxb12 + Mxb12T2 + Mxb12T3 + Mxb12T1 
        mxb13 = mxb13  + Mxb13T2 + Mxb13T3 + Mxb13T1
        mxl = mxl + MxlT2 + MxlT3 + MxlT1 

        myb21 = myb21 + Myb21T2 + Myb21T3 + Myb21T1 
        myb22 = myb22 + Myb22T2 + Myb22T3 + Myb22T1 
        myb23 = myb23  + Myb23T2 + Myb23T3 + Myb23T1
        myl = myl + MylT2 + MylT3 + MylT1 
        
        mzb31 = mzb31  + Mzb31T2 + Mzb31T3 + Mzb31T1
        mzb32 = mzb32  + Mzb32T2 + Mzb32T3 + Mzb32T1
        mzb33 = mzb33  + Mzb33T2 + Mzb33T3 + Mzb33T1
        mzL = mzL  + MzT2 + MzT3 + MzT1
    end
    
    
    au = astroconstants(0)[2];
    earthsundist = norm(rSunV);
    psrpcorr     = (au/earthsundist)^2.0
    mxb11 = mxb11*psrpcorr
    mxb12 = mxb12*psrpcorr
    mxb13 = mxb13*psrpcorr
    mxl = mxl*psrpcorr
    myb21 = myb21*psrpcorr
    myb22 = myb22*psrpcorr
    myb23 = myb23*psrpcorr
    myl = myl*psrpcorr        
    mzb31 = mzb31*psrpcorr
    mzb32 = mzb32*psrpcorr
    mzb33 = mzb33*psrpcorr
    mzL = mzL*psrpcorr

    out = zeros(6,1);

    HdotMcdeltaGdotOversdelta = mxb12+myb22+mzb32;
    hdotsdeltaG = (mxb11+myb21+mzb31);
    
    out[1] = mzL;

    out[2] = mxb13+myb23+mzb33;     

    out[3] = cdelta*out[2] + sdelta*HdotMcdeltaGdotOversdelta;

    out[4] = mxl*clA/GA/ssigma-myl*slA/GA/ssigma;
    
    out[6] = hdotsdeltaG/GA/sdelta;
    
    out[5] = -csigma*out[4] -cdelta*out[6];

    out[abs.(out).<1e-15] .=0.0;

   return out

end

function dragtorque_semianalytical_axisymm(satellite,vmt,csigma,ssigma,GA,cdelta,sdelta,clA,slA,chA,shA)
    
    VVT1 = satellite["constantstermsdragtorque"][1]
    VVT2 = satellite["constantstermsdragtorque"][2]
    VVT3 = satellite["constantstermsdragtorque"][3]
    VVT4 = satellite["constantstermsdragtorque"][4]
  
    mxl,myl,mzL,TBblock1,TBblock2,TBblock3 = getdragsaaveragedtermsattitudeeq_andoyer_cube(ssigma,csigma,cdelta,sdelta,clA,slA,chA,shA,vmt,VVT1,VVT2,VVT3,VVT4)
    out = zeros(6,1);
   
    out[1] = mzL;

    out[2] = TBblock3;     

    out[3] = cdelta*out[2] + sdelta*TBblock2;

    out[4] = mxl*clA/GA/ssigma-myl*slA/GA/ssigma;
    
    out[6] = TBblock1/GA/sdelta;
    
    out[5] = -csigma*out[4] -cdelta*out[6];

    out[abs.(out).<1e-15] .=0.0;

   return out

end

# semi-analyitical model 

function attitudeMeanEv_semianalytical_axisymm(du,u,p,t)

    # inputs
    planetsunparams       = p[1];
    inertialrefframeinfo  = p[2];
    satellite             = p[3];
    draginfo              = p[4];
    settings              = p[5];
    mjd2000_0             = p[6];
    
    IV = get(satellite,"MomentsOfInertia",[0;0;0]);
    A = IV[1]; C = IV[3];

      
    # unperturbed orbit evolution (equinoctial) 
    mu  = planetsunparams["muPlanet"];
    sma = u[7];
    P1  = u[8];
    P2  = u[9];
    Q1 =  u[10];
    Q2 =  u[11];
    eta = sqrt(1-P1^2-P2^2);
    du[7] = 0.0; #sma
    du[8] = 0.0  #P1
    du[9] = 0.0  #P2
    du[10] = 0.0  #Q1
    du[11] = 0.0  #Q2
    du[12] = sqrt(mu/sma^3); 

    # unperturbed attitude evolution
    LA  = u[1];
    GA  = u[2];
    
    du[1] = 0.0;
    du[2] = 0.0;
    du[3] = 0.0;
    du[4] = 0.0; #LA*(1/C-1/A);
    du[5] = GA/A;
    du[6] = 0.0;
   
    # perturbations for torque
    includeMagneticTorque  = settings["includeMagneticTorque"];
    includeSrpTorque       = settings["includeSrpTorque"];
    includeDragTorque      = settings["includeDragTorque"];

    # perturbations for orbit
    includeZonalHarmsAcc  = settings["includeZonalHarmsAcc"];
    includeThirdBodyAcc   = settings["includeThirdBodyAcc"];
    includeSunGravityAcc  = settings["includeSunGravityAcc"];
    includeSrpAcc         = settings["includeSrpAcc"];
    includeDragAcc         = settings["includeDragAcc"];

    if includeMagneticTorque || includeSrpTorque || includeDragTorque || includeZonalHarmsAcc || includeThirdBodyAcc || includeSunGravityAcc || includeSrpAcc || includeDragAcc 
        csigma = LA/GA;
        ssigma = sqrt(1-csigma^2);
        cdelta = u[3]/GA;
        sdelta = sqrt(1-cdelta^2);
        chA = cos(u[6]);
        shA = sin(u[6]);
        clA = cos(u[4]);
        slA = sin(u[4]);
       
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

        if  includeSrpTorque ||    includeSrpAcc || includeThirdBodyAcc || includeSunGravityAcc
            mjd2000 = mjd2000_0 + t/24.0/3600.0;
        end
    end

    if includeMagneticTorque
        muM = get(planetsunparams,"muM",0.0);
        IM = get(satellite,"intrinsicMagneticMoment",[0.0;0.0;0.0]);
        IMx = IM[1]; IMy = IM[2]; IMz = IM[3];
        cpa1 =  inertialrefframeinfo["cpa1"];
        cpa2 =  inertialrefframeinfo["cpa2"];
        cpa3 =  inertialrefframeinfo["cpa3"];
        
        mtcontr = magnetictorque_girf_axisym(csigma,ssigma,GA,cdelta,sdelta,clA,slA,chA,shA,cpa1,cpa2,cpa3,sma,eta,q11,q12,q13,q21,q22,q23,muM,IMx,IMy,IMz)
        du[1] = du[1] + mtcontr[1];
        du[2] = du[2] + mtcontr[2];
        du[3] = du[3] + mtcontr[3];
        du[4] = du[4] + mtcontr[4];
        du[5] = du[5] + mtcontr[5];
        du[6] = du[6] + mtcontr[6];   
    end

    if includeSrpTorque || includeSrpAcc || includeSunGravityAcc
        # sun ephemeris
        rSunV = -inertialrefframeinfo["ecliptic2inertial"]*celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000)

        if includeSrpTorque || includeSrpAcc
            # shadow
            includeEclipsesEffectsOnAttitude =  get(settings,"includeEclipsesEffectsOnAttitude",false);
            includeEclipsesEffectsOnOrbit    =  get(settings,"includeEclipsesEffectsOnOrbit",false);
            if includeEclipsesEffectsOnAttitude ||  includeEclipsesEffectsOnOrbit
                inshadow,TLin,TLout = eclipseinout(rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,planetsunparams["rPlanet"]);
                if inshadow
                    ELin  = truelong2ecclong(TLin,P1,P2);
                    ELout = truelong2ecclong(TLout,P1,P2);
                else
                    ELin = NaN;
                    ELout = NaN;
                end
            else
                inshadow = false;
                ELin = NaN;
                ELout = NaN;
            end
            if includeSrpTorque
                srpcontr = srptorque_semianalytical_axisymm(satellite,sma,P1,P2,q11,q12,q13,q21,q22,q23,rSunV,csigma,ssigma,GA,cdelta,sdelta,clA,slA,chA,shA,includeEclipsesEffectsOnAttitude,inshadow,ELin,ELout);
                du[1:6] = du[1:6] + srpcontr;
            end
            if includeSrpAcc
                srpacccontr = srpaveragedgauss_wrapper2_AV(sma,eta,P1,P2,Q1,Q2,GG,satellite,csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA,mu,rSunV, includeEclipsesEffectsOnOrbit,inshadow,ELin,ELout);
                du[7:12] = du[7:12] + srpacccontr;
            end
        end

        if  includeSunGravityAcc
            sungcontr = tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rSunV,mu,0.19891000000000E+31 * 6.67259e-20,5); 
            du[7:12] = du[7:12] + sungcontr;
        end

    end

    if includeDragTorque
        if includeDragTorque
            dragcontr =  dragtorque_semianalytical_axisymm(satellite,draginfo["averageoverM_attitude"],csigma,ssigma,GA,cdelta,sdelta,clA,slA,chA,shA);
            du[1:6] = du[1:6] + dragcontr; 
        end    

        if includeDragAcc
            dragacccontr = dragveragedgauss_cube(csigma,ssigma,clA,slA,cdelta,sdelta,chA,shA,sma,P1,P2,eta,Q1,Q2,GG,draginfo["averageoverM_orbit"],satellite,mu,draginfo["couplingstrategy"]);
            du[7:12] = du[7:12] + dragacccontr;
        end
    end

    if includeZonalHarmsAcc
        J2 = planetsunparams["zonalharmonicscoff"][1]
        J3 = planetsunparams["zonalharmonicscoff"][2]
        J4 = planetsunparams["zonalharmonicscoff"][3]
        J5 = planetsunparams["zonalharmonicscoff"][4]
        rPlanet = planetsunparams["rPlanet"]
        zharmcontr =  zonalharmonicsaveragedgauss(sma,eta,P1,P2,Q1,Q2,mu,rPlanet,J2,J3,J4,J5)
        du[7:12] = du[7:12] + zharmcontr;
    end

    if includeThirdBodyAcc 
        tbgcontr = zeros(6);
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
            tbgcontr = tbgcontr + tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rV_3b,mu,planetsunparams["perturbingbodiesgravparam"][jj],6)
        end
        du[7:12] = du[7:12] + tbgcontr;
    end

    return du;


end


##############################################################################################################################################################################
# AXISYMMETRIC BODIES (A=B=C), handling singularity sin(delta) = 0
##############################################################################################################################################################################

# perturbations
function magnetictorque_girf_axisym_andoyerlike(csigma,ssigma,JA2,JA3,cJA4,sJA4,JA6,JA7,cpa1,cpa2,cpa3,sma,eta,q11,q12,q13,q21,q22,q23,muM,IMx,IMy,IMz)
    out = zeros(7,1);
    if JA7==0.0 && JA6==0.0
        atanval = 0.0;
    else
        atanval = mod(atan(JA7,JA6),2.0*pi)
    end

    JA1dot_magnetic = 3*(((-q11*q13 - q21*q23)*JA3 + (q11*q12 + q21*q22)*JA6 - JA7*(-2/3 + q21^2 + q11^2))*cpa1 + ((-q12*q13 - q22*q23)*JA3 + (-2/3 + q12^2 + q22^2)*JA6 - JA7*(q11*q12 + q21*q22))*cpa2 + cpa3*((-q13^2 - q23^2 + 2/3)*JA3 + (q12*q13 + q22*q23)*JA6 - JA7*(q11*q13 + q21*q23)))*muM*ssigma*(IMx*cJA4 - IMy*sJA4)/(2*JA2*sma^3*eta^3);
    JA3dot_magnetic = -(3*muM*(((-2/3 + q21^2 + q11^2)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*JA6 + JA7*((q11*q12 + q21*q22)*cpa1 + (-2/3 + q12^2 + q22^2)*cpa2 + cpa3*(q12*q13 + q22*q23)))*((IMx*sJA4 + IMy*cJA4)*ssigma + IMz*csigma))/(2*JA2*sma^3*eta^3);
    JA4dot_magnetic = -3*((IMx*sJA4 + IMy*cJA4)*csigma - IMz*ssigma)*muM*(((q11*q13 + q21*q23)*JA3 + (-q11*q12 - q21*q22)*JA6 + JA7*(-2/3 + q21^2 + q11^2))*cpa1 + ((q12*q13 + q22*q23)*JA3 + (2/3 - q12^2 - q22^2)*JA6 + JA7*(q11*q12 + q21*q22))*cpa2 + cpa3*((q23^2 + q13^2 - 2/3)*JA3 + (-q12*q13 - q22*q23)*JA6 + JA7*(q11*q13 + q21*q23)))/(2*JA2^2*eta^3*sma^3*ssigma);
    JA5dot_magnetic = -csigma*JA4dot_magnetic - 3*((IMx*sJA4 + IMy*cJA4)*ssigma + IMz*csigma)*(((-2/3 + q21^2 + q11^2)*cpa1 + (q11*q12 + q21*q22)*cpa2 + cpa3*(q11*q13 + q21*q23))*JA6 + JA7*((q11*q12 + q21*q22)*cpa1 + (-2/3 + q12^2 + q22^2)*cpa2 + cpa3*(q12*q13 + q22*q23)))*atanval*muM/(2*JA2^2*eta^3*sma^3);
    JA6dot_magnetic = -3*((IMx*sJA4 + IMy*cJA4)*ssigma + IMz*csigma)*(((-q11^2 - q21^2 + 2/3)*cpa1 + (-q11*q13 - q21*q23)*cpa3 - (q11*q12 + q21*q22)*cpa2)*JA3 + ((q11*q13 + q21*q23)*cpa1 + (q12*q13 + q22*q23)*cpa2 + cpa3*(q23^2 + q13^2 - 2/3))*JA7)*muM/(2*JA2*sma^3*eta^3);
    JA7dot_magnetic = 3*((IMx*sJA4 + IMy*cJA4)*ssigma + IMz*csigma)*muM*(((q11*q12 + q21*q22)*cpa1 + (-2/3 + q12^2 + q22^2)*cpa2 + cpa3*(q12*q13 + q22*q23))*JA3 + JA6*((q11*q13 + q21*q23)*cpa1 + (q12*q13 + q22*q23)*cpa2 + cpa3*(q23^2 + q13^2 - 2/3)))/(2*JA2*sma^3*eta^3);
    

    out[1] = out[1] +  JA1dot_magnetic;
    out[3] = out[3] +  JA3dot_magnetic;
    out[4] = out[4] +  JA4dot_magnetic;
    out[5] = out[5] +  JA5dot_magnetic;
    out[6] = out[6] +  JA6dot_magnetic;
    out[7] = out[7] +  JA7dot_magnetic;
   
    # checks
    for jj = 1 : 7
        if abs( out[jj])<1e-15
            out[jj] = 0.0;
        end
    end

    return out;
end

function srptorque_semianalytical_axisymm_andoyerlike(satellite,sma,P1,P2,q11,q12,q13,q21,q22,q23,rSunV,csigma,ssigma,JA2,JA3,cJA4,sJA4,JA6,JA7,includeEclipsesEffects,inshadow,ELin,ELout)

    GA = JA2;
    cdelta = JA3/GA;
    if JA6==JA7 && JA6==0.0
        sdelta = 0.0;
        chA = 1.0;
        shA = 0.0;
    else
        sdelta = sqrt(JA6^2+JA7^2)/GA;
        chA  = JA6/sdelta/GA;
        shA  = JA7/sdelta/GA;
    end
    hA = mod(atan(shA,chA),2.0*pi)
    clA = cJA4;
    slA = sJA4;

    mxb11 = 0.0;     mxb12 = 0.0;     mxb13 = 0.0;     mxl = 0.0;
    myb21 = 0.0;     myb22 = 0.0;     myb23 = 0.0;     myl = 0.0;
    mzb31 = 0.0;     mzb32 = 0.0;     mzb33 = 0.0;     mzL = 0.0;      
 
    numberOfFacets = get(satellite,"numberOfFacets",0);
    facets = get(satellite,"facets",0.0);  

    vcu = getAveragedScSunUnitVector(rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,0,includeEclipsesEffects,inshadow,ELin,ELout); 
    cMxb11T1,cMxb12T1,cMxb13T1,cMxlT1,cMyb21T1,cMyb22T1,cMyb23T1,cMylT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzT1 = srpanalyticalaveragedcoeff_T1_andoyer_cube(ssigma,csigma,cdelta,sdelta,clA,slA,chA,shA,vcu);
    cMxb11T2,cMxb12T2,cMxb13T2,cMxlT2,cMyb21T2,cMyb22T2,cMyb23T2,cMylT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzT2 = srpanalyticalaveragedcoeff_T2_andoyer_cube(ssigma,csigma,cdelta,sdelta,clA,slA,chA,shA,vcu);
    cMxb11T3,cMxb12T3,cMxb13T3,cMxlT3,cMyb21T3,cMyb22T3,cMyb23T3,cMylT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzT3 = srpanalyticalaveragedcoeff_T3_andoyer_cube(ssigma,csigma,cdelta,sdelta,clA,slA,chA,shA,vcu);

    for kk = 1:numberOfFacets   
        
        VVT1 = facets[kk][5];
        VVT2 = facets[kk][6];
        VVT3 = facets[kk][7];
        
        Mxb11T1,Mxb12T1,Mxb13T1,MxlT1,Myb21T1,Myb22T1,Myb23T1,MylT1,Mzb31T1,Mzb32T1,Mzb33T1,MzT1 = srp_averagedterms(VVT1,cMxb11T1,cMxb12T1,cMxb13T1,cMxlT1,cMyb21T1,cMyb22T1,cMyb23T1,cMylT1,cMzb31T1,cMzb32T1,cMzb33T1,cMzT1);
        Mxb11T2,Mxb12T2,Mxb13T2,MxlT2,Myb21T2,Myb22T2,Myb23T2,MylT2,Mzb31T2,Mzb32T2,Mzb33T2,MzT2 = srp_averagedterms(VVT2,cMxb11T2,cMxb12T2,cMxb13T2,cMxlT2,cMyb21T2,cMyb22T2,cMyb23T2,cMylT2,cMzb31T2,cMzb32T2,cMzb33T2,cMzT2);
        Mxb11T3,Mxb12T3,Mxb13T3,MxlT3,Myb21T3,Myb22T3,Myb23T3,MylT3,Mzb31T3,Mzb32T3,Mzb33T3,MzT3 = srp_averagedterms(VVT3,cMxb11T3,cMxb12T3,cMxb13T3,cMxlT3,cMyb21T3,cMyb22T3,cMyb23T3,cMylT3,cMzb31T3,cMzb32T3,cMzb33T3,cMzT3);


        mxb11 = mxb11  + Mxb11T2 + Mxb11T3 + Mxb11T1
        mxb12 = mxb12 + Mxb12T2 + Mxb12T3 + Mxb12T1 
        mxb13 = mxb13  + Mxb13T2 + Mxb13T3 + Mxb13T1
        mxl = mxl + MxlT2 + MxlT3 + MxlT1 

        myb21 = myb21 + Myb21T2 + Myb21T3 + Myb21T1 
        myb22 = myb22 + Myb22T2 + Myb22T3 + Myb22T1 
        myb23 = myb23  + Myb23T2 + Myb23T3 + Myb23T1
        myl = myl + MylT2 + MylT3 + MylT1 
        
        mzb31 = mzb31  + Mzb31T2 + Mzb31T3 + Mzb31T1
        mzb32 = mzb32  + Mzb32T2 + Mzb32T3 + Mzb32T1
        mzb33 = mzb33  + Mzb33T2 + Mzb33T3 + Mzb33T1
        mzL = mzL  + MzT2 + MzT3 + MzT1
    end
   
    au = astroconstants(0)[2];
    earthsundist = norm(rSunV);
    psrpcorr     = (au/earthsundist)^2.0
    mxb11 = mxb11*psrpcorr
    mxb12 = mxb12*psrpcorr
    mxb13 = mxb13*psrpcorr
    mxl = mxl*psrpcorr
    myb21 = myb21*psrpcorr
    myb22 = myb22*psrpcorr
    myb23 = myb23*psrpcorr
    myl = myl*psrpcorr        
    mzb31 = mzb31*psrpcorr
    mzb32 = mzb32*psrpcorr
    mzb33 = mzb33*psrpcorr
    mzL = mzL*psrpcorr

    out = zeros(7,1);

    HdotMcdeltaGdotOversdelta = mxb12+myb22+mzb32;
    hdotsdeltaG = (mxb11+myb21+mzb31);
    
    out[1] = mzL;

    out[2] = mxb13+myb23+mzb33;     

    out[3] = cdelta*out[2] + sdelta*HdotMcdeltaGdotOversdelta;

    out[4] = mxl*clA/GA/ssigma-myl*slA/GA/ssigma;

    out[5] = -csigma*out[4] + hA*sdelta*HdotMcdeltaGdotOversdelta/GA;

    out[6] = -shA*hdotsdeltaG  + sdelta*chA*out[2] - cdelta*chA*HdotMcdeltaGdotOversdelta;
    
    out[7] =  chA*hdotsdeltaG  + sdelta*shA*out[2] - cdelta*shA*HdotMcdeltaGdotOversdelta;

    
    out[abs.(out).<1e-15].=0.0;

    return out

end

function dragtorque_semianalytical_axisymm_andoyerlike(satellite,vmt,csigma,ssigma,JA2,JA3,cJA4,sJA4,JA6,JA7)

    GA = JA2;
    cdelta = JA3/GA;
    if JA6==JA7 && JA6==0.0
        sdelta = 0.0;
        chA = 1.0;
        shA = 0.0;
    else
        sdelta = sqrt(JA6^2+JA7^2)/GA;
        chA  = JA6/sdelta/GA;
        shA  = JA7/sdelta/GA;
    end
    hA = mod(atan(shA,chA),2.0*pi)
    clA = cJA4;
    slA = sJA4;

        
    VVT1 = satellite["constantstermsdragtorque"][1]
    VVT2 = satellite["constantstermsdragtorque"][2]
    VVT3 = satellite["constantstermsdragtorque"][3]
    VVT4 = satellite["constantstermsdragtorque"][4]
  
    mxl,myl,mzL,TBblock1,TBblock2,TBblock3 = getdragsaaveragedtermsattitudeeq_andoyer_cube(ssigma,csigma,cdelta,sdelta,clA,slA,chA,shA,vmt,VVT1,VVT2,VVT3,VVT4)
   
    out = zeros(7,1);

    
    out[1] = mzL;

    out[2] = TBblock3;     

    out[3] = cdelta*out[2] + sdelta*TBblock2;

    out[4] = mxl*clA/GA/ssigma-myl*slA/GA/ssigma;

    out[5] = -csigma*out[4] + hA*sdelta*TBblock2/GA;

    out[6] = -shA*TBblock1  + sdelta*chA*out[2] - cdelta*chA*TBblock2;
    
    out[7] =  chA*TBblock1  + sdelta*shA*out[2] - cdelta*shA*TBblock2;

    
    out[abs.(out).<1e-15].=0.0;

    return out

end

# semi-analyitical model 
function attitudeMeanEv_semianalytical_axisymm_andoyerlike(du,u,p,t)
    # inputs
    
    planetsunparams        = p[1];
    inertialrefframeinfo   = p[2];
    satellite              = p[3];
    draginfo               = p[4];
    settings               = p[5];
    mjd2000_0              = p[6];

    IV = get(satellite,"MomentsOfInertia",[0;0;0]);
    A = IV[1]; 

    # unperturbed orbit evolution (equinoctial) 
    mu  = planetsunparams["muPlanet"];
    sma = u[8];
    P1  = u[9];
    P2  = u[10];
    Q1 =  u[11];
    Q2 =  u[12];
    eta = sqrt(1-P1^2-P2^2);
    du[8] = 0.0; #sma
    du[9] = 0.0  #P1
    du[10] = 0.0  #P2
    du[11] = 0.0  #Q1
    du[12] = 0.0  #Q2
    du[13] = sqrt(mu/sma^3);

    # unperturbed attitude evolution
    LA  = u[1];
    GA  = u[2];
    
    du[1] = 0.0;
    du[2] = 0.0;
    du[3] = 0.0;
    du[4] = 0.0; #LA*(1/C-1/A);
    du[5] = GA/A;
    du[6] = 0.0;
    du[7] = 0.0;

    # attitude perturbations
    # perturbations for torque
    includeMagneticTorque  = settings["includeMagneticTorque"];
    includeSrpTorque       = settings["includeSrpTorque"];
    includeDragTorque      = settings["includeDragTorque"];
   
    # perturbations for orbit
    includeZonalHarmsAcc  = settings["includeZonalHarmsAcc"];
    includeThirdBodyAcc   = settings["includeThirdBodyAcc"];
    includeSunGravityAcc  = settings["includeSunGravityAcc"];
    includeSrpAcc         = settings["includeSrpAcc"];
    includeDragAcc         = settings["includeDragAcc"];

    if includeMagneticTorque || includeSrpTorque || includeDragTorque || includeZonalHarmsAcc || includeThirdBodyAcc || includeSunGravityAcc || includeSrpAcc || includeDragAcc 
        csigma = LA/GA;
        ssigma = sqrt(1.0-csigma^2.0);
        clA = cos(u[4]);
        slA = sin(u[4]);
       
        GG = 1.0+Q1^2+Q2^2.0;
        q11 = (1.0-Q1^2.0+Q2^2.0)/GG;
        q12 = (2.0*Q1*Q2)/GG;
        q13 = (-2.0*Q1)/GG;
        q21 = (2.0*Q1*Q2)/GG;
        q22 = (1.0+Q1^2.0-Q2^2.0)/GG;
        q23 = (2.0*Q2)/GG;
        q31 = (2.0*Q1)/GG;
        q32 = (-2.0*Q2)/GG;
        q33 = (1.0-Q1^2.0-Q2^2.0)/GG;

        if  includeSrpTorque ||    includeSrpAcc || includeThirdBodyAcc || includeSunGravityAcc
            mjd2000 = mjd2000_0 + t/24.0/3600.0;
        end

    end

    if includeSrpAcc || includeDragAcc
        cdelta = u[3]/GA;
        if u[6]==u[7] && u[6]==0.0
            sdelta = 0.0
            chA = 1.0;
            shA = 0.0;
            gA = mod(u[5],2*pi);
        else
            sdelta = sqrt(u[7]^2+u[6]^2)/GA
            chA = u[6]/sdelta/GA
            shA = u[7]/sdelta/GA
        end
    end

    if includeMagneticTorque
        muM = get(planetsunparams,"muM",0.0);
        IM = get(satellite,"intrinsicMagneticMoment",[0.0;0.0;0.0]);
        IMx = IM[1]; IMy = IM[2]; IMz = IM[3];
        mtcontr = magnetictorque_girf_axisym_andoyerlike(csigma,ssigma,GA,u[3],clA,slA,u[6],u[7],inertialrefframeinfo["cpa1"],inertialrefframeinfo["cpa2"],inertialrefframeinfo["cpa3"],sma,eta,q11,q12,q13,q21,q22,q23,muM,IMx,IMy,IMz)
        du[1:7] = du[1:7] + mtcontr;           
    end

    if includeSrpTorque || includeSrpAcc || includeSunGravityAcc
        # sun ephemeris
        rSunV = -inertialrefframeinfo["ecliptic2inertial"]*celestialbodiesephemeris_position(planetsunparams["centralBodyIDX"],mjd2000)

        if includeSrpTorque || includeSrpAcc
            # shadow
            includeEclipsesEffectsOnAttitude =  get(settings,"includeEclipsesEffectsOnAttitude",false);
            includeEclipsesEffectsOnOrbit    =  get(settings,"includeEclipsesEffectsOnOrbit",false);
            if includeEclipsesEffectsOnAttitude || includeEclipsesEffectsOnOrbit
                inshadow,TLin,TLout = eclipseinout(rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,planetsunparams["rPlanet"]);
                if inshadow
                    ELin  = truelong2ecclong(TLin,P1,P2);
                    ELout = truelong2ecclong(TLout,P1,P2);
                else
                    ELin = NaN;
                    ELout = NaN;
                end
            else
                inshadow = false;
                ELin = NaN;
                ELout = NaN;
            end
            if includeSrpTorque
                srpcontr = srptorque_semianalytical_axisymm_andoyerlike(satellite,sma,P1,P2,q11,q12,q13,q21,q22,q23,rSunV,csigma,ssigma,GA,u[3],clA,slA,u[6],u[7],includeEclipsesEffectsOnAttitude,inshadow,ELin,ELout);
                du[1:7] = du[1:7] + srpcontr;
            end
            if includeSrpAcc
                srpacccontr = srpaveragedgauss_wrapper2_AV(sma,eta,P1,P2,Q1,Q2,GG,satellite,csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA,mu,rSunV,includeEclipsesEffectsOnOrbit,inshadow,ELin,ELout);
                du[8:13] = du[8:13] + srpacccontr;
            end
        end

        if  includeSunGravityAcc
            sungcontr = tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rSunV,mu,0.19891000000000E+31 * 6.67259e-20,5); 
            du[8:13] = du[8:13] + sungcontr;
        end

    end 
  
    if includeDragTorque     
            
        if includeDragTorque
            dragcontr = dragtorque_semianalytical_axisymm_andoyerlike(satellite,draginfo["averageoverM_attitude"],csigma,ssigma,GA,u[3],clA,slA,u[6],u[7]);
            du[1:7] = du[1:7] + dragcontr;
        end 

        if includeDragAcc
            dragacccontr = dragveragedgauss_cube(csigma,ssigma,clA,slA,cdelta,sdelta,chA,shA,sma,P1,P2,eta,Q1,Q2,GG,draginfo["averageoverM_orbit"],satellite,mu,draginfo["couplingstrategy"]);
            du[8:13] = du[8:13] + dragacccontr;
        end
    end
         
    if includeZonalHarmsAcc
        J2 = planetsunparams["zonalharmonicscoff"][1]
        J3 = planetsunparams["zonalharmonicscoff"][2]
        J4 = planetsunparams["zonalharmonicscoff"][3]
        J5 = planetsunparams["zonalharmonicscoff"][4]
        rPlanet = planetsunparams["rPlanet"]
        zharmcontr =  zonalharmonicsaveragedgauss(sma,eta,P1,P2,Q1,Q2,mu,rPlanet,J2,J3,J4,J5)
        du[8:13] = du[8:13] + zharmcontr;
    end

    if includeThirdBodyAcc 
        tbgcontr = zeros(6);
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
            tbgcontr = tbgcontr + tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rV_3b,mu,planetsunparams["perturbingbodiesgravparam"][jj],6)
        end
        du[8:13] = du[8:13] + tbgcontr;
    end
    
   
    return du;


end





# ################################ IGNORE !!!!!!!!!!!!!!!!!!!!!!!
# #################################
# function vectorfield(psig,skm,Jg,Jh,psil,psih,TL,k,IV,settings,inertialrefframeinfo,planetsunparams,numberOfFacets,facets,rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,EL1,EL2)
#     # skminterp,Jginterp,Jhinterp,psilinterp,psihinterp,TLinterp
   
#     A = IV[1]; B = IV[2]; C = IV[3];
#     k1r = sqrt(1+k);
#     smk = 1.0-skm
#     m   = smk/skm*k
#     KK  = ellF(pi/2.0,m);
#     PP  = ellP(-k,pi/2.0,m);

#     fun = zeros(6,1);
#     fun[1] = 0.0;
#     fun[2] = 0.0;
#     fun[3] = 0.0;
#     fun[4] = -pi/k1r*sqrt(skm)*Jg/2/A/C*(C-A)/KK;
#     fun[5] = Jg/A/C*(C-A)*(PP-smk*KK)/KK + Jg/A/C*(C*smk+A*skm);
#     fun[6] = 0.0;
    
#     includeGravityTorque  = get(settings,"includeGravityTorque",false);
#     includeMagneticTorque = get(settings,"includeMagneticTorque",false);
#     includeSrpTorque      = get(settings,"includeSrpTorque",false);
#     includeDragTorque     = get(settings,"includeDragTorque",false);
#     torquepert = zeros(6,1);
    
#     if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque
#         skmR = sqrt(skm);
#         smkR = sqrt(smk);
      
#         eF  = 2.0*KK*psil/pi;
#         sn = Elliptic.Jacobi.sn(eF,m);
#         cn = Elliptic.Jacobi.cn(eF,m);
#         dn = Elliptic.Jacobi.dn(eF,m);

#         lambda = mod(atan(sn,cn),2.0*pi);
#         eP  = ellP(-k,lambda,m) 
#         eE  = ellE(lambda,m);
#         EE  = ellE(pi/2.0,m)
        
#         g  = mod(psig + k1r/skmR*(eF/KK*PP-eP),2*pi);
#         cg = cos(g);
#         sg = sin(g);

#         cdelta = Jh/Jg;
#         sdelta = sqrt(1-cdelta^2);
        
#         cpsih = cos(psih);
#         spsih = sin(psih);

#         Ri2b,Rpsilpsig = i2brotmat_sadov(k,k1r,skmR,smkR,cn,sn,dn,cg,sg,cpsih,spsih,cdelta,sdelta);

#         rnorm = sma*(1-P1^2.0-P2^2.0)/(1+P1*sin(TL)+P2*cos(TL));

#         rUVI = [q11*cos(TL)+q21*sin(TL), q12*cos(TL)+q22*sin(TL), q13*cos(TL)+q23*sin(TL)];
#         rUV  = Ri2b*rUVI;

#         mu = planetsunparams["muPlanet"]

#         gtorque = [0.0;0.0;0.0];
#         if includeGravityTorque
#             gtorque[1] = 3.0*mu/(rnorm^3)*(C-B)*rUV[2]*rUV[3];
#             gtorque[2] = 3.0*mu/(rnorm^3)*(A-C)*rUV[3]*rUV[1];
#             gtorque[3] = 3.0*mu/(rnorm^3)*(B-A)*rUV[1]*rUV[2];
#         end

#         mtorque = [0.0;0.0;0.0];
#         if includeMagneticTorque
#             muM = planetsunparams["muM"];
#             Zeq = Ri2b*[inertialrefframeinfo["cpa1"];inertialrefframeinfo["cpa2"];inertialrefframeinfo["cpa3"]]
#             Hb  = (Zeq - 3.0*(Zeq[1]*rUV[1]+Zeq[2]*rUV[2]+Zeq[3]*rUV[3])*rUV)*muM/rnorm^3;
#             mtorque = LinearAlgebra.cross(IM,Hb);
#         end

#         srptorque = [0.0;0.0;0.0];
#         if includeSrpTorque 
#             # EL = truelong2ecclong(TL,P1,P2);
#             # ETilde = EL - (EL1+EL2)/2.0
#             # Esh    = (EL2-EL1)/2.0

#             # rPlanet = planetsunparams["rPlanet"];
#             # sc = (rSunV[1]*rV[1]+rSunV[2]*rV[2]+rSunV[3]*rV[3])/rsunnorm + sqrt(rnorm^2-rPlanet^2);
#             # shadowfun = 1/2*(1+tanh(1000*sc));
#             # shadowfun = 1.0-Esh/pi
#             # for kk=1:20
#             #     shadowfun = shadowfun - 2/kk/pi*sin(kk*Esh)*cos(kk*ETilde);
#             # end
#             shadowfun = 1.0;

#             uVect = Ri2b*rSunV; #Ri2b*(rSunV-rsnorm*rUVI);
#             ssd = norm(uVect); 
#             uVect = uVect/ssd;

#             for kk = 1:numberOfFacets 
#                 surfi = facets[kk][2];
#                 niv = facets[kk][4];
#                 rhoiv = facets[kk][3];
#                 niudot = niv[1]*uVect[1]+niv[2]*uVect[2]+niv[3]*uVect[3];
#                 gifun = max(0,niudot);
                
#                 cai = facets[kk][1][1];
#                 cdi = facets[kk][1][2];
#                 csi = facets[kk][1][3];
        
#                 srpforcekk_T1 = - (surfi*gifun*cai*uVect)*shadowfun;
#                 srpforcekk_T2 = - (surfi*gifun*cdi*niv)*shadowfun;
#                 srpforcekk_T3 = - (surfi*gifun*csi*niudot*niv)*shadowfun;

#                 srptorquekk_T1 = LinearAlgebra.cross(rhoiv,srpforcekk_T1);
#                 srptorquekk_T2 = LinearAlgebra.cross(rhoiv,srpforcekk_T2);
#                 srptorquekk_T3 = LinearAlgebra.cross(rhoiv,srpforcekk_T3);
#                 srptorquekk = srptorquekk_T1 + srptorquekk_T2 + srptorquekk_T3; 
#                 srptorque = srptorque + srptorquekk;
#             end        
#         end

#         exttorque = gtorque + mtorque + srptorque;
        
#         zn = eE-eF*EE/KK
#         T4 = (KK*m - PP*(m+k))/k;
#         mxb12 = exttorque[1]*Rpsilpsig[1,2];
#         myb22 = exttorque[2]*Rpsilpsig[2,2];
#         mzb32 = exttorque[3]*Rpsilpsig[3,2];
#         mxb13 = exttorque[1]*Rpsilpsig[1,3];        
#         myb23 = exttorque[2]*Rpsilpsig[2,3];       
#         mzb33 = exttorque[3]*Rpsilpsig[3,3];
#         mxdnsnMcnzn = exttorque[1]*(dn*sn-cn*zn);
#         mydncnPsnzn = exttorque[2]*(cn*dn+zn*sn)
#         mzdnznMmsncn = exttorque[3]*(dn*zn-m*sn*cn)
#         mxb11 = exttorque[1]*Rpsilpsig[1,1];
#         myb21 = exttorque[2]*Rpsilpsig[2,1];
#         mzb31 = exttorque[3]*Rpsilpsig[3,1];

#         torquepert  = zeros(6,1);

#         torquepert[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
    
#         torquepert[2] = mxb13+myb23+mzb33;

#         torquepert[3] = cdelta*torquepert[2] + sdelta*(mxb12+myb22+mzb32);

#         torquepert[4] = -pi/(2*Jg*KK*(1-m))*(mxdnsnMcnzn/smkR+
#         (1-m)*mydncnPsnzn/k1r/smkR+
#         mzdnznMmsncn/skmR);
    
#         torquepert[6] = (mxb11+myb21+mzb31)/Jg/sdelta;
            
#         torquepert[5] = -cdelta*torquepert[6] + 2/pi*sqrt(skm)*T4*k1r*torquepert[4] ;

#     end

#     fun = fun + torquepert;
#     return fun
# end

# function vectorfieldV2(psil,skm,Jg,Jhinterp,psiginterp,psihinterp,TLinterp,k,IV,settings,inertialrefframeinfo,planetsunparams,numberOfFacets,facets,rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,EL1,EL2)
#     # skminterp,Jginterp,Jhinterp,psilinterp,psihinterp,TLinterp
#    # skm  = skminterp(psil)
#    # Jg   = Jginterp(psil)
#    Jh   = Jhinterp(psil)
#     psig = mod(psiginterp(psil),2.0*pi)
#    psih = mod(psihinterp(psil),2.0*pi)
#    TL   = mod(TLinterp(psil),2.0*pi)
    

#     A = IV[1]; B = IV[2]; C = IV[3];
#     k1r = sqrt(1+k);
#     smk = 1.0-skm
#     m   = smk/skm*k
#     KK  = ellF(pi/2.0,m);
#     PP  = ellP(-k,pi/2.0,m);

#     fun = zeros(6,1);
#     fun[1] = 0.0;
#     fun[2] = 0.0;
#     fun[3] = 0.0;
#     fun[4] = -pi/k1r*sqrt(skm)*Jg/2/A/C*(C-A)/KK;
#     fun[5] = Jg/A/C*(C-A)*(PP-smk*KK)/KK + Jg/A/C*(C*smk+A*skm);
#     fun[6] = 0.0;
    
#     includeGravityTorque  = get(settings,"includeGravityTorque",false);
#     includeMagneticTorque = get(settings,"includeMagneticTorque",false);
#     includeSrpTorque      = get(settings,"includeSrpTorque",false);
#     includeDragTorque     = get(settings,"includeDragTorque",false);
#     torquepert = zeros(6,1);
    
#     if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque
#         skmR = sqrt(skm);
#         smkR = sqrt(smk);
      
#         eF  = 2.0*KK*mod(psil,2.0*pi)/pi;
#         lambda = mod(Elliptic.Jacobi.am(eF,m),2.0*pi);
#         sn = sin(lambda);
#         cn = cos(lambda);
#         dn = sqrt(1.0-m*sn+2.0);
        
#         eP  = ellP(-k,lambda,m) 
#         eE  = ellE(lambda,m);
#         EE  = ellE(pi/2.0,m)
        
#         g  = mod(psig + mod(k1r/skmR*(eF/KK*PP-eP),2.0*pi),2*pi);
#         cg = cos(g);
#         sg = sin(g);

#         cdelta = Jh/Jg;
#         sdelta = sqrt(1-cdelta^2);
        
#         cpsih = cos(psih);
#         spsih = sin(psih);

#         Ri2b,Rpsilpsig = i2brotmat_sadov(k,k1r,skmR,smkR,cn,sn,dn,cg,sg,cpsih,spsih,cdelta,sdelta);

#         rnorm = sma*(1-P1^2.0-P2^2.0)/(1+P1*sin(TL)+P2*cos(TL));
       
#         rUVI = [q11*cos(TL)+q21*sin(TL), q12*cos(TL)+q22*sin(TL), q13*cos(TL)+q23*sin(TL)];
#         rUV  = Ri2b*rUVI;
    

#         mu = planetsunparams["muPlanet"]

#         gtorque = [0.0;0.0;0.0];
#         if includeGravityTorque
#             gtorque[1] = 3.0*mu/(rnorm^3)*(C-B)*rUV[2]*rUV[3];
#             gtorque[2] = 3.0*mu/(rnorm^3)*(A-C)*rUV[3]*rUV[1];
#             gtorque[3] = 3.0*mu/(rnorm^3)*(B-A)*rUV[1]*rUV[2];
#         end

#         mtorque = [0.0;0.0;0.0];
#         if includeMagneticTorque
#             muM = planetsunparams["muM"];
#             Zeq = Ri2b*[inertialrefframeinfo["cpa1"];inertialrefframeinfo["cpa2"];inertialrefframeinfo["cpa3"]]
#             Hb  = (Zeq - 3.0*(Zeq[1]*rUV[1]+Zeq[2]*rUV[2]+Zeq[3]*rUV[3])*rUV)*muM/rnorm^3;
#             mtorque = LinearAlgebra.cross(IM,Hb);
#         end

#         srptorque = [0.0;0.0;0.0];
#         if includeSrpTorque 
#             EL = truelong2ecclong(TL,P1,P2);
#             # ETilde = EL - (EL1+EL2)/2.0
#             # Esh    = (EL2-EL1)/2.0

#             rPlanet = planetsunparams["rPlanet"];
#             sc = (rSunV[1]*rV[1]+rSunV[2]*rV[2]+rSunV[3]*rV[3])/rsunnorm + sqrt(rnorm^2-rPlanet^2);
#             shadowfun = 1/2*(1+tanh(1000*sc));
#             # shadowfun = 1.0-Esh/pi
#             # for kk=1:20
#             #     shadowfun = shadowfun - 2/kk/pi*sin(kk*Esh)*cos(kk*ETilde);
#             # end

#             uVect = Ri2b*(rSunV-rsnorm*rUVI);
#             ssd = norm(uVect); 
#             uVect = uVect/ssd;

#             for kk = 1:numberOfFacets 
#                 surfi = facets[kk][2];
#                 niv = facets[kk][4];
#                 rhoiv = facets[kk][3];
#                 niudot = niv[1]*uVect[1]+niv[2]*uVect[2]+niv[3]*uVect[3];
#                 gifun = max(0,niudot);
                
#                 cai = facets[kk][1][1];
#                 cdi = facets[kk][1][2];
#                 csi = facets[kk][1][3];
        
#                 srpforcekk_T1 = - (surfi*gifun*cai*uVect)*shadowfun;
#                 srpforcekk_T2 = - (surfi*gifun*cdi*niv)*shadowfun;
#                 srpforcekk_T3 = - (surfi*gifun*csi*niudot*niv)*shadowfun;

#                 srptorquekk_T1 = LinearAlgebra.cross(rhoiv,srpforcekk_T1);
#                 srptorquekk_T2 = LinearAlgebra.cross(rhoiv,srpforcekk_T2);
#                 srptorquekk_T3 = LinearAlgebra.cross(rhoiv,srpforcekk_T3);
#                 srptorquekk = srptorquekk_T1 + srptorquekk_T2 + srptorquekk_T3; 
#                 srptorque = srptorque + srptorquekk;
#             end        
#         end

#         exttorque = gtorque + mtorque + srptorque;
        
#         zn = eE-eF*EE/KK
#         T4 = (KK*m - PP*(m+k))/k;
#         mxb12 = exttorque[1]*Rpsilpsig[1,2];
#         myb22 = exttorque[2]*Rpsilpsig[2,2];
#         mzb32 = exttorque[3]*Rpsilpsig[3,2];
#         mxb13 = exttorque[1]*Rpsilpsig[1,3];        
#         myb23 = exttorque[2]*Rpsilpsig[2,3];       
#         mzb33 = exttorque[3]*Rpsilpsig[3,3];
#         mxdnsnMcnzn = exttorque[1]*(dn*sn-cn*zn);
#         mydncnPsnzn = exttorque[2]*(cn*dn+zn*sn)
#         mzdnznMmsncn = exttorque[3]*(dn*zn-m*sn*cn)
#         mxb11 = exttorque[1]*Rpsilpsig[1,1];
#         myb21 = exttorque[2]*Rpsilpsig[2,1];
#         mzb31 = exttorque[3]*Rpsilpsig[3,1];

#         torquepert  = zeros(6,1);

#         torquepert[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
    
#         torquepert[2] = mxb13+myb23+mzb33;

#         torquepert[3] = cdelta*torquepert[2] + sdelta*(mxb12+myb22+mzb32);

#         torquepert[4] = -pi/(2*Jg*KK*(1-m))*(mxdnsnMcnzn/smkR+
#         (1-m)*mydncnPsnzn/k1r/smkR+
#         mzdnznMmsncn/skmR);
    
#         torquepert[6] = (mxb11+myb21+mzb31)/Jg/sdelta;
            
#         torquepert[5] = -cdelta*torquepert[6] + 2/pi*sqrt(skm)*T4*k1r*torquepert[4] ;

#     end

#     fun = fun + torquepert;
#     return fun
# end

# function vectorfieldV2prova(psil,skminterp,Jginterp,Jhinterp,psig,psihinterp,TL,k,IV,settings,inertialrefframeinfo,planetsunparams,numberOfFacets,facets,rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,EL1,EL2)
#    skm  = skminterp(TL)
#     Jg   = Jginterp(TL)
#     Jh   = Jhinterp(TL)
#     psih = mod(psihinterp(TL),2.0*pi)
#     #    TL   = mod(TLinterp(psig),2.0*pi)
#    # psig = mod(psiginterp(psil),2.0*pi)

#     A = IV[1]; B = IV[2]; C = IV[3];
#     k1r = sqrt(1+k);
#     smk = 1.0-skm
#     m   = smk/skm*k
#     KK  = ellF(pi/2.0,m);
#     PP  = ellP(-k,pi/2.0,m);

#     fun = zeros(6,1);
#     fun[1] = 0.0;
#     fun[2] = 0.0;
#     fun[3] = 0.0;
#     fun[4] = -pi/k1r*sqrt(skm)*Jg/2/A/C*(C-A)/KK;
#     fun[5] = Jg/A/C*(C-A)*(PP-smk*KK)/KK + Jg/A/C*(C*smk+A*skm);
#     fun[6] = 0.0;
    
#     includeGravityTorque  = get(settings,"includeGravityTorque",false);
#     includeMagneticTorque = get(settings,"includeMagneticTorque",false);
#     includeSrpTorque      = get(settings,"includeSrpTorque",false);
#     includeDragTorque     = get(settings,"includeDragTorque",false);
#     torquepert = zeros(6,1);
    
#     if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque
#         skmR = sqrt(skm);
#         smkR = sqrt(smk);
      
#         eF  = 2.0*KK*mod(psil,2.0*pi)/pi;
#         lambda = mod(Elliptic.Jacobi.am(eF,m),2.0*pi)
#         sn = sin(lambda);
#         cn = cos(lambda);
#         dn = sqrt(1.0-m*sn+2.0);
        
#         eP  = ellP(-k,lambda,m) 
#         eE  = ellE(lambda,m);
#         EE  = ellE(pi/2.0,m)
        
#         g  = mod(psig + mod(k1r/skmR*(eF/KK*PP-eP),2.0*pi),2*pi);
#         cg = cos(g);
#         sg = sin(g);

#         cdelta = Jh/Jg;
#         sdelta = sqrt(1-cdelta^2);
        
#         cpsih = cos(psih);
#         spsih = sin(psih);

#         Ri2b,Rpsilpsig = i2brotmat_sadov(k,k1r,skmR,smkR,cn,sn,dn,cg,sg,cpsih,spsih,cdelta,sdelta);

#         rnorm = sma*(1-P1^2.0-P2^2.0)/(1+P1*sin(TL)+P2*cos(TL));
       
#         rUVI = [q11*cos(TL)+q21*sin(TL), q12*cos(TL)+q22*sin(TL), q13*cos(TL)+q23*sin(TL)];
#         rUV  = Ri2b*rUVI;
    

#         mu = planetsunparams["muPlanet"]

#         gtorque = [0.0;0.0;0.0];
#         if includeGravityTorque
#             gtorque[1] = 3.0*mu/(rnorm^3)*(C-B)*rUV[2]*rUV[3];
#             gtorque[2] = 3.0*mu/(rnorm^3)*(A-C)*rUV[3]*rUV[1];
#             gtorque[3] = 3.0*mu/(rnorm^3)*(B-A)*rUV[1]*rUV[2];
#         end

#         mtorque = [0.0;0.0;0.0];
#         if includeMagneticTorque
#             muM = planetsunparams["muM"];
#             Zeq = Ri2b*[inertialrefframeinfo["cpa1"];inertialrefframeinfo["cpa2"];inertialrefframeinfo["cpa3"]]
#             Hb  = (Zeq - 3.0*(Zeq[1]*rUV[1]+Zeq[2]*rUV[2]+Zeq[3]*rUV[3])*rUV)*muM/rnorm^3;
#             mtorque = LinearAlgebra.cross(IM,Hb);
#         end

#         srptorque = [0.0;0.0;0.0];
#         if includeSrpTorque 
#             EL = truelong2ecclong(TL,P1,P2);
#             # ETilde = EL - (EL1+EL2)/2.0
#             # Esh    = (EL2-EL1)/2.0

#             rPlanet = planetsunparams["rPlanet"];
#             sc = (rSunV[1]*rV[1]+rSunV[2]*rV[2]+rSunV[3]*rV[3])/rsunnorm + sqrt(rnorm^2-rPlanet^2);
#             shadowfun = 1/2*(1+tanh(1000*sc));
#             # shadowfun = 1.0-Esh/pi
#             # for kk=1:20
#             #     shadowfun = shadowfun - 2/kk/pi*sin(kk*Esh)*cos(kk*ETilde);
#             # end

#             uVect = Ri2b*(rSunV-rsnorm*rUVI);
#             ssd = norm(uVect); 
#             uVect = uVect/ssd;

#             for kk = 1:numberOfFacets 
#                 surfi = facets[kk][2];
#                 niv = facets[kk][4];
#                 rhoiv = facets[kk][3];
#                 niudot = niv[1]*uVect[1]+niv[2]*uVect[2]+niv[3]*uVect[3];
#                 gifun = max(0,niudot);
                
#                 cai = facets[kk][1][1];
#                 cdi = facets[kk][1][2];
#                 csi = facets[kk][1][3];
        
#                 srpforcekk_T1 = - (surfi*gifun*cai*uVect)*shadowfun;
#                 srpforcekk_T2 = - (surfi*gifun*cdi*niv)*shadowfun;
#                 srpforcekk_T3 = - (surfi*gifun*csi*niudot*niv)*shadowfun;

#                 srptorquekk_T1 = LinearAlgebra.cross(rhoiv,srpforcekk_T1);
#                 srptorquekk_T2 = LinearAlgebra.cross(rhoiv,srpforcekk_T2);
#                 srptorquekk_T3 = LinearAlgebra.cross(rhoiv,srpforcekk_T3);
#                 srptorquekk = srptorquekk_T1 + srptorquekk_T2 + srptorquekk_T3; 
#                 srptorque = srptorque + srptorquekk;
#             end        
#         end

#         exttorque = gtorque + mtorque + srptorque;
        
#         zn = eE-eF*EE/KK
#         T4 = (KK*m - PP*(m+k))/k;
#         mxb12 = exttorque[1]*Rpsilpsig[1,2];
#         myb22 = exttorque[2]*Rpsilpsig[2,2];
#         mzb32 = exttorque[3]*Rpsilpsig[3,2];
#         mxb13 = exttorque[1]*Rpsilpsig[1,3];        
#         myb23 = exttorque[2]*Rpsilpsig[2,3];       
#         mzb33 = exttorque[3]*Rpsilpsig[3,3];
#         mxdnsnMcnzn = exttorque[1]*(dn*sn-cn*zn);
#         mydncnPsnzn = exttorque[2]*(cn*dn+zn*sn)
#         mzdnznMmsncn = exttorque[3]*(dn*zn-m*sn*cn)
#         mxb11 = exttorque[1]*Rpsilpsig[1,1];
#         myb21 = exttorque[2]*Rpsilpsig[2,1];
#         mzb31 = exttorque[3]*Rpsilpsig[3,1];

#         torquepert  = zeros(6,1);

#         torquepert[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
    
#         torquepert[2] = mxb13+myb23+mzb33;

#         torquepert[3] = cdelta*torquepert[2] + sdelta*(mxb12+myb22+mzb32);

#         torquepert[4] = -pi/(2*Jg*KK*(1-m))*(mxdnsnMcnzn/smkR+
#         (1-m)*mydncnPsnzn/k1r/smkR+
#         mzdnznMmsncn/skmR);
    
#         torquepert[6] = (mxb11+myb21+mzb31)/Jg/sdelta;
            
#         torquepert[5] = -cdelta*torquepert[6] + 2/pi*sqrt(skm)*T4*k1r*torquepert[4] ;

#     end

#     fun = fun + torquepert;
#     return fun
# end

# function vectorfieldV3(t,skminterp,Jginterp,Jhinterp,psilinterp,psiginterp,psihinterp,TLinterp,k,IV,settings,inertialrefframeinfo,planetsunparams,numberOfFacets,facets,rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,EL1,EL2)
#     # skminterp,Jginterp,Jhinterp,psilinterp,psihinterp,TLinterp
#    skm  = skminterp(t)
#    Jg   = Jginterp(t)
#    Jh   = Jhinterp(t)
#    psil = mod(psilinterp(t),2.0*pi)
#    psig = mod(psiginterp(t),2.0*pi)
#    psih = mod(psihinterp(t),2.0*pi)
#    TL   = mod(TLinterp(t),2.0*pi)

#     A = IV[1]; B = IV[2]; C = IV[3];
#     k1r = sqrt(1+k);
#     smk = 1.0-skm
#     m   = smk/skm*k
#     KK  = ellF(pi/2.0,m);
#     PP  = ellP(-k,pi/2.0,m);

#     fun = zeros(6,1);
#     fun[1] = 0.0;
#     fun[2] = 0.0;
#     fun[3] = 0.0;
#     fun[4] = -pi/k1r*sqrt(skm)*Jg/2/A/C*(C-A)/KK;
#     fun[5] = Jg/A/C*(C-A)*(PP-smk*KK)/KK + Jg/A/C*(C*smk+A*skm);
#     fun[6] = 0.0;
    
#     includeGravityTorque  = get(settings,"includeGravityTorque",false);
#     includeMagneticTorque = get(settings,"includeMagneticTorque",false);
#     includeSrpTorque      = get(settings,"includeSrpTorque",false);
#     includeDragTorque     = get(settings,"includeDragTorque",false);
#     torquepert = zeros(6,1);
    
#     if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque
#         skmR = sqrt(skm);
#         smkR = sqrt(smk);
      
#         eF  = 2.0*KK*mod(psil,2.0*pi)/pi;
#         lambda = mod(Elliptic.Jacobi.am(eF,m),2.0*pi)
#         sn = sin(lambda)
#         cn = cos(lambda)
#         dn = sqrt(1.0-m*sn^2.0)
        
#         eP  = ellP(-k,lambda,m) 
#         eE  = ellE(lambda,m);
#         EE  = ellE(pi/2.0,m)
        
#         g  = mod(psig + mod(k1r/skmR*(eF/KK*PP-eP),2.0*pi),2*pi);
#         cg = cos(g);
#         sg = sin(g);

#         cdelta = Jh/Jg;
#         sdelta = sqrt(1-cdelta^2);
        
#         cpsih = cos(psih);
#         spsih = sin(psih);

#         Ri2b,Rpsilpsig = i2brotmat_sadov(k,k1r,skmR,smkR,cn,sn,dn,cg,sg,cpsih,spsih,cdelta,sdelta);
#         rnorm = sma*(1-P1^2.0-P2^2.0)/(1+P1*sin(TL)+P2*cos(TL));
#         rUVI = [q11*cos(TL)+q21*sin(TL), q12*cos(TL)+q22*sin(TL), q13*cos(TL)+q23*sin(TL)];
#         rUV  = Ri2b*rUVI;

#         mu = planetsunparams["muPlanet"]

#         gtorque = [0.0;0.0;0.0];
#         if includeGravityTorque
#             gtorque[1] = 3.0*mu/(rnorm^3)*(C-B)*rUV[2]*rUV[3];
#             gtorque[2] = 3.0*mu/(rnorm^3)*(A-C)*rUV[3]*rUV[1];
#             gtorque[3] = 3.0*mu/(rnorm^3)*(B-A)*rUV[1]*rUV[2];
#         end

#         mtorque = [0.0;0.0;0.0];
#         if includeMagneticTorque
#             muM = planetsunparams["muM"];
#             Zeq = Ri2b*[inertialrefframeinfo["cpa1"];inertialrefframeinfo["cpa2"];inertialrefframeinfo["cpa3"]]
#             Hb  = (Zeq - 3.0*(Zeq[1]*rUV[1]+Zeq[2]*rUV[2]+Zeq[3]*rUV[3])*rUV)*muM/rnorm^3;
#             mtorque = LinearAlgebra.cross(IM,Hb);
#         end

#         srptorque = [0.0;0.0;0.0];
#         if includeSrpTorque 
#             EL = truelong2ecclong(TL,P1,P2);
#             # ETilde = EL - (EL1+EL2)/2.0
#             # Esh    = (EL2-EL1)/2.0

#             rPlanet = planetsunparams["rPlanet"];
#             sc = (rSunV[1]*rV[1]+rSunV[2]*rV[2]+rSunV[3]*rV[3])/rsunnorm + sqrt(rnorm^2-rPlanet^2);
#             shadowfun = 1/2*(1+tanh(1000*sc));
#             # shadowfun = 1.0-Esh/pi
#             # for kk=1:20
#             #     shadowfun = shadowfun - 2/kk/pi*sin(kk*Esh)*cos(kk*ETilde);
#             # end

#             uVect = Ri2b*(rSunV-rsnorm*rUVI);
#             ssd = norm(uVect); 
#             uVect = uVect/ssd;

#             for kk = 1:numberOfFacets 
#                 surfi = facets[kk][2];
#                 niv = facets[kk][4];
#                 rhoiv = facets[kk][3];
#                 niudot = niv[1]*uVect[1]+niv[2]*uVect[2]+niv[3]*uVect[3];
#                 gifun = max(0,niudot);
                
#                 cai = facets[kk][1][1];
#                 cdi = facets[kk][1][2];
#                 csi = facets[kk][1][3];
        
#                 srpforcekk_T1 = - (surfi*gifun*cai*uVect)*shadowfun;
#                 srpforcekk_T2 = - (surfi*gifun*cdi*niv)*shadowfun;
#                 srpforcekk_T3 = - (surfi*gifun*csi*niudot*niv)*shadowfun;

#                 srptorquekk_T1 = LinearAlgebra.cross(rhoiv,srpforcekk_T1);
#                 srptorquekk_T2 = LinearAlgebra.cross(rhoiv,srpforcekk_T2);
#                 srptorquekk_T3 = LinearAlgebra.cross(rhoiv,srpforcekk_T3);
#                 srptorquekk = srptorquekk_T1 + srptorquekk_T2 + srptorquekk_T3; 
#                 srptorque = srptorque + srptorquekk;
#             end        
#         end

#         exttorque = gtorque + mtorque + srptorque;
        
#         zn = eE-eF*EE/KK
#         T4 = (KK*m - PP*(m+k))/k;
#         mxb11 = exttorque[1]*Rpsilpsig[1,1]; 
#         myb21 = exttorque[2]*Rpsilpsig[2,1];
#         mzb31 = exttorque[3]*Rpsilpsig[3,1];
#         mxb12 = exttorque[1]*Rpsilpsig[1,2];
#         myb22 = exttorque[2]*Rpsilpsig[2,2];
#         mzb32 = exttorque[3]*Rpsilpsig[3,2];
#         mxb13 = exttorque[1]*Rpsilpsig[1,3];        
#         myb23 = exttorque[2]*Rpsilpsig[2,3];       
#         mzb33 = exttorque[3]*Rpsilpsig[3,3];
#         mxdnsnMcnzn = exttorque[1]*(dn*sn-cn*zn);
#         mydncnPsnzn = exttorque[2]*(cn*dn+zn*sn)
#         mzdnznMmsncn = exttorque[3]*(dn*zn-m*sn*cn)


#         torquepert  = zeros(6,1);

#         torquepert[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
    
#         torquepert[2] = mxb13+myb23+mzb33;

#         torquepert[3] = cdelta*torquepert[2] + sdelta*(mxb12+myb22+mzb32);

#         torquepert[4] = -pi/(2*Jg*KK*(1-m))*(mxdnsnMcnzn/smkR+
#         (1-m)*mydncnPsnzn/k1r/smkR+
#         mzdnznMmsncn/skmR);
    
#         torquepert[6] = (mxb11+myb21+mzb31)/Jg/sdelta;
            
#         torquepert[5] = -cdelta*torquepert[6] + 2/pi*sqrt(skm)*T4*k1r*torquepert[4] ;

#     end

#     fun = fun + torquepert;
#     return fun
# end

# function vectorfieldV4(TL,skminterp,Jginterp,Jhinterp,psilinterp,psiginterp,psihinterp,k,IV,settings,inertialrefframeinfo,planetsunparams,numberOfFacets,facets,rSunV,sma,P1,P2,q11,q12,q13,q21,q22,q23,EL1,EL2)
#     # skminterp,Jginterp,Jhinterp,psilinterp,psihinterp,TLinterp
#    skm  = skminterp(TL)
#    Jg   = Jginterp(TL)
#    Jh   = Jhinterp(TL)
#    psil = mod(psilinterp(TL),2.0*pi)
#    psig = mod(psiginterp(TL),2.0*pi)
#    psih = mod(psihinterp(TL),2.0*pi)
 
#     A = IV[1]; B = IV[2]; C = IV[3];
#     k1r = sqrt(1+k);
#     smk = 1.0-skm
#     m   = smk/skm*k
#     KK  = ellF(pi/2.0,m);
#     PP  = ellP(-k,pi/2.0,m);

#     fun = zeros(6,1);
#     fun[1] = 0.0;
#     fun[2] = 0.0;
#     fun[3] = 0.0;
#     fun[4] = -pi/k1r*sqrt(skm)*Jg/2/A/C*(C-A)/KK;
#     fun[5] = Jg/A/C*(C-A)*(PP-smk*KK)/KK + Jg/A/C*(C*smk+A*skm);
#     fun[6] = 0.0;
    
#     includeGravityTorque  = get(settings,"includeGravityTorque",false);
#     includeMagneticTorque = get(settings,"includeMagneticTorque",false);
#     includeSrpTorque      = get(settings,"includeSrpTorque",false);
#     includeDragTorque     = get(settings,"includeDragTorque",false);
#     torquepert = zeros(6,1);
    
#     if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque
#         skmR = sqrt(skm);
#         smkR = sqrt(smk);
      
#         eF  = 2.0*KK*psil/pi;
#         sn = Elliptic.Jacobi.sn(eF,m);
#         cn = Elliptic.Jacobi.cn(eF,m);
#         dn = Elliptic.Jacobi.dn(eF,m);

#         lambda = mod(atan(sn,cn),2.0*pi);
#         eP  = ellP(-k,lambda,m) 
#         eE  = ellE(lambda,m);
#         EE  = ellE(pi/2.0,m)
        
#         g  = mod(psig + mod(k1r/skmR*(eF/KK*PP-eP),2.0*pi),2*pi);
#         cg = cos(g);
#         sg = sin(g);

#         cdelta = Jh/Jg;
#         sdelta = sqrt(1-cdelta^2);
        
#         cpsih = cos(psih);
#         spsih = sin(psih);

#         Ri2b,Rpsilpsig = i2brotmat_sadov(k,k1r,skmR,smkR,cn,sn,dn,cg,sg,cpsih,spsih,cdelta,sdelta);

#         rnorm = sma*(1-P1^2.0-P2^2.0)/(1+P1*sin(TL)+P2*cos(TL));
#         rUVI = [q11*cos(TL)+q21*sin(TL), q12*cos(TL)+q22*sin(TL), q13*cos(TL)+q23*sin(TL)];
#         rUV  = Ri2b*rUVI;

#         mu = planetsunparams["muPlanet"]

#         gtorque = [0.0;0.0;0.0];
#         if includeGravityTorque
#             gtorque[1] = 3.0*mu/(rnorm^3)*(C-B)*rUV[2]*rUV[3];
#             gtorque[2] = 3.0*mu/(rnorm^3)*(A-C)*rUV[3]*rUV[1];
#             gtorque[3] = 3.0*mu/(rnorm^3)*(B-A)*rUV[1]*rUV[2];
#         end

#         mtorque = [0.0;0.0;0.0];
#         if includeMagneticTorque
#             muM = planetsunparams["muM"];
#             Zeq = Ri2b*[inertialrefframeinfo["cpa1"];inertialrefframeinfo["cpa2"];inertialrefframeinfo["cpa3"]]
#             Hb  = (Zeq - 3.0*(Zeq[1]*rUV[1]+Zeq[2]*rUV[2]+Zeq[3]*rUV[3])*rUV)*muM/rnorm^3;
#             mtorque = LinearAlgebra.cross(IM,Hb);
#         end

#         srptorque = [0.0;0.0;0.0];
#         if includeSrpTorque 
#             EL = truelong2ecclong(TL,P1,P2);
#             # ETilde = EL - (EL1+EL2)/2.0
#             # Esh    = (EL2-EL1)/2.0

#             rPlanet = planetsunparams["rPlanet"];
#             sc = (rSunV[1]*rV[1]+rSunV[2]*rV[2]+rSunV[3]*rV[3])/rsunnorm + sqrt(rnorm^2-rPlanet^2);
#             shadowfun = 1/2*(1+tanh(1000*sc));
#             # shadowfun = 1.0-Esh/pi
#             # for kk=1:20
#             #     shadowfun = shadowfun - 2/kk/pi*sin(kk*Esh)*cos(kk*ETilde);
#             # end

#             uVect = Ri2b*(rSunV-rsnorm*rUVI);
#             ssd = norm(uVect); 
#             uVect = uVect/ssd;

#             for kk = 1:numberOfFacets 
#                 surfi = facets[kk][2];
#                 niv = facets[kk][4];
#                 rhoiv = facets[kk][3];
#                 niudot = niv[1]*uVect[1]+niv[2]*uVect[2]+niv[3]*uVect[3];
#                 gifun = max(0,niudot);
                
#                 cai = facets[kk][1][1];
#                 cdi = facets[kk][1][2];
#                 csi = facets[kk][1][3];
        
#                 srpforcekk_T1 = - (surfi*gifun*cai*uVect)*shadowfun;
#                 srpforcekk_T2 = - (surfi*gifun*cdi*niv)*shadowfun;
#                 srpforcekk_T3 = - (surfi*gifun*csi*niudot*niv)*shadowfun;

#                 srptorquekk_T1 = LinearAlgebra.cross(rhoiv,srpforcekk_T1);
#                 srptorquekk_T2 = LinearAlgebra.cross(rhoiv,srpforcekk_T2);
#                 srptorquekk_T3 = LinearAlgebra.cross(rhoiv,srpforcekk_T3);
#                 srptorquekk = srptorquekk_T1 + srptorquekk_T2 + srptorquekk_T3; 
#                 srptorque = srptorque + srptorquekk;
#             end        
#         end

#         exttorque = gtorque + mtorque + srptorque;
        
#         zn = eE-eF*EE/KK
#         T4 = (KK*m - PP*(m+k))/k;
#         mxb12 = exttorque[1]*Rpsilpsig[1,2];
#         myb22 = exttorque[2]*Rpsilpsig[2,2];
#         mzb32 = exttorque[3]*Rpsilpsig[3,2];
#         mxb13 = exttorque[1]*Rpsilpsig[1,3];        
#         myb23 = exttorque[2]*Rpsilpsig[2,3];       
#         mzb33 = exttorque[3]*Rpsilpsig[3,3];
#         mxdnsnMcnzn = exttorque[1]*(dn*sn-cn*zn);
#         mydncnPsnzn = exttorque[2]*(cn*dn+zn*sn)
#         mzdnznMmsncn = exttorque[3]*(dn*zn-m*sn*cn)
#         mxb11 = exttorque[1]*Rpsilpsig[1,1];
#         myb21 = exttorque[2]*Rpsilpsig[2,1];
#         mzb31 = exttorque[3]*Rpsilpsig[3,1];

#         torquepert  = zeros(6,1);

#         torquepert[1] = -2*mxb13/Jg*skm - 2*myb23*(1-m)/Jg/(1+k)*skm + 2*smk*mzb33/Jg;
    
#         torquepert[2] = mxb13+myb23+mzb33;

#         torquepert[3] = cdelta*torquepert[2] + sdelta*(mxb12+myb22+mzb32);

#         torquepert[4] = -pi/(2*Jg*KK*(1-m))*(mxdnsnMcnzn/smkR+
#         (1-m)*mydncnPsnzn/k1r/smkR+
#         mzdnznMmsncn/skmR);
    
#         torquepert[6] = (mxb11+myb21+mzb31)/Jg/sdelta;
            
#         torquepert[5] = -cdelta*torquepert[6] + 2/pi*sqrt(skm)*T4*k1r*torquepert[4] ;

#     end

#     fun = fun + torquepert;
#     return fun
# end