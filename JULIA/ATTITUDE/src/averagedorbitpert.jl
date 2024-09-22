############## zonal harmonics
"""
    zonalharmonicsaveragedgauss(sma,eta,P1,P2,Q1,Q2,mu,rPlanet,J2,J3,J4,J5)
    Returns the terms of the Gaussian Planetary equations due to zonal harmonics
    averaged with respect to the orbital mean anomaly

    INPUT  
    sma = semi-major axis 
    eta = sqrt(1-ec^2), ec = eccentricity
    P1  = ec*sin(om+OM), om = argument of the pericenter, OM = right ascension of the ascending node
    P2  = ec*cos(om+OM) 
    Q1  = tan(incl/2)sin(OM), incl = inclination
    Q2  = tan(incl/2)cos(OM)
    mu  = gravitational parameter of the central body
    rPlanet = mean radiu of the central body
    Ji, i = 2,3,4,5 = zonal harmonics of the planetRotation

    OUTPUT
    out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt], out = out(sma,P1,P2,Q1,Q2; mu, rPlanet, J2,J3,J4,J5)
    with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

"""
function zonalharmonicsaveragedgauss(sma,eta,P1,P2,Q1,Q2,mu,rPlanet,J2,J3,J4,J5)
    out = zeros(6);
    GG = 1+Q1^2+Q2^2.0
    smaOmuS = sqrt(sma/mu);

    if J2!=0.0
        sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN = getj2averagedterms(eta,P1,P2,Q1,Q2)
        kj2 = (3*mu*rPlanet^2*J2)/(2*sma^4*eta^8*GG^2)
        contr = zeros(6);
        contr[1] =  2/eta*sqrt(sma^3.0/mu)*kj2*(-P1*cTLfR + P2*sTLfR + phiLP1fT);
        contr[2] =  eta*smaOmuS*kj2*(-cTLfR + P1*phiLM1fT + phiLM1sTLfT + sTLfT - P2*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
        contr[3] =  eta*smaOmuS*kj2*( sTLfR + P2*phiLM1fT + phiLM1cTLfT + cTLfT + P1*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
        contr[4] =  0.5*eta*smaOmuS*kj2*GG*phiLM1sTLfN;
        contr[5] =  0.5*eta*smaOmuS*kj2*GG*phiLM1cTLfN;
        contr[6] =  -eta*smaOmuS*kj2*( (P1*sTLfR+P2*cTLfR +P1*phiLM1cTLfT -P2*phiLM1sTLfT + P1*cTLfT - P2*sTLfT)/(1+eta) + 2*eta*phiLM1fR + Q1*phiLM1cTLfN - Q2*phiLM1sTLfN); 

        for kk=1:6
            if abs(contr[kk])<1e-15
                contr[kk] = 0.0
            end
            out[kk] = out[kk] + contr[kk]
        end

    end

    if J3!=0.0
        sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN = getj3averagedterms(eta,P1,P2,GG,Q1,Q2)
        kj3 = (3*mu*J3*rPlanet^3)/(2*eta^10*sma^5*GG)
        contr = zeros(6);
        contr[1] =  2/eta*sqrt(sma^3.0/mu)*kj3*(-P1*cTLfR + P2*sTLfR + phiLP1fT);
        contr[2] =  eta*smaOmuS*kj3*(-cTLfR + P1*phiLM1fT + phiLM1sTLfT + sTLfT - P2*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
        contr[3] =  eta*smaOmuS*kj3*( sTLfR + P2*phiLM1fT + phiLM1cTLfT + cTLfT + P1*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
        contr[4] =  0.5*eta*smaOmuS*kj3*GG*phiLM1sTLfN;
        contr[5] =  0.5*eta*smaOmuS*kj3*GG*phiLM1cTLfN;
        contr[6] =  -eta*smaOmuS*kj3*( (P1*sTLfR+P2*cTLfR +P1*phiLM1cTLfT -P2*phiLM1sTLfT + P1*cTLfT - P2*sTLfT)/(1+eta) + 2*eta*phiLM1fR + Q1*phiLM1cTLfN - Q2*phiLM1sTLfN); 
        for kk=1:6
            if abs(contr[kk])<1e-15
                contr[kk] = 0.0
            end
            out[kk] = out[kk] + contr[kk]
        end
    end

    if J4!=0.0
        sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN = getj4averagedterms(eta,P1,P2,GG,Q1,Q2)
        kj4 = mu*J4*rPlanet^4/(eta^12*sma^6)
        contr = zeros(6);
        contr[1] =  2/eta*sqrt(sma^3.0/mu)*kj4*(-P1*cTLfR + P2*sTLfR + phiLP1fT);
        contr[2] =  eta*smaOmuS*kj4*(-cTLfR + P1*phiLM1fT + phiLM1sTLfT + sTLfT - P2*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
        contr[3] =  eta*smaOmuS*kj4*( sTLfR + P2*phiLM1fT + phiLM1cTLfT + cTLfT + P1*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
        contr[4] =  0.5*eta*smaOmuS*kj4*GG*phiLM1sTLfN;
        contr[5] =  0.5*eta*smaOmuS*kj4*GG*phiLM1cTLfN;
        contr[6] =  -eta*smaOmuS*kj4*( (P1*sTLfR+P2*cTLfR +P1*phiLM1cTLfT -P2*phiLM1sTLfT + P1*cTLfT - P2*sTLfT)/(1+eta) + 2*eta*phiLM1fR + Q1*phiLM1cTLfN - Q2*phiLM1sTLfN); 
        for kk=1:6
            if abs(contr[kk])<1e-15
                contr[kk] = 0.0
            end
            out[kk] = out[kk] + contr[kk]
        end
    end

    if J5!=0.0
        sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN = getj4averagedterms(eta,P1,P2,GG,Q1,Q2)
        kj5 = mu*J5*rPlanet^5/(2*eta^14*sma^7)
        contr = zeros(6);
        contr[1] =  2/eta*sqrt(sma^3.0/mu)*kj5*(-P1*cTLfR + P2*sTLfR + phiLP1fT);
        contr[2] =  eta*smaOmuS*kj5*(-cTLfR + P1*phiLM1fT + phiLM1sTLfT + sTLfT - P2*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
        contr[3] =  eta*smaOmuS*kj5*( sTLfR + P2*phiLM1fT + phiLM1cTLfT + cTLfT + P1*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
        contr[4] =  0.5*eta*smaOmuS*kj5*GG*phiLM1sTLfN;
        contr[5] =  0.5*eta*smaOmuS*kj5*GG*phiLM1cTLfN;
        contr[6] =  -eta*smaOmuS*kj5*( (P1*sTLfR+P2*cTLfR +P1*phiLM1cTLfT -P2*phiLM1sTLfT + P1*cTLfT - P2*sTLfT)/(1+eta) + 2*eta*phiLM1fR + Q1*phiLM1cTLfN - Q2*phiLM1sTLfN); 
        for kk=1:6
            if abs(contr[kk])<1e-15
                contr[kk] = 0.0
            end
            out[kk] = out[kk] + contr[kk]
        end
    end

    return out;
end

function getj2averagedterms(eta,P1,P2,Q1,Q2)
    
    # recurring terms
    t1 = Q1 ^ 2
    t2 = t1 ^ 2
    t3 = Q2 ^ 2
    t7 = t3 ^ 2
    t11 = P2 * Q1
    t12 = Q2 * t11
    t15 = eta ^ 2
    t16 = t15 * eta
    t19 = t1 * P2
    t23 = P1 * Q1
    t24 = Q2 * t23
    t27 = t3 * P2  
    t55 = t1 + t3 - 1

    sTLfR = -t16 * (P1 * (t2 + t1 * (2 * t3 - 1) + t7 - 7 * t3 + 1) + 6 * t12)
    cTLfR = -(t2 * P2 + t7 * P2 + 2 * t3 * t19 + P2 - 7 * t19 + 6 * t24 - t27) * t16
    phiLM1fR = -(2 * t3 * t1 - 4 * t1 + t2 - 4 * t3 + t7 + 1) * t16
    phiLP1fT = -6 * (P1 * Q2 - t11) * (P2 * Q2 + t23) * t16
    phiLM1fT = 0
    phiLM1sTLfT = P2 * (t1 * t16 - t3 * t16) - 2 * Q1 * Q2 * P1 * t16
    phiLM1cTLfT = (P1 * (t1 - t3) + 2 * t12) * t16
    sTLfT = -4 * (t24 - t19 / 2 + t27 / 2) * t16
    cTLfT = 2 * phiLM1cTLfT
    phiLM1sTLfN = 2 * t55 * Q2 * t16
    phiLM1cTLfN = -2 * t55 * Q1 * t16
   
    
    return  sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN;
end

function getj3averagedterms(eta,P1,P2,GG,Q1,Q2)

    t1 = P1 ^ 2
    t3 = P2 ^ 2
    t6 = Q2 ^ 2
    t7 = t6 * Q2
    t9 = P1 * P2
    t10 = t6 * Q1
    t11 = t10 * t9
    t13 = GG ^ 2
    t14 = Q1 ^ 2
    t15 = 10/3 * t14
    t16 = t13 - t15
    t25 = P1 * Q1
    t30 = eta ^ 2
    t31 = t30 * eta
    t33 = 0.1e1 / t13
    t37 = Q1 * t1 * t13
    t38 = 3 * t37
    t40 = P2 * Q2
    t41 = t40 * P1 * t13
    t42 = 6 * t41
    t44 = Q1 * t3 * t13
    t46 = t14 * Q1
    t47 = t46 * t1
    t48 = 10 * t47
    t50 = t6 * Q1 * t1
    t51 = 30 * t50
    t53 = Q2 * t14 * t9
    t54 = 60 * t53
    t55 = t7 * t9
    t56 = 20 * t55
    t57 = t46 * t3
    t60 = t6 * Q1 * t3
    t61 = 30 * t60
    t63 = 4 * Q1 * t13
    t64 = 20 * t46
    t65 = 20 * t10
    t66 = t38 - t42 + 9 * t44 - t48 - t51 + t54 + t56 - 50 * t57 - t61 + t63 - t64 - t65
    t68 = t25 + t40
    t76 = 10 * t14
    t80 = 4/3 * t13
    t81 = 20/3 * t14
    t86 = 5 * t14
    t87 = 5 * t6
    t89 = (t13 - t86 - t87) * t31
    t96 = t38 + 2 * t41 + t44 - t48 - t51 + 20 * t53 - t56 - 10 * t57 + 10 * t60 + t63 - t64 - t65
    t99 = 10 * t1
    t100 = 10 * t3
    t103 = 20 * t11
    t105 = t1 * (t13 + t76)
    t109 = t3 * (3 * t13 - 30 * t14)
    t117 = 2 * (t13 - t76) * Q1 * t9
    t127 = 9 * t37 + t42 + 3 * t44 - 30 * t47 - 90 * t50 + t54 - 60 * t55 - 30 * t57 + t61 + t63 - t64 - t65
    t140 = P2 * Q1
    t144 = t14 + t6 - 1
    t146 = t33 * t31

    sTLfR = -9 * t33 * t31 * (t7 * (-50/9 * t1 - 10/9 * t3 - 20/9) + 20/3 * t11 + Q2 * (t1 * t16 + t3 * (t13 / 3 - t15) + 4/9 * t13 - 20/9 * t14) - 2/3 * P2 * t16 * t25)
    phiLM1fR = -8 * t33 * (P1 * Q2 - t140) * t89
    cTLfR = t33 * t66 * t31
    phiLP1fT = 3 * t33 * (t1 * (t13 - t15 - 10 * t6) + 40/3 * Q1 * Q2 * t9 + t3 * (t13 - t76 - 10/3 * t6) + t80 - t81 - 20/3 * t6) * t31 * t68
    phiLM1fT = 2 * t33 * t68 * t89
    phiLM1sTLfT = t33 * t96 * t31 / 4
    phiLM1cTLfT = t33 * t31 * (t7 * (-t99 - t100 - 20) + t103 + Q2 * (t105 + t109 + 4 * t13 - 20 * t14) + t117) / 4
    sTLfT = t33 * t127 * t31 / 4
    cTLfT = 3/4 * t33 * t31 * (t7 * (-t99 - t100 - 20/3) + t103 + Q2 * (t105 + t109 + t80 - t81) + t117)
    phiLM1sTLfN = -t146 * t144 * (P1 * (t13 - t86 - 15 * t6) + 10 * Q2 * t140)
    phiLM1cTLfN = -t146 * (P2 * (t13 - 15 * t14 - t87) + 10 * Q2 * t25) * t144

    return  sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN;

end

function getj4averagedterms(eta,P1,P2,GG,Q1,Q2)

    t1 = eta ^ 2
    t2 = t1 * eta
    t3 = Q1 ^ 2
    t4 = t3 ^ 2
    t6 = GG ^ 2
    t8 = Q2 ^ 2
    t9 = 350/3 * t8
    t12 = t8 * t6
    t14 = t6 ^ 2
    t15 = t8 ^ 2
    t18 = P1 ^ 2
    t19 = t18 * P1
    t23 = t6 - 7/2 * t3 - 35/6 * t8
    t25 = P2 * Q2
    t29 = P2 ^ 2
    t30 = 175/3 * t29
    t44 = (t29 + 2) * t6
    t45 = t8 * t44
    t51 = 35/6 * t29
    t54 = 7/2 * t29
    t63 = 0.1e1 / t14
    t72 = t29 * P2
    t74 = Q1 * Q2
    t77 = t6 - 35/6 * t3 - 7/2 * t8
    t82 = 175/3 * t18
    t89 = (t18 + 2) * t6
    t96 = (t18 + 2/3) * t6
    t103 = 7/2 * t18
    t106 = 35/6 * t18
    t108 = t8 * (-t106 - 28/3)
    t111 = P1 * Q2
    t117 = P2 * Q1
    t119 = P1 * Q1
    t121 = (t119 + t25) * (t111 - t117)
    t123 = P1 * P2
    t127 = 2 * t6
    t134 = 14/3 * t3
    t135 = 14/3 * t8
    t148 = t72 * (35/24 * t4 + t3 * (-t6 / 4 - 7/4 * t8) + t12 / 4 - 7/8 * t15)
    t149 = t3 * Q1
    t152 = t8 * Q2
    t156 = t29 * (-7/2 * Q2 * t149 * P1 + 7/2 * t152 * t119)
    t157 = 21/8 * t18
    t161 = 21/4 * t8 * t18
    t188 = t19 * (-7/2 * t4 + t3 * (t6 - 7 * t8) - t12 + 35/6 * t15)
    t194 = t18 * (14 * Q2 * t149 * P2 - 14 * t152 * t117)
    t195 = 35/2 * t29
    t199 = 21 * t8 * t29
    t203 = 21/2 * t29
    t249 = t29 * t6
    t250 = 3 * t249
    t274 = (t3 + t8 - 1) * t2
    t282 = 14/3 * t18
    t325 = t18 * t6
    t367 = -560 * t152 * Q1 * t123 - 560 * Q2 * t149 * t123 + 420 * t8 * t3 * t29 + 70 * t15 * t29 + 350 * t4 * t29 - 40 * t3 * t6 + 280 * t8 * t3 - 40 * t12 + 2 * t14 + 140 * t15 + 140 * t4

    sTLfR = 45/16 * t63 * (t19 * (35/3 * t4 + t3 * (-20/3 * t6 + t9) - 100/3 * t12 + t14 + 1225/9 * t15) + 40 * t18 * t25 * Q1 * t23 + P1 * (t4 * (t30 + 280/9) + t3 * (t8 * (210 * t29 + 560/3) - 20 * (t29 + 2/3) * t6) + t15 * (t30 + 1400/9) - 20 * t45 + (t29 + 4/3) * t14) + 40/3 * t25 * (t3 * (-t51 - 28/3) + t8 * (-t54 - 28/3) + t44) * Q1) * t2
    cTLfR = 45/16 * t63 * (t72 * (1225/9 * t4 + t3 * (-100/3 * t6 + t9) - 20/3 * t12 + t14 + 35/3 * t15) + 40 * t29 * t77 * P1 * t74 + P2 * (t4 * (1400/9 + t82) + t3 * (t8 * (560/3 + 210 * t18) - 20 * t89) + t15 * (t82 + 280/9) - 20 * t8 * t96 + (t18 + 4/3) * t14) + 40/3 * t111 * (t3 * (-t103 - 28/3) + t108 + t89) * Q1) * t2
    phiLM1fR = 15/16 * t63 * (120 * Q2 * t117 * P1 * t6 + 420 * t8 * t3 * t18 + 3 * t18 * t14 + 3 * t29 * t14 + 350 * t15 * t18 + 70 * t4 * t18 - 90 * t3 * t249 - 30 * t8 * t249 - 30 * t3 * t325 - 90 * t8 * t325 + t367) * t2
    phiLP1fT = 75/2 * t63 * (t18 * t23 + 14/3 * t74 * t123 + t29 * t77 + t127 - 28/3 * t3 - 28/3 * t8) * t2 * t121
    phiLM1fT = 45/2 * t63 * (t6 - t134 - t135) * t2 * t121
    phiLM1sTLfT = 15/2 * t63 * t2 * (t148 + t156 + P2 * (t4 * (t157 + 7) + t3 * (t161 - 3/4 * t89) + 3/4 * t8 * (t108 + t89)) + t111 * Q1 * (t3 * (-t103 - 14) + t8 * (-t106 - 14) + (t18 + 3) * t6))
    phiLM1cTLfT = -15/8 * t63 * t2 * (t188 + t194 + P1 * (t4 * (-t195 - 28) + t3 * (t199 + 3 * t44) + t15 * (t203 + 28) - 3 * t45) + 4 * (t3 * (-t51 - 14) + t8 * (-t54 - 14) + (t29 + 3) * t6) * Q2 * t117)
    sTLfT = 30 * t63 * t2 * (t148 + t156 + P2 * (t4 * (t157 + 7/3) + t3 * (t161 - 3/4 * t96) + 3/4 * t8 * (t8 * (-t106 - 28/9) + t96)) + t111 * Q1 * (t3 * (-t103 - 14/3) + t8 * (-t106 - 14/3) + (t18 + 1) * t6))
    cTLfT = -15/2 * t63 * t2 * (t188 + t194 + P1 * (t4 * (-t195 - 28/3) + t3 * (t250 + t199 + t127) + t15 * (t203 + 28/3) + t8 * (-t250 - t127)) + 4 * t25 * (t3 * (-t51 - 14/3) + t8 * (-t54 - 14/3) + (t29 + 1) * t6) * Q1)
    phiLM1sTLfN = -135/8 * t63 * (t152 * (-70/9 * t18 - 14/9 * t29 - 28/9) + 28/3 * t8 * Q1 * t123 + Q2 * (t3 * (-t282 - 14/3 * t29 - 28/9) + (t18 + t29 / 3 + 4/9) * t6) - 2/3 * (t6 - t134) * P1 * t117) * t274
    phiLM1cTLfN = 45/8 * t63 * (t149 * (-t282 - 70/3 * t29 - 28/3) + 28 * Q2 * t3 * t123 + Q1 * (t8 * (-14 * t18 - 14 * t29 - 28/3) + (t18 + 3 * t29 + 4/3) * t6) - 2 * t111 * (t6 - t135) * P2) * t274

   
    return  sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN;

end

function getj5averagedterms(eta,P1,P2,GG,Q1,Q2)
 
    t1 = P1 ^ 2
    t2 = t1 ^ 2
    t4 = P2 ^ 2
    t8 = t4 ^ 2
    t12 = Q2 ^ 2
    t13 = t12 ^ 2
    t14 = t13 * Q2
    t19 = P2 * Q1
    t20 = t13 * t19
    t23 = GG ^ 2
    t25 = Q1 ^ 2
    t33 = 112/3 * t23
    t34 = 168 * t25
    t45 = 112/25 * t23
    t48 = t12 * Q2
    t50 = P1 * Q1
    t51 = 18/5 * t25
    t64 = t23 ^ 2
    t65 = t25 ^ 2
    t66 = 63/5 * t65
    t67 = t25 * t23
    t68 = 7 * t67
    t69 = t64 + t66 - t68
    t71 = 6/5 * t64
    t76 = 12/5 * t64
    t78 = 112/5 * t67
    t81 = t64 / 5
    t82 = 147/5 * t65
    t85 = 4/5 * t64
    t86 = 84 * t65
    t90 = 8/25 * t64
    t98 = 2 * t64
    t106 = eta ^ 2
    t107 = t106 * eta
    t110 = 0.1e1 / t64 / GG
    t112 = 63/5 * t2
    t119 = t65 * Q1
    t121 = P1 * Q2
    t128 = 7 * t23
    t129 = 126 * t12
    t130 = -t128 + t129
    t131 = t2 * t130
    t136 = 504 * t12
    t147 = 672/5 * t12
    t148 = 112/5 * t23
    t150 = t25 * Q1
    t152 = 6 * t12
    t158 = 8/3 * t23
    t159 = 12 * t12
    t162 = P2 * Q2
    t163 = t25 * t162
    t166 = t12 * t23
    t167 = 35 * t166
    t168 = 147 * t13
    t169 = t64 - t167 + t168
    t170 = t2 * t169
    t171 = 6 * t64
    t176 = 4 * t64
    t177 = 112 * t166
    t178 = 420 * t13
    t181 = 5 * t64
    t182 = 63 * t13
    t185 = 12 * t64
    t186 = 252 * t13
    t189 = 336/5 * t13
    t190 = 112/5 * t166
    t191 = 8/5 * t64
    t198 = 7 * t166
    t199 = 63/5 * t13
    t212 = t50 + t162
    t218 = t1 * P1
    t228 = 126 * t13
    t239 = 12 * t25
    t245 = 147 * t65
    t246 = 35 * t23
    t252 = 112 * t23
    t259 = 336/5 * t65
    t262 = t2 * (t25 * t130 - t167 + t168 + t64 + t66) + 56 * t218 * t19 * Q2 * (t23 - t51 - t152) + t1 * (t4 * (126 * t65 + t25 * (2268/5 * t12 - 42 * t23) + t98 - 42 * t166 + t228) + t86 + t25 * (-t33 + t136) + t176 - t177 + t178) + 56 * P1 * t19 * (t4 * (t23 - 6 * t25 - 18/5 * t12) + t158 - t239 - t159) * Q2 + t8 * (t245 + t25 * (-t246 + t129) + t64 - t198 + t199) + t4 * (420 * t65 + t25 * (-t252 + t136) + t176 - 112/3 * t166 + 84 * t13) + t259 + t25 * (t147 - t148) + t189 - t190 + t191
    t264 = t110 * t107
    t268 = 105 * t4
    t271 = P1 * P2
    t273 = Q2 * t150 * t271
    t285 = 9/2 * t12
    t291 = 105 * t1
    t306 = 378/5 * t4
    t309 = 147/5 * t8
    t310 = 252 * t4
    t313 = 5/9 * t4
    t322 = t4 * (-126/5 * t23 + 756/5 * t12)
    t323 = 336/5 * t23
    t329 = t8 * (-t128 - 42 * t12)
    t337 = t1 * (t23 - 30/7 * t12)
    t341 = t4 * (t23 / 7 + 18/7 * t12)
    t350 = t4 * (t71 - t228 + 42/5 * t166)
    t359 = t8 * (49/5 * t166 + t81 - 189/5 * t13)
    t360 = 336/5 * t166
    t366 = t1 * t169
    t369 = t4 * (t64 - 21 * t166 + t182)
    t370 = 168 * t166
    t380 = 147 * t2
    t381 = 378 * t4
    t384 = 63 * t8
    t388 = 9/5 * t4
    t395 = t2 * (-t246 - 210 * t25)
    t399 = t4 * (-126 * t23 + 756 * t25)
    t400 = 336 * t23
    t406 = t8 * (-t246 + 630 * t25)
    t415 = t1 * (t23 + 18 * t25)
    t418 = t4 * (t128 - 30 * t25)
    t423 = t12 * t19
    t429 = t2 * (t64 + 49 * t67 - 189 * t65)
    t431 = 630 * t65
    t433 = t4 * (t171 + 42 * t67 - t431)
    t434 = 336 * t67
    t441 = t8 * (t181 - 175 * t67 + 735 * t65)
    t448 = 112 * t67
    t455 = t1 * (t64 - 21 * t67 + 63 * t65)
    t458 = t4 * (t64 - 35 * t67 + t245)
    t459 = 168 * t67
    t552 = t25 + t12 - 1
    t555 = 210 * t12
    t581 = (t4 + 2) * t23
    t623 = (t1 + 2) * t23
    t656 = t1 * t23
    t665 = t4 * t23
    t693 = -504 * t48 * Q1 * t271 + 630 * t12 * t25 * t4 + 1008 * t12 * t25 + 189 * t13 * t4 + 441 * t65 * t4 + 504 * t13 + t185 - 504 * t273 - t370 - t459 + 504 * t65
    
    sTLfR = 1125/16 * t110 * t107 * (t14 * (1323/25 * t2 + t1 * (882/25 * t4 + 588/5) + 336/25 + 63/25 * t8 + 84/5 * t4) - 588/5 * t20 * (t1 + 3/7 * t4 + 10/7) * P1 + t48 * (t2 * (-49/3 * t23 + 294/5 * t25) + t1 * (t4 * (-14 * t23 + 756/5 * t25) - t33 + t34) + t8 * (-7/5 * t23 + 126/5 * t25) + t4 * (-112/15 * t23 + 504/5 * t25) - t45 + 672/25 * t25) + 28 * t12 * (t1 * (t23 - t51) + t4 * (3/5 * t23 - t51) + 8/5 * t23 - 36/5 * t25) * P2 * t50 + Q2 * (t2 * t69 + t1 * (t4 * (t71 + 378/5 * t65 - 126/5 * t67) + t76 + 252/5 * t65 - t78) + t8 * (t81 + t82 - t68) + t4 * (t85 + t86 - t78) - 112/25 * t67 + t90 + 336/25 * t65) - 4/5 * (t1 * t69 + t4 * (t64 - 35/3 * t67 + t82) + t98 - 56/3 * t67 + 42 * t65) * P2 * t50)
    cTLfR = -225/16 * t110 * t107 * (t119 * (t112 + t1 * (84 + 882/5 * t4) + 1323/5 * t8 + 588 * t4 + 336/5) - 252 * t65 * (t1 + 7/3 * t4 + 10/3) * P2 * t121 + t150 * (t131 + t1 * (t4 * (-70 * t23 + 756 * t12) - t33 + t136) + t8 * (-245/3 * t23 + 294 * t12) + t4 * (-560/3 * t23 + 840 * t12) + t147 - t148) + 84 * t163 * (t1 * (t23 - t152) + t4 * (5/3 * t23 - t152) + t158 - t159) * P1 + Q1 * (t170 + t1 * (t4 * (t171 - 126 * t166 + 378 * t13) + t176 - t177 + t178) + t8 * (t181 - t167 + t182) + t4 * (t185 - t177 + t186) + t189 - t190 + t191) - 4 * t162 * (t1 * (t64 - 35/3 * t166 + 147/5 * t13) + t4 * (t64 - t198 + t199) + t98 - 56/3 * t166 + 42 * t13) * P1)
    phiLM1fR = 15/2 * t110 * (112 * Q2 * t19 * P1 * t23 + 630 * t12 * t25 * t1 + 441 * t13 * t1 + 9 * t1 * t64 + 189 * t65 * t1 - 140 * t12 * t656 - 84 * t12 * t665 - 84 * t25 * t656 - 140 * t25 * t665 + 9 * t4 * t64 + t693) * (t121 - t19) * t107
    phiLP1fT = -225/16 * t264 * t262 * t212
    phiLM1fT = -45/4 * t264 * (t65 * (21 * t1 + t268 + 56) - 168 * t273 + t25 * (t12 * (126 * t1 + 126 * t4 + 112) - 28/3 * (t1 + 3 * t4 + 2) * t23) + 112/3 * Q1 * P2 * (t23 - t285) * t121 + t13 * (t291 + 21 * t4 + 56) - 28 * t12 * (t1 + t4 / 3 + 2/3) * t23 + (t1 + t4 + 4/3) * t64) * t212
    phiLM1sTLfT = -75/32 * t110 * t107 * (t119 * (t112 + t1 * (756/5 + t306) + t309 + t310 + 336/5) - 756/5 * t65 * t162 * (t1 + t313 + 14/3) * P1 + t150 * (t131 + t1 * (t322 - t323 + 4536/5 * t12) + t329 + t4 * (-t323 - 504/5 * t12) + t147 - t148) + 196/5 * t163 * (t337 + t341 + 24/7 * t23 - 36/7 * t12) * P1 + Q1 * (t170 + t1 * (t350 - 1008/5 * t166 + 36/5 * t64 + 756 * t13) + t359 + t4 * (t360 + t76 - 1764/5 * t13) + t189 - t190 + t191) + 4/5 * t162 * (t366 + t369 + t171 - t370 + 630 * t13) * P1)
    phiLM1cTLfT = -15/32 * t110 * t107 * (t14 * (t380 + t1 * (t381 + 1260) + t384 + 756 * t4 + 336) - 420 * t20 * (t1 + t388 + 42/5) * P1 + t48 * (t395 + t1 * (t399 - t400 - 504 * t25) + t406 + t4 * (-t400 + 4536 * t25) - t252 + 672 * t25) + 28 * t423 * (t415 + t418 + 24 * t23 - 36 * t25) * P1 + Q2 * (t429 + t1 * (t433 + t185 + t434 - 1764 * t65) + t441 + t4 * (36 * t64 - 1008 * t67 + 3780 * t65) + 8 * t64 - t448 + 336 * t65) + 4 * t19 * (t455 + t458 + t171 - t459 + t431) * P1)
    sTLfT = -375/32 * t110 * (t119 * (t112 + t1 * (252/5 + t306) + t309 + 84 * t4 + 336/25) - 756/5 * t65 * P2 * (t1 + t313 + 14/9) * t121 + t150 * (t131 + t1 * (t322 - t148 + 1512/5 * t12) + t329 + t4 * (-t148 - 168/5 * t12) + 672/25 * t12 - t45) + 196/5 * t25 * P2 * (t337 + t341 + 8/7 * t23 - 12/7 * t12) * t121 + Q1 * (t170 + t1 * (t350 - t360 + t76 + t186) + t359 + t4 * (t190 + t85 - 588/5 * t13) + 336/25 * t13 - 112/25 * t166 + t90) + 4/5 * t162 * (t366 + t369 + t98 - 56 * t166 + 210 * t13) * P1) * t107
    cTLfT = -75/32 * t110 * t107 * (t14 * (t380 + t1 * (t381 + 420) + 336/5 + t384 + t310) - 420 * t20 * (t1 + t388 + 14/5) * P1 + t48 * (t395 + t1 * (t399 - t252 - t34) + t406 + t4 * (-t252 + 1512 * t25) - t148 + 672/5 * t25) + 28 * t423 * (t415 + t418 + 8 * t23 - t239) * P1 + Q2 * (t429 + t1 * (t433 + t176 + t448 - 588 * t65) + t441 + t4 * (t185 - t434 + 1260 * t65) - t78 + t191 + t259) + 4 * t19 * (t455 + t458 + t98 - 56 * t67 + 210 * t65) * P1)
    phiLM1sTLfN = 45/8 * t264 * (t218 * (21 * t65 + t25 * (-28/3 * t23 + t555) - 140/3 * t166 + t64 + 245 * t13) + 56 * t1 * t19 * (t23 - 9/2 * t25 - 15/2 * t12) * Q2 + P1 * (t65 * (t268 + 56) + t25 * (t12 * (t381 + 336) - 28 * (t4 + 2/3) * t23) + t13 * (t268 + 280) - 28 * t12 * t581 + (t4 + 4/3) * t64) + 56/3 * t19 * Q2 * (t25 * (-15/2 * t4 - 12) + t12 * (-9/2 * t4 - 12) + t581)) * t552
    phiLM1cTLfN = 45/8 * t264 * (t4 * P2 * (245 * t65 + t25 * (-140/3 * t23 + t555) + t64 - 28/3 * t166 + 21 * t13) + 56 * t4 * Q1 * Q2 * (t23 - 15/2 * t25 - t285) * P1 + P2 * (t65 * (t291 + 280) + t25 * (t12 * (378 * t1 + 336) - 28 * t623) + t13 * (t291 + 56) - 28 * t12 * (t1 + 2/3) * t23 + (t1 + 4/3) * t64) + 56/3 * Q1 * (t25 * (-9/2 * t1 - 12) + t12 * (-15/2 * t1 - 12) + t623) * t121) * t552

    
    return sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN;
end

############### srp
# for triaxial and axysimmetric (A=B!=C)
"""
    srpaveragedgauss_wrapper1_AV(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    Returns the terms of the Gaussian Planetary equations due to the solar radiation pressure
    averaged with respect to the orbital mean anomaly 

    Function suitable only for triaxial satellites or axisymmetric satellite with one moment of inertia different from the other two.  

    This function if just an interface. Currently, it is an inteface  for srpaveragedgauss_wrapper1_AV_OURMODEL, which assumes a model of the light pressure force coherent
    with the model implemented for the attitude dynamics. If the user wants to use the simplified model model in the semi-analytical
    propagator, assuming the srp force aligned with the Earth-Sun direction,  it is sufficient to change this function: instead of calling   srpaveragedgauss_wrapper1_AV_OURMODEL 
    it has to call srpaveragedgauss_wrapper1_AV_Udirection (the inputs and outputs remain unchanged).
    
    INPUT  
    sma = semi-major axis
    eta = sqrt(1-ec^2), ec = eccentricity
    P1  = ec*sin(om+OM), om = argument of the pericenter, OM = right ascension of the ascending node
    P2  = ec*cos(om+OM) 
    Q1  = tan(incl/2)sin(OM), incl = inclination
    Q2  = tan(incl/2)cos(OM)
    GG  = 1+Q1^2+Q2^2
    satellite = Dict with the characteristics of the satellite **
    k,k1r,smk,skm,,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih, = Sadov variables related quantities
    mu  = gravitational parameter of the central body
    rSunV = position vector of the sun
    includeEclipsesEffects = boolean : if true the eclipses effects are considered, if false the body is considered in sun-light (even if it is not)
    passageinshadowoccurs = boolean: if true part of the orbit is in shadow; if false the orbit is in sun-light (if includeEclipsesEffects, set passageinshadowoccurs = false)
    Elinshadow  = value of the eccentric longitude at the entrance of the shadow region (if passageinshadowoccurs, set  Elinshadow = NaN)
    Eloutshadow = value of the eccentric longitude at the exit from the shadow region (if passageinshadowoccurs, set  Eloutshadow = NaN)

    **
    satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD);
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

    
    ***
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
    to GV and the xy plane:
    k1r = sqrt(1+k)
    smk = m/(k+m)
    skm = k/(m+k)
    T1 = ((m-1.0)*KK+EE)/m/KK;
    T2 = 2.0*(KK-EE)/KK/m^2.0 + (EE-2.0*KK)/m/KK;
    T3 = EE/KK
    KK,EE = complete elliptic integrals of first and second kind with characteristics m
    cdelta = cos(delta), sdelta = sin(delta), delta = inclination angle between 
    the XY plane and the plane perpendicular to GV
    cpsih = cos(psih), spsih = sin(psih), psih = angle between X axis and N
                
    OUTPUT
    out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt]
    with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

"""
function srpaveragedgauss_wrapper1_AV(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    out = srpaveragedgauss_wrapper1_AV_OURMODEL(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow);
    # out = srpaveragedgauss_wrapper1_AV_Udirection(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow);
  
    return out;
end

"""
    the srp acceleration is assumed in the opposite direction of the planet-sun unit vector
"""
function srpaveragedgauss_wrapper1_AV_Udirection(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
 
    PsrpCRAOm = srpPsrpCRAOm_averaged_psilpsig_V2(satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,rSunV)
    au = astroconstants(0)[2];
    earthsundist = norm(rSunV);
    psrpcorr     = (au/earthsundist)^2.0
    PsrpCRAOm =  PsrpCRAOm*psrpcorr
    out = srpaveragedgauss(sma,eta,P1,P2,Q1,Q2,GG,mu,PsrpCRAOm,-rSunV/norm(rSunV),includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow);
    return out;

end

"""
   model of srp acceleration coherent with the one used for attitude. 
"""
function srpaveragedgauss_wrapper1_AV_OURMODEL(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    Asrp = srpAcc_averaged_psilpsig(satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,rSunV/norm(rSunV))
    PsrpCRAOm = norm(Asrp);
    accdirection = Asrp/PsrpCRAOm;
    au = astroconstants(0)[2];
    earthsundist = norm(rSunV);
    psrpcorr     = (au/earthsundist)^2.0
    PsrpCRAOm =  PsrpCRAOm*psrpcorr

    out = srpaveragedgauss(sma,eta,P1,P2,Q1,Q2,GG,mu,PsrpCRAOm,accdirection,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow);
    return out;

end

# for axysimmetric (A=B=C)
"""
    srpaveragedgauss_wrapper1_AV(sma,eta,P1,P2,Q1,Q2,GG,satellite,k,k1r,smk,skm,m,T1,T2,T3,KK,cdelta,sdelta,cpsih,spsih,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    Returns the terms of the Gaussian Planetary equations due to the solar radiation pressure
    averaged with respect to the orbital mean anomaly 

    Function suitable only for axisymmetric satellites with equal principal moments of inertia. 

    This function if just an interface. Currently, it is an inteface  for srpaveragedgauss_wrapper2_AV_OURMODEL, which assumes a model of the light pressure force coherent
    with the model implemented for the attitude dynamics. If the user wants to use the simplified model model in the semi-analytical
    propagator, assuming the srp force aligned with the Earth-Sun direction,  it is sufficient to change this function: instead of calling   srpaveragedgauss_wrapper1_AV_OURMODEL 
    it has to call srpaveragedgauss_wrapper2_AV_Udirection (the inputs and outputs remain unchanged).
    
    INPUT  
    sma = semi-major axis
    eta = sqrt(1-ec^2), ec = eccentricity
    P1  = ec*sin(om+OM), om = argument of the pericenter, OM = right ascension of the ascending node
    P2  = ec*cos(om+OM) 
    Q1  = tan(incl/2)sin(OM), incl = inclination
    Q2  = tan(incl/2)cos(OM)
    GG  = 1+Q1^2+Q2^2
    satellite = Dict with the characteristics of the satellite **
    csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA = Andoyer-Serret variables related quantities ***
    mu  = gravitational parameter of the central body
    rSunV = position vector of the sun
    includeEclipsesEffects = boolean : if true the eclipses effects are considered, if false the body is considered in sun-light (even if it is not)
    passageinshadowoccurs = boolean: if true part of the orbit is in shadow; if false the orbit is in sun-light (if includeEclipsesEffects, set passageinshadowoccurs = false)
    Elinshadow  = value of the eccentric longitude at the entrance of the shadow region (if passageinshadowoccurs, set  Elinshadow = NaN)
    Eloutshadow = value of the eccentric longitude at the exit from the shadow region (if passageinshadowoccurs, set  Eloutshadow = NaN)

    **
    satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD);
    IV = [A,B,C] : principal moment of inertia
    IM = intrinsic magnetic moment of the satellite [A m^2]
    n  = number of facets in which the satellite surface is divided
    facets = [facet1...facetn]
        where  
            facetk = [facet_coeff,facet_area,facet_Vinfo,facet_nv]
            facet_coeff = [cai,cdi,csi]
            facet_area  = area of the facet [m^2]
            facet_Vinfo  = [rhoiv,vv1,vv2,vv3], rhoiv = centre of mass to facet centroid vector [m], vvi = verteces of the facet [m], i=1..3
            facet_nv    = normal unit vector  

            with 
                cai = (1-rf*sp)
                cdi = (2/3*(1-sp)*rf+2/3*(1-rf))
                csi = 2*rf*sp*psrp
                with rf reflectivity and sp specular coefficient (see Benson&Sheers2021)
    CD: aerodynamic coefficient

    ***
    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    Let GV=diag([A,B,C])omV be angular momentum of the body 
    (omV = [p,q,r]^T is the angular velocity) 
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular
    to GV and the xy plane.
    csigma=cos(sigma), ssigma=sin(sigma), sigma = inclination angle between
            the plane perpendicular to GV and the xy plane
    cdelta=cos(delta), sdelta=sin(delta), delta = inclination angle between 
            the XY plane and the plane perpendicular to GV.
    chA=cos(hA),shA=sin(hA), hA = angle between X axis and N
    clA=cos(lA),slA=sin(lA), lA = angle between N'' and x axis
   

    OUTPUT
    out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt]
    with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

"""
function srpaveragedgauss_wrapper2_AV(sma,eta,P1,P2,Q1,Q2,GG,satellite,csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    out = srpaveragedgauss_wrapper2_AV_OURMODEL(sma,eta,P1,P2,Q1,Q2,GG,satellite,csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow);
end

"""
    the srp acceleration is assumed in the opposite direction of the planet-sun unit vector
"""
function srpaveragedgauss_wrapper2_AV_Udirection(sma,eta,P1,P2,Q1,Q2,GG,satellite,csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    
    PsrpCRAOm = srpPsrpCRAOm_averaged_axisym(satellite,csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA,rSunV);

    out = srpaveragedgauss(sma,eta,P1,P2,Q1,Q2,GG,mu,PsrpCRAOm,-rSunV/norm(rSunV),includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow);
    return out;

end

"""
   model of srp acceleration coherent with the one used for attitude. 
"""
function srpaveragedgauss_wrapper2_AV_OURMODEL(sma,eta,P1,P2,Q1,Q2,GG,satellite,csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    
    Asrp = srpAcc_averaged_psilpsig_axisym(satellite,csigma,ssigma,cdelta,sdelta,clA,slA,chA,shA,rSunV/norm(rSunV));
    
    PsrpCRAOm = norm(Asrp);
    accdirection = Asrp/PsrpCRAOm;

    out = srpaveragedgauss(sma,eta,P1,P2,Q1,Q2,GG,mu,PsrpCRAOm,accdirection,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow);
    return out;

end

# functions
"""
    srpaveragedgauss(sma,eta,P1,P2,Q1,Q2,GG,mu,PsrpCRAOm,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    Returns the terms of the Gaussian Planetary equations due to the solar radiation pressure
    averaged with respect to the orbital mean anomaly

    INPUT  
    sma = semi-major axis 
    P1  = ec*sin(om+OM), om = argument of the pericenter, OM = right ascension of the ascending node
    P2  = ec*cos(om+OM) 
    eta = sqrt(1-P1^2-P2^2)
    Q1  = tan(incl/2)sin(OM), incl = inclination
    Q2  = tan(incl/2)cos(OM)
    GG  = 1+Q1^2+Q2^2
    mu  = gravitational parameter of the central body
    PsrpCRAOm = Psrp*CR*AOM, Psrp = solar radiation pressure, CR = reflectivity coefficients, AOM = area over mass radiation
    rSunV = position vector of the sun
    includeEclipsesEffects = boolean : if true the eclipses effects are considered, if false the body is considered in sun-light (even if it is not)
    passageinshadowoccurs = boolean: if true part of the orbit is in shadow; if false the orbit is in sun-light (if includeEclipsesEffects, set passageinshadowoccurs = false)
    Elinshadow  = value of the eccentric longitude at the entrance of the shadow region (if passageinshadowoccurs, set  Elinshadow = NaN)
    Eloutshadow = value of the eccentric longitude at the exit from the shadow region (if passageinshadowoccurs, set  Eloutshadow = NaN)

    OUTPUT
    out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt], out = out(sma,P1,P2,Q1,Q2; mu, rSunV/norm(rSunV), PsrpCRAOm)
    with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

"""
function srpaveragedgauss(sma,eta,P1,P2,Q1,Q2,GG,mu,PsrpCRAOm,forcedirection,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)

    srpaccUV = forcedirection;
    if includeEclipsesEffects && passageinshadowoccurs
        IntegralsV = getdefintegrals_0_EL_forconstantaccelerationininertialframeaveraged(sma,P1,P2,ELinshadow) - getdefintegrals_0_EL_forconstantaccelerationininertialframeaveraged(sma,P1,P2,ELoutofshadow);
        if ELinshadow<ELoutofshadow
            IntegralsV = IntegralsV + getdefintegrals_0_EL_forconstantaccelerationininertialframeaveraged(sma,P1,P2,2.0*pi);
        end
        IntegralsV = IntegralsV /2.0/pi;         
    else
        IntegralsV = getdefintegrals_0_EL_forconstantaccelerationininertialframeaveraged(sma,P1,P2,2.0*pi)/2.0/pi;
    end
    IntegralsV[abs.(IntegralsV).<1e-15].=0.0;
    # println([ELinshadow,ELoutofshadow,sma,P1,P2],IntegralsV)

    sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN = getconstantaccelerationininertialframeaveragedterms(srpaccUV[1],srpaccUV[2],srpaccUV[3],IntegralsV,sma,P1,P2,eta,Q1,Q2,GG);
    out = zeros(6);
    ksrp = PsrpCRAOm;
    smaOmuS = sqrt(sma/mu);
    out[1] =  2/eta*sqrt(sma^3.0/mu)*ksrp*(-P1*cTLfR + P2*sTLfR + phiLP1fT);
    out[2] =  eta*smaOmuS*ksrp*(-cTLfR + P1*phiLM1fT + phiLM1sTLfT + sTLfT - P2*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
    out[3] =  eta*smaOmuS*ksrp*( sTLfR + P2*phiLM1fT + phiLM1cTLfT + cTLfT + P1*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
    out[4] =  0.5*eta*smaOmuS*ksrp*GG*phiLM1sTLfN;
    out[5] =  0.5*eta*smaOmuS*ksrp*GG*phiLM1cTLfN;
    out[6] =  -eta*smaOmuS*ksrp*( (P1*sTLfR+P2*cTLfR +P1*phiLM1cTLfT -P2*phiLM1sTLfT + P1*cTLfT - P2*sTLfT)/(1+eta) + 2*eta*phiLM1fR + Q1*phiLM1cTLfN - Q2*phiLM1sTLfN); 
    out[abs.(out).<1e-15].=0.0;
    return out;

end

function getdefintegrals_0_EL_forconstantaccelerationininertialframeaveraged(sma,P1,P2,EL)

    Integrals = zeros(12);
    Integrals[1] = EL
    Integrals[2] = -cos(EL)+1
    Integrals[3] = sin(EL)
    Integrals[4] = sin(EL) ^ 2 / 2
    Integrals[5] = sin(EL) * cos(EL) / 2 + EL / 2

    if EL<=pi
        ANGLE = EL
    else
        ANGLE = pi;
    end

    if ANGLE == 0.0
        ATN0 = -atan(P1/sqrt(1-P1^2-P2^2));
        ATN1 = 0.0
        ATN2 = -ATN0;
    elseif ANGLE == pi
        ATN0 = pi/2.0;
        ATN1 = -pi/2.0;
        ATN2 = -pi/2.0
    else
        ATN0 = atan((tan(ANGLE/2)*P2 + tan(ANGLE/2) - P1)/sqrt(-P1^2 - P2^2 + 1))
        ATN1 = atan((cos(ANGLE) - 1)/sin(ANGLE))
        ATN2 = atan((P2*cos(ANGLE) + sin(ANGLE)*P1 + cos(EL) - P2 - 1)/(sqrt(-P1^2 - P2^2 + 1)*sin(ANGLE)));
    end

    Integrals0180_0 = zeros(7)
    Integrals0180_0[1] = -2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * atan(P1 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2)) / sma
    Integrals0180_0[2] = -2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * (-P2 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * log(1 - P2) / 2 + atan(P1 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2)) * P1) / sma / (P1 ^ 2 + P2 ^ 2)
    Integrals0180_0[3] = -(P1 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * log(1 - P2) + 2 * atan(P1 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2)) * P2) * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2)
    Integrals0180_0[4] = 2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * (P1 * P2 * (P1 ^ 2 + P2 ^ 2 - 2) * atan(P1 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2)) + sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * ((-P1 ^ 2 / 2 + P2 ^ 2 / 2) * log(1 - P2) + P1 ^ 2 * P2 + P2 ^ 3)) / sma / (P1 ^ 2 + P2 ^ 2) ^ 2
    Integrals0180_0[5] = -2 * ((P1 ^ 4 + (P2 ^ 2 - 1) * P1 ^ 2 + P2 ^ 2) * atan(P1 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2)) + P1 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * (P1 ^ 2 + P2 ^ 2 + P2 * log(1 - P2))) * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2) ^ 2
    Integrals0180_0[6] = -2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * (P1 * (P1 ^ 4 + (-P2 ^ 2 - 1) * P1 ^ 2 - 2 * P2 ^ 4 + 3 * P2 ^ 2) * atan(P1 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2)) - sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * ((P1 ^ 4 + (P2 ^ 2 - 3) * P1 ^ 2 + P2 ^ 2) * P2 * log(1 - P2) - 2 * P1 ^ 4 + 2 * P2 ^ 4) / 2) / (P1 ^ 2 + P2 ^ 2) ^ 3 / sma
    Integrals0180_0[7] = -((6 * P1 ^ 4 * P2 + (6 * P2 ^ 3 - 6 * P2) * P1 ^ 2 + 2 * P2 ^ 3) * atan(P1 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2)) + sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * ((P1 ^ 4 + (P2 ^ 2 - 1) * P1 ^ 2 + 3 * P2 ^ 2) * log(1 - P2) + 4 * P1 ^ 2 * P2 + 4 * P2 ^ 3) * P1) * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / (P1 ^ 2 + P2 ^ 2) ^ 3 / sma

    Integrals0180 = zeros(7)
    Integrals0180[1] = 2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ATN0 / sma
    Integrals0180[2] = (2 * P1 * ATN1 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) + P2 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * log(1 - sin(ANGLE) * P1 - P2 * cos(ANGLE)) - 2 * ATN2 * P1) * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2)
    Integrals0180[3] = (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2) * (-P1 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * log(1 - sin(ANGLE) * P1 - P2 * cos(ANGLE)) + 2 * P2 * ATN1 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) - 2 * ATN2 * P2)
    Integrals0180[4] = -(-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * (-2 * P1 * P2 * (P1 ^ 2 + P2 ^ 2 - 2) * ATN2 + ((P1 ^ 2 - P2 ^ 2) * log(1 - sin(ANGLE) * P1 - P2 * cos(ANGLE)) - 4 * ATN1 * P1 * P2 + (P1 ^ 2 + P2 ^ 2) * (sin(ANGLE) * P1 - P2 * cos(ANGLE) - P2)) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1)) / sma / (P1 ^ 2 + P2 ^ 2) ^ 2
    Integrals0180[5] = -2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ((P1 ^ 4 + (P2 ^ 2 - 1) * P1 ^ 2 + P2 ^ 2) * ATN2 + (log(1 - sin(ANGLE) * P1 - P2 * cos(ANGLE)) * P1 * P2 + (P1 ^ 2 - P2 ^ 2) * ATN1 + (P1 ^ 2 + P2 ^ 2) * (cos(ANGLE) * P1 + P2 * sin(ANGLE) + P1) / 2) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1)) / sma / (P1 ^ 2 + P2 ^ 2) ^ 2
    Integrals0180[6] = (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * (-2 * P1 * (P1 ^ 4 + (-P2 ^ 2 - 1) * P1 ^ 2 - 2 * P2 ^ 4 + 3 * P2 ^ 2) * ATN2 + ((P1 ^ 4 + (P2 ^ 2 - 3) * P1 ^ 2 + P2 ^ 2) * P2 * log(1 - sin(ANGLE) * P1 - P2 * cos(ANGLE)) + P1 * (P1 ^ 4 - P2 ^ 4 - 2 * P1 ^ 2 + 6 * P2 ^ 2) * ATN1 - (P1 ^ 2 + P2 ^ 2) * ((-P1 ^ 2 * P2 - P2 ^ 3) * cos(ANGLE) ^ 2 + ((P1 ^ 3 + P1 * P2 ^ 2) * sin(ANGLE) + 2 * P1 ^ 2 - 2 * P2 ^ 2) * cos(ANGLE) + 4 * sin(ANGLE) * P1 * P2 + (P2 + 2) * P1 ^ 2 + P2 ^ 2 * (P2 - 2)) / 2) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1)) / (P1 ^ 2 + P2 ^ 2) ^ 3 / sma
    Integrals0180[7] = -(-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ((6 * P1 ^ 4 * P2 + (6 * P2 ^ 3 - 6 * P2) * P1 ^ 2 + 2 * P2 ^ 3) * ATN2 + ((P1 ^ 5 + (P2 ^ 2 - 1) * P1 ^ 3 + 3 * P1 * P2 ^ 2) * log(1 - sin(ANGLE) * P1 - P2 * cos(ANGLE)) + (-3 * P1 ^ 4 * P2 + (-4 * P2 ^ 3 + 6 * P2) * P1 ^ 2 - P2 ^ 5 - 2 * P2 ^ 3) * ATN1 + (P1 ^ 2 + P2 ^ 2) * ((P1 ^ 3 + P1 * P2 ^ 2) * cos(ANGLE) ^ 2 + ((P1 ^ 2 + P2 ^ 2) * sin(ANGLE) + 4 * P1) * P2 * cos(ANGLE) + (-2 * P1 ^ 2 + 2 * P2 ^ 2) * sin(ANGLE) - P1 * (P1 ^ 2 + P2 ^ 2 - 4 * P2)) / 2) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1)) / (P1 ^ 2 + P2 ^ 2) ^ 3 / sma
    
    Integrals[6:12]  = Integrals0180 - Integrals0180_0; 


    if EL > pi
        ANGLE = 2.0*pi-EL;
        if ANGLE == 0.0
            ATN1 = 0.0;
            ATN3 = atan(P1/sqrt(1-P1^2-P2^2))
            ATN4 = ATN3
        elseif ANGLE == pi
            ATN1 = -pi/2.0;
            ATN3 = pi/2.0;
            ATN4 = pi/2.0
        else
            ATN1 = atan((cos(ANGLE) - 1)/sin(ANGLE))
            ATN3 = atan((sin(ANGLE)*P1 - (cos(ANGLE) - 1)*(P2 + 1))/(sqrt(-P1^2 - P2^2 + 1)*sin(ANGLE)))
            ATN4 = atan((tan(ANGLE/2)*P2 + tan(ANGLE/2) + P1)/sqrt(-P1^2 - P2^2 + 1));
        end

        Integrals180360_180 = zeros(7)
        Integrals180360_180[1] = (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * pi / sma
        Integrals180360_180[2] = (-P1 * pi * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) - P2 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * log(P2 + 1) + pi * P1) * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2)
        Integrals180360_180[3] = (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2) * (sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * P1 * log(P2 + 1) - P2 * pi * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) + pi * P2)
        Integrals180360_180[4] = -(-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * (((-P1 ^ 2 + P2 ^ 2) * log(P2 + 1) + 2 * pi * P1 * P2) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) + P1 * P2 * (P1 ^ 2 + P2 ^ 2 - 2) * pi) / (P1 ^ 2 + P2 ^ 2) ^ 2 / sma
        Integrals180360_180[5] = (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ((2 * log(P2 + 1) * P1 * P2 + pi * P1 ^ 2 - pi * P2 ^ 2) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) + (P1 ^ 4 + (P2 ^ 2 - 1) * P1 ^ 2 + P2 ^ 2) * pi) / (P1 ^ 2 + P2 ^ 2) ^ 2 / sma
        Integrals180360_180[6] = -(-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ((2 * (P1 ^ 4 + (P2 ^ 2 - 3) * P1 ^ 2 + P2 ^ 2) * P2 * log(P2 + 1) + P1 * (P1 ^ 4 - P2 ^ 4 - 2 * P1 ^ 2 + 6 * P2 ^ 2) * pi) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) - 2 * P1 * (P1 ^ 4 + (-P2 ^ 2 - 1) * P1 ^ 2 - 2 * P2 ^ 4 + 3 * P2 ^ 2) * pi) / (P1 ^ 2 + P2 ^ 2) ^ 3 / sma / 2
        Integrals180360_180[7] = -3/2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ((-2/3 * P1 * ((P1 ^ 2 + 3) * P2 ^ 2 + P1 ^ 4 - P1 ^ 2) * log(P2 + 1) + (P2 ^ 4 / 3 + (4/3 * P1 ^ 2 + 2/3) * P2 ^ 2 + P1 ^ 4 - 2 * P1 ^ 2) * P2 * pi) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) - 2 * ((P1 ^ 2 + 1/3) * P2 ^ 2 + P1 ^ 4 - P1 ^ 2) * P2 * pi) / (P1 ^ 2 + P2 ^ 2) ^ 3 / sma
        
        Integrals180360 = zeros(7)
        Integrals180360[1] = 2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ATN4 / sma
        Integrals180360[2] = (2 * P1 * ATN1 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) - P2 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * log(sin(ANGLE) * P1 - P2 * cos(ANGLE) + 1) + 2 * ATN3 * P1) * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2)
        Integrals180360[3] = (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2) * (sqrt(-P1 ^ 2 - P2 ^ 2 + 1) * P1 * log(sin(ANGLE) * P1 - P2 * cos(ANGLE) + 1) + 2 * P2 * ATN1 * sqrt(-P1 ^ 2 - P2 ^ 2 + 1) + 2 * ATN3 * P2)
        Integrals180360[4] = (-2 * P1 * P2 * (P1 ^ 2 + P2 ^ 2 - 2) * ATN3 + ((P1 ^ 2 - P2 ^ 2) * log(sin(ANGLE) * P1 - P2 * cos(ANGLE) + 1) + 4 * ATN1 * P1 * P2 - (P1 ^ 2 + P2 ^ 2) * (sin(ANGLE) * P1 + P2 * cos(ANGLE) + P2)) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1)) * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) / sma / (P1 ^ 2 + P2 ^ 2) ^ 2
        Integrals180360[5] = 2 * (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ((P1 ^ 4 + (P2 ^ 2 - 1) * P1 ^ 2 + P2 ^ 2) * ATN3 + (log(sin(ANGLE) * P1 - P2 * cos(ANGLE) + 1) * P1 * P2 + (-P1 ^ 2 + P2 ^ 2) * ATN1 + (P1 ^ 2 + P2 ^ 2) * (cos(ANGLE) * P1 - P2 * sin(ANGLE) + P1) / 2) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1)) / sma / (P1 ^ 2 + P2 ^ 2) ^ 2
        Integrals180360[6] = -(-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * (-2 * P1 * (P1 ^ 4 + (-P2 ^ 2 - 1) * P1 ^ 2 - 2 * P2 ^ 4 + 3 * P2 ^ 2) * ATN3 + ((P1 ^ 4 + (P2 ^ 2 - 3) * P1 ^ 2 + P2 ^ 2) * P2 * log(sin(ANGLE) * P1 - P2 * cos(ANGLE) + 1) - P1 * (P1 ^ 4 - P2 ^ 4 - 2 * P1 ^ 2 + 6 * P2 ^ 2) * ATN1 + ((P1 ^ 2 * P2 + P2 ^ 3) * cos(ANGLE) ^ 2 + ((P1 ^ 3 + P1 * P2 ^ 2) * sin(ANGLE) - 2 * P1 ^ 2 + 2 * P2 ^ 2) * cos(ANGLE) + 4 * sin(ANGLE) * P1 * P2 + (-P2 - 2) * P1 ^ 2 - P2 ^ 3 + 2 * P2 ^ 2) * (P1 ^ 2 + P2 ^ 2) / 2) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1)) / (P1 ^ 2 + P2 ^ 2) ^ 3 / sma
        Integrals180360[7] = (-P1 ^ 2 - P2 ^ 2 + 1) ^ (-1/2) * ((6 * P1 ^ 4 * P2 + (6 * P2 ^ 3 - 6 * P2) * P1 ^ 2 + 2 * P2 ^ 3) * ATN3 + ((P1 ^ 5 + (P2 ^ 2 - 1) * P1 ^ 3 + 3 * P1 * P2 ^ 2) * log(sin(ANGLE) * P1 - P2 * cos(ANGLE) + 1) + (3 * P1 ^ 4 * P2 + (4 * P2 ^ 3 - 6 * P2) * P1 ^ 2 + P2 ^ 5 + 2 * P2 ^ 3) * ATN1 + (P1 ^ 2 + P2 ^ 2) * ((P1 ^ 3 + P1 * P2 ^ 2) * cos(ANGLE) ^ 2 - ((P1 ^ 2 + P2 ^ 2) * sin(ANGLE) - 4 * P1) * P2 * cos(ANGLE) + (2 * P1 ^ 2 - 2 * P2 ^ 2) * sin(ANGLE) - P1 * (P1 ^ 2 + P2 ^ 2 - 4 * P2)) / 2) * sqrt(-P1 ^ 2 - P2 ^ 2 + 1)) / (P1 ^ 2 + P2 ^ 2) ^ 3 / sma
        
        Integrals[6:12] = Integrals[6:12] + Integrals180360_180 -  Integrals180360;
    end
    
    return Integrals;

end

function getconstantaccelerationininertialframeaveragedterms(fx,fy,fz,IntegralsV,sma,P1,P2,eta,Q1,Q2,GG)

    if size(IntegralsV)[1]!=12
        Base.error("wrong input: something is wrong");
    end
    
    # init
    sTLfR = zeros(12);
    cTLfR = zeros(12);
    phiLM1fR = zeros(12);
    phiLP1fT = zeros(12);
    phiLM1fT = zeros(12);
    phiLM1sTLfT = zeros(12);
    phiLM1cTLfT = zeros(12);
    sTLfT = zeros(12);
    cTLfT = zeros(12);
    phiLM1sTLfN = zeros(12);
    phiLM1cTLfN = zeros(12);
    
    # recurring combinations of terms 
    
    t1 = P1 ^ 2
    t2 = fy * t1
    t3 = P1 * fx
    t4 = P2 * t3
    t6 = Q2 ^ 2
    t8 = Q1 * fx
    t9 = t8 + fz
    t10 = 2 * t9
    t11 = t1 * t10
    t12 = P2 * Q1
    t13 = fy * P1
    t14 = t13 * t12
    t16 = 2 * t8
    t17 = 2 * fz
    t20 = Q1 ^ 2
    t21 = t20 + 1
    t22 = t21 * fy
    t23 = t1 * t22
    t24 = fx * t20
    t26 = 2 * Q1 * fz
    t27 = t24 + t26 - fx
    t29 = P1 * t27 * P2
    t31 = eta ^ 2
    t35 = P2 ^ 2
    t36 = t35 * fy
    t44 = -1 + P2
    t45 = P2 + 1
    t46 = t45 * t44
    t51 = fy * t20
    t61 = t35 - 2
    t62 = P1 * t61
    t64 = t44 ^ 2
    t66 = t45 ^ 2
    t74 = t66 * t64
    t79 = t61 * P2
    t80 = t27 * P1
    t85 = 0.1e1 / GG
    t86 = 1 + eta
    t87 = t86 ^ 2
    t88 = 1 / t87
    t89 = t88 * t85
    sTLfR[6] = t89 * sma * (t31 * (t6 * (-t2 + t4 - fy) + Q2 * (t11 + 2 * t14 + t16 + t17) + t23 - t29 + t22) + eta * (t6 * (-2 * t2 + 3 * t4 + 2 * t36 - 2 * fy) + Q2 * (4 * t1 * t9 - 4 * t9 * t46 + 6 * t14) + t1 * (2 * t51 + 2 * fy) - 3 * t29 - 2 * t46 * t22) + t6 * (-t62 * P2 * fx - t66 * t64 * fy - t2) + Q2 * (-2 * P1 * t61 * fy * t12 + 2 * t9 * t74 + t11) + t23 + t80 * t79 + t74 * t22)
    t90 = Q2 * fy
    t91 = t90 - fz
    t92 = 2 * t91
    t94 = t6 + 1
    t95 = t94 * fx
    t96 = Q1 * t92 - t24 + t95
    t97 = P2 * t35
    t98 = t97 * t96
    t99 = Q2 * Q1
    t102 = fy * t6
    t104 = 2 * Q2 * fz
    t105 = 2 * fx * t99 + fy - t102 + t104 + t51
    t106 = t105 * P1
    t107 = t35 * t106
    t109 = -t92
    t110 = Q1 * t109
    t111 = fx * t6
    t112 = t24 + t110 - t111 - fx
    t113 = t1 + eta + 1
    t114 = t113 * t112
    t116 = t86 * P1
    t117 = t105 * t116
    t121 = 1 / t86
    t122 = t85 * t121
    sTLfR[7] = t122 * (P2 * t114 + 2 * t107 - 2 * t117 + t98) * sma
    t123 = P2 * t13
    t124 = fx * t35
    t145 = -t111 - fx
    t148 = (t102 - t104 - fy) * P1
    t154 = t85 * t121 * sma
    sTLfR[8] = -2 * t154 * P1 * (t20 * (t123 - t124 - fx * t31 / 2 - eta * fx / 2) + Q1 * (t31 * t91 + eta * t91 + 2 * (P1 * Q2 * fx + P2 * t91) * P2) + t31 * t95 / 2 + eta * t95 / 2 - (P2 * t145 + t148) * P2)
    t157 = t35 ^ 2
    t158 = t157 * t96
    t159 = t97 * t106
    t160 = 2 * t159
    t161 = t1 - t31 + 1
    t163 = t35 * t112 * t161
    t164 = t105 * P2
    t168 = t87 * t112 * eta
    sTLfR[9] = -t89 * sma * (-2 * t164 * t116 + t158 + t160 + t163 + t168)
    t174 = -Q2 * t10 - fy + t102 - t51
    t175 = t157 * t174
    t178 = -2 * fy * t99 - fx - t111 + t24 + t26
    t179 = t178 * P1
    t180 = t97 * t179
    t181 = 2 * t180
    t183 = Q2 * t10 - t102 + t22
    t184 = 2 * eta
    t185 = t1 + t184 + 2
    t188 = eta - 1
    t189 = t188 * P1
    t190 = t86 * P2
    t191 = t178 * t190
    t192 = t191 * t189
    sTLfR[10] = t89 * (t35 * t185 * t183 - t87 * t183 + t175 - t181 - t192) * sma
    t196 = P2 * sma
    t197 = t1 + t87
    t200 = eta + 2
    t201 = t200 * P1
    t202 = t105 * t86
    cTLfR[6] = -t88 * t85 * (P2 * t112 * t197 - t202 * t201 + t107) * t196
    t208 = -t1 + t35 - eta - 1
    t214 = -t208
    cTLfR[7] = -t154 * P2 * (t6 * (fy * t208 + 2 * t4) + Q2 * (2 * t9 * t214 + 4 * t14) + fy * t214 * t21 - 2 * t29)
    t225 = P2 * t96
    t227 = t225 + t106 / 2
    cTLfR[8] = -2 * t122 * (t31 * t227 + eta * t227 + (t225 + t106) * t35) * sma
    t237 = t35 * t161 * t183
    t238 = P1 * eta
    t241 = t183 * eta
    t242 = t87 * t241
    cTLfR[9] = t89 * (-2 * t191 * t238 + t175 - t181 + t237 + t242) * sma
    t245 = 2 * t31
    t246 = t1 - t245 - t184
    t249 = t105 * t190
    t250 = t249 * t189
    cTLfR[10] = t89 * sma * (t35 * t112 * t246 - t87 * t112 * t31 + t158 + t160 + t250)
    t256 = P1 * t1
    t258 = eta + 3
    t259 = t258 * P2
    t263 = t35 - 3/2 * eta - 3/2
    t269 = 1 / t31
    t270 = t269 * t122
    phiLM1fR[6] = -t270 * sma * (-2 * P1 * t183 * t263 - t1 * t178 * t259 + t256 * t86 * t183 - t191)
    t272 = t178 * P2
    t273 = t256 * t272
    t275 = t35 - 3 * eta - 3
    t278 = eta + 3/2
    t279 = P1 * t278
    t282 = t35 - eta - 1
    phiLM1fR[7] = -t270 * sma * (t1 * t275 * t183 + t282 * t183 + 2 * t279 * t272 + t273)
    t287 = t256 * t164
    t288 = t31 / 3
    t289 = eta / 3
    t290 = t35 + t288 + t289
    t295 = P1 * (t35 - t184 - 5/2)
    t299 = t31 / 2
    t300 = eta / 2
    t301 = t35 * t278 + t299 + t300
    phiLM1fR[8] = t270 * sma * (-3 * t1 * t290 * t112 - 2 * t112 * t301 - 2 * t295 * t164 + t287)
    t306 = t97 * t174
    t307 = t258 * P1
    t308 = t35 * t178
    t311 = t1 * t200 + eta + 1
    phiLM1fR[9] = -2 * t270 * sma * (P2 * t311 * t183 - t178 * t86 * t238 - t308 * t307 + t306)
    t324 = eta + 5
    t326 = t35 * t324 - t184 - 2
    t330 = t35 * t258 + t184 + t245
    phiLM1fR[10] = t270 * sma * (-t1 * t258 * t112 * P2 - P1 * t326 * t105 + t112 * t330 * P2 + t256 * t202)
    t335 = 3 * t180
    t336 = t1 + t289 + 1/3
    t340 = t246 * P1
    phiLM1fR[11] = t270 * sma * (-t86 * t1 * t183 + 3 * t35 * t336 * t183 + t272 * t340 + t175 - t335)
    t346 = t157 * t112
    t347 = 3 * t159
    t348 = t1 - t288 - t289
    t352 = t185 * P1
    t354 = eta * t1
    t355 = t86 * t112
    phiLM1fR[12] = -t270 * (-3 * t35 * t112 * t348 + t164 * t352 - t355 * t354 + t346 - t347) * sma
    t360 = sma * t31
    t361 = P1 * t112
    t362 = t361 + t164
    phiLP1fT[6] = -t85 * t362 * t360
    t371 = P1 * P2 * Q2 * fx
    t372 = 2 * t371
    t378 = P2 * t148
    phiLP1fT[7] = t85 * t121 * (t20 * (t86 * fx + t123 - t124) + Q1 * (t35 * t92 - 2 * t91 * t86 + t372) + t35 * t95 - t378 - t86 * t95) * t360
    t383 = t31 * t183
    t384 = P2 * fy
    t389 = P1 * Q1 * fy
    t393 = P2 * t22
    phiLP1fT[8] = t85 * t121 * (t383 + t241 + (t6 * (-t3 - t384) + Q2 * (P2 * t10 - 2 * t389) + t393 + t80) * P2) * t360
    t399 = t256 * t355
    t400 = t1 * t105
    t401 = t400 * t259
    t402 = t112 * t263
    phiLM1fT[6] = -t270 * (-2 * P1 * t402 + t249 + t399 + t401) * sma
    t412 = t112 * t282
    phiLM1fT[7] = t270 * sma * (-t1 * t275 * t112 + 2 * t279 * t164 + t287 - t412)
    phiLM1fT[8] = t270 * sma * (3 * t1 * t290 * t183 + 2 * t183 * t301 - 2 * t295 * t272 + t273)
    t424 = t35 * t105
    phiLM1fT[9] = -2 * t270 * sma * (P2 * t112 * t311 + t202 * t238 + t424 * t307 + t98)
    phiLM1fT[10] = t270 * sma * (-P1 * t326 * t112 - t330 * t164 + t399 + t401)
    phiLM1fT[11] = -t270 * sma * (t86 * t112 * t1 - 3 * t35 * t336 * t112 + t164 * t340 + t346 - t347)
    t451 = t86 * t1
    phiLM1fT[12] = -t270 * sma * (3 * t35 * t183 * t348 + t451 * t241 + t272 * t352 + t175 - t335)
    t457 = t256 * t105 * t190
    t461 = 5 * eta
    t465 = t282 ^ 2
    t466 = t112 * t465
    t469 = t269 * t89
    phiLM1sTLfT[6] = t469 * sma * (t457 - 2 * t451 * t402 - P1 * (t35 * t200 - t245 - t461 - 3) * t164 + t466)
    t470 = t87 * t112
    t472 = 4 * eta
    t473 = t35 - t31 - t472 - 3
    phiLM1sTLfT[7] = -t469 * sma * (P1 * t275 * t412 - t400 * t473 * P2 - t105 * t282 * t190 + t256 * t470)
    t483 = P2 * t157
    t484 = t483 * t112
    t485 = t157 * t106
    t490 = (t1 + t461 + 5) * P1
    t496 = t105 * t87
    t497 = t496 * t238
    phiLM1sTLfT[8] = -t469 * sma * (t484 - 2 * t485 - 2 * t97 * t114 + t424 * t490 + t190 * t112 * (t1 * t324 + eta + 1) + 2 * t497)
    t502 = t1 + 1
    t503 = t502 * t183
    t504 = eta * t31
    t511 = -t9
    t516 = t1 + 2/3
    t527 = 2 * t35 * t179
    t528 = t1 + 1/2
    phiLM1sTLfT[9] = 2 * t469 * (t504 * t503 / 2 + t31 * (P2 * t179 + t503) + eta * (t157 * (t102 / 2 + Q2 * t511 - t22 / 2) - t180 + 3/2 * t35 * t516 * t183 + t272 * (t1 + 3) * P1 + t503 / 2) + (t306 - t527 + 2 * P2 * t528 * t183 + t178 * (t1 + 2) * P1) * P2) * sma
    t547 = eta + 5/3
    phiLM1sTLfT[10] = -t469 * sma * (t457 - 4 * t1 * (t35 * (eta + 5/4) - t87 / 2) * t112 - 3 * P1 * (t35 * t547 + t86 * t200 * t188 / 3) * t164 + t466)
    t559 = t483 * t174
    t560 = t157 * t179
    t561 = 3 * t560
    t562 = t1 - t288 + 1/3
    phiLM1sTLfT[11] = -t469 * sma * (t559 - t561 + 3 * t97 * t183 * t562 + t308 * (t1 + t472 + 4) * P1 + t190 * (t1 * t188 + eta + t31) * t183 - t178 * t87 * P1)
    t578 = 3 * t485
    t579 = 2/3 * eta
    t585 = (t1 - t245 + 2) * P1
    phiLM1sTLfT[12] = t469 * (t484 - t578 - 3 * t97 * t112 * (t1 + t579 + 2/3) + t424 * t585 + 2 * t190 * t112 * (t1 + t300 + 1/2) + t497) * sma
    t598 = 2 * t35
    phiLM1cTLfT[6] = t469 * sma * P2 * (t399 + 2 * t1 * t278 * t164 - P1 * (-t245 + eta * (t35 - 5) + t598 - 3) * t112 + t105 * t87 * P2)
    phiLM1cTLfT[7] = -t469 * (-t1 * t112 * t473 + t117 * t259 - t86 * t412 + t287) * t196
    t626 = t178 * t87 * t238
    phiLM1cTLfT[8] = -t469 * (-2 * t560 + 3 * t97 * (t1 + t288 + 4/3 * eta + 1) * t183 + t308 * t490 + 2 * P2 * t113 * t86 * t241 + 2 * t626) * sma
    t631 = P2 * t106
    t635 = -t94 * fx + t110 + t24
    t636 = t502 * t635
    t637 = t636 / 3
    phiLM1cTLfT[9] = 3 * t469 * sma * (t504 * (2/3 * t631 + t637) + t31 * (2 * t631 + 2/3 * t636) + eta * (t157 * (-t24 / 3 + 2/3 * Q1 * t91 + t95 / 3) + 4/3 * t159 + t35 * t635 * t516 + 4/3 * t631 + t637) + 4/3 * (t35 * (-t24 / 2 + Q1 * t91 + t95 / 2) + 3/2 * t631 + t635 * t528) * t35)
    t687 = t87 * t383
    phiLM1cTLfT[10] = -t469 * sma * (-2 * t157 * t278 * t183 - 3 * t97 * t547 * t179 + 2 * t35 * (-t504 - 3 * t31 + eta * (t1 - 2) + 3/2 * t1) * t183 + P2 * t178 * (t1 - t31 - eta + 2) * t116 - t687)
    phiLM1cTLfT[11] = t469 * sma * (t484 - t578 - 3 * t97 * t112 * t562 + t424 * (t1 - 4 * t31 - t472) * P1 - t190 * t112 * (eta * t502 - t1 + t31) - t496 * t31 * P1)
    phiLM1cTLfT[12] = t469 * sma * (t559 - t561 + 3 * t97 * (t1 - 2/3 * t31 - t579) * t183 + t308 * t585 + 2 * P2 * t86 * (t1 - t299 - t300) * t241 + t626)
    t719 = fx * t1
    t722 = t1 * t109
    t723 = 2 * t90
    t726 = t1 * t145
    sTLfT[6] = t89 * sma * (t31 * (t20 * (t719 + t123 + fx) + Q1 * (t722 + t372 - t723 + t17) + t726 - t378 - t111 - fx) + eta * (t20 * (2 * t719 + 3 * t123 - 2 * t124 + 2 * fx) + Q1 * (-4 * t1 * t91 + 4 * t91 * t46 + 6 * t371) + 2 * t1 * t145 - 3 * t378 + 2 * t46 * t95) + t20 * (t66 * t64 * fx - t62 * t384 + t719) + Q1 * (-2 * P1 * t61 * fx * P2 * Q2 - 2 * t91 * t74 + t722) + t726 + t148 * t79 - t74 * t95)
    sTLfT[7] = -t122 * sma * (P2 * t113 * t183 + 2 * t178 * t116 + t306 - t527)
    sTLfT[8] = -2 * t154 * P1 * (t6 * (-t4 - t36 - fy * t31 / 2 - eta * fy / 2) + Q2 * (t31 * t9 + eta * t9 - 2 * (P2 * t511 + t389) * P2) + t31 * t22 / 2 + eta * t22 / 2 + (t393 + t80) * P2)
    sTLfT[9] = t89 * (2 * t272 * t116 + t175 - t181 + t237 + t242) * sma
    sTLfT[10] = t89 * sma * (t35 * t112 * t185 + t158 + t160 + t250 - t470)
    cTLfT[6] = t88 * t85 * sma * (P2 * t105 * t197 + t355 * t201 - t35 * t361) * P2
    cTLfT[7] = -t154 * (t20 * (fx * t214 + 2 * t123) + Q1 * (-2 * t91 * t214 + 4 * t371) - fx * t214 * t94 - 2 * t378) * P2
    t831 = t164 + t361 / 2
    cTLfT[8] = -2 * t122 * sma * (eta * t831 + t31 * t831 + t362 * t35)
    cTLfT[9] = t89 * sma * (2 * t249 * t238 + t158 + t160 + t163 + t168)
    cTLfT[10] = -t89 * sma * (t35 * t246 * t183 + t175 - t181 - t192 - t687)
    t850 = fz * (t20 + t6 - 1) - t16 + t723
    phiLM1sTLfN[1] = -t270 * P1 * (t35 - t184 - 2) * t850
    phiLM1sTLfN[2] = -t270 * (t354 + t1 - t35 + eta + 1) * t850
    phiLM1sTLfN[3] = -t121 * t269 * t85 * P1 * P2 * t200 * t850
    phiLM1sTLfN[4] = t270 * P2 * t214 * t850
    t868 = fz * t20 + fz * t6 - fz - t16 + t723
    phiLM1sTLfN[5] = t270 * P1 * (t598 - eta - 1) * t868
    phiLM1cTLfN[1] = t270 * P2 * t113 * t850
    phiLM1cTLfN[2] = phiLM1sTLfN[3]
    phiLM1cTLfN[3] = -t270 * (t31 + eta * (t35 + 1) + t598) * t850
    phiLM1cTLfN[4] = t270 * P1 * (t598 + t31 + eta) * t868
    phiLM1cTLfN[5] = -t270 * P2 * (t1 - t35 - t31 - eta) * t850


    #
    sTLfR       = LinearAlgebra.dot(sTLfR,IntegralsV)
    cTLfR       = LinearAlgebra.dot(cTLfR,IntegralsV)
    phiLM1fR    = LinearAlgebra.dot(phiLM1fR,IntegralsV)
    phiLP1fT    = LinearAlgebra.dot(phiLP1fT,IntegralsV)
    phiLM1fT    = LinearAlgebra.dot(phiLM1fT,IntegralsV)
    phiLM1sTLfT = LinearAlgebra.dot(phiLM1sTLfT,IntegralsV)
    phiLM1cTLfT = LinearAlgebra.dot(phiLM1cTLfT,IntegralsV)
    sTLfT       = LinearAlgebra.dot(sTLfT,IntegralsV)
    cTLfT       = LinearAlgebra.dot(cTLfT,IntegralsV)
    phiLM1sTLfN = LinearAlgebra.dot(phiLM1sTLfN,IntegralsV)
    phiLM1cTLfN = LinearAlgebra.dot(phiLM1cTLfN,IntegralsV)

    return sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN
end

############### 3bp
"""
    tbaveragedgauss_wrapper1(sma,P1,P2,Q1,Q2,rV_3b,muPlanet,mu3b,deg)
    Returns the terms of the Gaussian Planetary equations due to gravity influence
    of a third body averaged with respect to the orbital mean anomaly

    INPUT  
    sma = semi-major axis 
    P1  = ec*sin(om+OM), om = argument of the pericenter, OM = right ascension of the ascending node
    P2  = ec*cos(om+OM) 
    Q1  = tan(incl/2)sin(OM), incl = inclination
    Q2  = tan(incl/2)cos(OM)
    rV_3b = position vector of the third body (in the inertial ref frame)
    muPlanet  = gravitational parameter of the central body
    mu3b      = gravitational parameter of the third body
    deg = integer between 2 and 5 --> the equations of motion are expanded in a series of the small parameter r/|rV_3b| with r the 
        distance between the satellite and the central body. 

    OUTPUT
    out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt], out = out(sma,P1,P2,Q1,Q2; muPlanet, rV_3b, mu3b)
    with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

"""
function tbaveragedgauss_wrapper1(sma,P1,P2,Q1,Q2,rV_3b,muPlanet,mu3b,deg)
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
    out = tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rV_3b,muPlanet,mu3b,deg);
    return out;
end

"""
    tbaveragedgauss_wrapper1(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rV_3b,muPlanet,mu3b,deg)
    Returns the terms of the Gaussian Planetary equations due to gravity influence
    of a third body averaged with respect to the orbital mean anomaly

    INPUT  
    sma = semi-major axis 
    P1  = ec*sin(om+OM), om = argument of the pericenter, OM = right ascension of the ascending node
    P2  = ec*cos(om+OM) 
    eta = sqrt(1-P1^2-P2^2)
    Q1  = tan(incl/2)sin(OM), incl = inclination
    Q2  = tan(incl/2)cos(OM)
    GG  = 1+Q1^2+Q2^2
    q11 = (1.0-Q1^2.0+Q2^2.0)/GG;
    q12 = (2.0*Q1*Q2)/GG;
    q13 = (-2.0*Q1)/GG;
    q21 = (2.0*Q1*Q2)/GG;
    q22 = (1.0+Q1^2.0-Q2^2.0)/GG;
    q23 = (2.0*Q2)/GG;
    q31 = (2.0*Q1)/GG;
    q32 = (-2.0*Q2)/GG;
    q33 = (1.0-Q1^2.0-Q2^2.0)/GG;
    rV_3b = position vector of the third body (in the inertial ref frame)
    muPlanet  = gravitational parameter of the central body
    mu3b      = gravitational parameter of the third body
    deg = integer between 2 and 5 --> the equations of motion are expanded in a series of the small parameter r/|rV_3b| with r the 
        distance between the satellite and the central body. 

    OUTPUT
    out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt], out = out(sma,P1,P2,Q1,Q2; muPlanet, rV_3b, mu3b)
    with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

"""
function tbaveragedgauss(sma,P1,P2,eta,Q1,Q2,GG,q11,q12,q13,q21,q22,q23,q31,q32,q33,rV_3b,muPlanet,mu3b,deg)
    r3 = norm(rV_3b);
    xuv3I = rV_3b[1]/r3;
    yuv3I = rV_3b[2]/r3;
    zuv3I = rV_3b[3]/r3;

    xuv3 = q11*xuv3I + q12*yuv3I + q13*zuv3I
    yuv3 = q21*xuv3I + q22*yuv3I + q23*zuv3I
    zuv3 = q31*xuv3I + q32*yuv3I + q33*zuv3I
    
    n = sqrt(muPlanet/sma^3.0);
    
    out = getaveragedterms_3bp_deg2(n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3b);

    if deg>=3
        out = out + getaveragedterms_3bp_deg3(sma,n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3b);
    end

    if deg>=4
        out = out + getaveragedterms_3bp_deg4(sma,n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3b);
    end

    if deg>=5
        out = out + getaveragedterms_3bp_deg5(sma,n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3b);
    end

    if deg>=6
        out = out + getaveragedterms_3bp_deg6(sma,n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3b);
    end

    out[abs.(out).<1e-15].=0.0;
    # println("out = ",out)
    return out;
end

function getaveragedterms_3bp_deg2(n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3)

    dsmadt = 0.0;
    t1 = P1 ^ 2
    t2 = P2 * t1
    t4 = Q1 * xuv3 * zuv3
    t7 = Q2 * yuv3
    t8 = zuv3 * t7
    t11 = P2 ^ 2
    t12 = t11 * P1
    t14 = yuv3 * Q1 * zuv3
    t17 = xuv3 * Q2
    t18 = zuv3 * t17
    t21 = P2 * t11
    t23 = xuv3 * zuv3
    t27 = yuv3 * zuv3
    t30 = P1 * t1
    t34 = xuv3 ^ 2
    t37 = yuv3 ^ 2
    t40 = xuv3 * yuv3
    t63 = 15 * P1 * xuv3 * yuv3 - 3 * t23 * P2 * Q1 + 3 * t27 * P2 * Q2 - 12 * t23 * Q1 * t21 - 3 * t27 * Q2 * t21 - 15 * yuv3 * xuv3 * t30 + 12 * t34 * P2 - 3 * t37 * P2 - 15 * t14 * t12 + 15 * t18 * t12 - 15 * t40 * t12 - 12 * t34 * t2 + 3 * t37 * t2 + 3 * t4 * t2 + 12 * t8 * t2 - 12 * t34 * t21 + 3 * t37 * t21 - 3 * P2 + 3 * t2 + 3 * t21
    t64 = 0.1e1 / n
    t66 = 0.1e1 / eta
    t67 = r3 ^ 2
    t69 = 0.1e1 / r3 / t67
    t70 = t69 * t66
    dP1dt = t70 * t64 * t63 / 2
    t115 = 3 * t23 * P1 * Q1 - 3 * t27 * P1 * Q2 - 15 * P2 * xuv3 * yuv3 - 3 * t23 * Q1 * t30 - 12 * t27 * Q2 * t30 + 15 * yuv3 * xuv3 * t21 + 3 * t34 * P1 - 12 * t37 * P1 - 3 * t34 * t12 + 12 * t37 * t12 + 12 * t4 * t12 + 3 * t8 * t12 + 15 * t14 * t2 - 15 * t18 * t2 + 15 * t40 * t2 - 3 * t34 * t30 + 12 * t37 * t30 + 3 * P1 - 3 * t12 - 3 * t30
    dP2dt = t64 * t69 * t66 * t115 / 2
    t119 = Q1 ^ 2
    t120 = Q1 * t119
    t121 = t120 * t1
    t123 = t119 * t1
    t124 = Q2 ^ 2
    t125 = yuv3 * t124
    t128 = Q1 * t1
    t129 = Q2 * t124
    t130 = xuv3 * t129
    t132 = t124 ^ 2
    t136 = P1 * P2
    t137 = Q2 * t120
    t138 = yuv3 * t137
    t141 = t124 * t119
    t142 = xuv3 * t141
    t145 = t129 * Q1
    t146 = yuv3 * t145
    t152 = t120 * t11
    t155 = t119 * t11
    t157 = Q1 * t11
    t167 = 8 * yuv3 * t124 * t1 + 4 * yuv3 * t132 * t1 - yuv3 * t132 * t11 + 5 * xuv3 * t132 * t136 + t17 * t121 + 4 * t125 * t123 - t125 * t155 + t130 * t128 + 2 * t17 * t128 - 4 * t130 * t157 - 5 * t138 * t136 + 5 * t142 * t136 - 5 * t146 * t136 - 4 * t17 * t152
    t168 = Q2 * Q1
    t169 = yuv3 * t168
    t172 = xuv3 * t124
    t180 = xuv3 * t137
    t181 = yuv3 * t141
    t182 = xuv3 * t145
    t189 = xuv3 * t168
    t192 = -2 * yuv3 * t124 * t11 + 4 * yuv3 * t1 - yuv3 * t11 + yuv3 * t132 - 10 * t169 * t136 + 10 * t172 * t136 + 5 * xuv3 * t136 - 8 * t17 * t157 + 2 * t125 - t180 + t181 - t182 - 2 * t189 + yuv3
    t195 = 0.1e1 / GG
    t197 = t64 * t70
    dQ1dt = 3/4 * t197 * t195 * zuv3 * (t167 + t192)
    t199 = t119 ^ 2
    t205 = yuv3 * t129
    t228 = xuv3 * t199 * t1 - 4 * xuv3 * t199 * t11 - 5 * yuv3 * t199 * t136 + 4 * t7 * t121 + t172 * t123 + 2 * xuv3 * t123 + 4 * t205 * t128 + 8 * t7 * t128 + 5 * t180 * t136 - 5 * t181 * t136 + 5 * t182 * t136 - t7 * t152 - 4 * t172 * t155 - t205 * t157
    t247 = -10 * yuv3 * t119 * t136 + xuv3 * t1 - 4 * xuv3 * t11 - 2 * xuv3 * t119 + 10 * t189 * t136 - 5 * yuv3 * t136 - 8 * xuv3 * t155 - 2 * t7 * t157 - xuv3 * t199 + t138 - t142 + t146 + 2 * t169 - xuv3
    dQ2dt = -3/4 * t197 * t195 * zuv3 * (t228 + t247)
    t253 = xuv3 * eta
    t254 = zuv3 * t253
    t257 = Q2 * t1
    t259 = eta * yuv3 * zuv3
    t270 = Q2 * t11
    t277 = eta * t1
    t293 = eta * t11
    t298 = -15 * t259 * Q1 * t136 + 15 * t254 * Q2 * t136 - 30 * yuv3 * t253 * t136 + 3 * t23 * t128 + 3 * t254 * t128 - 15 * t14 * t136 + 15 * t18 * t136 - 12 * t23 * t157 - 12 * t254 * t157 + 12 * t259 * t257 + 12 * t27 * t257 - 3 * t259 * t270 - 3 * t27 * t270 + 3 * t34 * t277 - 12 * t37 * t277 - 12 * t34 * t293 + 3 * t37 * t293
    t326 = -3 * t23 * eta * Q1 + 3 * t27 * eta * Q2 - 3 * t34 * eta - 3 * t37 * eta + 3 * t34 * t1 + 3 * t37 * t1 + 3 * t34 * t11 + 3 * t37 * t11 + 2 * eta - 2 * t1 - 2 * t11 + 3 * t277 + 3 * t293 - 3 * t34 - 3 * t37 - 3 * t4 + 3 * t8 + 2
    dMLdt = 1 / (1 + eta) * t70 * t64 * (t298 + t326) / 2

    out = mu3*[dsmadt,dP1dt,dP2dt,dQ1dt,dQ2dt,dMLdt]
    return out;
        
end

function getaveragedterms_3bp_deg3(sma,n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3)

    dsmadt = 0.0;
    t1 = P1 ^ 2
    t2 = P1 * t1
    t3 = P2 * t2
    t5 = yuv3 * xuv3
    t6 = zuv3 * t5
    t9 = xuv3 ^ 2
    t10 = t9 * Q2
    t11 = zuv3 * t10
    t14 = yuv3 ^ 2
    t15 = t14 * Q2
    t16 = zuv3 * t15
    t19 = P2 ^ 2
    t20 = t19 * t1
    t22 = zuv3 * t9 * Q1
    t26 = zuv3 * t14 * Q1
    t32 = P2 * t19
    t33 = t32 * P1
    t41 = t19 ^ 2
    t42 = Q1 * t41
    t43 = zuv3 * t9
    t46 = zuv3 * t14
    t52 = 10 * t6 * Q1 * t3 - 60 * t6 * Q1 * t33 + 60 * t6 * Q2 * t20 - 10 * t6 * Q2 * t41 - 5 * t11 * t3 + 30 * t11 * t33 + 20 * t16 * t3 - 15 * t16 * t33 + 15 * t22 * t20 - 30 * t26 * t20 - 20 * t43 * t42 + 5 * t46 * t42
    t53 = t1 ^ 2
    t54 = xuv3 * t9
    t57 = xuv3 * t53
    t60 = t9 * yuv3
    t63 = yuv3 * t14
    t68 = t14 * xuv3
    t77 = xuv3 * t41
    t80 = zuv3 * Q2
    t83 = Q1 * zuv3
    t88 = -30 * t14 * t57 + 15 * t14 * t77 - 15 * t54 * t20 - 15 * t68 * t20 + 3 * t83 * t20 - 60 * t60 * t3 + 10 * t63 * t3 - 3 * t80 * t3 - 60 * t60 * t33 + 10 * t63 * t33 - 3 * t80 * t33 - 20 * t54 * t41 + 5 * t54 * t53
    t90 = P2 * P1
    t100 = Q1 * t19
    t115 = xuv3 * t1
    t120 = -10 * t6 * Q1 * t90 + 10 * t6 * Q2 * t19 - 10 * t54 * t1 - 15 * t43 * t100 - 5 * t46 * t100 + 5 * t11 * t90 + 25 * t14 * t115 + 15 * t16 * t90 + 12 * xuv3 * t20 + 6 * yuv3 * t3 + 6 * yuv3 * t33 + 3 * zuv3 * t42 + 3 * t57
    t128 = xuv3 * t19
    t141 = 4 * zuv3 * t100 - 20 * t14 * t128 + 15 * t54 * t19 + 60 * t60 * t90 - 10 * t63 * t90 - 4 * t80 * t90 - 6 * yuv3 * t90 + t115 - 5 * t128 + 5 * t54 + 5 * t68 + 9 * t77 - 4 * xuv3
    t147 = 0.1e1 / n / eta
    t148 = r3 ^ 2
    t149 = t148 ^ 2
    t150 = 0.1e1 / t149
    t151 = t150 * t147
    dP1dt = -15/16 * t151 * sma * (t52 + t88 + t120 + t141)
    t157 = Q2 * t53
    t183 = -60 * t6 * Q1 * t20 + 10 * t6 * Q1 * t53 + 60 * t6 * Q2 * t3 - 10 * t6 * Q2 * t33 + 30 * t11 * t20 - 5 * t43 * t157 + 20 * t46 * t157 - 15 * t16 * t20 + 15 * t22 * t3 - 20 * t22 * t33 - 30 * t26 * t3 + 5 * t26 * t33
    t212 = -30 * yuv3 * t9 * t41 + 15 * yuv3 * t9 * t53 - 3 * zuv3 * t157 - 15 * t60 * t20 - 15 * t63 * t20 - 3 * t80 * t20 + 10 * t54 * t3 - 60 * t68 * t3 + 3 * t83 * t3 + 10 * t54 * t33 - 60 * t68 * t33 + 5 * t63 * t41 - 20 * t63 * t53
    t217 = Q2 * t1
    t244 = -10 * t6 * Q1 * t1 + 10 * t6 * Q2 * t90 - 20 * yuv3 * t9 * t1 + 15 * t63 * t1 + 12 * yuv3 * t20 + 5 * t43 * t217 + 15 * t46 * t217 - 15 * t22 * t90 - 5 * t26 * t90 + 6 * xuv3 * t3 + 3 * t83 * t33 + 6 * xuv3 * t33 + 9 * yuv3 * t53
    t268 = 25 * yuv3 * t9 * t19 - 5 * yuv3 * t1 - 10 * t63 * t19 + yuv3 * t19 - 4 * zuv3 * t217 + 3 * yuv3 * t41 - 10 * t54 * t90 + 60 * t68 * t90 + 4 * t83 * t90 - 6 * xuv3 * t90 + 5 * t60 + 5 * t63 - 4 * yuv3
    dP2dt = 15/16 * t151 * sma * (t183 + t212 + t244 + t268)
    t273 = Q1 ^ 2
    t274 = Q1 * t273
    t275 = t274 * t2
    t277 = Q2 * xuv3 * yuv3
    t280 = t273 * t2
    t281 = Q2 ^ 2
    t282 = t9 * t281
    t285 = t14 * t281
    t288 = Q1 * t2
    t289 = Q2 * t281
    t291 = yuv3 * xuv3 * t289
    t294 = t281 ^ 2
    t295 = t294 * t2
    t300 = P2 * t1
    t301 = Q2 * t274
    t302 = t9 * t301
    t305 = t14 * t301
    t308 = t273 * t300
    t310 = yuv3 * xuv3 * t281
    t313 = t289 * Q1
    t314 = t9 * t313
    t318 = t14 * t313
    t325 = t19 * P1
    t329 = t281 * t273
    t330 = t9 * t329
    t333 = t14 * t329
    t336 = Q1 * t325
    t345 = t274 * t32
    t350 = t273 * t32
    t353 = 60 * yuv3 * xuv3 * t294 * t300 - 15 * t14 * t294 * t325 - 60 * t277 * t274 * t325 + 30 * t9 * t294 * t325 - 20 * t10 * t345 + 5 * t15 * t345 - 60 * t291 * t336 - 30 * t318 * t300 - 10 * t310 * t350 + 30 * t330 * t325 - 15 * t333 * t325
    t355 = Q1 * t32
    t356 = t9 * t289
    t359 = t14 * t289
    t370 = t281 * t2
    t380 = Q2 * Q1
    t381 = t9 * t380
    t384 = t14 * t380
    t399 = t274 * P1
    t402 = t273 * P1
    t407 = -120 * t277 * t336 - 10 * t277 * t399 + 60 * t282 * t325 + 5 * t282 * t402 - 30 * t285 * t325 + 15 * t285 * t402 - 3 * t294 * t325 + 120 * t310 * t300 + 30 * t381 * t300 - 60 * t384 * t300 - 3 * t329 * t325
    t410 = P1 * Q1
    t413 = t294 * P1
    t429 = t274 * P2
    t435 = t273 * P2
    t438 = P2 * Q1
    t447 = t9 * t2
    t459 = 10 * t5 * t294 * P2 + 20 * t14 * t2 - 6 * t281 * t325 + 6 * t380 * t300 + 60 * t5 * t300 + 10 * t310 * t435 + 30 * t9 * t325 - 15 * t356 * t438 - 5 * t359 * t438 - 6 * t370 - 5 * t447
    t468 = t281 * P1
    t475 = xuv3 * t32
    t493 = t9 * P1
    t499 = P2 * xuv3
    t503 = 20 * t5 * t281 * P2 + 15 * t14 * P1 + 8 * Q2 * t438 - 30 * t10 * t438 - 10 * t15 * t438 + 10 * yuv3 * t499 - 4 * P1 - 3 * t2 - 3 * t325 - 8 * t468 + 5 * t493
    t511 = 0.1e1 / GG * t150 * t147
    dQ1dt = -15/32 * t511 * sma * (-10 * t5 * t294 * t32 - 20 * t5 * t281 * t32 + t503 + t459 + 3 * t301 * t300 + 3 * t313 * t300 - 20 * t356 * t355 + 5 * t359 * t355 - 3 * t281 * t280 + 20 * t277 * t288 + 10 * t291 * t288 - 5 * t9 * t295 + 20 * t14 * t295 + 15 * t302 * t300 - 30 * t305 * t300 + 60 * t310 * t308 + 15 * t314 * t300 + 10 * t277 * t275 - 5 * t282 * t280 + 20 * t285 * t280 + t353 - 20 * t277 * t410 + 10 * t9 * t468 + 30 * t14 * t468 + 6 * Q2 * t355 - 10 * t291 * t410 + 5 * t9 * t413 + 15 * t14 * t413 + 3 * Q2 * t345 + 3 * t289 * t355 - 40 * t10 * t355 + 10 * t15 * t355 - 15 * t10 * t429 - 5 * t15 * t429 - 10 * t9 * t370 + 40 * t14 * t370 - 10 * yuv3 * t475 + 4 * Q2 * t429 + 4 * t289 * t438 - 15 * t14 * t325 - 4 * t281 * t402 - 3 * t295 + t407 - 4 * t413) * zuv3
    t514 = t273 ^ 2
    t542 = Q1 * t300
    t560 = t514 * t32
    t567 = -60 * yuv3 * xuv3 * t514 * t325 - 60 * t310 * t273 * t325 + 5 * t14 * t560 - 10 * t277 * t345 + 60 * t291 * t542 - 30 * t333 * t300 + 30 * t302 * t325 - 15 * t305 * t325 + 30 * t314 * t325 - 15 * t318 * t325 - 20 * t9 * t560
    t617 = -120 * yuv3 * xuv3 * t273 * t325 - 10 * t5 * t514 * P1 - 60 * t14 * t273 * t300 + 30 * t9 * t273 * t300 + 5 * t10 * t399 + 15 * t15 * t399 + 120 * t277 * t542 - 3 * t301 * t325 - 3 * t313 * t325 + 60 * t381 * t325 - 30 * t384 * t325
    t635 = t514 * P2
    t663 = 10 * yuv3 * xuv3 * t2 - 6 * Q2 * t288 - 30 * t14 * t300 + 10 * t277 * t429 - 15 * t282 * t435 - 5 * t285 * t435 + 10 * t291 * t438 + 15 * t9 * t300 - 6 * t380 * t325 - 60 * t5 * t325 + 6 * t308
    t703 = -10 * P1 * xuv3 * yuv3 - 5 * t14 * P2 - 15 * t9 * P2 - 8 * Q2 * t410 - 10 * t14 * t435 + 20 * t277 * t438 - 30 * t9 * t435 + 4 * P2 + 3 * t300 + 3 * t32 + 8 * t435
    dQ2dt = 15/32 * t511 * sma * (10 * t5 * t514 * t2 + 15 * t9 * t514 * t300 - 30 * t14 * t514 * t300 + 60 * t277 * t274 * t300 + t703 + 20 * t5 * t280 - 3 * t289 * t288 - 10 * t10 * t288 + 15 * t330 * t300 + 20 * t359 * t288 - 5 * t356 * t288 + 10 * t310 * t280 + 20 * t15 * t275 - 5 * t10 * t275 - 20 * t9 * t32 + 5 * t14 * t32 + 4 * t281 * t435 - 4 * Q2 * t399 + 5 * t356 * t410 + 15 * t359 * t410 + 3 * t281 * t350 - 40 * t9 * t350 + 10 * t14 * t350 - 20 * t277 * t355 - 15 * t9 * t635 - 5 * t14 * t635 - 10 * t310 * t402 + 40 * t15 * t288 + 3 * t514 * t300 + 3 * t329 * t300 - 20 * t282 * t350 + 5 * t285 * t350 - 10 * t291 * t355 - 3 * Q2 * t275 - 20 * t5 * t402 - 4 * t289 * t410 + 10 * t10 * t410 + 30 * t15 * t410 + 6 * t350 + 3 * t560 + t617 + t567 + 4 * t635 + t663) * zuv3
    t713 = Q2 * t2
    t714 = t9 * eta
    t715 = zuv3 * t714
    t719 = zuv3 * t14 * eta
    t726 = Q2 * t300
    t727 = eta * xuv3
    t729 = yuv3 * zuv3 * t727
    t734 = Q2 * t325
    t744 = Q2 * t32
    t754 = eta * t2
    t768 = 15 * t54 * eta * t300 - 10 * t6 * eta * t744 + 15 * t22 * t300 - 30 * t26 * t300 + 10 * t6 * t288 + 5 * t719 * t355 - 5 * t43 * t713 + 20 * t46 * t713 + 60 * t6 * t726 + 15 * t60 * t754 - 20 * t63 * t754
    t791 = eta * t32
    t797 = eta * zuv3
    t823 = 3 * Q1 * eta * zuv3 * t300 - 3 * Q2 * eta * zuv3 * t325 - 10 * t6 * eta * t410 + 10 * t63 * t2 + 10 * t54 * t300 + 10 * t68 * t300 + 10 * t60 * t325 + 10 * t63 * t325 + 10 * yuv3 * t447 + 15 * t68 * t791 - 3 * t797 * t713
    t826 = P1 * Q2
    t841 = Q2 * P2
    t865 = P1 * eta
    t874 = 9 * eta * yuv3 * t325 + 9 * t727 * t300 + 3 * t83 * t300 - 3 * t80 * t325 + 3 * zuv3 * t355 - 10 * t6 * t410 + 5 * t43 * t826 + 15 * t46 * t826 - 15 * t60 * t865 - 15 * t63 * t865 + 9 * xuv3 * t791
    t882 = eta * P2
    t918 = -10 * t63 * P1 + 8 * P1 * yuv3 - 10 * t54 * P2 - 10 * t14 * t499 + 4 * t797 * t438 + 4 * zuv3 * t438 - 4 * zuv3 * t826 + 12 * yuv3 * t865 + 12 * xuv3 * t882 - 8 * t475 + 8 * t499
    dMLdt = -15/16 * t151 / (1 + eta) * sma * (10 * t6 * eta * t288 - 90 * t14 * t727 * t300 - 90 * yuv3 * t714 * t325 + 15 * t63 * eta * t325 + 10 * t6 * eta * t841 + t918 + 9 * yuv3 * t754 - 15 * t715 * t438 - 5 * t719 * t438 - 3 * zuv3 * t713 - 15 * t68 * t882 - 8 * yuv3 * t2 - 8 * xuv3 * t300 - 8 * yuv3 * t325 - 4 * t797 * t826 - 10 * yuv3 * t493 - 15 * t43 * t438 + 5 * t715 * t826 + 15 * t719 * t826 + 3 * t797 * t355 + 10 * t54 * t32 + 10 * t14 * t475 - 5 * t46 * t438 + 10 * t6 * t841 - 15 * t54 * t882 + t823 + 20 * t719 * t713 + 15 * t715 * t542 - 30 * t719 * t542 - 20 * t715 * t355 + 60 * t729 * t726 - 60 * t729 * t336 - 5 * t715 * t713 + 30 * t715 * t734 - 15 * t719 * t734 - 15 * t16 * t325 - 20 * t54 * t791 - 60 * t6 * t336 + 30 * t11 * t325 - 20 * t43 * t355 + 5 * t46 * t355 - 10 * t6 * t744 + t874 + t768)
    
    out = mu3*[dsmadt,dP1dt,dP2dt,dQ1dt,dQ2dt,dMLdt]
    return out;
        
end

function getaveragedterms_3bp_deg4(sma,n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3)

    dsmadt = 0.0;
    t1 = P1 ^ 2
    t2 = t1 ^ 2
    t3 = P2 * t2
    t5 = xuv3 ^ 2
    t6 = xuv3 * t5
    t7 = t6 * Q1 * zuv3
    t10 = Q1 * t3
    t11 = yuv3 ^ 2
    t12 = t11 * xuv3
    t13 = zuv3 * t12
    t16 = Q2 * t3
    t17 = yuv3 * t5
    t18 = zuv3 * t17
    t21 = yuv3 * t11
    t22 = t21 * Q2
    t23 = zuv3 * t22
    t26 = P1 * t1
    t27 = P2 ^ 2
    t28 = t27 * t26
    t29 = Q1 * t28
    t33 = zuv3 * t21 * Q1
    t36 = t6 * Q2
    t37 = zuv3 * t36
    t40 = Q2 * t28
    t43 = P2 * t27
    t44 = t43 * t1
    t47 = Q1 * t44
    t50 = Q2 * t44
    t55 = t27 ^ 2
    t56 = t55 * P1
    t57 = Q1 * t56
    t60 = 42 * t13 * t10 + 294 * t13 * t40 - 357 * t13 * t47 - 42 * t18 * t16 + 147 * t18 * t29 + 357 * t18 * t50 - 294 * t18 * t57 + 56 * t23 * t3 - 84 * t23 * t44 - 98 * t33 * t28 - 49 * t37 * t28 - 7 * t7 * t3 + 84 * t7 * t44
    t65 = Q2 * t56
    t68 = P2 * t55
    t69 = Q1 * t68
    t70 = zuv3 * t6
    t75 = Q2 * t68
    t78 = zuv3 * t21
    t81 = P1 * t2
    t85 = xuv3 * t81
    t88 = t5 ^ 2
    t91 = t5 * t11
    t94 = t11 ^ 2
    t97 = yuv3 * t6
    t100 = t21 * xuv3
    t103 = 49 * yuv3 * t6 * t81 + 49 * t100 * t28 - 147 * t13 * t65 + 42 * t13 * t69 - 42 * t18 * t75 - 98 * t21 * t85 - 245 * t97 * t28 + 42 * t88 * t3 - 357 * t91 * t3 + 42 * t94 * t3 + 49 * t33 * t56 + 98 * t37 * t56 - 56 * t70 * t69 + 7 * t78 * t75
    t117 = t5 * t68
    t123 = Q1 * xuv3 * zuv3
    t126 = Q2 * yuv3
    t127 = zuv3 * t126
    t131 = yuv3 * Q1 * zuv3
    t134 = xuv3 * Q2
    t135 = zuv3 * t134
    t140 = 147 * t100 * t56 + 84 * t11 * t117 - 3 * t123 * t3 + 15 * t123 * t44 - 18 * t127 * t3 + 21 * t131 * t28 - 21 * t135 * t28 - 14 * t88 * t44 - 273 * t91 * t44 + 35 * t94 * t44 - 294 * t97 * t56 - 56 * t88 * t68 - 7 * t94 * t68
    t143 = P2 * t1
    t146 = Q1 * t143
    t149 = Q2 * t143
    t158 = t27 * P1
    t159 = Q1 * t158
    t166 = Q2 * t158
    t169 = xuv3 * zuv3
    t172 = yuv3 * zuv3
    t175 = Q1 * t43
    t178 = -15 * t127 * t44 - 35 * t13 * t146 + 147 * t13 * t166 + 21 * t131 * t56 - 21 * t135 * t56 + 84 * t23 * t143 + 14 * t7 * t143 + 35 * t18 * t149 - 49 * t33 * t158 + 49 * t37 * t158 - 147 * t18 * t159 + 18 * t169 * t69 + 3 * t172 * t75 - 84 * t70 * t175
    t183 = Q2 * t43
    t194 = xuv3 * yuv3
    t200 = xuv3 * t26
    t211 = -98 * yuv3 * t6 * t26 + 15 * t11 * t3 + 9 * t11 * t44 - 35 * t13 * t175 - 84 * t88 * t143 + 322 * t91 * t143 + 35 * t18 * t183 - 14 * t78 * t183 + 84 * t194 * t28 + 49 * t21 * t200 + 15 * t5 * t3 + 51 * t5 * t44 + 21 * yuv3 * t85
    t225 = t5 * t43
    t239 = -196 * t100 * t158 - 119 * t11 * t225 - 6 * t11 * t68 - t123 * t143 - 41 * t127 * t143 + 42 * t131 * t158 - 42 * t135 * t158 - 35 * t94 * t143 + 245 * t97 * t158 + 41 * t169 * t175 + 63 * t194 * t56 + 14 * t88 * t43 + 14 * t94 * t43 + 36 * t117
    t242 = P2 * Q1
    t247 = P2 * Q2
    t265 = xuv3 * P1
    t268 = 49 * yuv3 * t6 * P1 - 16 * t11 * t143 - 7 * t13 * t242 + 26 * t5 * t143 - 21 * t194 * t158 + t172 * t183 + 7 * t18 * t247 + 21 * yuv3 * t200 + 49 * t21 * t265 - 7 * t70 * t242 + 7 * t78 * t247 - 3 * t3 - 6 * t44
    t275 = t5 * P2
    t289 = t11 * P2 + 42 * t88 * P2 - 7 * t94 * P2 + 35 * t11 * t275 + 5 * t11 * t43 + 4 * t169 * t242 - 4 * t172 * t247 - 42 * yuv3 * t265 + 4 * P2 - t143 + 5 * t225 - 41 * t275 - t43 - 3 * t68
    t293 = sma ^ 2
    t295 = r3 ^ 2
    t296 = t295 ^ 2
    t300 = 0.1e1 / n / r3 / t296
    t301 = 0.1e1 / eta
    t302 = t301 * t300
    dP1dt = 15/16 * t302 * t293 * (t60 + t103 + t140 + t178 + t211 + t239 + t268 + t289)
    t304 = Q1 * t81
    t309 = Q2 * t81
    t332 = -147 * t18 * t10 - 294 * t13 * t16 + 357 * t13 * t29 - 42 * t13 * t304 + 42 * t18 * t309 - 357 * t18 * t40 + 294 * t18 * t47 + 84 * t23 * t28 - 84 * t7 * t28 + 98 * t33 * t3 + 49 * t37 * t3 + 7 * t70 * t304 - 56 * t78 * t309
    t349 = t5 * t81
    t362 = 294 * t100 * t3 - 84 * t11 * t349 + 147 * t13 * t50 - 42 * t13 * t57 + 42 * t18 * t65 - 7 * t23 * t56 - 35 * t88 * t28 + 273 * t91 * t28 - 147 * t97 * t3 - 49 * t33 * t44 - 98 * t37 * t44 + 56 * t7 * t56 + 7 * t88 * t81 + 56 * t94 * t81
    t379 = xuv3 * t68
    t392 = 98 * yuv3 * t6 * t68 + 245 * t100 * t44 - 15 * t123 * t28 - 21 * t131 * t3 + 21 * t135 * t3 + 3 * t169 * t304 + 18 * t172 * t309 - 49 * t21 * t379 + 14 * t94 * t28 - 49 * t97 * t44 - 42 * t88 * t56 + 357 * t91 * t56 - 42 * t94 * t56
    t395 = Q1 * t26
    t400 = Q2 * t26
    t423 = -18 * t123 * t56 + 15 * t127 * t28 - 3 * t127 * t56 - 147 * t13 * t149 + 35 * t13 * t395 - 21 * t131 * t44 + 21 * t135 * t44 + 49 * t33 * t143 - 49 * t37 * t143 + 147 * t18 * t146 + 84 * t7 * t158 - 35 * t18 * t400 - 14 * t70 * t395 - 84 * t78 * t400
    t443 = t5 * t26
    t452 = -51 * t11 * t28 + 119 * t11 * t443 - 36 * t11 * t81 + 35 * t13 * t159 + 196 * t97 * t143 + 14 * t23 * t158 - 35 * t18 * t166 - 63 * t194 * t3 - 84 * t194 * t44 - 14 * t88 * t26 - 14 * t94 * t26 - 9 * t5 * t28 + 6 * t349
    t470 = xuv3 * t43
    t482 = -49 * yuv3 * t6 * t43 - 245 * t100 * t143 - 15 * t11 * t56 - 41 * t123 * t158 - 42 * t131 * t143 + 42 * t135 * t143 + 35 * t88 * t158 - 322 * t91 * t158 + 84 * t94 * t158 + t169 * t395 + 41 * t172 * t400 + 98 * t21 * t470 - 21 * yuv3 * t379 - 15 * t5 * t56
    t485 = P1 * Q1
    t490 = P1 * Q2
    t507 = -26 * t11 * t158 - 5 * t11 * t26 - t127 * t158 + 7 * t13 * t485 + 21 * t194 * t143 + 16 * t5 * t158 - 7 * t18 * t490 + 7 * t70 * t485 - 7 * t78 * t490 + 6 * t28 - 5 * t443 + 3 * t56 + 3 * t81
    t510 = t5 * P1
    t520 = xuv3 * P2
    t532 = -49 * yuv3 * t6 * P2 + 41 * t11 * P1 + 7 * t88 * P1 - 42 * t94 * P1 - 35 * t11 * t510 - 4 * t169 * t485 + 4 * t172 * t490 - 49 * t21 * t520 - 21 * yuv3 * t470 + 42 * yuv3 * t520 - 4 * P1 + t158 + t26 - t510
    dP2dt = 15/16 * t302 * t293 * (t332 + t362 + t392 + t423 + t452 + t482 + t507 + t532)
    t538 = Q1 ^ 2
    t539 = Q1 * t538
    t540 = t539 * t2
    t543 = t11 * t134
    t546 = t538 * t2
    t547 = Q2 ^ 2
    t549 = yuv3 * t5 * t547
    t552 = t21 * t547
    t555 = Q1 * t2
    t556 = Q2 * t547
    t557 = t6 * t556
    t560 = xuv3 * t556
    t561 = t11 * t560
    t564 = t547 ^ 2
    t565 = t564 * t2
    t570 = P2 * t26
    t571 = t539 * t570
    t572 = t5 * Q2
    t573 = yuv3 * t572
    t576 = Q2 * t539
    t577 = t21 * t576
    t580 = t547 * t538
    t581 = t6 * t580
    t584 = 42 * t17 * t565 - 56 * t21 * t565 + 7 * t36 * t540 - 42 * t543 * t540 + 42 * t549 * t546 - 56 * t552 * t546 + 7 * t557 * t555 - 42 * t561 * t555 + 98 * t577 * t570 + 49 * t581 * t570 - 147 * t573 * t571
    t585 = t538 * t570
    t586 = xuv3 * t547
    t587 = t11 * t586
    t590 = Q1 * t570
    t592 = yuv3 * t5 * t556
    t595 = t556 * Q1
    t596 = t21 * t595
    t599 = t6 * t564
    t602 = xuv3 * t564
    t603 = t11 * t602
    t606 = t27 * t1
    t607 = t6 * t576
    t610 = t539 * t606
    t613 = t538 * t606
    t616 = t21 * t580
    t619 = t6 * t595
    t622 = Q1 * t606
    t626 = yuv3 * t5 * t564
    t629 = 357 * t543 * t610 - 357 * t549 * t613 + 357 * t561 * t622 + 98 * t596 * t570 + 49 * t599 * t570 - 294 * t603 * t570 - 294 * t587 * t585 - 147 * t592 * t590 - 84 * t607 * t606 + 84 * t616 * t606 - 84 * t619 * t606 - 357 * t626 * t606
    t631 = t21 * t564
    t634 = t43 * P1
    t635 = t539 * t634
    t642 = t538 * t634
    t645 = Q1 * t634
    t654 = t539 * t55
    t659 = t538 * t55
    t662 = 56 * t36 * t654 - 42 * t543 * t654 + 42 * t549 * t659 + 294 * t573 * t635 - 49 * t577 * t634 - 98 * t581 * t634 + 147 * t587 * t642 + 294 * t592 * t645 - 49 * t596 * t634 - 98 * t599 * t634 + 147 * t603 * t634 + 84 * t631 * t606
    t665 = Q1 * t55
    t670 = t564 * t55
    t677 = yuv3 * t547
    t688 = t547 * t2
    t691 = 3 * t134 * t540 + 42 * t17 * t670 + 84 * t17 * t688 - 7 * t21 * t670 + 14 * t36 * t555 - 84 * t543 * t555 + 18 * t677 * t546 - 7 * t552 * t659 + 3 * t560 * t555 + 56 * t557 * t665 - 42 * t561 * t665 + 18 * yuv3 * t565
    t696 = yuv3 * t576
    t699 = xuv3 * t580
    t702 = yuv3 * t595
    t707 = Q2 * Q1
    t708 = t21 * t707
    t713 = t6 * t547
    t718 = xuv3 * t576
    t721 = yuv3 * t580
    t724 = -112 * t21 * t688 - 588 * t587 * t570 + 21 * t602 * t570 - 21 * t696 * t570 + 21 * t699 * t570 - 21 * t702 * t570 + 196 * t708 * t570 + 98 * t713 * t570 - 294 * t573 * t590 - 15 * t718 * t606 + 15 * t721 * t606
    t725 = xuv3 * t595
    t728 = t6 * t707
    t733 = yuv3 * t564
    t740 = t539 * t1
    t745 = t538 * t1
    t750 = Q1 * t1
    t755 = -14 * t36 * t740 + 714 * t543 * t622 + 35 * t543 * t740 - 714 * t549 * t606 - 35 * t549 * t745 + 168 * t552 * t606 - 84 * t552 * t745 - 14 * t557 * t750 + 35 * t561 * t750 - 15 * t725 * t606 - 168 * t728 * t606 + 15 * t733 * t606
    t757 = t564 * t1
    t778 = P2 * P1
    t779 = t539 * t778
    t784 = -35 * t17 * t757 - 84 * t21 * t757 + 588 * t573 * t645 + 147 * t573 * t779 + 49 * t577 * t778 + 294 * t587 * t634 + 21 * t602 * t634 - 21 * t696 * t634 + 21 * t699 * t634 - 21 * t702 * t634 - 98 * t708 * t634 - 196 * t713 * t634
    t787 = t538 * t778
    t790 = Q1 * t778
    t811 = -18 * t134 * t654 + 112 * t36 * t665 - 84 * t543 * t665 - 18 * t560 * t665 - 49 * t581 * t778 - 147 * t587 * t787 + 147 * t592 * t790 + 49 * t596 * t778 - 49 * t599 * t778 - 147 * t603 * t778 - 3 * t677 * t659 - 3 * yuv3 * t670
    t815 = t547 * t55
    t820 = t539 * t27
    t825 = t538 * t27
    t830 = Q1 * t27
    t835 = t564 * t27
    t842 = 6 * t134 * t555 + 84 * t17 * t815 - 35 * t17 * t835 - 14 * t21 * t815 + 14 * t21 * t835 + 84 * t36 * t820 + 35 * t543 * t820 - 35 * t549 * t825 + 14 * t552 * t825 + 84 * t557 * t830 + 35 * t561 * t830
    t845 = t5 * t2
    t850 = yuv3 * t707
    t859 = xuv3 * t707
    t869 = -294 * t12 * t570 + t134 * t740 - 357 * t17 * t606 - 56 * t21 * t2 + 84 * t21 * t606 + 42 * t586 * t570 + 49 * t6 * t570 - 42 * t850 * t570 + 30 * t677 * t606 - 30 * t859 * t606 + 36 * yuv3 * t688 + 42 * yuv3 * t845
    t880 = t547 * t1
    t895 = 147 * t12 * t634 - 70 * t17 * t880 - 168 * t21 * t880 - 28 * t36 * t750 + 70 * t543 * t750 + t560 * t750 + 42 * t586 * t634 - 98 * t6 * t634 - 42 * t850 * t634 + 41 * t677 * t745 - 42 * t696 * t778 + 41 * yuv3 * t757
    t914 = t5 * t55
    t921 = -36 * t134 * t665 - 41 * t134 * t820 - 7 * t21 * t55 + 294 * t573 * t790 - 294 * t587 * t778 + 42 * t602 * t778 + 42 * t699 * t778 - 42 * t702 * t778 + 98 * t708 * t778 - 98 * t713 * t778 - 6 * yuv3 * t815 + 42 * yuv3 * t914
    t932 = t547 * t27
    t944 = 7 * t12 * t576 - 7 * t17 * t580 - 70 * t17 * t932 + 28 * t21 * t932 + 168 * t36 * t830 + 70 * t543 * t830 - 41 * t560 * t830 - t677 * t825 - yuv3 * t835 + 7 * t607 - 7 * t616 + 7 * t619
    t959 = t5 * t1
    t968 = -84 * t21 * t1 + 7 * t12 * t595 + 2 * t134 * t750 + 18 * yuv3 * t2 + 21 * xuv3 * t570 + 15 * yuv3 * t606 + 21 * xuv3 * t634 - 84 * t850 * t778 + 82 * yuv3 * t880 - 35 * yuv3 * t959 - 7 * t626 - 7 * t631
    t982 = t5 * t27
    t991 = -147 * t12 * t778 - 82 * t134 * t830 + 14 * t21 * t27 - 3 * yuv3 * t55 + 84 * t586 * t778 - 49 * t6 * t778 - 2 * yuv3 * t932 - 35 * yuv3 * t982 - 4 * t718 + 4 * t721 - 4 * t725 + 14 * t728
    t1007 = 41 * yuv3 * t1 + 14 * t12 * t707 - yuv3 * t27 + 42 * xuv3 * t778 - 7 * t17 - 7 * t21 - 14 * t549 - 14 * t552 + 8 * t677 + 4 * t733 - 8 * t859 + 4 * yuv3
    t1016 = t300 * t301 / GG
    dQ1dt = -15/32 * t1016 * t293 * zuv3 * (t584 + t629 + t662 + t691 + t724 + t755 + t784 + t811 + t842 + t869 + t895 + t921 + t944 + t968 + t991 + t1007)
    t1019 = t538 ^ 2
    t1020 = t1019 * t2
    t1035 = t21 * t556
    t1039 = yuv3 * t5 * t1019
    t1042 = t21 * t1019
    t1047 = -42 * t12 * t1020 + 7 * t6 * t1020 - 56 * t1035 * t555 - 147 * t1039 * t570 + 98 * t1042 * t570 - 56 * t22 * t540 + 42 * t573 * t540 - 42 * t587 * t546 + 7 * t713 * t546 + 42 * t592 * t555 + 49 * t607 * t570
    t1058 = t6 * t1019
    t1061 = xuv3 * t1019
    t1062 = t11 * t1061
    t1075 = -84 * t1058 * t606 + 357 * t1062 * t606 - 294 * t543 * t571 - 147 * t549 * t585 - 294 * t561 * t590 + 98 * t616 * t570 + 49 * t619 * t570 - 357 * t573 * t610 + 84 * t577 * t606 - 84 * t581 * t606 + 357 * t587 * t613 - 357 * t592 * t622
    t1095 = t1019 * t55
    t1102 = 294 * t1039 * t634 - 49 * t1042 * t634 - 42 * t12 * t1095 + 56 * t6 * t1095 + 147 * t543 * t635 + 294 * t549 * t642 + 147 * t561 * t645 + 42 * t573 * t654 + 84 * t596 * t606 - 98 * t607 * t634 - 49 * t616 * t634 - 98 * t619 * t634
    t1123 = yuv3 * t556
    t1128 = 3 * xuv3 * t1020 - 7 * t1035 * t665 + 18 * t1123 * t555 - 84 * t12 * t546 + 18 * t126 * t540 - 7 * t22 * t654 + 3 * t586 * t546 + 14 * t6 * t546 + 84 * t573 * t555 - 42 * t587 * t659 + 42 * t592 * t665 + 56 * t713 * t659
    t1133 = yuv3 * t1019
    t1141 = yuv3 * t5 * t538
    t1144 = t21 * t538
    t1157 = -15 * t1061 * t606 - 21 * t1133 * t570 - 294 * t1141 * t570 + 196 * t1144 * t570 - 112 * t22 * t555 - 588 * t543 * t590 + 21 * t718 * t570 - 21 * t721 * t570 + 21 * t725 * t570 + 98 * t728 * t570 + 15 * t696 * t606
    t1160 = t6 * t538
    t1163 = xuv3 * t538
    t1164 = t11 * t1163
    t1173 = t1019 * t1
    t1186 = -168 * t1160 * t606 + 714 * t1164 * t606 + 35 * t12 * t1173 - 14 * t6 * t1173 - 84 * t22 * t740 - 714 * t573 * t622 - 35 * t573 * t740 + 35 * t587 * t745 - 15 * t699 * t606 + 15 * t702 * t606 + 168 * t708 * t606 - 14 * t713 * t745
    t1212 = -84 * t1035 * t750 + 147 * t1039 * t778 + 49 * t1042 * t778 - 21 * t1133 * t634 + 588 * t1141 * t634 - 98 * t1144 * t634 + 294 * t543 * t645 - 35 * t592 * t750 + 21 * t718 * t634 - 21 * t721 * t634 + 21 * t725 * t634 - 196 * t728 * t634
    t1237 = -18 * xuv3 * t1095 - 3 * t1123 * t665 - 84 * t12 * t659 - 3 * t126 * t654 - 147 * t543 * t779 + 147 * t549 * t787 - 147 * t561 * t790 - 18 * t586 * t659 + 112 * t6 * t659 - 49 * t607 * t778 + 49 * t616 * t778 - 49 * t619 * t778
    t1245 = t1019 * t27
    t1264 = 14 * t1035 * t830 + 35 * t12 * t1245 + 84 * t6 * t1245 - 14 * t22 * t665 + 14 * t22 * t820 + 6 * xuv3 * t546 + 84 * t573 * t665 - 35 * t573 * t820 + 35 * t587 * t825 - 35 * t592 * t830 + 84 * t713 * t825
    t1269 = xuv3 * t2
    t1272 = yuv3 * t538
    t1290 = -42 * t11 * t1269 - 30 * t1163 * t606 + xuv3 * t1173 + 357 * t12 * t606 + 36 * t126 * t555 - 42 * t1272 * t570 - 147 * t17 * t570 + 7 * t6 * t2 + 98 * t21 * t570 + 42 * t859 * t570 - 84 * t6 * t606 + 30 * t850 * t606
    t1315 = 41 * t1123 * t750 - 42 * t1133 * t778 + 70 * t12 * t745 + 41 * t126 * t740 - 42 * t1272 * t634 + 294 * t17 * t634 - 49 * t21 * t634 - 168 * t22 * t750 - 70 * t573 * t750 + t586 * t745 - 28 * t6 * t745 + 42 * t859 * t634
    t1336 = xuv3 * t55
    t1341 = -42 * t11 * t1336 + 294 * t1141 * t778 + 98 * t1144 * t778 - 41 * xuv3 * t1245 - 6 * t126 * t665 - 294 * t543 * t790 + 56 * t6 * t55 - 36 * xuv3 * t659 + 42 * t718 * t778 - 42 * t721 * t778 + 42 * t725 * t778 - 98 * t728 * t778
    t1362 = -t1123 * t830 + 70 * t12 * t825 - t126 * t820 - 7 * t17 * t576 + 28 * t22 * t830 - 70 * t573 * t830 - 41 * t586 * t825 + 168 * t6 * t825 + 7 * t1058 + 7 * t1062 - 7 * t577 + 7 * t581
    t1379 = xuv3 * t1
    t1386 = -14 * t6 * t1 + 35 * t11 * t1379 + 7 * t12 * t580 + 82 * t126 * t750 - 84 * t1272 * t778 - 7 * t17 * t595 - 21 * yuv3 * t570 - 15 * xuv3 * t606 - 21 * yuv3 * t634 + 2 * xuv3 * t745 + 3 * t1269 - 7 * t596
    t1401 = xuv3 * t27
    t1408 = 35 * t11 * t1401 - 2 * t126 * t830 + 147 * t17 * t778 + 49 * t21 * t778 + 84 * t6 * t27 + 84 * t859 * t778 - 82 * xuv3 * t825 - 4 * t1061 + 14 * t1160 - 18 * t1336 + 4 * t696 - 4 * t699
    t1422 = -14 * t17 * t707 - 42 * yuv3 * t778 - 8 * t1163 + 14 * t1164 + 7 * t12 + t1379 - 41 * t1401 + 7 * t6 + 4 * t702 - 14 * t708 + 8 * t850 - 4 * xuv3
    dQ2dt = 15/32 * t1016 * t293 * zuv3 * (t1047 + t1075 + t1102 + t1128 + t1157 + t1186 + t1212 + t1237 + t1264 + t1290 + t1315 + t1341 + t1362 + t1386 + t1408 + t1422)
    t1430 = t6 * eta
    t1431 = zuv3 * t1430
    t1437 = Q2 * t2
    t1442 = zuv3 * t21 * eta
    t1445 = t5 * eta
    t1446 = t172 * t1445
    t1451 = Q2 * t570
    t1454 = eta * xuv3
    t1456 = zuv3 * t11 * t1454
    t1463 = Q2 * t606
    t1468 = -210 * t13 * eta * t555 + 210 * t18 * eta * t1437 + 245 * t1431 * t1451 + 35 * t1431 * t555 - 420 * t1431 * t622 - 280 * t1442 * t1437 + 420 * t1442 * t1463 + 490 * t1442 * t590 - 1785 * t1446 * t1463 - 735 * t1446 * t590 - 1470 * t1456 * t1451 + 1785 * t1456 * t622
    t1473 = Q2 * t634
    t1483 = Q2 * t55
    t1497 = -210 * t13 * eta * t665 + 210 * t18 * eta * t1483 - 210 * t13 * t555 - 490 * t1431 * t1473 + 280 * t1431 * t665 + 210 * t18 * t1437 - 280 * t78 * t1437 - 35 * t1442 * t1483 - 245 * t1442 * t645 + 1470 * t1446 * t645 + 735 * t1456 * t1473 + 35 * t70 * t555
    t1499 = eta * t2
    t1514 = yuv3 * t1430
    t1517 = t21 * t1454
    t1526 = -1470 * t13 * t1451 + 1785 * t13 * t622 - 1785 * t18 * t1463 + 35 * t88 * t1499 - 420 * t91 * t1499 + 280 * t94 * t1499 - 980 * t1514 * t570 + 1960 * t1517 * t570 - 735 * t18 * t590 + 490 * t33 * t570 + 245 * t37 * t570 - 420 * t7 * t606
    t1529 = t88 * eta
    t1532 = t11 * t1445
    t1535 = t94 * eta
    t1556 = 735 * t13 * t1473 - 210 * t13 * t665 + 210 * t18 * t1483 + 1960 * t1514 * t634 - 980 * t1517 * t634 - 420 * t1529 * t606 + 3570 * t1532 * t606 - 420 * t1535 * t606 + 1470 * t18 * t645 + 420 * t23 * t606 - 245 * t33 * t634 - 490 * t37 * t634 + 280 * t70 * t665
    t1561 = eta * t55
    t1568 = zuv3 * t1454
    t1572 = eta * yuv3 * zuv3
    t1587 = -175 * t11 * t845 + 90 * t1572 * t1437 + 105 * t1568 * t1451 - 35 * t78 * t1483 + 280 * t88 * t1561 - 420 * t91 * t1561 + 35 * t94 * t1561 + 15 * t1568 * t555 - 105 * t1572 * t590 + 35 * t88 * t2 - 210 * t94 * t2 - 490 * t97 * t570
    t1605 = Q2 * t1
    t1617 = 175 * t13 * eta * t750 - 175 * t18 * eta * t1605 - 490 * t100 * t570 - 70 * t1431 * t750 - 420 * t1442 * t1605 + 75 * t1572 * t1463 + 105 * t1568 * t1473 - 75 * t1568 * t622 - 105 * t1572 * t645 - 175 * t88 * t606 - 350 * t91 * t606 - 175 * t94 * t606 - 490 * t97 * t634
    t1625 = Q2 * t778
    t1645 = 175 * t13 * eta * t830 - 490 * t100 * t634 - 175 * t11 * t914 - 245 * t1431 * t1625 + 420 * t1431 * t830 + 245 * t1442 * t790 + 735 * t1446 * t790 - 735 * t1456 * t1625 - 15 * t1572 * t1483 - 90 * t1568 * t665 - 210 * t88 * t55 + 35 * t94 * t55
    t1646 = Q2 * t27
    t1664 = yuv3 * t1454
    t1673 = t11 * eta
    t1676 = -175 * t18 * eta * t1646 - 180 * t11 * t1499 - 75 * t123 * t606 + 75 * t127 * t606 - 105 * t131 * t570 + 105 * t135 * t570 + 90 * t172 * t1437 + 70 * t1442 * t1646 - 150 * t1445 * t606 + 30 * t5 * t1499 - 420 * t1664 * t570 - 150 * t1673 * t606 + 15 * t169 * t555
    t1688 = eta * t1
    t1705 = 175 * t13 * t750 - 105 * t131 * t634 + 105 * t135 * t634 - 175 * t18 * t1605 - 420 * t78 * t1605 - 420 * t1664 * t634 - 70 * t88 * t1688 + 350 * t91 * t1688 + 420 * t94 * t1688 + 735 * t18 * t790 + 245 * t33 * t778 - 70 * t70 * t750
    t1730 = 30 * t11 * t1561 - 735 * t13 * t1625 + 175 * t13 * t830 - 15 * t172 * t1483 + 980 * t1514 * t778 + 980 * t1517 * t778 - 180 * t5 * t1561 - 175 * t18 * t1646 + 70 * t78 * t1646 - 90 * t169 * t665 - 245 * t37 * t778 + 420 * t70 * t830
    t1732 = eta * t27
    t1756 = -70 * t88 * t1 + 205 * t11 * t2 + 200 * t11 * t606 + 105 * t11 * t959 + 5 * t1568 * t750 + 205 * t1572 * t1605 + 420 * t88 * t1732 + 350 * t91 * t1732 - 70 * t94 * t1732 + 420 * t194 * t570 + 200 * t5 * t606 - 5 * t845
    t1782 = 175 * t94 * t1 + 490 * t100 * t778 - 5 * t11 * t55 + 105 * t11 * t982 + 210 * t1568 * t1625 - 205 * t1568 * t830 - 5 * t1572 * t1646 - 210 * t1572 * t790 + 420 * t194 * t634 + 175 * t88 * t27 - 70 * t94 * t27 + 490 * t97 * t778 + 205 * t914
    t1789 = Q1 * eta
    t1792 = Q2 * eta
    t1812 = 35 * zuv3 * eta * Q1 * t6 + 30 * eta * t606 - 410 * t11 * t1688 + 35 * t13 * t1789 - 210 * t131 * t778 + 210 * t135 * t778 + 205 * t172 * t1605 + 10 * t5 * t1688 + 5 * t169 * t750 - 35 * t18 * t1792 - 35 * t78 * t1792 + 15 * t1499
    t1834 = 35 * t169 * Q1 * t11 + 10 * t11 * t1732 - 5 * t172 * t1646 - 840 * t1664 * t778 - 205 * t169 * t830 - 35 * t572 * t172 - 410 * t5 * t1732 + 35 * t1529 + 70 * t1532 + 35 * t1535 + 15 * t1561 - 35 * t23 + 35 * t7
    t1853 = -165 * t11 * t1 + 45 * t11 * t27 - 20 * t169 * t1789 + 20 * t172 * t1792 - 420 * t194 * t778 - 20 * t2 - 20 * t55 - 40 * t606 + 35 * t88 + 70 * t91 + 45 * t959 - 165 * t982
    t1866 = 35 * t94 + 40 * t1688 + 40 * t1732 - 20 * t123 + 20 * t127 - 40 * t1445 - 40 * t1673 + 12 * t1 + 12 * t27 - 40 * t5 - 40 * t11 + 8 * eta + 8
    dMLdt = -3/16 * t302 / (1 + eta) * t293 * (t1468 + t1497 + t1526 + t1556 + t1587 + t1617 + t1645 + t1676 + t1705 + t1730 + t1756 + t1782 + t1812 + t1834 + t1853 + t1866)
    
    out = mu3*[dsmadt,dP1dt,dP2dt,dQ1dt,dQ2dt,dMLdt]
    return out;
        
end

function getaveragedterms_3bp_deg5(sma,n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3)

    dsmadt = 0.0;
    t1 = P1 ^ 2
    t2 = t1 ^ 2
    t3 = P1 * t2
    t4 = P2 * t3
    t5 = Q1 * t4
    t6 = xuv3 ^ 2
    t7 = xuv3 * t6
    t8 = yuv3 * t7
    t9 = zuv3 * t8
    t12 = yuv3 ^ 2
    t13 = yuv3 * t12
    t14 = t13 * xuv3
    t15 = zuv3 * t14
    t18 = t6 ^ 2
    t19 = t18 * Q2
    t20 = zuv3 * t19
    t23 = Q2 * t4
    t24 = t12 * t6
    t25 = zuv3 * t24
    t28 = t12 ^ 2
    t29 = t28 * Q2
    t30 = zuv3 * t29
    t33 = P2 ^ 2
    t34 = t33 * t2
    t36 = zuv3 * t18 * Q1
    t39 = Q1 * t34
    t43 = zuv3 * t28 * Q1
    t46 = Q2 * t34
    t51 = P1 * t1
    t52 = P2 * t33
    t53 = t52 * t51
    t54 = Q1 * t53
    t59 = -1344 * t15 * t46 - 168 * t15 * t5 + 1932 * t15 * t54 - 21 * t20 * t4 + 252 * t25 * t23 - 1008 * t25 * t39 - 168 * t30 * t4 + 105 * t36 * t34 + 336 * t43 * t34 + 672 * t9 * t46 + 84 * t9 * t5 - 1344 * t9 * t54
    t62 = Q2 * t53
    t67 = t33 ^ 2
    t68 = t67 * t1
    t71 = Q1 * t68
    t76 = Q2 * t68
    t81 = P2 * t67
    t82 = t81 * P1
    t83 = Q1 * t82
    t90 = Q2 * t82
    t95 = 1344 * t15 * t76 - 672 * t15 * t83 + 336 * t20 * t53 - 336 * t20 * t82 - 2898 * t25 * t62 + 2898 * t25 * t71 + 1008 * t25 * t90 + 420 * t30 * t53 - 105 * t30 * t82 - 420 * t36 * t68 - 336 * t43 * t68 - 1932 * t9 * t76 + 1344 * t9 * t83
    t97 = t33 * t67
    t98 = Q1 * t97
    t99 = zuv3 * t18
    t104 = zuv3 * t28
    t107 = Q2 * t97
    t112 = t1 * t2
    t113 = xuv3 * t18
    t116 = t7 * t112
    t119 = xuv3 * t112
    t122 = t18 * yuv3
    t125 = t6 * t13
    t128 = yuv3 * t28
    t133 = 21 * t104 * t98 - 84 * t15 * t107 + 168 * t9 * t107 + 21 * t113 * t112 - 231 * t113 * t34 - 336 * t12 * t116 + 336 * t28 * t119 - 672 * t122 * t4 + 1932 * t125 * t4 - 168 * t128 * t4 - 252 * t25 * t98 + 168 * t99 * t98
    t134 = t7 * t12
    t137 = t28 * xuv3
    t160 = t7 * t97
    t163 = -84 * t113 * t68 + 168 * t113 * t97 - 420 * t12 * t160 + 672 * t122 * t53 + 1344 * t122 * t82 + 588 * t125 * t53 - 1344 * t125 * t82 - 84 * t128 * t53 + 84 * t128 * t82 + 2562 * t134 * t34 + 2478 * t134 * t68 - 672 * t137 * t34 - 903 * t137 * t68
    t166 = xuv3 * t97
    t169 = yuv3 * xuv3
    t170 = zuv3 * t169
    t173 = t6 * Q2
    t174 = zuv3 * t173
    t177 = t12 * Q2
    t178 = zuv3 * t177
    t182 = zuv3 * t6 * Q1
    t186 = zuv3 * t12 * Q1
    t197 = P2 * t51
    t198 = Q1 * t197
    t203 = 84 * t15 * t198 + 105 * t28 * t166 + 224 * t170 * t46 + 28 * t170 * t5 - 196 * t170 * t54 - 14 * t174 * t4 + 98 * t174 * t53 + 84 * t178 * t4 + 42 * t178 * t53 + 42 * t182 * t34 - 112 * t186 * t34 - 168 * t9 * t198
    t206 = Q2 * t197
    t217 = t33 * t1
    t220 = Q1 * t217
    t225 = Q2 * t217
    t234 = -1344 * t15 * t225 + 196 * t170 * t76 - 224 * t170 * t83 + 112 * t174 * t82 - 42 * t182 * t68 - 98 * t186 * t68 + 42 * t20 * t197 - 420 * t30 * t197 - 126 * t25 * t206 - 210 * t36 * t217 + 336 * t43 * t217 + 882 * t25 * t220 - 588 * t9 * t225
    t238 = t52 * P1
    t239 = Q1 * t238
    t246 = Q2 * t238
    t251 = zuv3 * t6
    t254 = zuv3 * t12
    t259 = Q1 * t67
    t266 = Q2 * t67
    t269 = -42 * t104 * t259 - 28 * t170 * t107 + 588 * t15 * t239 - 42 * t178 * t82 - 336 * t20 * t238 + 210 * t30 * t238 + 1344 * t9 * t239 - 882 * t25 * t246 + 126 * t25 * t259 - 84 * t251 * t98 + 14 * t254 * t98 + 420 * t99 * t259 - 84 * t9 * t266
    t275 = t6 * yuv3
    t282 = xuv3 * t12
    t298 = 630 * t12 * t7 * t2 - 63 * t113 * t2 - 112 * t12 * t119 + 1344 * t122 * t197 - 1344 * t125 * t197 - 28 * t13 * t4 + 28 * t13 * t53 + 168 * t15 * t266 - 196 * t275 * t4 - 644 * t275 * t53 - 406 * t282 * t34 - 28 * t7 * t34 + 14 * t116
    t325 = 462 * t113 * t217 - 672 * t122 * t238 + 1932 * t125 * t238 + 84 * t128 * t197 - 168 * t128 * t238 + 56 * t13 * t82 - 2478 * t134 * t217 + 1218 * t137 * t217 - 448 * t275 * t82 - 224 * t282 * t68 - 182 * t7 * t68 - 140 * t160
    t330 = t7 * t67
    t333 = xuv3 * t67
    t336 = zuv3 * Q2
    t339 = zuv3 * Q1
    t356 = 84 * t113 * t67 + 70 * t12 * t166 + 546 * t12 * t330 + 28 * t170 * t198 - 14 * t174 * t197 + 294 * t178 * t197 + 42 * t182 * t217 - 322 * t186 * t217 - 231 * t28 * t333 - 5 * t336 * t4 - 10 * t336 * t53 + 5 * t339 * t34 + 10 * t339 * t68
    t368 = P2 * P1
    t369 = Q1 * t368
    t376 = Q2 * t368
    t385 = 84 * t15 * t369 + 644 * t170 * t225 - 644 * t170 * t239 + 322 * t174 * t238 - 42 * t178 * t238 - 21 * t20 * t368 - 126 * t25 * t376 - 294 * t251 * t259 - 105 * t30 * t368 - 5 * t336 * t82 + 84 * t9 * t369 + 5 * zuv3 * t98
    t390 = Q1 * t33
    t397 = Q2 * t33
    t407 = xuv3 * t2
    t414 = 21 * t104 * t390 - 210 * t12 * t407 - 84 * t15 * t397 - 28 * t170 * t266 - 448 * t275 * t197 + 126 * t25 * t390 + 14 * t254 * t259 + 35 * xuv3 * t34 + 105 * t99 * t390 - 84 * t9 * t397 + 20 * yuv3 * t4 + 40 * yuv3 * t53 + 5 * t119
    t427 = t7 * t1
    t430 = xuv3 * t1
    t443 = 63 * t113 * t1 - 252 * t12 * t427 - 672 * t122 * t368 - 588 * t125 * t368 + 56 * t13 * t197 - 28 * t13 * t238 + 14 * t282 * t217 - 238 * t7 * t217 - 196 * t275 * t238 - 315 * t28 * t430 + 55 * xuv3 * t68 + 20 * yuv3 * t82
    t452 = t7 * t33
    t455 = xuv3 * t33
    t470 = -231 * t113 * t33 - 28 * t12 * t333 - 84 * t12 * t452 + 84 * t128 * t368 - 56 * t170 * t369 + 28 * t174 * t368 + 84 * t178 * t368 - 20 * t336 * t197 + 20 * t339 * t217 - 20 * t336 * t238 + 147 * t28 * t455 + 25 * t166 - 154 * t330
    t495 = 294 * t12 * t430 - 28 * t13 * t368 + 56 * t170 * t397 + 20 * yuv3 * t197 + 50 * xuv3 * t217 + 20 * yuv3 * t238 - 84 * t251 * t390 - 28 * t254 * t390 + 20 * zuv3 * t259 + 644 * t275 * t368 + 35 * t333 + 15 * t407 - 42 * t427
    t513 = -70 * t12 * t455 - 8 * t336 * t368 - 40 * yuv3 * t368 + 8 * zuv3 * t390 - 21 * t113 - 42 * t134 - 21 * t137 + 28 * t282 - 12 * t430 + 266 * t452 - 52 * t455 + 28 * t7 - 8 * xuv3
    t518 = sma ^ 2
    t519 = sma * t518
    t521 = 0.1e1 / eta
    t522 = r3 ^ 2
    t523 = t522 ^ 2
    t525 = 0.1e1 / t522 / t523
    t527 = 0.1e1 / n
    t528 = t527 * t525 * t521
    dP1dt = 105/128 * t528 * t519 * (t59 + t95 + t133 + t163 + t203 + t234 + t269 + t298 + t325 + t356 + t385 + t414 + t443 + t470 + t495 + t513)
    t530 = Q1 * t112
    t535 = Q2 * t112
    t556 = -168 * t104 * t535 - 1344 * t15 * t23 + 1932 * t15 * t39 - 168 * t15 * t530 + 672 * t9 * t23 - 1008 * t25 * t5 + 252 * t25 * t535 + 105 * t36 * t4 - 1344 * t9 * t39 + 336 * t43 * t4 + 84 * t9 * t530 - 21 * t99 * t535
    t583 = 1344 * t15 * t62 - 672 * t15 * t71 + 336 * t20 * t34 - 336 * t20 * t68 - 2898 * t25 * t46 + 2898 * t25 * t54 + 1008 * t25 * t76 + 420 * t30 * t34 - 105 * t30 * t68 - 420 * t36 * t53 - 336 * t43 * t53 - 1932 * t9 * t62 + 1344 * t9 * t71
    t598 = t6 * t112
    t611 = 105 * yuv3 * t18 * t112 + 168 * t128 * t112 + 84 * t113 * t4 - 903 * t122 * t34 - 420 * t13 * t598 - 1344 * t134 * t4 + 1344 * t137 * t4 - 84 * t15 * t90 - 252 * t25 * t83 + 168 * t36 * t82 + 21 * t43 * t82 + 168 * t9 * t90
    t637 = t6 * t97
    t640 = 336 * yuv3 * t18 * t97 - 84 * t113 * t53 - 168 * t113 * t82 - 672 * t122 * t68 + 2478 * t125 * t34 + 2562 * t125 * t68 - 84 * t128 * t34 - 231 * t128 * t68 - 336 * t13 * t637 + 588 * t134 * t53 + 1932 * t134 * t82 + 672 * t137 * t53 - 672 * t137 * t82
    t663 = Q1 * t2
    t668 = 21 * t128 * t97 + 84 * t15 * t663 + 224 * t170 * t23 - 196 * t170 * t39 + 28 * t170 * t530 + 98 * t174 * t34 + 42 * t178 * t34 + 42 * t182 * t4 - 112 * t186 * t4 - 14 * t251 * t535 + 84 * t254 * t535 - 168 * t9 * t663
    t669 = Q2 * t2
    t696 = -420 * t104 * t669 - 1344 * t15 * t206 + 196 * t170 * t62 - 224 * t170 * t71 + 112 * t174 * t68 - 42 * t182 * t53 - 98 * t186 * t53 - 210 * t36 * t197 + 336 * t43 * t197 + 882 * t25 * t198 - 588 * t9 * t206 - 126 * t25 * t669 + 42 * t99 * t669
    t724 = 588 * t15 * t220 - 28 * t170 * t90 - 42 * t178 * t68 - 84 * t182 * t82 + 14 * t186 * t82 - 336 * t20 * t217 + 210 * t30 * t217 + 1344 * t9 * t220 - 882 * t25 * t225 + 420 * t36 * t238 - 42 * t43 * t238 + 126 * t25 * t239 - 84 * t9 * t246
    t742 = t6 * t2
    t753 = -231 * yuv3 * t18 * t2 - 140 * t13 * t112 - 168 * t113 * t197 + 84 * t128 * t2 - 182 * t13 * t34 + 546 * t13 * t742 + 168 * t15 * t246 - 224 * t275 * t34 - 448 * t282 * t4 - 644 * t282 * t53 + 56 * t7 * t4 + 28 * t7 * t53 + 70 * yuv3 * t598
    t781 = 84 * t113 * t238 + 1218 * t122 * t217 - 2478 * t125 * t217 + 462 * t128 * t217 - 28 * t13 * t68 + 1932 * t134 * t197 - 1344 * t134 * t238 - 672 * t137 * t197 + 1344 * t137 * t238 - 406 * t275 * t68 - 196 * t282 * t82 - 28 * t7 * t82
    t786 = t6 * t67
    t809 = -63 * t128 * t67 + 630 * t13 * t786 + 14 * t13 * t97 + 28 * t170 * t663 + 42 * t182 * t197 - 322 * t186 * t197 - 14 * t251 * t669 + 294 * t254 * t669 - 10 * t336 * t34 + 5 * t339 * t4 + 10 * t339 * t53 - 5 * zuv3 * t535 - 112 * yuv3 * t637
    t821 = Q1 * t1
    t826 = Q2 * t1
    t837 = -105 * t104 * t826 + 84 * t15 * t821 + 644 * t170 * t206 - 644 * t170 * t220 + 322 * t174 * t217 - 42 * t178 * t217 - 294 * t182 * t238 - 126 * t25 * t826 - 5 * t336 * t68 + 5 * t339 * t82 + 84 * t9 * t821 - 21 * t99 * t826
    t864 = 25 * yuv3 * t112 - 154 * t13 * t2 - 84 * t15 * t376 - 28 * t170 * t246 + 14 * t186 * t238 + 126 * t25 * t369 + 55 * yuv3 * t34 + 105 * t36 * t368 + 21 * t43 * t368 - 84 * t9 * t376 + 20 * xuv3 * t4 + 40 * xuv3 * t53 - 28 * yuv3 * t742
    t880 = t6 * t1
    t893 = 147 * yuv3 * t18 * t1 - 231 * t128 * t1 + 84 * t113 * t368 - 238 * t13 * t217 - 84 * t13 * t880 - 196 * t282 * t197 - 28 * t7 * t197 + 14 * t275 * t217 - 448 * t282 * t238 + 56 * t7 * t238 + 35 * yuv3 * t68 + 20 * xuv3 * t82
    t905 = t6 * t33
    t922 = -315 * yuv3 * t18 * t33 + 63 * t128 * t33 - 252 * t13 * t905 - 588 * t134 * t368 - 672 * t137 * t368 - 56 * t170 * t821 + 20 * t339 * t197 - 20 * t336 * t217 + 28 * t251 * t826 + 84 * t254 * t826 - 20 * zuv3 * t669 - 210 * yuv3 * t786 + 5 * yuv3 * t97
    t950 = 266 * t13 * t1 + 56 * t170 * t376 - 84 * t182 * t368 - 28 * t186 * t368 + 20 * xuv3 * t197 + 35 * yuv3 * t2 + 50 * yuv3 * t217 + 20 * t339 * t238 + 20 * xuv3 * t238 + 644 * t282 * t368 - 28 * t7 * t368 + 15 * yuv3 * t67 - 70 * yuv3 * t880
    t971 = -52 * yuv3 * t1 - 42 * t13 * t33 - 12 * yuv3 * t33 + 8 * t339 * t368 - 40 * xuv3 * t368 - 8 * zuv3 * t826 + 294 * yuv3 * t905 - 21 * t122 - 42 * t125 - 21 * t128 + 28 * t13 + 28 * t275 - 8 * yuv3
    dP2dt = -105/128 * t528 * t519 * (t556 + t583 + t611 + t640 + t668 + t696 + t724 + t753 + t781 + t809 + t837 + t864 + t893 + t922 + t950 + t971)
    t979 = Q1 ^ 2
    t980 = Q1 * t979
    t981 = t980 * t3
    t983 = yuv3 * t7 * Q2
    t986 = xuv3 * Q2
    t987 = t13 * t986
    t990 = t979 * t3
    t991 = Q2 ^ 2
    t992 = t18 * t991
    t995 = t6 * t991
    t996 = t12 * t995
    t999 = t28 * t991
    t1002 = Q1 * t3
    t1003 = Q2 * t991
    t1005 = yuv3 * t7 * t1003
    t1008 = xuv3 * t1003
    t1009 = t13 * t1008
    t1012 = t991 ^ 2
    t1013 = t1012 * t3
    t1020 = P2 * t2
    t1021 = Q2 * t980
    t1022 = t18 * t1021
    t1025 = 84 * t1005 * t1002 - 168 * t1009 * t1002 - 21 * t18 * t1013 + 252 * t24 * t1013 - 168 * t28 * t1013 + 105 * t1022 * t1020 + 84 * t983 * t981 - 168 * t987 * t981 - 21 * t992 * t990 + 252 * t996 * t990 - 168 * t999 * t990
    t1026 = t980 * t1020
    t1027 = t12 * t173
    t1030 = t28 * t1021
    t1033 = t979 * t1020
    t1035 = yuv3 * t7 * t991
    t1038 = xuv3 * t991
    t1039 = t13 * t1038
    t1042 = t1003 * Q1
    t1043 = t18 * t1042
    t1046 = Q1 * t1020
    t1047 = t6 * t1003
    t1048 = t12 * t1047
    t1051 = t28 * t1042
    t1055 = yuv3 * t7 * t1012
    t1058 = xuv3 * t1012
    t1059 = t13 * t1058
    t1062 = t33 * t51
    t1063 = t980 * t1062
    t1068 = t991 * t979
    t1069 = t18 * t1068
    t1072 = 336 * t1030 * t1020 + 105 * t1043 * t1020 + 336 * t1051 * t1020 + 672 * t1055 * t1020 - 1344 * t1059 * t1020 - 1008 * t1027 * t1026 + 672 * t1035 * t1033 - 1344 * t1039 * t1033 - 1008 * t1048 * t1046 + 336 * t1069 * t1062 - 1344 * t983 * t1063 + 1932 * t987 * t1063
    t1074 = t979 * t1062
    t1077 = t28 * t1068
    t1080 = Q1 * t1062
    t1085 = t18 * t1012
    t1088 = t6 * t1012
    t1089 = t12 * t1088
    t1092 = t28 * t1012
    t1095 = t52 * t1
    t1098 = t980 * t1095
    t1103 = t979 * t1095
    t1108 = -1344 * t1005 * t1080 + 1932 * t1009 * t1080 - 420 * t1022 * t1095 + 2898 * t1027 * t1098 - 336 * t1030 * t1095 - 1932 * t1035 * t1103 + 1344 * t1039 * t1103 + 420 * t1077 * t1062 + 336 * t1085 * t1062 - 2898 * t1089 * t1062 + 420 * t1092 * t1062 - 2898 * t996 * t1074
    t1111 = Q1 * t1095
    t1120 = t67 * P1
    t1121 = t980 * t1120
    t1128 = t979 * t1120
    t1133 = Q1 * t1120
    t1138 = 1344 * t1005 * t1133 - 672 * t1009 * t1133 - 420 * t1043 * t1095 + 2898 * t1048 * t1111 - 336 * t1051 * t1095 - 1932 * t1055 * t1095 + 1344 * t1059 * t1095 - 336 * t1069 * t1120 - 105 * t1077 * t1120 + 1344 * t983 * t1121 - 672 * t987 * t1121 + 1008 * t996 * t1128
    t1147 = t980 * t81
    t1154 = t979 * t81
    t1159 = Q1 * t81
    t1160 = t18 * t1003
    t1165 = t28 * t1003
    t1168 = -252 * t1027 * t1147 + 168 * t1035 * t1154 - 84 * t1039 * t1154 - 252 * t1048 * t1159 - 336 * t1085 * t1120 + 1008 * t1089 * t1120 - 105 * t1092 * t1120 + 168 * t19 * t1147 + 21 * t29 * t1147 + 168 * t1160 * t1159 + 21 * t1165 * t1159
    t1169 = t1012 * t81
    t1174 = yuv3 * t986
    t1179 = t12 * t991
    t1182 = yuv3 * t1008
    t1193 = t991 * t3
    t1198 = 28 * t1182 * t1002 + 168 * t983 * t1002 - 336 * t987 * t1002 + 84 * t12 * t1013 - 14 * t6 * t1013 - 84 * t14 * t1169 + 168 * t8 * t1169 + 28 * t1174 * t981 + 84 * t1179 * t990 - 42 * t18 * t1193 + 504 * t24 * t1193 - 14 * t995 * t990
    t1202 = t6 * t1021
    t1205 = t12 * t1021
    t1208 = yuv3 * t1038
    t1211 = t6 * t1042
    t1214 = t12 * t1042
    t1217 = Q2 * Q1
    t1218 = t18 * t1217
    t1223 = t28 * t1217
    t1226 = yuv3 * t1058
    t1233 = 1344 * t1035 * t1020 - 2688 * t1039 * t1020 + 42 * t1202 * t1020 - 112 * t1205 * t1020 + 42 * t1211 * t1020 - 112 * t1214 * t1020 + 210 * t1218 * t1020 + 672 * t1223 * t1020 + 224 * t1226 * t1020 - 2016 * t1027 * t1046 + 224 * t1208 * t1033 - 336 * t28 * t1193
    t1236 = t6 * t1068
    t1239 = t12 * t1068
    t1250 = t12 * t1012
    t1259 = t980 * t51
    t1262 = 98 * t1088 * t1062 + 98 * t1236 * t1062 + 42 * t1239 * t1062 + 42 * t1250 * t1062 + 672 * t992 * t1062 - 5796 * t996 * t1062 + 840 * t999 * t1062 - 196 * t1174 * t1063 - 196 * t1182 * t1080 - 2688 * t983 * t1080 + 3864 * t987 * t1080 - 168 * t983 * t1259
    t1268 = t979 * t51
    t1275 = Q1 * t51
    t1280 = t1012 * t51
    t1291 = -168 * t1005 * t1275 + 84 * t1009 * t1275 - 42 * t1202 * t1095 - 98 * t1205 * t1095 + 84 * t987 * t1259 + 42 * t992 * t1268 - 126 * t996 * t1268 - 420 * t999 * t1268 + 42 * t18 * t1280 - 126 * t24 * t1280 - 420 * t28 * t1280
    t1310 = P2 * t1
    t1313 = t980 * t1310
    t1318 = -210 * t1022 * t1310 + 5796 * t1027 * t1111 + 882 * t1027 * t1313 + 336 * t1030 * t1310 - 3864 * t1035 * t1095 + 2688 * t1039 * t1095 - 42 * t1211 * t1095 - 98 * t1214 * t1095 - 840 * t1218 * t1095 - 672 * t1223 * t1095 + 196 * t1226 * t1095 + 196 * t1208 * t1103
    t1320 = t979 * t1310
    t1327 = Q1 * t1310
    t1346 = -588 * t1035 * t1320 - 1344 * t1039 * t1320 - 210 * t1043 * t1310 + 882 * t1048 * t1327 + 336 * t1051 * t1310 - 588 * t1055 * t1310 - 1344 * t1059 * t1310 + 112 * t1236 * t1120 - 42 * t1239 * t1120 - 224 * t1174 * t1121 - 224 * t1182 * t1133 + 2688 * t983 * t1133
    t1359 = t33 * P1
    t1360 = t980 * t1359
    t1367 = t979 * t1359
    t1372 = Q1 * t1359
    t1375 = 1344 * t1005 * t1372 - 336 * t1069 * t1359 + 210 * t1077 * t1359 + 112 * t1088 * t1120 - 42 * t1250 * t1120 - 672 * t992 * t1120 + 2016 * t996 * t1120 - 210 * t999 * t1120 - 1344 * t987 * t1133 + 1344 * t983 * t1360 + 588 * t987 * t1360 - 882 * t996 * t1367
    t1394 = t12 * t1003
    t1403 = 588 * t1009 * t1372 - 504 * t1027 * t1159 - 84 * t1047 * t1159 - 336 * t1085 * t1359 - 882 * t1089 * t1359 + 210 * t1092 * t1359 - 84 * t173 * t1147 + 14 * t177 * t1147 - 28 * t1208 * t1154 + 14 * t1394 * t1159 + 336 * t19 * t1159 + 42 * t29 * t1159
    t1406 = t991 * t81
    t1411 = t980 * t52
    t1418 = t979 * t52
    t1423 = Q1 * t52
    t1430 = t1012 * t52
    t1433 = 126 * t1027 * t1411 - 84 * t1035 * t1418 + 168 * t1039 * t1418 + 126 * t1048 * t1423 + 420 * t1160 * t1423 - 42 * t1165 * t1423 - 28 * t169 * t1169 - 168 * t14 * t1406 + 336 * t8 * t1406 + 420 * t19 * t1411 - 42 * t29 * t1411 - 84 * t8 * t1430
    t1446 = t18 * t3
    t1448 = t6 * t3
    t1457 = t6 * t1217
    t1460 = 56 * t1174 * t1002 + 5 * t1021 * t1020 + 5 * t1042 * t1020 + 84 * t1457 * t1020 + 168 * t12 * t1193 - 28 * t6 * t1193 + 252 * t12 * t1448 + 168 * t14 * t1430 - 168 * t28 * t3 - 5 * t991 * t990 - 5 * t1013 - 21 * t1446
    t1461 = t12 * t1217
    t1486 = -10 * t1012 * t1062 + 448 * t1208 * t1020 - 1344 * t14 * t1020 - 224 * t1461 * t1020 + 672 * t8 * t1020 - 10 * t1068 * t1062 + 84 * t1179 * t1062 + 336 * t18 * t1062 - 2898 * t24 * t1062 + 420 * t28 * t1062 + 196 * t995 * t1062 - 392 * t1174 * t1080
    t1507 = t991 * t51
    t1514 = 28 * t1174 * t1259 + 294 * t1179 * t1268 + 28 * t1182 * t1275 + 294 * t12 * t1280 - 14 * t995 * t1268 - 336 * t983 * t1275 + 168 * t987 * t1275 - 14 * t6 * t1280 + 84 * t18 * t1507 - 252 * t24 * t1507 - 840 * t28 * t1507
    t1539 = 10 * t1021 * t1095 + 10 * t1042 * t1095 + 392 * t1208 * t1095 + 1344 * t14 * t1095 - 84 * t1457 * t1095 - 196 * t1461 * t1095 - 1932 * t8 * t1095 + 42 * t1202 * t1310 - 322 * t1205 * t1310 + 644 * t1208 * t1320 + 42 * t1211 * t1310 - 322 * t1214 * t1310
    t1565 = -5 * t1012 * t1120 + 1764 * t1027 * t1327 - 1176 * t1035 * t1310 - 2688 * t1039 * t1310 - 5 * t1068 * t1120 - 84 * t1179 * t1120 - 336 * t18 * t1120 + 224 * t995 * t1120 - 448 * t1174 * t1133 - 420 * t1218 * t1310 + 672 * t1223 * t1310 + 644 * t1226 * t1310
    t1590 = 322 * t1088 * t1359 + 1008 * t24 * t1120 - 105 * t28 * t1120 - 644 * t1174 * t1360 - 644 * t1182 * t1372 + 322 * t1236 * t1359 - 42 * t1239 * t1359 - 42 * t1250 * t1359 - 672 * t992 * t1359 - 1764 * t996 * t1359 + 2688 * t983 * t1372 + 1176 * t987 * t1372
    t1595 = t980 * P1
    t1600 = t979 * P1
    t1607 = P1 * Q1
    t1612 = t1012 * P1
    t1619 = 84 * t1005 * t1607 + 84 * t1009 * t1607 + 420 * t999 * t1359 + 84 * t983 * t1595 + 84 * t987 * t1595 - 21 * t992 * t1600 - 126 * t996 * t1600 - 105 * t999 * t1600 - 21 * t18 * t1612 - 126 * t24 * t1612 - 105 * t28 * t1612
    t1630 = t7 * t81
    t1633 = xuv3 * t81
    t1646 = 5 * Q2 * t1147 + 5 * t1003 * t1159 - 294 * t1047 * t1423 - 168 * t173 * t1159 + 28 * t177 * t1159 - 28 * t1208 * t1418 - 84 * t13 * t1633 + 14 * t1394 * t1423 - 56 * t169 * t1406 - 294 * t173 * t1411 + 14 * t177 * t1411 + 168 * yuv3 * t1630
    t1656 = t991 * t52
    t1661 = t980 * P2
    t1668 = t979 * P2
    t1673 = P2 * Q1
    t1676 = 252 * t1027 * t1423 + 126 * t1027 * t1661 - 84 * t1035 * t1668 - 84 * t1039 * t1668 + 105 * t1160 * t1673 + 336 * t14 * t1656 + 840 * t19 * t1423 - 84 * t29 * t1423 - 28 * t169 * t1430 - 168 * t8 * t1656 + 105 * t19 * t1661 + 21 * t29 * t1661
    t1681 = t1012 * P2
    t1700 = 10 * t1217 * t1020 + 224 * t169 * t1020 + 126 * t1048 * t1673 + 42 * t12 * t1062 + 98 * t6 * t1062 - 20 * t991 * t1062 + 21 * t1165 * t1673 + 84 * t12 * t3 - 84 * t14 * t1681 - 84 * t8 * t1681 - 10 * t1193 - 14 * t1448
    t1713 = t18 * t51
    t1715 = t6 * t51
    t1726 = 20 * t1021 * t1310 + 20 * t1217 * t1095 + 196 * t169 * t1095 + 56 * t1174 * t1275 + 588 * t12 * t1507 - 126 * t12 * t1715 - 20 * t991 * t1268 - 28 * t6 * t1507 - 420 * t28 * t51 - 20 * t1280 + 42 * t1713
    t1751 = -20 * t1012 * t1359 + 20 * t1042 * t1310 - 20 * t1068 * t1359 - 42 * t12 * t1120 + 112 * t6 * t1120 - 10 * t991 * t1120 - 1288 * t1174 * t1372 + 1288 * t1208 * t1310 - 1344 * t14 * t1310 + 84 * t1457 * t1310 - 644 * t1461 * t1310 - 588 * t8 * t1310
    t1777 = -56 * t1174 * t1595 - 84 * t1179 * t1359 + 84 * t1179 * t1600 - 56 * t1182 * t1607 - 336 * t18 * t1359 - 882 * t24 * t1359 + 210 * t28 * t1359 + 644 * t995 * t1359 + 28 * t995 * t1600 + 168 * t983 * t1607 + 168 * t987 * t1607 + 28 * t6 * t1612
    t1780 = t991 * P1
    t1801 = t7 * t52
    t1804 = 10 * Q2 * t1159 + 20 * Q2 * t1411 + 20 * t1003 * t1423 + 84 * t12 * t1612 - 588 * t173 * t1423 + 28 * t177 * t1423 - 28 * yuv3 * t1633 - 56 * t169 * t1656 - 42 * t18 * t1780 - 252 * t24 * t1780 - 210 * t28 * t1780 - 84 * yuv3 * t1801
    t1807 = xuv3 * t52
    t1828 = t991 * P2
    t1833 = 252 * t1027 * t1673 - 84 * t1047 * t1673 + 56 * t1208 * t1668 + 168 * t13 * t1807 - 28 * t1394 * t1673 - 168 * t14 * t1828 - 84 * t173 * t1661 - 28 * t177 * t1661 + 210 * t19 * t1673 + 42 * t29 * t1673 + 56 * t169 * t1681 - 168 * t8 * t1828
    t1853 = -42 * t12 * t1359 + 294 * t12 * t51 + 40 * t1217 * t1310 + 644 * t169 * t1310 + 322 * t6 * t1359 - 40 * t991 * t1359 - 8 * t991 * t1600 - 10 * t1062 - 5 * t1120 - 40 * t1507 - 14 * t1715 - 5 * t3
    t1862 = t18 * P1
    t1864 = t6 * P1
    t1879 = -105 * t28 * P1 + 40 * Q2 * t1423 + 8 * Q2 * t1661 + 8 * t1003 * t1673 - 112 * t1174 * t1607 + 168 * t12 * t1780 - 126 * t12 * t1864 - 168 * t173 * t1673 + 56 * t6 * t1780 - 28 * yuv3 * t1807 - 8 * t1612 - 21 * t1862
    t1884 = t7 * P2
    t1887 = xuv3 * P2
    t1901 = 84 * t12 * P1 + 16 * Q2 * t1673 - 84 * t13 * t1887 - 56 * t177 * t1673 + 112 * t169 * t1828 - 84 * yuv3 * t1884 + 56 * yuv3 * t1887 - 8 * P1 - 20 * t1359 - 16 * t1780 + 28 * t1864 - 20 * t51
    t1912 = t525 * t527 / GG * t521
    dQ1dt = 105/256 * t1912 * t519 * (t1346 + t1318 + t1291 + t1262 + t1233 + t1198 + t1168 + t1108 + t1853 + t1833 + t1804 + t1777 + t1751 + t1726 + t1700 + t1676 + t1646 + t1619 + t1590 + t1565 + t1539 + t1514 + t1486 + t1460 + t1138 + t1433 + t1403 + t1375 + t1072 + t1901 + t1879 + t1025) * zuv3
    t1914 = t979 ^ 2
    t1915 = t1914 * t3
    t1936 = t18 * t1914
    t1939 = 252 * t1048 * t1002 - 21 * t1160 * t1002 - 168 * t1165 * t1002 + 105 * t1936 * t1020 + 252 * t1027 * t981 + 84 * t1035 * t990 - 168 * t1039 * t990 - 168 * t14 * t1915 - 21 * t19 * t981 + 84 * t8 * t1915 - 168 * t29 * t981
    t1940 = t6 * t1914
    t1941 = t12 * t1940
    t1944 = t28 * t1914
    t1962 = yuv3 * t7 * t1914
    t1965 = xuv3 * t1914
    t1966 = t13 * t1965
    t1971 = 672 * t1005 * t1046 - 1344 * t1009 * t1046 + 105 * t1069 * t1020 + 336 * t1077 * t1020 - 1008 * t1941 * t1020 + 336 * t1944 * t1020 + 336 * t1022 * t1062 + 672 * t983 * t1026 - 1344 * t987 * t1026 - 1008 * t996 * t1033 - 1344 * t1962 * t1062 + 1932 * t1966 * t1062
    t1997 = -2898 * t1027 * t1063 + 420 * t1030 * t1062 - 1344 * t1035 * t1074 + 1932 * t1039 * t1074 + 336 * t1043 * t1062 - 2898 * t1048 * t1080 + 420 * t1051 * t1062 - 420 * t1936 * t1095 + 2898 * t1941 * t1095 - 336 * t1944 * t1095 - 1932 * t983 * t1098 + 1344 * t987 * t1098
    t2022 = -1932 * t1005 * t1111 + 1344 * t1009 * t1111 - 336 * t1022 * t1120 + 1008 * t1027 * t1121 - 105 * t1030 * t1120 + 1344 * t1035 * t1128 - 672 * t1039 * t1128 - 420 * t1069 * t1095 - 336 * t1077 * t1095 + 2898 * t996 * t1103 + 1344 * t1962 * t1120 - 672 * t1966 * t1120
    t2031 = t1914 * t81
    t2048 = -336 * t1043 * t1120 + 1008 * t1048 * t1133 - 105 * t1051 * t1120 + 168 * t983 * t1147 - 84 * t987 * t1147 + 168 * t992 * t1154 - 252 * t996 * t1154 + 21 * t999 * t1154 + 168 * t18 * t2031 - 252 * t24 * t2031 + 21 * t28 * t2031
    t2073 = 504 * t1027 * t1002 - 14 * t1047 * t1002 + 84 * t1394 * t1002 - 42 * t19 * t1002 + 168 * t1005 * t1159 - 84 * t1009 * t1159 + 28 * t1208 * t990 - 336 * t14 * t990 + 28 * t169 * t1915 - 14 * t173 * t981 + 84 * t177 * t981 + 168 * t8 * t990
    t2079 = t12 * t1914
    t2088 = t18 * t979
    t2091 = t6 * t979
    t2092 = t12 * t2091
    t2095 = t28 * t979
    t2104 = -336 * t29 * t1002 + 42 * t1236 * t1020 - 112 * t1239 * t1020 + 42 * t1940 * t1020 - 112 * t2079 * t1020 + 210 * t2088 * t1020 - 2016 * t2092 * t1020 + 672 * t2095 * t1020 + 224 * t1174 * t1026 + 224 * t1182 * t1046 + 1344 * t983 * t1046 - 2688 * t987 * t1046
    t2105 = yuv3 * t1965
    t2115 = yuv3 * t7 * t979
    t2118 = xuv3 * t979
    t2119 = t13 * t2118
    t2132 = t1914 * t51
    t2135 = -5796 * t1027 * t1080 + 98 * t1202 * t1062 + 42 * t1205 * t1062 + 98 * t1211 * t1062 + 42 * t1214 * t1062 + 672 * t1218 * t1062 + 840 * t1223 * t1062 - 196 * t2105 * t1062 - 2688 * t2115 * t1062 + 3864 * t2119 * t1062 - 196 * t1208 * t1074 - 168 * t8 * t2132
    t2161 = -126 * t1027 * t1259 - 168 * t1035 * t1268 + 84 * t1039 * t1268 - 126 * t1048 * t1275 - 42 * t1940 * t1095 - 98 * t2079 * t1095 + 42 * t1160 * t1275 - 420 * t1165 * t1275 + 42 * t19 * t1259 - 420 * t29 * t1259 + 84 * t14 * t2132
    t2186 = -42 * t1236 * t1095 - 98 * t1239 * t1095 - 840 * t2088 * t1095 + 5796 * t2092 * t1095 - 672 * t2095 * t1095 + 196 * t1174 * t1098 + 196 * t1182 * t1111 - 3864 * t983 * t1111 + 2688 * t987 * t1111 - 210 * t1936 * t1310 + 882 * t1941 * t1310 + 336 * t1944 * t1310
    t2212 = -588 * t1005 * t1327 - 1344 * t1009 * t1327 - 210 * t1069 * t1310 + 336 * t1077 * t1310 + 112 * t1202 * t1120 - 42 * t1205 * t1120 - 224 * t2105 * t1120 + 2688 * t2115 * t1120 - 224 * t1208 * t1128 - 588 * t983 * t1313 - 1344 * t987 * t1313 + 882 * t996 * t1320
    t2237 = -336 * t1022 * t1359 + 2016 * t1027 * t1133 - 882 * t1027 * t1360 + 210 * t1030 * t1359 + 1344 * t1035 * t1367 + 112 * t1211 * t1120 - 42 * t1214 * t1120 - 672 * t1218 * t1120 - 210 * t1223 * t1120 - 1344 * t2119 * t1120 + 1344 * t1962 * t1359 + 588 * t1966 * t1359
    t2264 = 588 * t1039 * t1367 - 336 * t1043 * t1359 - 882 * t1048 * t1372 + 210 * t1051 * t1359 - 28 * t1174 * t1147 + 14 * t1179 * t1154 + 336 * t18 * t1154 - 504 * t24 * t1154 + 42 * t28 * t1154 - 84 * t995 * t1154 + 14 * t12 * t2031 - 84 * t6 * t2031
    t2271 = t1914 * t52
    t2290 = -84 * t1005 * t1423 - 28 * t1182 * t1159 + 336 * t983 * t1159 - 168 * t987 * t1159 - 84 * t983 * t1411 + 168 * t987 * t1411 + 420 * t992 * t1418 + 126 * t996 * t1418 - 42 * t999 * t1418 + 420 * t18 * t2271 + 126 * t24 * t2271 - 42 * t28 * t2271
    t2307 = xuv3 * t3
    t2316 = t12 * t979
    t2319 = 84 * yuv3 * t7 * t3 - 5 * Q2 * t981 - 5 * t1003 * t1002 - 28 * t173 * t1002 + 168 * t177 * t1002 + 168 * t1009 * t1423 + 5 * t1068 * t1020 + 5 * t1914 * t1020 + 84 * t2091 * t1020 - 224 * t2316 * t1020 - 168 * t13 * t2307 + 56 * t169 * t990
    t2330 = yuv3 * t2118
    t2345 = 105 * t18 * t1020 - 1008 * t24 * t1020 + 336 * t28 * t1020 - 10 * t1021 * t1062 - 10 * t1042 * t1062 + 448 * t1174 * t1046 + 1932 * t14 * t1062 + 196 * t1457 * t1062 + 84 * t1461 * t1062 - 392 * t2330 * t1062 - 1344 * t8 * t1062 + 28 * t169 * t2132
    t2372 = -252 * t1027 * t1275 - 14 * t1047 * t1275 + 10 * t1914 * t1095 + 28 * t1208 * t1268 - 14 * t173 * t1259 + 294 * t177 * t1259 + 168 * t14 * t1268 - 336 * t8 * t1268 + 294 * t1394 * t1275 + 84 * t19 * t1275 - 840 * t29 * t1275
    t2397 = 10 * t1068 * t1095 - 420 * t18 * t1095 - 84 * t2091 * t1095 - 196 * t2316 * t1095 + 2898 * t24 * t1095 - 336 * t28 * t1095 + 392 * t1174 * t1111 + 644 * t1174 * t1313 + 42 * t1236 * t1310 - 322 * t1239 * t1310 + 42 * t1940 * t1310 - 322 * t2079 * t1310
    t2423 = -5 * t1021 * t1120 - 5 * t1042 * t1120 + 224 * t1457 * t1120 - 84 * t1461 * t1120 - 448 * t2330 * t1120 + 1344 * t8 * t1120 + 644 * t1182 * t1327 - 420 * t2088 * t1310 + 1764 * t2092 * t1310 + 672 * t2095 * t1310 - 1176 * t983 * t1327 - 2688 * t987 * t1327
    t2448 = -1764 * t1027 * t1372 - 672 * t14 * t1120 + 322 * t1202 * t1359 - 42 * t1205 * t1359 - 644 * t1208 * t1367 + 322 * t1211 * t1359 - 42 * t1214 * t1359 - 672 * t1218 * t1359 + 420 * t1223 * t1359 - 644 * t2105 * t1359 + 2688 * t2115 * t1359 + 1176 * t2119 * t1359
    t2451 = t1914 * P1
    t2473 = -126 * t1027 * t1595 + 84 * t1035 * t1600 + 84 * t1039 * t1600 - 126 * t1048 * t1607 - 21 * t1160 * t1607 - 105 * t1165 * t1607 + 84 * t14 * t2451 - 21 * t19 * t1595 - 105 * t29 * t1595 + 84 * t8 * t2451 + 5 * t2031
    t2484 = t6 * t81
    t2499 = 28 * t12 * t1154 - 168 * t6 * t1154 + 5 * t991 * t1154 - 56 * t1174 * t1159 - 28 * t1174 * t1411 + 14 * t1179 * t1418 + 14 * t12 * t2271 - 252 * t12 * t2484 - 294 * t995 * t1418 + 168 * t18 * t81 - 294 * t6 * t2271 + 21 * t28 * t81
    t2513 = t1914 * P2
    t2526 = -28 * t1182 * t1423 + 840 * t18 * t1418 + 252 * t24 * t1418 - 84 * t28 * t1418 - 168 * t983 * t1423 + 336 * t987 * t1423 - 84 * t983 * t1661 - 84 * t987 * t1661 + 105 * t992 * t1668 + 105 * t18 * t2513 + 126 * t24 * t2513 + 21 * t28 * t2513
    t2550 = -10 * Q2 * t1002 - 20 * Q2 * t1259 - 84 * t1005 * t1673 - 84 * t1009 * t1673 - 112 * t12 * t1020 + 42 * t6 * t1020 - 20 * t1217 * t1062 - 196 * t169 * t1062 + 126 * t996 * t1668 + 21 * t999 * t1668 + 28 * yuv3 * t2307 + 10 * t1033
    t2565 = xuv3 * t51
    t2577 = -168 * yuv3 * t7 * t51 - 20 * t1003 * t1275 + 20 * t1068 * t1310 - 98 * t12 * t1095 - 42 * t6 * t1095 + 56 * t169 * t1268 - 28 * t173 * t1275 + 588 * t177 * t1275 + 84 * t13 * t2565 + 20 * t1914 * t1310 + 20 * t1103
    t2602 = -20 * t1021 * t1359 - 20 * t1042 * t1359 - 10 * t1217 * t1120 - 224 * t169 * t1120 + 1288 * t1174 * t1327 - 210 * t18 * t1310 + 84 * t2091 * t1310 - 644 * t2316 * t1310 + 882 * t24 * t1310 + 336 * t28 * t1310 + 644 * t1457 * t1359 - 1288 * t2330 * t1359
    t2628 = 28 * t1047 * t1607 - 56 * t1208 * t1600 + 588 * t14 * t1359 - 84 * t1461 * t1359 + 1344 * t8 * t1359 + 84 * t1394 * t1607 + 168 * t14 * t1600 + 28 * t173 * t1595 + 84 * t177 * t1595 + 168 * t8 * t1600 - 42 * t19 * t1607 - 56 * t169 * t2451
    t2648 = t6 * t52
    t2651 = -252 * t1027 * t1607 - 56 * t1174 * t1423 + 28 * t12 * t1418 + 126 * t12 * t2648 + 14 * t12 * t81 - 588 * t6 * t1418 + 20 * t991 * t1418 - 210 * t29 * t1607 + 420 * t18 * t52 + 10 * t1154 + 20 * t2271 - 84 * t2484
    t2678 = 56 * t1174 * t1661 - 28 * t1179 * t1668 + 56 * t1182 * t1673 - 28 * t12 * t2513 + 210 * t18 * t1668 + 252 * t24 * t1668 + 42 * t28 * t1668 - 84 * t995 * t1668 - 168 * t983 * t1673 - 168 * t987 * t1673 - 84 * t6 * t2513 - 42 * t28 * t52
    t2700 = -40 * Q2 * t1275 - 8 * Q2 * t1595 - 8 * t1003 * t1607 - 322 * t12 * t1310 - 40 * t1217 * t1359 + 42 * t6 * t1310 - 644 * t169 * t1359 - 112 * t169 * t1600 + 28 * yuv3 * t2565 + 5 * t1020 + 10 * t1095 + 40 * t1320
    t2709 = xuv3 * P1
    t2724 = 84 * yuv3 * t7 * P1 - 56 * t12 * t1668 + 14 * t12 * t52 + 84 * t13 * t2709 + 56 * t173 * t1607 + 168 * t177 * t1607 - 168 * t6 * t1668 + 8 * t991 * t1668 + 40 * t1418 + 8 * t2513 - 294 * t2648 + 5 * t81
    t2729 = t6 * P2
    t2745 = -28 * t12 * P2 + 105 * t18 * P2 + 21 * t28 * P2 - 16 * Q2 * t1607 + 112 * t1174 * t1673 + 126 * t12 * t2729 - 56 * yuv3 * t2709 + 8 * P2 + 20 * t1310 + 16 * t1668 - 84 * t2729 + 20 * t52
    dQ2dt = -105/256 * t1912 * t519 * zuv3 * (t2161 + t2550 + t2745 + t2700 + t1997 + t2651 + t2526 + t2473 + t2264 + t2448 + t2678 + t2602 + t2186 + t2212 + t2345 + t2577 + t2135 + t2397 + t1971 + t1939 + t2290 + t2319 + t2237 + t2048 + t2499 + t2724 + t2022 + t2073 + t2104 + t2372 + t2423 + t2628)
    t2755 = eta * t1002
    t2760 = Q2 * t3
    t2761 = t18 * eta
    t2762 = zuv3 * t2761
    t2769 = zuv3 * t28 * eta
    t2774 = t6 * eta
    t2775 = t254 * t2774
    t2780 = Q2 * t1020
    t2781 = t7 * eta
    t2782 = yuv3 * zuv3
    t2783 = t2782 * t2781
    t2786 = eta * xuv3
    t2788 = zuv3 * t13 * t2786
    t2793 = 252 * t25 * eta * t2760 + 105 * t2762 * t1046 + 336 * t2769 * t1046 - 1008 * t2775 * t1046 - 1344 * t2783 * t1080 - 168 * t15 * t2755 + 84 * t9 * t2755 - 21 * t2762 * t2760 - 168 * t2769 * t2760 + 672 * t2783 * t2780 - 1344 * t2788 * t2780
    t2796 = Q2 * t1062
    t2809 = Q2 * t1095
    t2818 = 1932 * t2788 * t1080 - 420 * t2762 * t1111 - 336 * t2769 * t1111 + 2898 * t2775 * t1111 + 1344 * t2783 * t1133 - 672 * t2788 * t1133 + 336 * t2762 * t2796 + 420 * t2769 * t2796 - 2898 * t2775 * t2796 - 1932 * t2783 * t2809 + 1344 * t2788 * t2809
    t2820 = Q2 * t1120
    t2834 = Q2 * t81
    t2835 = eta * t2834
    t2846 = -252 * t25 * eta * t1159 - 168 * t15 * t1002 + 84 * t9 * t1002 + 168 * t2762 * t1159 + 21 * t2769 * t1159 - 84 * t15 * t2835 - 21 * t99 * t2760 - 336 * t2762 * t2820 - 105 * t2769 * t2820 + 1008 * t2775 * t2820 + 168 * t9 * t2835
    t2851 = eta * t3
    t2868 = t113 * eta
    t2871 = t12 * t2781
    t2874 = 105 * t2868 * t1020 - 1680 * t2871 * t1020 + 105 * t36 * t1020 + 336 * t43 * t1020 - 168 * t104 * t2760 - 1008 * t25 * t1046 + 105 * t122 * t2851 - 420 * t125 * t2851 + 168 * t128 * t2851 - 1344 * t15 * t2780 + 252 * t25 * t2760 + 672 * t9 * t2780
    t2877 = t28 * t2786
    t2890 = yuv3 * t2761
    t2893 = t13 * t2774
    t2896 = t128 * eta
    t2903 = 1680 * t2877 * t1020 + 336 * t20 * t1062 - 1680 * t2890 * t1062 + 4830 * t2893 * t1062 - 420 * t2896 * t1062 + 420 * t30 * t1062 + 1932 * t15 * t1080 - 1344 * t9 * t1080 - 420 * t36 * t1095 + 2898 * t25 * t1111 - 2898 * t25 * t2796
    t2926 = -420 * t2868 * t1095 + 4830 * t2871 * t1095 - 1680 * t2877 * t1095 - 336 * t43 * t1095 - 336 * t20 * t1120 - 105 * t30 * t1120 - 672 * t15 * t1133 + 1344 * t9 * t1133 + 1344 * t15 * t2809 + 1008 * t25 * t2820 - 1932 * t9 * t2809
    t2944 = eta * t81
    t2951 = 21 * t104 * t1159 + 1680 * t2890 * t1120 - 1680 * t2893 * t1120 + 105 * t2896 * t1120 + 168 * t113 * t2944 - 252 * t25 * t1159 + 168 * t99 * t1159 - 420 * t134 * t2944 + 105 * t137 * t2944 - 84 * t15 * t2834 + 168 * t9 * t2834
    t2954 = zuv3 * t2774
    t2958 = zuv3 * t12 * eta
    t2971 = t2782 * t2786
    t2980 = 84 * t113 * t1020 - 588 * t134 * t1020 - 672 * t137 * t1020 + 42 * t2954 * t1046 - 112 * t2958 * t1046 - 168 * t128 * t3 - 84 * t13 * t1448 + 84 * yuv3 * t1446 + 28 * t170 * t2755 - 14 * t2954 * t2760 + 84 * t2958 * t2760 + 224 * t2971 * t2780
    t2996 = eta * t1275
    t3001 = Q2 * t51
    t3009 = -126 * t25 * eta * t3001 - 588 * t122 * t1062 - 672 * t125 * t1062 - 84 * t128 * t1062 - 196 * t2971 * t1080 + 84 * t15 * t2996 + 42 * t2762 * t3001 - 420 * t2769 * t3001 + 98 * t2954 * t2796 + 42 * t2958 * t2796 - 168 * t9 * t2996
    t3028 = Q2 * t1310
    t3033 = -84 * t113 * t1095 - 672 * t134 * t1095 - 588 * t137 * t1095 - 42 * t2954 * t1111 - 98 * t2958 * t1111 - 210 * t2762 * t1327 + 336 * t2769 * t1327 + 882 * t2775 * t1327 - 588 * t2783 * t3028 - 1344 * t2788 * t3028 + 196 * t2971 * t2809
    t3051 = Q2 * t1359
    t3058 = -672 * t122 * t1120 - 588 * t125 * t1120 + 84 * t128 * t1120 - 224 * t2971 * t1133 + 1344 * t2783 * t1372 + 588 * t2788 * t1372 - 336 * t2762 * t3051 + 210 * t2769 * t3051 - 882 * t2775 * t3051 + 112 * t2954 * t2820 - 42 * t2958 * t2820
    t3078 = Q2 * t52
    t3079 = eta * t3078
    t3086 = 126 * t25 * eta * t1423 + 28 * t170 * t1002 - 168 * t113 * t81 - 84 * t2954 * t1159 + 14 * t2958 * t1159 - 84 * t12 * t1630 + 420 * t2762 * t1423 - 42 * t2769 * t1423 + 168 * t15 * t3079 + 84 * t28 * t1633 - 28 * t170 * t2835 - 84 * t9 * t3079
    t3105 = t12 * t2786
    t3112 = 42 * t182 * t1020 - 112 * t186 * t1020 + 70 * t2781 * t1020 - 560 * t3105 * t1020 + 98 * t174 * t1062 - 196 * t170 * t1080 - 140 * t13 * t2851 + 224 * t170 * t2780 - 14 * t251 * t2760 + 84 * t254 * t2760 + 70 * t275 * t2851
    t3115 = yuv3 * t2774
    t3118 = t13 * eta
    t3131 = eta * t51
    t3138 = -420 * t104 * t3001 + 42 * t178 * t1062 - 490 * t3115 * t1062 - 70 * t3118 * t1062 - 210 * t122 * t3131 + 210 * t125 * t3131 + 84 * t15 * t1275 - 168 * t9 * t1275 + 420 * t128 * t3131 - 126 * t25 * t3001 + 42 * t99 * t3001
    t3162 = -42 * t182 * t1095 - 98 * t186 * t1095 - 70 * t2781 * t1095 - 490 * t3105 * t1095 - 210 * t2868 * t1310 - 210 * t36 * t1310 + 336 * t43 * t1310 + 882 * t25 * t1327 - 1344 * t15 * t3028 + 196 * t170 * t2809 - 588 * t9 * t3028
    t3187 = 112 * t174 * t1120 - 42 * t178 * t1120 - 560 * t3115 * t1120 + 70 * t3118 * t1120 - 224 * t170 * t1133 + 1470 * t2871 * t1310 + 1680 * t2877 * t1310 - 336 * t20 * t1359 + 210 * t30 * t1359 + 588 * t15 * t1372 + 1344 * t9 * t1372 - 882 * t25 * t3051
    t3214 = -42 * t104 * t1423 - 84 * t251 * t1159 + 14 * t254 * t1159 + 1680 * t2890 * t1359 + 1470 * t2893 * t1359 - 210 * t2896 * t1359 + 126 * t25 * t1423 + 420 * t99 * t1423 - 28 * t170 * t2834 + 70 * t282 * t2944 - 140 * t7 * t2944
    t3219 = eta * t52
    t3226 = eta * zuv3
    t3234 = Q1 * eta * zuv3
    t3241 = 644 * t282 * t1020 + 5 * t3234 * t1020 - 28 * t7 * t1020 + 420 * t113 * t3219 + 196 * t13 * t3 + 210 * t134 * t3219 - 210 * t137 * t3219 - 28 * yuv3 * t1448 + 168 * t15 * t3078 - 5 * t3226 * t2760 - 84 * t9 * t3078
    t3244 = Q2 * eta * zuv3
    t3267 = 168 * t13 * t1062 + 616 * t275 * t1062 - 10 * t3244 * t1062 + 10 * t3234 * t1095 + 168 * t7 * t1095 + 84 * t128 * t51 - 84 * t13 * t1715 + 28 * t170 * t2996 - 168 * yuv3 * t1713 - 14 * t2954 * t3001 + 294 * t2958 * t3001
    t3292 = 616 * t282 * t1095 - 28 * t13 * t1120 + 644 * t275 * t1120 - 5 * t3244 * t1120 - 168 * t113 * t1310 + 420 * t134 * t1310 + 588 * t137 * t1310 + 42 * t2954 * t1327 - 322 * t2958 * t1327 - 644 * t2971 * t1372 + 322 * t2954 * t3051 + 644 * t2971 * t3028
    t3303 = eta * t1607
    t3308 = P1 * Q2
    t3319 = -126 * t25 * eta * t3308 + 5 * t3226 * t1159 + 588 * t122 * t1359 + 420 * t125 * t1359 - 168 * t128 * t1359 + 84 * t15 * t3303 - 21 * t2762 * t3308 - 105 * t2769 * t3308 - 42 * t2958 * t3051 + 84 * t9 * t3303 + 196 * t1630
    t3341 = Q2 * P2
    t3342 = eta * t3341
    t3345 = 126 * t25 * eta * t1673 + 84 * t113 * t52 - 28 * t12 * t1633 - 84 * t12 * t1801 - 294 * t2954 * t1423 + 14 * t2958 * t1423 + 105 * t2762 * t1673 + 21 * t2769 * t1673 - 28 * t170 * t3079 - 168 * t28 * t1807 - 84 * t9 * t3342
    t3359 = eta * yuv3
    t3370 = 25 * t2786 * t1020 + 5 * t339 * t1020 + 50 * t3359 * t1062 - 10 * t336 * t1062 + 28 * t170 * t1275 - 84 * t15 * t3342 - 14 * t251 * t3001 + 294 * t254 * t3001 + 70 * t275 * t3131 - 5 * zuv3 * t2760 + 25 * yuv3 * t2851
    t3395 = 50 * t2786 * t1095 + 10 * t339 * t1095 + 25 * t3359 * t1120 - 5 * t336 * t1120 - 490 * t13 * t3131 + 42 * t182 * t1310 - 322 * t186 * t1310 + 70 * t2781 * t1310 - 1610 * t3105 * t1310 + 322 * t174 * t1359 - 644 * t170 * t1372 + 644 * t170 * t3028
    t3415 = P1 * eta
    t3422 = -105 * t104 * t3308 + 105 * t122 * t3415 + 210 * t125 * t3415 + 105 * t128 * t3415 - 42 * t178 * t1359 - 1610 * t3115 * t1359 + 70 * t3118 * t1359 + 84 * t15 * t1607 + 84 * t9 * t1607 - 126 * t25 * t3308 - 21 * t99 * t3308
    t3445 = 21 * t104 * t1673 + 5 * zuv3 * t1159 - 294 * t251 * t1423 + 14 * t254 * t1423 + 126 * t25 * t1673 + 105 * t99 * t1673 - 28 * t170 * t3078 + 70 * t282 * t3219 + 25 * xuv3 * t2944 - 490 * t7 * t3219 - 84 * t9 * t3341
    t3449 = eta * P2
    t3470 = -40 * xuv3 * t1020 - 80 * yuv3 * t1062 - 80 * xuv3 * t1095 + 105 * t113 * t3449 - 84 * t13 * t51 + 210 * t134 * t3449 + 105 * t137 * t3449 - 84 * t15 * t3341 + 140 * yuv3 * t1715 - 40 * yuv3 * t3 - 20 * t3226 * t3001
    t3495 = -40 * yuv3 * t1120 + 140 * t13 * t1359 + 168 * t13 * t1864 - 532 * t282 * t1310 + 20 * t3234 * t1310 + 140 * t7 * t1310 - 532 * t275 * t1359 - 20 * t3244 * t1359 - 56 * t170 * t3303 + 84 * yuv3 * t1862 + 28 * t2954 * t3308 + 84 * t2958 * t3308
    t3518 = 84 * t128 * P1 + 84 * t113 * P2 + 140 * t12 * t1807 + 168 * t12 * t1884 + 20 * t3226 * t1423 - 84 * t2954 * t1673 - 28 * t2958 * t1673 + 56 * t170 * t3342 + 84 * t28 * t1887 - 40 * t1633 - 84 * t1801
    t3541 = -140 * t13 * t3415 + 100 * t2786 * t1310 + 20 * t339 * t1310 + 100 * t3359 * t1359 - 20 * t336 * t1359 - 56 * t170 * t1607 + 28 * t251 * t3308 + 84 * t254 * t3308 - 140 * t275 * t3415 - 20 * zuv3 * t3001 + 100 * yuv3 * t3131
    t3565 = 8 * xuv3 * t1310 + 8 * yuv3 * t1359 + 20 * zuv3 * t1423 - 84 * t251 * t1673 - 28 * t254 * t1673 + 56 * t170 * t3341 - 140 * t282 * t3449 + 100 * xuv3 * t3219 - 8 * t3226 * t3308 - 140 * t7 * t3449 + 8 * yuv3 * t51
    t3587 = -112 * t13 * P1 + 32 * P1 * yuv3 - 112 * t12 * t1887 + 8 * t3226 * t1673 + 8 * zuv3 * t1673 - 112 * yuv3 * t1864 - 8 * zuv3 * t3308 + 40 * yuv3 * t3415 + 40 * xuv3 * t3449 + 8 * t1807 - 112 * t1884 + 32 * t1887
    dMLdt = 105/128 * t521 / (1 + eta) * t525 * t527 * t519 * (t3086 + t2903 + t3033 + t3267 + t3187 + t3112 + t2951 + t3241 + t2846 + t3009 + t2818 + t3587 + t3470 + t3319 + t3292 + t2874 + t2793 + t3445 + t3541 + t3138 + t3058 + t3214 + t3345 + t2926 + t3162 + t2980 + t3422 + t3395 + t3565 + t3370 + t3495 + t3518)

    out = mu3*[dsmadt,dP1dt,dP2dt,dQ1dt,dQ2dt,dMLdt]
    return out;
        
end

function getaveragedterms_3bp_deg6(sma,n,P1,P2,eta,Q1,Q2,GG,r3,xuv3,yuv3,zuv3,mu3)

    dsmadt = 0.0;
    t1 = sma ^ 2
    t2 = t1 ^ 2
    t3 = P1 ^ 2
    t4 = t3 ^ 2
    t5 = t4 * t3
    t6 = P2 * t5
    t7 = xuv3 ^ 2
    t8 = t7 ^ 2
    t9 = t8 * xuv3
    t11 = zuv3 * t9 * Q1
    t14 = Q1 * t6
    t15 = t7 * xuv3
    t16 = yuv3 ^ 2
    t17 = t16 * t15
    t18 = zuv3 * t17
    t21 = t16 ^ 2
    t22 = t21 * xuv3
    t23 = zuv3 * t22
    t26 = t6 * Q2
    t27 = yuv3 * t8
    t28 = zuv3 * t27
    t31 = t16 * yuv3
    t32 = t31 * t7
    t33 = zuv3 * t32
    t36 = t21 * yuv3
    t37 = t36 * Q2
    t38 = zuv3 * t37
    t41 = t4 * P1
    t42 = P2 ^ 2
    t43 = t42 * t41
    t44 = Q1 * t43
    t50 = zuv3 * t36 * Q1
    t53 = t9 * Q2
    t54 = zuv3 * t53
    t58 = Q2 * t43
    t63 = t42 * P2
    t64 = t63 * t4
    t67 = Q1 * t64
    t72 = Q2 * t64
    t79 = t3 * P1
    t80 = t42 ^ 2
    t81 = t80 * t79
    t82 = Q1 * t81
    t89 = 3960 * t11 * t64 + 23760 * t18 * t58 - 54450 * t18 * t67 - 23760 * t23 * t58 + 39600 * t23 * t67 + 27225 * t28 * t72 - 39600 * t28 * t82 - 79200 * t33 * t72 + 94050 * t33 * t82 + 7920 * t38 * t64 - 7920 * t50 * t81
    t93 = Q2 * t81
    t98 = t80 * P2
    t99 = t98 * t3
    t102 = Q1 * t99
    t107 = Q2 * t99
    t114 = t80 * t42
    t115 = t114 * P1
    t116 = Q1 * t115
    t126 = Q2 * t115
    t131 = t80 * t63
    t132 = Q1 * t131
    t133 = zuv3 * t9
    t140 = Q2 * t131
    t145 = zuv3 * t36
    t148 = 1485 * t50 * t115 - 4752 * t54 * t115 - 23760 * t33 * t116 + 23760 * t18 * t126 - 7425 * t23 * t126 + 2112 * t133 * t132 - 5280 * t18 * t132 + 1320 * t23 * t132 + 165 * t145 * t140 + 2640 * t28 * t140 - 2640 * t33 * t140
    t151 = t4 * t79
    t155 = t15 * t151
    t158 = xuv3 * t151
    t161 = t8 * t7
    t164 = t16 * t8
    t167 = t21 * t7
    t170 = t21 * t16
    t173 = yuv3 * t9
    t176 = t31 * t15
    t179 = t36 * xuv3
    t205 = 23760 * t173 * t115 - 39600 * t176 * t115 - 3960 * t161 * t64 - 3168 * t161 * t99 + 51975 * t164 * t64 + 71280 * t164 * t99 - 14850 * t167 * t64 - 50490 * t167 * t99 + 2475 * t170 * t99 + 54450 * t176 * t81 - 16335 * t179 * t81
    t211 = t8 * t131
    t214 = t7 * t131
    t220 = t15 * Q1 * zuv3
    t223 = t16 * xuv3
    t224 = zuv3 * t223
    t227 = yuv3 * t7
    t228 = zuv3 * t227
    t231 = zuv3 * t31
    t232 = Q2 * t231
    t239 = zuv3 * t31 * Q1
    t242 = t15 * Q2
    t243 = zuv3 * t242
    t254 = P2 * t4
    t257 = Q1 * t254
    t260 = Q2 * t254
    t267 = 495 * t11 * t254 - 4950 * t18 * t257 + 1350 * t220 * t64 + 6480 * t224 * t58 - 6750 * t224 * t67 + 6750 * t228 * t72 - 4050 * t228 * t82 - 2160 * t239 * t43 - 810 * t243 * t43 - 7920 * t38 * t254 + 2475 * t28 * t260
    t277 = t42 * t79
    t278 = Q1 * t277
    t287 = Q2 * t277
    t299 = t63 * t3
    t302 = Q1 * t299
    t307 = Q2 * t299
    t320 = 6750 * t228 * t107 - 7920 * t11 * t299 + 810 * t239 * t115 + 2160 * t243 * t115 - 6480 * t228 * t116 + 49500 * t18 * t302 + 24750 * t23 * t302 - 1350 * t232 * t99 - 24750 * t28 * t307 + 7920 * t38 * t299 - 49500 * t33 * t307
    t324 = t80 * P1
    t325 = Q1 * t324
    t334 = Q2 * t324
    t339 = zuv3 * t15
    t349 = Q1 * t98
    t354 = Q2 * t98
    t365 = t16 * t7
    t368 = yuv3 * t15
    t371 = t31 * xuv3
    t374 = 7920 * t133 * t349 + 90 * t231 * t140 - 495 * t145 * t354 + 810 * yuv3 * t155 - 2160 * t31 * t158 - 2475 * t23 * t349 + 4950 * t33 * t354 - 6750 * t365 * t6 - 3240 * t43 * t368 - 6210 * t371 * t43 + 675 * t8 * t6
    t380 = t15 * t41
    t383 = xuv3 * t41
    t423 = -10800 * t368 * t115 + 7920 * t161 * t299 - 54450 * t164 * t299 + 74250 * t167 * t299 - 4950 * t170 * t299 + 44550 * t173 * t277 - 74250 * t176 * t277 + 22770 * t179 * t277 + 1215 * t21 * t99 - 11340 * t365 * t99 - 2160 * t8 * t99
    t438 = t8 * t98
    t441 = t7 * t98
    t447 = Q1 * xuv3
    t448 = zuv3 * t447
    t451 = yuv3 * Q2
    t452 = zuv3 * t451
    t456 = Q1 * yuv3 * zuv3
    t459 = Q2 * xuv3
    t460 = zuv3 * t459
    t477 = 1350 * t224 * t257 - 1350 * t228 * t260 + 7200 * t232 * t254 + 225 * t456 * t43 - 225 * t460 * t43 - 25 * t448 * t6 + 150 * t448 * t64 - 200 * t452 * t6 - 375 * t452 * t64 + 450 * t456 * t81 - 450 * t460 * t81
    t503 = P2 * t3
    t506 = Q1 * t503
    t511 = Q2 * t503
    t528 = -495 * t11 * t503 + 225 * t456 * t115 - 225 * t460 * t115 + 1980 * t18 * t506 - 25650 * t228 * t325 + 2475 * t23 * t506 + 1350 * t239 * t324 + 8550 * t243 * t324 - 990 * t28 * t511 - 4950 * t33 * t511 - 3960 * t38 * t503
    t532 = t42 * P1
    t533 = Q1 * t532
    t542 = Q2 * t532
    t547 = xuv3 * zuv3
    t550 = yuv3 * zuv3
    t560 = Q1 * t63
    t567 = Q2 * t63
    t580 = 3960 * t133 * t560 + 495 * t145 * t567 + 225 * yuv3 * t158 + 375 * t16 * t6 + 4950 * t18 * t560 + 1350 * t224 * t349 - 1350 * t228 * t354 + 990 * t23 * t560 - 2475 * t28 * t567 - 1980 * t33 * t567 + 150 * t7 * t6
    t583 = yuv3 * xuv3
    t610 = t15 * t79
    t613 = xuv3 * t79
    t630 = 4455 * yuv3 * t9 * t79 + 225 * t16 * t99 + 3960 * t161 * t503 - 22275 * t164 * t503 - 23760 * t167 * t503 - 450 * t371 * t277 - 10800 * t365 * t299 - 5850 * t8 * t299 - 1980 * t31 * t610 - 6435 * t36 * t613 + 1350 * t7 * t99
    t655 = t8 * t63
    t658 = t7 * t63
    t675 = 2475 * t16 * t655 - 3960 * t161 * t63 - 495 * t170 * t63 + 5940 * t21 * t658 + 135 * t21 * t98 - 75 * t448 * t254 - 1125 * t452 * t254 + 1200 * t456 * t277 - 1200 * t460 * t277 + 1050 * t448 * t299 - 1050 * t452 * t299
    t712 = P2 * Q1
    t719 = Q2 * P2
    t724 = 165 * t133 * t712 + 330 * t18 * t712 - 1890 * t224 * t560 + 1890 * t228 * t567 + 165 * t23 * t712 - 270 * t231 * t567 - 165 * t28 * t719 - 330 * t33 * t719 - 4050 * t339 * t560 + 1125 * t547 * t349 + 75 * t550 * t354
    t763 = t15 * P1
    t766 = xuv3 * P1
    t769 = -1485 * yuv3 * t9 * P1 + 600 * t16 * t299 - 1215 * t21 * t503 + 2550 * t7 * t299 - 2970 * t31 * t763 + 2475 * t583 * t324 - 1485 * t36 * t766 + 25110 * t365 * t503 + 23490 * t368 * t532 - 6210 * t371 * t532 - 3375 * t8 * t503
    t783 = t8 * P2
    t786 = t7 * P2
    t813 = -180 * t224 * t712 + 180 * t228 * t719 + 180 * t231 * t719 - 180 * t339 * t712 + 60 * t448 * t503 - 780 * t452 * t503 + 720 * t456 * t532 - 720 * t460 * t532 + 780 * t547 * t560 - 60 * t550 * t567 - 75 * t254
    t850 = -60 * t16 * P2 - 135 * t21 * P2 + 1890 * t16 * t786 + 40 * t547 * t712 - 40 * t550 * t719 - 720 * yuv3 * t766 + 40 * P2 + 60 * t503 + 60 * t63 + 2025 * t783 - 780 * t786
    t857 = 0.1e1 / eta
    t858 = r3 ^ 2
    t860 = t858 ^ 2
    t862 = 0.1e1 / t860 / t858 / r3
    t863 = t862 * t857
    t864 = 0.1e1 / n
    t865 = t864 * t863
    dP1dt = -21/128 * t865 * (2475 * t583 * t81 - 20250 * t368 * t277 + 2475 * t170 * t503 + 1125 * t583 * t115 - 14850 * t368 * t324 - 22275 * t173 * t532 - 11880 * t176 * t532 + 10395 * t179 * t532 - 75 * t16 * t131 + 540 * t16 * t441 + 270 * t220 * t503 - 1890 * t224 * t506 + 1890 * t228 * t511 + 4050 * t232 * t503 + 1200 * t456 * t324 - 1200 * t460 * t324 - 6480 * t228 * t533 - 2160 * t239 * t532 + 2160 * t243 * t532 + 6480 * t224 * t542 - 165 * t145 * t719 + 975 * yuv3 * t383 + 900 * t7 * t254 + 675 * t16 * t254 + 3450 * t583 * t277 - 3510 * yuv3 * t610 + 6390 * t31 * t613 - 75 * t16 * t98 - 4590 * t16 * t658 + 135 * t21 * t63 - 1320 * t161 * P2 - 2475 * t16 * t783 - 990 * t21 * t786 + 165 * t170 * P2 - 25 * t131 + t320 - 25 * t6 + t374 - 2160 * t211 + 600 * t214 - 75 * t98 - 75 * t99 + t148 - 1320 * t28 * t26 + 5280 * t33 * t26 - 2112 * t38 * t6 - 165 * t11 * t6 + t89 + t528 - 4455 * yuv3 * t9 * t41 + 1485 * yuv3 * t9 * t151 + t267 + t769 - 150 * t299 + t724 - 480 * yuv3 * t613 - 270 * t7 * t503 - 990 * t16 * t503 - 2880 * t583 * t532 + 2160 * yuv3 * t763 + 2160 * t31 * t766 + 210 * t16 * t63 + 2640 * t18 * t14 - 2640 * t23 * t14 + 7425 * t28 * t44 - 23760 * t33 * t44 + 4752 * t50 * t43 - 1485 * t54 * t43 + 7920 * t54 * t81 - 94050 * t18 * t93 + 39600 * t23 * t93 - 7920 * t11 * t99 + 79200 * t18 * t102 - 27225 * t23 * t102 - 39600 * t28 * t107 + 54450 * t33 * t107 - 3960 * t38 * t99 + 23760 * t28 * t116 - 7920 * t31 * t155 + 4752 * t36 * t158 + 1320 * t161 * t6 - 27225 * t164 * t6 + 39600 * t167 * t6 - 2640 * t170 * t6 - 22275 * t173 * t43 + 86130 * t176 * t43 - 19008 * t179 * t43 + 7425 * t179 * t115 + 2112 * t161 * t131 - 7920 * t16 * t211 + 3960 * t21 * t214 - 165 * t170 * t131 - 90 * t220 * t6 + 720 * t224 * t14 - 720 * t228 * t26 - 1470 * t658 + t580 + 2160 * t16 * t214 - 135 * t21 * t131 + 3168 * t161 * t98 + 7920 * t16 * t438 - 8910 * t21 * t441 + 495 * t170 * t98 + 4050 * t228 * t278 - 8550 * t239 * t277 - 1350 * t243 * t277 + 25650 * t224 * t287 + 375 * t448 * t99 - 150 * t452 * t99 + 2700 * t220 * t299 - 27000 * t224 * t302 + 27000 * t228 * t307 - 2700 * t232 * t299 - 4050 * t224 * t334 + 7425 * t28 * t533 + 8910 * t33 * t533 + 1485 * t50 * t532 - 1485 * t54 * t532 - 8910 * t18 * t542 - 7425 * t23 * t542 + 200 * t547 * t132 + 25 * t550 * t140 - 7200 * t339 * t349 + 1575 * t583 * t43 + 540 * yuv3 * t380 - 6390 * t31 * t383 + 900 * t7 * t64 + 675 * t16 * t64 + 675 * t8 * t254 - 20250 * t365 * t254 + 1350 * t21 * t254 + 5175 * t655 + 1650 * t441 + t477 + t675 - 75 * t64 + t813 + t850 - 5040 * t438 + t423 + t205 + t630 + 1440 * t232 * t6 + 2430 * t228 * t44 - 1350 * t239 * t81 + 1350 * t243 * t81 + 4050 * t224 * t93 - 14850 * t28 * t278 + 14850 * t33 * t278 + 7920 * t50 * t277 + 2970 * t54 * t277 - 14850 * t18 * t287 - 39600 * t23 * t287 - 6750 * t224 * t102 - 2430 * t224 * t126 + 39600 * t28 * t325 + 14850 * t33 * t325 - 2970 * t50 * t324 - 7920 * t54 * t324 - 14850 * t18 * t334 + 14850 * t23 * t334 - 1440 * t339 * t132 + 720 * t224 * t132 - 720 * t228 * t140 + 12870 * t31 * t380 + 3168 * t36 * t383 + 675 * t8 * t64 - 20250 * t365 * t64 + 1350 * t21 * t64 - 3960 * t161 * t254 + 51975 * t164 * t254 - 14850 * t167 * t254 - 14850 * t368 * t81 + 4050 * t371 * t115 + 54450 * t176 * t324 - 16335 * t179 * t324) * t2
    t868 = Q1 * t151
    t875 = Q2 * t151
    t913 = -3960 * t11 * t43 - 23760 * t18 * t26 + 54450 * t18 * t44 + 23760 * t23 * t26 - 39600 * t23 * t44 - 27225 * t28 * t58 + 39600 * t28 * t67 + 79200 * t33 * t58 - 94050 * t67 * t33 - 7920 * t38 * t43 + 7920 * t50 * t64
    t958 = 23760 * t33 * t102 - 23760 * t18 * t107 + 7425 * t23 * t107 - 2112 * t11 * t115 - 165 * t38 * t115 + 5280 * t18 * t116 - 1320 * t23 * t116 - 2640 * t28 * t126 + 2640 * t33 * t126 - 1485 * t50 * t99 + 4752 * t54 * t99
    t963 = t8 * t151
    t966 = t7 * t151
    t1006 = 2640 * t161 * t115 - 39600 * t164 * t115 + 14850 * t164 * t81 - 51975 * t167 * t81 + 3168 * t170 * t43 + 3960 * t170 * t81 + 16335 * t173 * t64 + 19008 * t173 * t99 - 54450 * t176 * t64 - 86130 * t176 * t99 + 22275 * t179 * t99
    t1015 = t15 * t131
    t1018 = xuv3 * t131
    t1044 = Q1 * t41
    t1049 = Q2 * t41
    t1056 = -495 * t133 * t1044 + 4950 * t18 * t1044 + 7920 * t145 * t1049 - 2475 * t28 * t1049 - 1350 * t220 * t43 - 6480 * t224 * t26 + 6750 * t224 * t44 - 6750 * t228 * t58 + 4050 * t228 * t67 + 2160 * t239 * t6 + 810 * t243 * t6
    t1103 = 6480 * t228 * t102 + 7920 * t11 * t277 - 49500 * t18 * t278 - 6750 * t228 * t93 - 24750 * t23 * t278 + 1350 * t232 * t81 - 810 * t239 * t99 - 2160 * t243 * t99 - 7920 * t38 * t277 + 24750 * t28 * t287 + 49500 * t33 * t287
    t1147 = -7920 * t11 * t324 - 90 * t232 * t115 + 2160 * t21 * t151 - 2160 * t16 * t966 + 2475 * t23 * t325 + 495 * t38 * t324 - 4950 * t33 * t334 - 4050 * t368 * t6 + 10800 * t371 * t6 - 1215 * t8 * t43 + 135 * t963
    t1156 = t8 * t41
    t1159 = t7 * t41
    t1195 = 4950 * t161 * t277 - 74250 * t164 * t277 + 54450 * t167 * t277 - 7920 * t170 * t277 - 22770 * t173 * t299 + 74250 * t176 * t299 - 44550 * t179 * t299 - 675 * t21 * t81 + 20250 * t365 * t81 + 6210 * t368 * t99 + 3240 * t371 * t99
    t1214 = t15 * t98
    t1217 = xuv3 * t98
    t1243 = -1350 * t224 * t1044 + 1350 * t228 * t1049 - 7200 * t231 * t1049 - 150 * t448 * t43 + 375 * t452 * t43 - 225 * t456 * t6 - 450 * t456 * t64 + 225 * t460 * t6 + 450 * t460 * t64 + 25 * t547 * t868 + 200 * t550 * t875
    t1269 = Q1 * t79
    t1276 = Q2 * t79
    t1293 = 495 * t133 * t1269 - 1980 * t18 * t1269 - 2475 * t23 * t1269 + 3960 * t145 * t1276 + 990 * t28 * t1276 + 4950 * t33 * t1276 + 25650 * t228 * t302 - 1350 * t239 * t299 - 8550 * t243 * t299 - 225 * t456 * t99 + 225 * t460 * t99
    t1337 = -3960 * t11 * t532 - 600 * t16 * t151 - 4950 * t18 * t533 - 1350 * t224 * t325 + 1350 * t228 * t334 - 990 * t23 * t533 + 2475 * t28 * t542 + 1980 * t33 * t542 - 495 * t38 * t532 - 1125 * t583 * t6 + 75 * t966
    t1364 = t8 * t79
    t1367 = t7 * t79
    t1384 = -5940 * t16 * t1364 - 2475 * t21 * t1367 + 495 * t161 * t79 + 3960 * t170 * t79 - 10395 * t173 * t503 + 11880 * t176 * t503 + 22275 * t179 * t503 + 5850 * t21 * t277 + 450 * t368 * t299 + 20250 * t371 * t299 - 1575 * t583 * t99
    t1414 = t15 * t63
    t1417 = xuv3 * t63
    t1432 = 6435 * yuv3 * t9 * t63 + 75 * t547 * t1044 + 1125 * t550 * t1049 + 6390 * yuv3 * t1214 - 540 * t31 * t1217 + 1980 * t31 * t1414 - 4455 * t36 * t1417 - 1200 * t456 * t254 + 1200 * t460 * t254 - 1050 * t448 * t277 + 1050 * t452 * t277
    t1469 = P1 * Q1
    t1476 = P1 * Q2
    t1481 = -165 * t133 * t1469 - 330 * t18 * t1469 - 165 * t23 * t1469 + 165 * t28 * t1476 + 330 * t33 * t1476 + 4050 * t220 * t532 + 1890 * t224 * t533 - 1890 * t228 * t542 + 270 * t232 * t532 - 1125 * t448 * t324 - 75 * t452 * t324
    t1520 = 4590 * t16 * t1367 - 900 * t16 * t324 + 3375 * t21 * t532 - 5175 * t21 * t79 - 3450 * t583 * t299 - 675 * t7 * t324 - 25110 * t365 * t532 + 6210 * t368 * t503 - 23490 * t371 * t503 + 1215 * t8 * t532 + 25 * t115
    t1525 = t8 * P1
    t1528 = t7 * P1
    t1542 = t15 * P2
    t1545 = P2 * xuv3
    t1570 = -60 * t547 * t1269 + 780 * t550 * t1276 + 180 * t224 * t1469 + 180 * t339 * t1469 - 180 * t228 * t1476 - 180 * t231 * t1476 - 780 * t448 * t532 + 60 * t452 * t532 - 720 * t456 * t503 + 720 * t460 * t503 + 75 * t41
    t1607 = 780 * t16 * P1 + 480 * yuv3 * t1417 - 40 * t547 * t1469 + 40 * t550 * t1476 - 2160 * yuv3 * t1542 - 2160 * t31 * t1545 + 720 * yuv3 * t1545 - 40 * P1 + 60 * t1528 - 60 * t532 - 60 * t79
    dP2dt = -21/128 * t865 * (-4050 * t228 * t257 + 8550 * t239 * t254 + 1350 * t243 * t254 - 25650 * t224 * t260 - 375 * t448 * t81 + 150 * t452 * t81 - 2700 * t220 * t277 + 27000 * t224 * t278 - 27000 * t228 * t287 + 2700 * t232 * t277 + 75 * t81 + t958 + t1243 - 675 * t7 * t81 + t1520 + 75 * t43 + t1006 + t1570 + t1103 - 135 * t1156 + 75 * t1159 - 135 * t1364 - 210 * t1367 - 4752 * yuv3 * t9 * t131 - 3168 * yuv3 * t9 * t98 + 1485 * yuv3 * t9 * P2 + t1293 + t1337 + 135 * t1525 - 1350 * t243 * t64 - 4050 * t224 * t72 + 14850 * t28 * t257 - 14850 * t33 * t257 - 7920 * t50 * t254 - 2970 * t54 * t254 + 14850 * t18 * t260 + 39600 * t23 * t260 + 6750 * t224 * t82 + 2430 * t224 * t107 - 39600 * t28 * t302 - 14850 * t33 * t302 + 2970 * t50 * t299 + 7920 * t54 * t299 + 14850 * t18 * t307 - 14850 * t23 * t307 + 1440 * t220 * t115 - 720 * t224 * t116 + 720 * t228 * t126 + 11340 * t365 * t43 + 2160 * t21 * t43 - 495 * t161 * t41 + 8910 * t16 * t1156 - 7920 * t21 * t1159 - 3168 * t170 * t41 + 14850 * t371 * t64 + 16335 * t173 * t254 - 54450 * t176 * t254 - 1350 * t8 * t81 + 6750 * t365 * t115 - 675 * t21 * t115 + 14850 * t164 * t324 - 51975 * t167 * t324 + 3960 * t170 * t324 + 2160 * yuv3 * t1015 - 810 * t31 * t1018 - 12870 * t31 * t1214 + 4455 * t36 * t1217 + 1470 * t16 * t79 + 2880 * t583 * t503 + 990 * t7 * t532 + 270 * t16 * t532 - 1890 * t16 * t1528 - 2025 * t21 * P1 - 975 * yuv3 * t1217 - 6390 * yuv3 * t1414 + 3510 * t31 * t1417 + 2970 * t31 * t1542 + 1485 * t36 * t1545 - 1650 * t16 * t41 - 2475 * t583 * t254 - 600 * t7 * t277 - 2550 * t16 * t277 - 165 * t161 * P1 + 990 * t16 * t1525 + 2475 * t21 * t1528 + 1320 * t170 * P1 + 165 * t145 * t1476 - 270 * t339 * t1269 + 1890 * t224 * t1269 - 1890 * t228 * t1276 - 4050 * t231 * t1276 - 1200 * t456 * t299 + 1200 * t460 * t299 + 6480 * t228 * t506 + 2160 * t239 * t503 + 25 * t151 + t1607 + t1481 + t1384 - 7425 * t28 * t14 + 23760 * t33 * t14 - 4752 * t50 * t6 + 1485 * t54 * t6 - 7920 * t54 * t64 + 94050 * t18 * t72 - 39600 * t23 * t72 + 7920 * t11 * t81 - 79200 * t18 * t82 + 27225 * t23 * t82 + 39600 * t28 * t93 - 54450 * t33 * t93 + 3960 * t38 * t81 - 23760 * t28 * t102 + 165 * t161 * t151 - 3960 * t16 * t963 + 7920 * t21 * t966 - 2112 * t170 * t151 - 7425 * t173 * t6 + 39600 * t176 * t6 - 23760 * t179 * t6 - 2475 * t161 * t43 + 50490 * t164 * t43 - 71280 * t167 * t43 + 27225 * t167 * t115 - 1320 * t170 * t115 + 7920 * t31 * t1015 - 1485 * t36 * t1018 + 90 * t339 * t868 - 720 * t224 * t868 + 720 * t228 * t875 - 1440 * t231 * t875 - 2430 * t228 * t14 + 1350 * t239 * t64 + t913 + 150 * t277 + 165 * t133 * t868 - 2640 * t18 * t868 + 2640 * t23 * t868 + 1320 * t28 * t875 - 5280 * t33 * t875 + 2112 * t145 * t875 + t1147 - 2160 * t243 * t503 - 6480 * t224 * t511 - 375 * t7 * t115 - 150 * t16 * t115 - 1350 * t8 * t324 + 20250 * t365 * t324 - 675 * t21 * t324 - 2475 * t161 * t532 + 23760 * t164 * t532 + 22275 * t167 * t532 - 3960 * t170 * t532 - 225 * yuv3 * t1018 - 225 * t7 * t43 - 1350 * t16 * t43 - 540 * t16 * t1159 + 5040 * t21 * t41 - 2475 * t583 * t64 + 14850 * t371 * t254 - 900 * t16 * t81 + 10800 * t365 * t277 - 1485 * t50 * t503 + 1485 * t54 * t503 + 8910 * t18 * t511 + 7425 * t23 * t511 - 200 * t448 * t115 - 25 * t452 * t115 + 7200 * t220 * t324 + 4050 * t224 * t307 - 7425 * t28 * t506 - 8910 * t33 * t506 + t1432 + 75 * t324 + t1056 + t1195) * t2
    t1616 = Q1 ^ 2
    t1617 = t1616 * Q1
    t1618 = t1617 * t5
    t1621 = t16 * t242
    t1624 = t21 * t459
    t1627 = t1616 * t5
    t1628 = Q2 ^ 2
    t1630 = yuv3 * t8 * t1628
    t1633 = t7 * t1628
    t1634 = t31 * t1633
    t1637 = t36 * t1628
    t1640 = Q1 * t5
    t1641 = t1628 * Q2
    t1642 = t9 * t1641
    t1645 = t15 * t1641
    t1646 = t16 * t1645
    t1649 = xuv3 * t1641
    t1650 = t21 * t1649
    t1653 = t1628 ^ 2
    t1654 = t1653 * t5
    t1661 = P2 * t41
    t1662 = t1617 * t1661
    t1663 = t8 * Q2
    t1664 = yuv3 * t1663
    t1667 = t7 * Q2
    t1668 = t31 * t1667
    t1671 = Q2 * t1617
    t1672 = t36 * t1671
    t1675 = t1628 * t1616
    t1676 = t9 * t1675
    t1679 = t1616 * t1661
    t1680 = t15 * t1628
    t1681 = t16 * t1680
    t1684 = xuv3 * t1628
    t1685 = t21 * t1684
    t1688 = Q1 * t1661
    t1690 = yuv3 * t8 * t1641
    t1693 = t7 * t1641
    t1694 = t31 * t1693
    t1697 = -2640 * t1621 * t1618 + 2640 * t1624 * t1618 + 165 * t53 * t1618 + 1320 * t1630 * t1627 - 5280 * t1634 * t1627 + 2112 * t1637 * t1627 + 165 * t1642 * t1640 - 2640 * t1646 * t1640 + 2640 * t1650 * t1640 + 1320 * t27 * t1654 - 5280 * t32 * t1654 + 2112 * t36 * t1654 - 4752 * t1672 * t1661 + 1485 * t1676 * t1661 - 7425 * t1664 * t1662 + 23760 * t1668 * t1662 - 23760 * t1681 * t1679 + 23760 * t1685 * t1679 - 7425 * t1690 * t1688 + 23760 * t1694 * t1688
    t1698 = t1641 * Q1
    t1699 = t36 * t1698
    t1702 = t9 * t1653
    t1705 = t15 * t1653
    t1706 = t16 * t1705
    t1709 = xuv3 * t1653
    t1710 = t21 * t1709
    t1713 = t42 * t4
    t1714 = t9 * t1671
    t1717 = t1617 * t1713
    t1722 = t1616 * t1713
    t1727 = t36 * t1675
    t1730 = t9 * t1698
    t1733 = Q1 * t1713
    t1739 = yuv3 * t8 * t1653
    t1742 = t7 * t1653
    t1743 = t31 * t1742
    t1746 = t36 * t1653
    t1749 = t63 * t79
    t1750 = t1617 * t1749
    t1759 = 54450 * t1621 * t1717 - 39600 * t1624 * t1717 - 27225 * t1630 * t1722 + 79200 * t1634 * t1722 + 54450 * t1646 * t1733 - 39600 * t1650 * t1733 - 4752 * t1699 * t1661 + 1485 * t1702 * t1661 - 23760 * t1706 * t1661 + 23760 * t1710 * t1661 + 39600 * t1664 * t1750 - 94050 * t1668 * t1750 + 7920 * t1672 * t1749 - 7920 * t1676 * t1749 - 3960 * t1714 * t1713 - 7920 * t1727 * t1713 - 3960 * t1730 * t1713 - 27225 * t1739 * t1713 + 79200 * t1743 * t1713 - 7920 * t1746 * t1713
    t1761 = t1616 * t1749
    t1766 = Q1 * t1749
    t1779 = t80 * t3
    t1782 = t1617 * t1779
    t1787 = t1616 * t1779
    t1796 = Q1 * t1779
    t1807 = -79200 * t1621 * t1782 + 27225 * t1624 * t1782 + 39600 * t1630 * t1787 - 54450 * t1634 * t1787 - 79200 * t1646 * t1796 + 27225 * t1650 * t1796 + 94050 * t1681 * t1761 - 39600 * t1685 * t1761 + 39600 * t1690 * t1766 - 94050 * t1694 * t1766 + 7920 * t1699 * t1749 - 7920 * t1702 * t1749 + 94050 * t1706 * t1749 - 39600 * t1710 * t1749 + 7920 * t1714 * t1779 + 3960 * t1727 * t1779 + 7920 * t1730 * t1779 + 39600 * t1739 * t1779 - 54450 * t1743 * t1779 + 3960 * t1746 * t1779
    t1808 = t98 * P1
    t1809 = t1617 * t1808
    t1818 = t1616 * t1808
    t1823 = Q1 * t1808
    t1836 = t1617 * t114
    t1843 = t1616 * t114
    t1850 = Q1 * t114
    t1855 = 5280 * t1621 * t1836 - 1320 * t1624 * t1836 - 2640 * t1630 * t1843 + 2640 * t1634 * t1843 - 165 * t1637 * t1843 - 2112 * t1642 * t1850 + 5280 * t1646 * t1850 - 23760 * t1664 * t1809 + 23760 * t1668 * t1809 - 1485 * t1672 * t1808 + 4752 * t1676 * t1808 - 23760 * t1681 * t1818 + 7425 * t1685 * t1818 - 23760 * t1690 * t1823 + 23760 * t1694 * t1823 - 1485 * t1699 * t1808 + 4752 * t1702 * t1808 - 23760 * t1706 * t1808 + 7425 * t1710 * t1808 - 2112 * t53 * t1836
    t1860 = t1653 * t114
    t1869 = t16 * t459
    t1872 = yuv3 * t1633
    t1875 = t31 * t1628
    t1880 = t16 * t1649
    t1893 = t1628 * t5
    t1900 = yuv3 * t1667
    t1903 = t31 * t1671
    t1906 = -720 * t1869 * t1618 + 90 * t242 * t1618 - 5280 * t1621 * t1640 + 5280 * t1624 * t1640 + 720 * t1872 * t1627 - 1440 * t1875 * t1627 + 90 * t1645 * t1640 - 720 * t1880 * t1640 + 330 * t53 * t1640 - 1320 * t1650 * t1850 + 720 * t227 * t1654 - 1440 * t31 * t1654 + 2160 * t1903 * t1661 - 2430 * t1900 * t1662 - 2640 * t27 * t1860 + 2640 * t32 * t1860 - 165 * t36 * t1860 + 2640 * t27 * t1893 - 10560 * t32 * t1893 + 4224 * t36 * t1893
    t1907 = t15 * t1675
    t1910 = t16 * t1684
    t1913 = yuv3 * t1693
    t1916 = t31 * t1698
    t1923 = Q2 * Q1
    t1924 = t36 * t1923
    t1929 = t16 * t1709
    t1932 = t9 * t1628
    t1939 = t15 * t1671
    t1946 = t15 * t1698
    t1951 = t9 * t1923
    t1958 = 108900 * t1621 * t1733 - 79200 * t1624 * t1733 - 47520 * t1681 * t1661 + 47520 * t1685 * t1661 + 810 * t1705 * t1661 + 810 * t1907 * t1661 + 2160 * t1916 * t1661 - 9504 * t1924 * t1661 - 6480 * t1929 * t1661 + 2970 * t1932 * t1661 - 14850 * t1664 * t1688 + 47520 * t1668 * t1688 - 6480 * t1910 * t1679 - 2430 * t1913 * t1688 - 1350 * t1939 * t1713 - 1350 * t1946 * t1713 - 7920 * t1951 * t1713 + 6750 * t1869 * t1717 - 6750 * t1872 * t1722 + 6750 * t1880 * t1733
    t1960 = yuv3 * t1742
    t1969 = t1617 * t4
    t1974 = t1616 * t4
    t1979 = Q1 * t4
    t1984 = t1653 * t4
    t2005 = 4950 * t1621 * t1969 - 54450 * t1630 * t1713 - 2475 * t1630 * t1974 + 158400 * t1634 * t1713 - 15840 * t1637 * t1713 + 7920 * t1637 * t1974 - 495 * t1642 * t1979 + 4950 * t1646 * t1979 + 79200 * t1664 * t1766 - 188100 * t1668 * t1766 - 6750 * t1960 * t1713 + 1350 * t1903 * t1749 - 1350 * t1907 * t1749 + 1350 * t1916 * t1749 + 4050 * t1900 * t1750 - 4050 * t1910 * t1761 + 4050 * t1913 * t1766 - 495 * t53 * t1969 - 2475 * t27 * t1984 + 7920 * t36 * t1984
    t2018 = P2 * t79
    t2019 = t1617 * t2018
    t2029 = t1616 * t2018
    t2034 = Q1 * t2018
    t2051 = t31 * t1675
    t2054 = 14850 * t1681 * t2029 + 39600 * t1685 * t2029 + 14850 * t1690 * t2034 - 14850 * t1694 * t2034 - 7920 * t1699 * t2018 - 2970 * t1702 * t2018 + 14850 * t1706 * t2018 + 39600 * t1710 * t2018 + 1350 * t2051 * t1779 + 6750 * t1869 * t1782 - 6750 * t1872 * t1787
    t2069 = t31 * t1653
    t2078 = t42 * t3
    t2081 = t1617 * t2078
    t2086 = t1616 * t2078
    t2095 = Q1 * t2078
    t2104 = -158400 * t1621 * t1796 - 49500 * t1621 * t2081 + 54450 * t1624 * t1796 - 24750 * t1624 * t2081 + 79200 * t1630 * t1779 + 24750 * t1630 * t2086 - 108900 * t1634 * t1779 + 49500 * t1634 * t2086 + 7920 * t1637 * t1779 - 49500 * t1646 * t2095 - 24750 * t1650 * t2095 + 7920 * t1714 * t2078 - 7920 * t1727 * t2078 + 7920 * t1730 * t2078 + 24750 * t1739 * t2078 + 49500 * t1743 * t2078 + 15840 * t1951 * t1779 - 6750 * t1960 * t1779 + 1350 * t2069 * t1779 + 6750 * t1880 * t1796
    t2135 = t63 * P1
    t2136 = t1617 * t2135
    t2145 = t1616 * t2135
    t2148 = -47520 * t1664 * t1823 - 39600 * t1664 * t2136 + 47520 * t1668 * t1823 - 14850 * t1668 * t2136 + 2970 * t1672 * t2135 + 7920 * t1676 * t2135 - 47520 * t1681 * t1808 + 14850 * t1681 * t2145 + 14850 * t1685 * t1808 - 2160 * t1705 * t1808 - 7920 * t1746 * t2078 - 810 * t1903 * t1808 - 2160 * t1907 * t1808 - 810 * t1916 * t1808 - 2970 * t1924 * t1808 + 2430 * t1929 * t1808 + 9504 * t1932 * t1808 + 6480 * t1900 * t1809 + 2430 * t1910 * t1818 + 6480 * t1913 * t1823
    t2152 = Q1 * t2135
    t2187 = t1628 * t114
    t2192 = 10560 * t1621 * t1850 - 2640 * t1624 * t1850 + 1440 * t1645 * t1850 - 14850 * t1685 * t2145 - 39600 * t1690 * t2152 - 14850 * t1694 * t2152 + 2970 * t1699 * t2135 + 7920 * t1702 * t2135 + 14850 * t1706 * t2135 - 14850 * t1710 * t2135 - 720 * t1869 * t1836 + 1440 * t242 * t1836 + 720 * t1872 * t1843 - 90 * t1875 * t1843 - 720 * t1880 * t1850 - 4224 * t53 * t1850 + 720 * t227 * t1860 - 90 * t31 * t1860 - 5280 * t27 * t2187 + 5280 * t32 * t2187
    t2195 = t1617 * t80
    t2200 = t1616 * t80
    t2205 = Q1 * t80
    t2210 = t1653 * t80
    t2217 = yuv3 * t1628
    t2232 = t8 * t5
    t2235 = t7 * t5
    t2240 = 25 * t459 * t1618 + 2475 * t1624 * t2195 + 200 * t2217 * t1627 - 4950 * t1634 * t2200 + 495 * t1637 * t2200 + 25 * t1649 * t1640 - 1440 * t1869 * t1640 + 180 * t242 * t1640 - 7920 * t1642 * t2205 + 2475 * t1650 * t2205 + 200 * yuv3 * t1654 + 1440 * t227 * t1893 - 2880 * t31 * t1893 - 330 * t36 * t2187 - 7920 * t53 * t2195 - 4950 * t32 * t2210 + 495 * t36 * t2210 + 1320 * yuv3 * t2232 - 5280 * t31 * t2235 + 2112 * t36 * t5
    t2243 = yuv3 * t1671
    t2246 = xuv3 * t1675
    t2249 = yuv3 * t1698
    t2254 = t31 * t1923
    t2269 = xuv3 * t1671
    t2272 = yuv3 * t1675
    t2275 = xuv3 * t1698
    t2278 = t15 * t1923
    t2283 = yuv3 * t1653
    t2292 = 1620 * t1680 * t1661 - 23760 * t17 * t1661 + 225 * t1709 * t1661 - 12960 * t1910 * t1661 + 23760 * t22 * t1661 - 225 * t2243 * t1661 + 225 * t2246 * t1661 - 225 * t2249 * t1661 + 4320 * t2254 * t1661 + 1485 * t9 * t1661 - 4860 * t1900 * t1688 - 13500 * t1872 * t1713 - 150 * t2269 * t1713 + 375 * t2272 * t1713 - 150 * t2275 * t1713 - 2700 * t2278 * t1713 + 375 * t2283 * t1713 - 27225 * t27 * t1713 + 79200 * t32 * t1713 + 13500 * t1869 * t1733
    t2311 = t1628 * t4
    t2334 = 9900 * t1621 * t1979 - 2700 * t1680 * t1749 + 450 * t1709 * t1749 - 7920 * t36 * t1713 - 8100 * t1910 * t1749 - 450 * t2243 * t1749 + 450 * t2246 * t1749 - 450 * t2249 * t1749 + 2700 * t2254 * t1749 - 7920 * t9 * t1749 + 8100 * t1900 * t1766 - 1350 * t1869 * t1969 + 1350 * t1872 * t1974 - 7200 * t1875 * t1974 - 1350 * t1880 * t1979 - 990 * t53 * t1979 + 1350 * t227 * t1984 - 7200 * t31 * t1984 - 4950 * t27 * t2311 + 15840 * t36 * t2311
    t2376 = 29700 * t1664 * t2034 - 29700 * t1668 * t2034 + 29700 * t1681 * t2018 + 79200 * t1685 * t2018 + 94050 * t17 * t1749 + 1350 * t1705 * t2018 - 39600 * t22 * t1749 - 375 * t2269 * t1779 + 150 * t2272 * t1779 - 375 * t2275 * t1779 + 13500 * t1869 * t1796 - 4050 * t1900 * t2019 + 8550 * t1903 * t2018 + 1350 * t1907 * t2018 - 25650 * t1910 * t2029 - 4050 * t1913 * t2034 + 8550 * t1916 * t2018 - 15840 * t1924 * t2018 - 25650 * t1929 * t2018 - 5940 * t1932 * t2018
    t2418 = t1617 * t3
    t2421 = -99000 * t1621 * t2095 - 49500 * t1624 * t2095 + 49500 * t1630 * t2078 + 99000 * t1634 * t2078 - 15840 * t1637 * t2078 + 27000 * t1880 * t2095 - 2700 * t1946 * t2078 + 15840 * t1951 * t2078 - 27000 * t1960 * t2078 + 2700 * t2069 * t2078 + 495 * t53 * t2418
    t2431 = t1616 * t3
    t2438 = Q1 * t3
    t2445 = t1653 * t3
    t2470 = -1980 * t1621 * t2418 - 2475 * t1624 * t2418 + 990 * t1630 * t2431 + 4950 * t1634 * t2431 + 3960 * t1637 * t2431 + 495 * t1642 * t2438 - 1980 * t1646 * t2438 - 2475 * t1650 * t2438 - 4320 * t1680 * t1808 + 225 * t1709 * t1808 + 4860 * t1910 * t1808 - 225 * t2243 * t1808 + 225 * t2246 * t1808 - 225 * t2249 * t1808 - 1620 * t2254 * t1808 + 4752 * t9 * t1808 + 12960 * t1900 * t1823 + 990 * t27 * t2445 + 4950 * t32 * t2445 + 3960 * t36 * t2445
    t2503 = P2 * P1
    t2504 = t1617 * t2503
    t2513 = -79200 * t1664 * t2152 - 7425 * t1664 * t2504 - 29700 * t1668 * t2152 - 8910 * t1668 * t2504 - 1485 * t1672 * t2503 + 1485 * t1676 * t2503 + 29700 * t1681 * t2135 - 29700 * t1685 * t2135 - 23760 * t17 * t1808 - 8550 * t1705 * t2135 + 7425 * t22 * t1808 + 25650 * t1900 * t2136 - 1350 * t1903 * t2135 - 8550 * t1907 * t2135 + 4050 * t1910 * t2145 + 25650 * t1913 * t2152 - 1350 * t1916 * t2135 + 5940 * t1924 * t2135 + 4050 * t1929 * t2135 + 15840 * t1932 * t2135
    t2515 = t1616 * t2503
    t2520 = Q1 * t2503
    t2549 = t8 * t114
    t2552 = t7 * t114
    t2559 = -165 * t36 * t114 - 200 * t1649 * t1850 + 8910 * t1681 * t2515 + 7425 * t1685 * t2515 - 7425 * t1690 * t2520 - 8910 * t1694 * t2520 - 1485 * t1699 * t2503 + 1485 * t1702 * t2503 + 8910 * t1706 * t2503 + 7425 * t1710 * t2503 - 200 * t459 * t1836 - 25 * t2217 * t1843 - 1440 * t1869 * t1850 + 2880 * t242 * t1850 - 25 * yuv3 * t1860 + 1440 * t227 * t2187 - 180 * t31 * t2187 + 7200 * t242 * t2195 - 2640 * yuv3 * t2549 + 2640 * t31 * t2552
    t2574 = t1628 * t80
    t2579 = t1617 * t42
    t2586 = t1616 * t42
    t2593 = Q1 * t42
    t2600 = t1653 * t42
    t2605 = -4950 * t1621 * t2579 + 4950 * t1624 * t2205 - 990 * t1624 * t2579 + 2475 * t1630 * t2586 + 1980 * t1634 * t2586 - 495 * t1637 * t2586 - 3960 * t1642 * t2593 + 7200 * t1645 * t2205 - 4950 * t1646 * t2593 - 990 * t1650 * t2593 - 1350 * t1869 * t2195 + 1350 * t1872 * t2200 - 1350 * t1880 * t2205 - 15840 * t53 * t2205 + 1350 * t227 * t2210 - 9900 * t32 * t2574 + 990 * t36 * t2574 - 3960 * t53 * t2579 + 2475 * t27 * t2600 + 1980 * t32 * t2600
    t2618 = yuv3 * t1923
    t2627 = xuv3 * t1923
    t2648 = t8 * t4
    t2651 = 810 * t15 * t1661 + 50 * t459 * t1640 + 75 * t1649 * t1979 + 450 * t1684 * t1661 - 6480 * t223 * t1661 - 450 * t2618 * t1661 + 750 * t2217 * t1713 - 6750 * t227 * t1713 - 300 * t2627 * t1713 - 2700 * t1869 * t1979 + 400 * yuv3 * t1893 + 75 * t459 * t1969 + 1125 * t2217 * t1974 + 1125 * yuv3 * t1984 + 720 * yuv3 * t2235 + 2700 * t227 * t2311 - 14400 * t31 * t2311 - 495 * t36 * t2600 - 2475 * yuv3 * t2648 - 1440 * t31 * t5
    t2692 = -1350 * t15 * t1749 + 2700 * t1680 * t2018 + 900 * t1684 * t1749 + 14850 * t17 * t2018 + 1200 * t1709 * t2018 - 4050 * t223 * t1749 - 900 * t2618 * t1749 + 300 * t2217 * t1779 - 6750 * t227 * t1779 - 750 * t2627 * t1779 + 1350 * t31 * t1779 - 8100 * t1900 * t2034 - 51300 * t1910 * t2018 + 39600 * t22 * t2018 - 1200 * t2243 * t2018 + 1200 * t2246 * t2018 - 1200 * t2249 * t2018 + 17100 * t2254 * t2018 - 2970 * t9 * t2018 + 7920 * t36 * t4
    t2734 = -3960 * t1621 * t2438 - 4950 * t1624 * t2438 - 270 * t1645 * t2438 + 54000 * t1869 * t2095 + 1890 * t1869 * t2418 - 54000 * t1872 * t2078 - 1890 * t1872 * t2431 + 5400 * t1875 * t2078 - 4050 * t1875 * t2431 + 1890 * t1880 * t2438 - 1050 * t2269 * t2078 + 1050 * t2272 * t2078 - 1050 * t2275 * t2078 - 5400 * t2278 * t2078 + 1050 * t2283 * t2078 + 24750 * t27 * t2078 + 49500 * t32 * t2078 - 7920 * t36 * t2078 - 270 * t242 * t2418 + 990 * t53 * t2438
    t2739 = t1628 * t3
    t2779 = -17100 * t1680 * t2135 + 14850 * t17 * t2135 + 1200 * t1709 * t2135 + 51300 * t1900 * t2152 + 6480 * t1900 * t2504 + 8100 * t1910 * t2135 - 14850 * t22 * t2135 + 1200 * t2246 * t2135 - 1200 * t2249 * t2135 - 2700 * t2254 * t2135 + 7920 * t9 * t2135
    t2824 = -90 * t31 * t114 - 1125 * t1649 * t2205 - 14850 * t1664 * t2520 - 17820 * t1668 * t2520 + 17820 * t1681 * t2503 + 14850 * t1685 * t2503 - 2160 * t1705 * t2503 - 400 * t459 * t1850 + 2160 * t1903 * t2503 - 2160 * t1907 * t2503 - 6480 * t1910 * t2515 + 6480 * t1913 * t2520 + 2160 * t1916 * t2503 - 2970 * t1924 * t2503 - 6480 * t1929 * t2503 + 2970 * t1932 * t2503 - 50 * yuv3 * t2187 - 1125 * t459 * t2195 - 75 * t2217 * t2200 + 720 * yuv3 * t2552
    t2833 = t7 * t80
    t2860 = t1628 * t42
    t2867 = -9900 * t1621 * t2593 - 1980 * t1624 * t2593 + 4050 * t1645 * t2593 - 2700 * t1869 * t2205 + 1890 * t1869 * t2579 - 1890 * t1872 * t2586 + 270 * t1875 * t2586 + 1890 * t1880 * t2593 + 14400 * t242 * t2205 - 75 * yuv3 * t2210 + 2700 * t227 * t2574 - 1890 * t227 * t2600 + 4050 * t242 * t2579 - 7920 * t53 * t2593 + 270 * t31 * t2600 + 4950 * t27 * t2860 - 4950 * t31 * t2833 + 3960 * t32 * t2860 - 990 * t36 * t2860 + 495 * t36 * t80
    t2897 = t7 * t4
    t2904 = 225 * xuv3 * t1661 - 330 * t17 * t1671 - 165 * t22 * t1671 + 165 * t27 * t1675 + 330 * t32 * t1675 - 330 * t17 * t1698 - 165 * t22 * t1698 + 375 * yuv3 * t1713 + 450 * xuv3 * t1749 + 150 * t459 * t1979 + 2250 * yuv3 * t2311 + 1350 * yuv3 * t2897 - 7200 * t31 * t4 + 200 * yuv3 * t5 - 165 * t1714 + 165 * t1727 - 165 * t1730 + 165 * t1739 + 330 * t1743 + 165 * t1746
    t2939 = t8 * t3
    t2942 = t7 * t3
    t2947 = 1350 * t15 * t2018 - 60 * t1649 * t2438 + 2400 * t1684 * t2018 + 150 * yuv3 * t1779 + 3780 * t1869 * t2438 - 25650 * t223 * t2018 - 2400 * t2618 * t2018 + 2100 * t2217 * t2078 - 27000 * t227 * t2078 - 2100 * t2627 * t2078 + 2700 * t31 * t2078 + 780 * t2217 * t2431 - 3780 * t227 * t2739 - 60 * t459 * t2418 - 540 * t242 * t2438 + 780 * yuv3 * t2445 - 8100 * t31 * t2739 + 990 * yuv3 * t2939 + 4950 * t31 * t2942 + 3960 * t36 * t3
    t2990 = -25 * yuv3 * t114 - 8550 * t15 * t2135 - 4320 * t1680 * t2503 + 2400 * t1684 * t2135 + 8910 * t17 * t2503 + 720 * t1709 * t2503 + 225 * xuv3 * t1808 + 12960 * t1900 * t2520 - 12960 * t1910 * t2503 + 4050 * t223 * t2135 - 2400 * t2618 * t2135 + 7425 * t22 * t2503 - 2250 * t459 * t2205 - 720 * t2243 * t2503 + 720 * t2246 * t2503 - 720 * t2249 * t2503 + 4320 * t2254 * t2503 + 1485 * t9 * t2503 - 150 * yuv3 * t2574 + 1350 * yuv3 * t2833
    t3007 = t8 * t42
    t3010 = t7 * t42
    t3029 = -780 * t1649 * t2593 + 180 * t223 * t1671 - 180 * t227 * t1675 + 180 * t223 * t1698 - 660 * t17 * t1923 + 3780 * t1869 * t2593 - 330 * t22 * t1923 + 60 * t2217 * t2586 - 3780 * t227 * t2860 + 8100 * t242 * t2593 - 780 * t459 * t2579 + 60 * yuv3 * t2600 + 540 * t31 * t2860 + 2475 * yuv3 * t3007 + 1980 * t31 * t3010 - 495 * t36 * t42 + 180 * t1939 + 180 * t1946 - 330 * t1951 - 180 * t2051
    t3066 = -2160 * t15 * t2503 + 1440 * t1684 * t2503 + 1200 * xuv3 * t2018 + 1050 * yuv3 * t2078 + 1200 * xuv3 * t2135 - 6480 * t223 * t2503 - 120 * t459 * t2438 - 1440 * t2618 * t2503 - 1560 * t459 * t2593 + 1560 * yuv3 * t2739 + 120 * yuv3 * t2860 - 1890 * yuv3 * t2942 - 4050 * t31 * t3 + 1125 * yuv3 * t4 - 75 * yuv3 * t80 + 330 * t1630 + 660 * t1634 + 330 * t1637 - 180 * t1960 - 180 * t2069
    t3095 = 720 * xuv3 * t2503 + 780 * yuv3 * t3 + 60 * yuv3 * t42 + 80 * t2217 - 180 * t227 - 80 * t2627 + 165 * t27 - 180 * t31 + 330 * t32 + 165 * t36 + 40 * yuv3
    t3106 = t863 * t864 / GG
    dQ1dt = 21/256 * t3106 * t2 * zuv3 * (t2104 + t2421 + t2376 + t2334 + t2292 + 40 * t2283 + 360 * t2278 - 40 * t2275 + 40 * t2272 - 40 * t2269 + t2240 + t2192 + t2148 + t1697 + t3095 + t2904 + t1958 + t2470 + t2867 + t2824 + t2779 + t2734 + t1807 + t2692 + t2651 + t2605 + t2559 + t2513 + t2054 + t2947 + t3066 + t2990 + t2005 + t3029 + t1759 + t1906 - 360 * t1875 + 150 * t2283 * t1779 - 13500 * t1872 * t1779 + 2700 * t1875 * t1779 + 39600 * t27 * t1779 - 54450 * t32 * t1779 + 3960 * t36 * t1779 - 2700 * t1939 * t2078 + 27000 * t1869 * t2081 - 27000 * t1872 * t2086 + 2700 * t2051 * t2078 - 1200 * t2243 * t2135 + 15840 * t1924 * t1749 - 1350 * t1705 * t1749 - 4050 * t1929 * t1749 - 15840 * t1932 * t1749 + 188100 * t1681 * t1749 - 79200 * t1685 * t1749 + 14850 * t1664 * t2019 - 14850 * t1668 * t2019 - 7920 * t1672 * t2018 - 2970 * t1676 * t2018 - 1890 * t227 * t2445 - 4050 * t31 * t2445 + 1980 * t27 * t2739 + 9900 * t32 * t2739 + 7920 * t36 * t2739 - 450 * t2618 * t1808 + 450 * t1684 * t1808 - 2160 * t15 * t1808 + 2430 * t223 * t1808 - 1890 * yuv3 * t3010 + 270 * t31 * t42 + 360 * t223 * t1923 + t1855 - 360 * t1872)
    t3108 = t1616 ^ 2
    t3109 = t3108 * t5
    t3132 = t36 * t1641
    t3136 = yuv3 * t8 * t3108
    t3139 = t7 * t3108
    t3140 = t31 * t3139
    t3143 = t36 * t3108
    t3156 = 1320 * t1664 * t1618 - 5280 * t1668 * t1618 + 2112 * t37 * t1618 - 23760 * t1621 * t1662 + 23760 * t1624 * t1662 - 2640 * t1681 * t1627 + 2640 * t1685 * t1627 + 165 * t1932 * t1627 - 7425 * t1630 * t1679 + 23760 * t1634 * t1679 + 1320 * t1690 * t1640 - 5280 * t1694 * t1640 + 2112 * t3132 * t1640 + 1485 * t1714 * t1661 - 7425 * t3136 * t1661 + 23760 * t3140 * t1661 - 4752 * t3143 * t1661 - 2640 * t17 * t3109 + 2640 * t22 * t3109 + 165 * t9 * t3109
    t3165 = t9 * t3108
    t3168 = t15 * t3108
    t3169 = t16 * t3168
    t3172 = xuv3 * t3108
    t3173 = t21 * t3172
    t3202 = -23760 * t1646 * t1688 + 23760 * t1650 * t1688 - 4752 * t1727 * t1661 + 1485 * t1730 * t1661 - 27225 * t1664 * t1717 + 79200 * t1668 * t1717 - 7920 * t1672 * t1713 - 3960 * t1676 * t1713 + 54450 * t1681 * t1722 - 39600 * t1685 * t1722 - 27225 * t1690 * t1733 + 79200 * t1694 * t1733 - 7920 * t1699 * t1713 - 3960 * t3165 * t1713 + 54450 * t3169 * t1713 - 39600 * t3173 * t1713 - 7920 * t1714 * t1749 + 39600 * t3136 * t1749 - 94050 * t3140 * t1749 + 7920 * t3143 * t1749
    t3244 = 94050 * t1621 * t1750 - 39600 * t1624 * t1750 + 39600 * t1630 * t1761 - 94050 * t1634 * t1761 + 94050 * t1646 * t1766 - 39600 * t1650 * t1766 + 39600 * t1664 * t1782 - 54450 * t1668 * t1782 + 3960 * t1672 * t1779 + 7920 * t1676 * t1779 - 79200 * t1681 * t1787 + 27225 * t1685 * t1787 + 39600 * t1690 * t1796 - 54450 * t1694 * t1796 + 3960 * t1699 * t1779 + 7920 * t1727 * t1749 - 7920 * t1730 * t1749 + 7920 * t3165 * t1779 - 79200 * t3169 * t1779 + 27225 * t3173 * t1779
    t3269 = t3108 * t114
    t3286 = -23760 * t1621 * t1809 + 7425 * t1624 * t1809 - 23760 * t1630 * t1818 + 23760 * t1634 * t1818 - 23760 * t1646 * t1823 + 7425 * t1650 * t1823 - 2640 * t1664 * t1836 + 2640 * t1668 * t1836 + 5280 * t1681 * t1843 + 5280 * t17 * t3269 + 4752 * t1714 * t1808 - 1485 * t1727 * t1808 + 4752 * t1730 * t1808 - 23760 * t3136 * t1808 + 23760 * t3140 * t1808 - 1485 * t3143 * t1808 - 165 * t37 * t1836 - 2112 * t1932 * t1843 - 1320 * t22 * t3269 - 2112 * t9 * t3269
    t3303 = t31 * Q2
    t3318 = t31 * t1641
    t3327 = yuv3 * t3139
    t3330 = t31 * t3108
    t3333 = 90 * t15 * t3109 + 720 * t1900 * t1618 - 1440 * t3303 * t1618 + 90 * t1680 * t1627 - 5280 * t17 * t1627 - 720 * t1910 * t1627 + 5280 * t22 * t1627 + 330 * t9 * t1627 + 2640 * t1664 * t1640 - 10560 * t1668 * t1640 + 720 * t1913 * t1640 - 1440 * t3318 * t1640 + 4224 * t37 * t1640 - 2430 * t3327 * t1661 + 2160 * t3330 * t1661 - 1320 * t1685 * t1843 - 2640 * t1690 * t1850 + 2640 * t1694 * t1850 - 165 * t3132 * t1850 - 720 * t223 * t3109
    t3343 = yuv3 * t8 * t1616
    t3346 = t7 * t1616
    t3347 = t31 * t3346
    t3350 = t36 * t1616
    t3365 = t16 * t3172
    t3374 = t9 * t1616
    t3377 = t15 * t1616
    t3378 = t16 * t3377
    t3381 = xuv3 * t1616
    t3382 = t21 * t3381
    t3385 = -47520 * t1621 * t1688 + 47520 * t1624 * t1688 + 810 * t1939 * t1661 + 810 * t1946 * t1661 + 2970 * t1951 * t1661 + 2160 * t2051 * t1661 - 14850 * t3343 * t1661 + 47520 * t3347 * t1661 - 9504 * t3350 * t1661 - 6480 * t1869 * t1662 - 2430 * t1872 * t1679 - 6480 * t1880 * t1688 - 1350 * t1907 * t1713 - 1350 * t3168 * t1713 + 6750 * t3365 * t1713 - 7920 * t3374 * t1713 + 108900 * t3378 * t1713 - 79200 * t3382 * t1713 - 6750 * t1900 * t1717 + 6750 * t1910 * t1722
    t3395 = t3108 * t4
    t3428 = -54450 * t1664 * t1733 - 2475 * t1664 * t1969 + 158400 * t1668 * t1733 + 4950 * t1681 * t1974 - 2475 * t1690 * t1979 + 4950 * t17 * t3395 - 15840 * t1924 * t1713 - 6750 * t1913 * t1733 - 1350 * t1939 * t1749 + 1350 * t2051 * t1749 + 4050 * t3327 * t1749 + 1350 * t3330 * t1749 + 79200 * t3343 * t1749 - 188100 * t3347 * t1749 - 4050 * t1869 * t1750 + 4050 * t1872 * t1761 - 495 * t1932 * t1974 + 7920 * t37 * t1969 + 7920 * t3132 * t1979 - 495 * t9 * t3395
    t3472 = 14850 * t1621 * t2019 + 39600 * t1624 * t2019 + 14850 * t1630 * t2029 - 14850 * t1634 * t2029 + 14850 * t1646 * t2034 + 39600 * t1650 * t2034 - 7920 * t1727 * t2018 - 2970 * t1730 * t2018 + 1350 * t1903 * t1779 + 6750 * t3365 * t1779 - 6750 * t1900 * t1782
    t3517 = 79200 * t1664 * t1796 + 24750 * t1664 * t2081 - 108900 * t1668 * t1796 + 49500 * t1668 * t2081 - 7920 * t1672 * t2078 + 7920 * t1676 * t2078 - 49500 * t1681 * t2086 - 24750 * t1685 * t2086 + 24750 * t1690 * t2095 + 49500 * t1694 * t2095 + 1350 * t1916 * t1779 + 7920 * t1924 * t1779 + 15840 * t3374 * t1779 - 158400 * t3378 * t1779 + 54450 * t3382 * t1779 + 6750 * t1910 * t1787 - 6750 * t1913 * t1796 + 7920 * t3165 * t2078 - 49500 * t3169 * t2078 - 24750 * t3173 * t2078
    t3558 = -47520 * t1621 * t1823 + 14850 * t1621 * t2136 + 14850 * t1624 * t1823 - 7920 * t1699 * t2078 + 7920 * t1714 * t2135 - 2160 * t1939 * t1808 - 2160 * t1946 * t1808 + 9504 * t1951 * t1808 - 810 * t2051 * t1808 + 6480 * t3327 * t1808 - 810 * t3330 * t1808 - 47520 * t3343 * t1808 + 47520 * t3347 * t1808 - 2970 * t3350 * t1808 + 2430 * t1869 * t1809 + 6480 * t1872 * t1818 + 2430 * t1880 * t1823 - 39600 * t3136 * t2135 - 14850 * t3140 * t2135 + 2970 * t3143 * t2135
    t3600 = 1440 * t15 * t3269 - 14850 * t1624 * t2136 - 39600 * t1630 * t2145 - 14850 * t1634 * t2145 + 14850 * t1646 * t2152 - 14850 * t1650 * t2152 - 5280 * t1664 * t1850 + 5280 * t1668 * t1850 + 1440 * t1680 * t1843 + 10560 * t17 * t1843 + 2970 * t1727 * t2135 + 7920 * t1730 * t2135 + 720 * t1900 * t1836 - 90 * t3303 * t1836 - 720 * t1910 * t1843 - 2640 * t22 * t1843 - 4224 * t9 * t1843 + 720 * t1913 * t1850 - 90 * t3318 * t1850 - 720 * t223 * t3269
    t3603 = t3108 * t80
    t3630 = yuv3 * t1641
    t3639 = t15 * t5
    t3642 = xuv3 * t5
    t3645 = 180 * t15 * t1627 - 2640 * t16 * t3639 + 200 * t451 * t1618 + 25 * t1684 * t1627 - 1440 * t223 * t1627 + 1440 * t1900 * t1640 - 2880 * t3303 * t1640 + 200 * t3630 * t1640 - 4950 * t1668 * t2195 + 2475 * t1685 * t2200 - 4950 * t1694 * t2205 - 330 * t37 * t1850 - 7920 * t1932 * t2200 + 2640 * t21 * t3642 + 495 * t37 * t2195 + 2475 * t22 * t3603 + 495 * t3132 * t2205 + 25 * xuv3 * t3109 - 7920 * t9 * t3603 + 165 * t9 * t5
    t3648 = yuv3 * t3108
    t3655 = yuv3 * t3346
    t3658 = t31 * t1616
    t3681 = t16 * t3381
    t3692 = 225 * t2269 * t1661 - 225 * t2272 * t1661 + 225 * t2275 * t1661 + 1620 * t2278 * t1661 - 7425 * t27 * t1661 + 23760 * t32 * t1661 - 4752 * t36 * t1661 - 225 * t3648 * t1661 - 4860 * t3655 * t1661 + 4320 * t3658 * t1661 - 12960 * t1869 * t1688 + 54450 * t17 * t1713 + 375 * t2243 * t1713 - 150 * t2246 * t1713 + 375 * t2249 * t1713 - 150 * t3172 * t1713 - 2700 * t3377 * t1713 + 13500 * t3681 * t1713 - 3960 * t9 * t1713 - 13500 * t1900 * t1733
    t3733 = -4950 * t1664 * t1979 + 9900 * t17 * t1974 - 39600 * t22 * t1713 + 450 * t2269 * t1749 - 450 * t2272 * t1749 + 450 * t2275 * t1749 - 2700 * t2278 * t1749 + 39600 * t27 * t1749 - 450 * t3648 * t1749 + 8100 * t3655 * t1749 + 2700 * t3658 * t1749 - 8100 * t1869 * t1766 + 1350 * t1900 * t1969 - 1350 * t1910 * t1974 + 1350 * t1913 * t1979 - 7200 * t3303 * t1969 - 990 * t9 * t1974 - 7200 * t3318 * t1979 + 15840 * t37 * t1979 - 1350 * t223 * t3395
    t3775 = 29700 * t1621 * t2034 + 79200 * t1624 * t2034 - 94050 * t32 * t1749 + 7920 * t36 * t1749 + 150 * t2243 * t1779 - 375 * t2246 * t1779 - 375 * t3172 * t1779 + 13500 * t3681 * t1779 - 25650 * t1869 * t2019 - 4050 * t1872 * t2029 - 25650 * t1880 * t2034 + 1350 * t1939 * t2018 + 1350 * t1946 * t2018 - 5940 * t1951 * t2018 + 8550 * t2051 * t2018 - 4050 * t3327 * t2018 + 8550 * t3330 * t2018 + 29700 * t3343 * t2018 - 29700 * t3347 * t2018 - 15840 * t3350 * t2018
    t3817 = t3108 * t3
    t3820 = 49500 * t1664 * t2095 + 99000 * t1668 * t2095 - 2700 * t1907 * t2078 + 27000 * t1910 * t2086 - 27000 * t1913 * t2095 + 2700 * t1916 * t2078 - 15840 * t1924 * t2078 + 15840 * t3374 * t2078 - 99000 * t3378 * t2078 - 49500 * t3382 * t2078 + 495 * t9 * t3817
    t3866 = 990 * t1664 * t2418 + 4950 * t1668 * t2418 - 1980 * t1681 * t2431 - 2475 * t1685 * t2431 + 990 * t1690 * t2438 + 4950 * t1694 * t2438 - 1980 * t17 * t3817 + 225 * t2269 * t1808 - 225 * t2272 * t1808 + 225 * t2275 * t1808 - 4320 * t2278 * t1808 - 23760 * t27 * t1808 - 225 * t3648 * t1808 + 12960 * t3655 * t1808 - 1620 * t3658 * t1808 + 4860 * t1869 * t1823 + 495 * t1932 * t2431 - 2475 * t22 * t3817 + 3960 * t37 * t2418 + 3960 * t3132 * t2438
    t3907 = 29700 * t1621 * t2152 - 29700 * t1624 * t2152 + 1485 * t1714 * t2503 + 23760 * t32 * t1808 - 1485 * t36 * t1808 + 4050 * t1869 * t2136 + 25650 * t1872 * t2145 + 4050 * t1880 * t2152 - 8550 * t1939 * t2135 - 8550 * t1946 * t2135 + 15840 * t1951 * t2135 - 1350 * t2051 * t2135 + 25650 * t3327 * t2135 - 1350 * t3330 * t2135 - 79200 * t3343 * t2135 - 29700 * t3347 * t2135 + 5940 * t3350 * t2135 - 7425 * t3136 * t2503 - 8910 * t3140 * t2503 - 1485 * t3143 * t2503
    t3943 = t15 * t114
    t3946 = xuv3 * t114
    t3951 = -2112 * t9 * t114 + 2880 * t15 * t1843 + 7200 * t15 * t3603 + 5280 * t16 * t3943 + 8910 * t1621 * t2504 + 7425 * t1624 * t2504 - 7425 * t1630 * t2515 - 8910 * t1634 * t2515 + 8910 * t1646 * t2520 + 7425 * t1650 * t2520 - 200 * t1684 * t1843 - 1485 * t1727 * t2503 + 1485 * t1730 * t2503 - 25 * t451 * t1836 - 1440 * t223 * t1843 + 1440 * t1900 * t1850 - 180 * t3303 * t1850 - 25 * t3630 * t1850 - 1320 * t21 * t3946 - 200 * xuv3 * t3269
    t3970 = t3108 * t42
    t3993 = 2475 * t1664 * t2579 - 9900 * t1668 * t2205 + 1980 * t1668 * t2579 + 7200 * t1680 * t2200 - 4950 * t1681 * t2586 - 990 * t1685 * t2586 + 2475 * t1690 * t2593 + 1980 * t1694 * t2593 - 4950 * t17 * t3970 + 1350 * t1900 * t2195 - 1350 * t1910 * t2200 + 1350 * t1913 * t2205 - 3960 * t1932 * t2586 + 4950 * t22 * t2200 - 990 * t22 * t3970 - 15840 * t9 * t2200 + 990 * t37 * t2205 - 1350 * t223 * t3603 - 495 * t37 * t2579 - 3960 * t9 * t3970
    t4005 = yuv3 * t1616
    t4036 = -1350 * t15 * t1713 - 720 * t16 * t3642 + 50 * xuv3 * t1627 + 400 * t451 * t1640 - 2430 * t227 * t1661 + 450 * t2627 * t1661 + 2160 * t31 * t1661 - 450 * t4005 * t1661 + 75 * t1684 * t1974 + 6750 * t223 * t1713 + 750 * t2618 * t1713 - 300 * t3381 * t1713 + 2700 * t1900 * t1979 + 1125 * t451 * t1969 - 2700 * t223 * t1974 - 14400 * t3303 * t1979 + 1125 * t3630 * t1979 - 495 * t3132 * t2593 + 75 * xuv3 * t3395 + 90 * t3639
    t4078 = 4950 * t16 * t15 * t4 + 4050 * t227 * t1749 + 900 * t2627 * t1749 + 1350 * t31 * t1749 - 900 * t4005 * t1749 + 6750 * t223 * t1779 + 300 * t2618 * t1779 - 750 * t3381 * t1779 - 51300 * t1869 * t2034 + 1200 * t2269 * t2018 - 1200 * t2272 * t2018 + 1200 * t2275 * t2018 + 2700 * t2278 * t2018 + 14850 * t27 * t2018 - 14850 * t32 * t2018 - 7920 * t36 * t2018 - 1200 * t3648 * t2018 - 8100 * t3655 * t2018 + 17100 * t3658 * t2018 - 495 * t9 * t4
    t4120 = -270 * t15 * t3817 - 270 * t1680 * t2431 - 49500 * t17 * t2078 - 3960 * t17 * t2431 - 54000 * t1900 * t2095 - 1890 * t1900 * t2418 + 1890 * t1910 * t2431 - 24750 * t22 * t2078 + 1050 * t2243 * t2078 - 1050 * t2246 * t2078 + 1050 * t2249 * t2078 + 5400 * t2254 * t2078 - 1050 * t3172 * t2078 - 5400 * t3377 * t2078 + 54000 * t3681 * t2078 + 7920 * t9 * t2078 - 4950 * t22 * t2431 + 1890 * t223 * t3817 - 4050 * t3303 * t2418 + 990 * t9 * t2431
    t4164 = 8100 * t1869 * t2152 + 1200 * t2269 * t2135 - 1200 * t2272 * t2135 + 1200 * t2275 * t2135 - 17100 * t2278 * t2135 - 39600 * t27 * t2135 - 14850 * t32 * t2135 + 2970 * t36 * t2135 + 51300 * t3655 * t2135 - 2700 * t3658 * t2135 + 6480 * t3327 * t2503
    t4208 = -720 * t16 * t3946 + 17820 * t1621 * t2520 + 14850 * t1624 * t2520 - 1125 * t1684 * t2200 - 400 * xuv3 * t1843 - 50 * t451 * t1850 - 6480 * t1869 * t2504 + 6480 * t1872 * t2515 - 6480 * t1880 * t2520 - 2160 * t1939 * t2503 - 2160 * t1946 * t2503 + 2970 * t1951 * t2503 + 2160 * t2051 * t2503 - 75 * t451 * t2195 + 2160 * t3330 * t2503 - 14850 * t3343 * t2503 - 17820 * t3347 * t2503 - 2970 * t3350 * t2503 - 1125 * xuv3 * t3603 + 1440 * t3943
    t4219 = xuv3 * t80
    t4250 = 14400 * t15 * t2200 + 4050 * t15 * t3970 + 4950 * t1664 * t2593 + 3960 * t1668 * t2593 + 4050 * t1680 * t2586 - 9900 * t17 * t2586 + 2700 * t1900 * t2205 - 1890 * t1900 * t2579 + 1890 * t1910 * t2586 - 1890 * t1913 * t2593 + 2475 * t21 * t4219 - 1980 * t22 * t2586 - 2700 * t223 * t2200 - 75 * t3630 * t2205 + 1890 * t223 * t3970 + 270 * t3303 * t2579 - 7920 * t9 * t2586 + 270 * t3318 * t2593 - 990 * t37 * t2593 - 7920 * t9 * t80
    t4279 = xuv3 * t4
    t4286 = -1350 * t16 * t4279 - 225 * yuv3 * t1661 + 165 * t27 * t1671 + 330 * t32 * t1671 - 330 * t17 * t1675 - 165 * t22 * t1675 + 165 * t27 * t1698 + 330 * t32 * t1698 - 150 * xuv3 * t1713 - 450 * yuv3 * t1749 + 150 * xuv3 * t1974 + 2250 * t451 * t1979 - 2400 * t4005 * t2018 + 165 * t1672 - 165 * t1676 + 165 * t1699 - 165 * t3165 - 330 * t3169 - 165 * t3173 + 25 * t3642
    t4321 = t15 * t3
    t4324 = xuv3 * t3
    t4329 = -2700 * t15 * t2078 - 540 * t15 * t2431 - 1980 * t16 * t4321 - 60 * t1684 * t2431 - 375 * xuv3 * t1779 - 225 * yuv3 * t1808 - 3780 * t1900 * t2438 - 4050 * t227 * t2018 + 2400 * t2627 * t2018 + 8550 * t31 * t2018 + 27000 * t223 * t2078 + 2100 * t2618 * t2078 - 2100 * t3381 * t2078 - 2475 * t21 * t4324 + 3780 * t223 * t2431 + 780 * t451 * t2418 - 8100 * t3303 * t2438 + 780 * t3630 * t2438 + 495 * t9 * t3 - 60 * xuv3 * t3817
    t4371 = 7200 * t15 * t80 - 1350 * t16 * t4219 - 12960 * t1869 * t2520 + 25650 * t227 * t2135 + 2400 * t2627 * t2135 - 1350 * t31 * t2135 - 2400 * t4005 * t2135 - 2250 * xuv3 * t2200 - 150 * t451 * t2205 + 720 * t2269 * t2503 - 720 * t2272 * t2503 + 720 * t2275 * t2503 - 4320 * t2278 * t2503 - 7425 * t27 * t2503 - 8910 * t32 * t2503 - 1485 * t36 * t2503 - 720 * t3648 * t2503 + 12960 * t3655 * t2503 + 4320 * t3658 * t2503 - 200 * t3946
    t4390 = t15 * t42
    t4393 = xuv3 * t42
    t4407 = 8100 * t15 * t2586 - 4950 * t16 * t4390 - 180 * t227 * t1671 + 180 * t223 * t1675 - 780 * t1684 * t2586 - 3780 * t1900 * t2593 - 990 * t21 * t4393 + 3780 * t223 * t2586 + 60 * t451 * t2579 + 540 * t3303 * t2593 + 60 * t3630 * t2593 - 780 * xuv3 * t3970 - 3960 * t9 * t42 - 180 * t1903 + 180 * t1907 + 180 * t3168 + 180 * t3365 - 330 * t3374 - 660 * t3378 - 330 * t3382
    t4444 = 1890 * t16 * t4324 - 180 * t227 * t1698 + 330 * t27 * t1923 + 660 * t32 * t1923 - 1200 * yuv3 * t2018 - 1050 * xuv3 * t2078 - 1200 * yuv3 * t2135 + 6480 * t227 * t2503 - 120 * xuv3 * t2431 + 1560 * t451 * t2438 + 1440 * t2627 * t2503 + 2160 * t31 * t2503 - 1440 * t4005 * t2503 - 1560 * xuv3 * t2586 + 120 * t451 * t2593 - 180 * t1916 + 330 * t1924 - 1125 * t4219 + 75 * t4279 - 270 * t4321
    t4470 = -720 * yuv3 * t2503 + 180 * t15 - 330 * t17 - 165 * t22 + 180 * t223 + 80 * t2618 - 80 * t3381 - 60 * t4324 - 780 * t4393 - 165 * t9 - 40 * xuv3
    dQ2dt = -21/256 * t3106 * zuv3 * t2 * (t4078 - 360 * t2254 + 40 * t2249 - 40 * t2246 + 40 * t2243 + t3600 + t4444 + t3993 + 4050 * t4390 + t4120 + t4036 + t3385 + 7920 * t9 * t1779 - 79200 * t17 * t1779 + 27225 * t22 * t1779 - 2700 * t3168 * t2078 + 27000 * t3365 * t2078 - 27000 * t1900 * t2081 + 2700 * t1903 * t2078 - 1890 * t1913 * t2438 - 4050 * t3318 * t2438 + 1980 * t1664 * t2438 + 9900 * t1668 * t2438 + 7920 * t37 * t2438 - 450 * t4005 * t1808 + 450 * t2627 * t1808 + 6480 * t227 * t1808 - 810 * t31 * t1808 - 1200 * t3648 * t2135 + 1890 * t16 * t4393 - 360 * t227 * t1923 + t3156 + t3558 + t3775 + t3820 + t3472 - 40 * t3172 + 360 * t3377 + t3951 + t4286 + t3333 + 15840 * t3350 * t1749 - 1350 * t1946 * t1749 - 4050 * t1880 * t1766 - 15840 * t1951 * t1749 + 188100 * t1621 * t1766 - 79200 * t1624 * t1766 + 14850 * t3136 * t2018 - 14850 * t3140 * t2018 - 7920 * t3143 * t2018 - 2970 * t1714 * t2018 + 150 * t2249 * t1779 - 13500 * t1900 * t1796 + 2700 * t2254 * t1779 + t4164 + t4208 + t3645 + t3866 + t4371 + t3733 + t3692 + 360 * t3681 + t4250 + t3244 + t3517 + t3907 + t4329 + t3428 + t3202 + t4407 + t3286 + t4470)
    t4481 = t9 * eta
    t4482 = zuv3 * t4481
    t4485 = eta * t1640
    t4490 = Q2 * t5
    t4491 = eta * t4490
    t4497 = zuv3 * t36 * eta
    t4500 = t8 * eta
    t4501 = t550 * t4500
    t4504 = t7 * eta
    t4505 = t231 * t4504
    t4510 = Q2 * t1661
    t4513 = t15 * eta
    t4514 = zuv3 * t16
    t4515 = t4514 * t4513
    t4518 = eta * xuv3
    t4519 = zuv3 * t21
    t4520 = t4519 * t4518
    t4529 = Q2 * t1713
    t4538 = 1155 * t4482 * t1640 - 33264 * t4497 * t1688 - 51975 * t4501 * t1688 + 166320 * t4505 * t1688 - 27720 * t4482 * t1733 + 381150 * t4515 * t1733 - 277200 * t4520 * t1733 + 277200 * t4501 * t1766 - 18480 * t18 * t4485 + 18480 * t23 * t4485 + 9240 * t28 * t4491 - 36960 * t33 * t4491 + 10395 * t4482 * t4510 + 14784 * t4497 * t4490 - 55440 * t4497 * t4529 - 190575 * t4501 * t4529 + 554400 * t4505 * t4529 - 166320 * t4515 * t4510 + 166320 * t4520 * t4510
    t4543 = Q2 * t1749
    t4556 = Q2 * t1779
    t4569 = Q2 * t1808
    t4578 = eta * t1850
    t4583 = 55440 * t4497 * t1766 - 658350 * t4505 * t1766 + 55440 * t4482 * t1796 - 554400 * t4515 * t1796 + 190575 * t4520 * t1796 + 36960 * t18 * t4578 - 10395 * t4497 * t1823 - 166320 * t4501 * t1823 + 166320 * t4505 * t1823 - 14784 * t4482 * t1850 - 9240 * t23 * t4578 - 55440 * t4482 * t4543 + 33264 * t4482 * t4569 + 27720 * t4497 * t4556 + 277200 * t4501 * t4556 - 381150 * t4505 * t4556 + 658350 * t4515 * t4543 - 166320 * t4515 * t4569 - 277200 * t4520 * t4543 + 51975 * t4520 * t4569
    t4585 = Q2 * t114
    t4586 = eta * t4585
    t4605 = eta * t5
    t4626 = 1155 * t133 * t1640 + 14784 * t145 * t4490 + 1155 * t161 * t4605 - 27720 * t164 * t4605 - 18480 * t18 * t1640 + 18480 * t23 * t1640 - 33264 * t50 * t1661 + 10395 * t54 * t1661 + 55440 * t167 * t4605 - 51975 * t28 * t1688 + 166320 * t33 * t1688 - 14784 * t170 * t4605 - 166320 * t18 * t4510 + 166320 * t23 * t4510 + 9240 * t28 * t4490 - 18480 * t28 * t4586 - 36960 * t33 * t4490 + 18480 * t33 * t4586 - 1155 * t4497 * t4585
    t4627 = yuv3 * t4481
    t4630 = t31 * t4513
    t4633 = t36 * t4518
    t4648 = t161 * eta
    t4651 = t16 * t4500
    t4654 = t21 * t4504
    t4657 = t170 * eta
    t4674 = -27720 * t11 * t1713 - 62370 * t4627 * t1661 + 332640 * t4630 * t1661 - 199584 * t4633 * t1661 - 55440 * t38 * t1713 - 27720 * t4648 * t1713 + 571725 * t4651 * t1713 - 831600 * t4654 * t1713 + 55440 * t4657 * t1713 + 381150 * t18 * t1733 - 277200 * t23 * t1733 + 332640 * t4627 * t1749 + 55440 * t50 * t1749 - 55440 * t54 * t1749 + 277200 * t28 * t1766 - 658350 * t33 * t1766 + 658350 * t18 * t4543 - 277200 * t23 * t4543 - 190575 * t28 * t4529 + 554400 * t33 * t4529
    t4715 = 55440 * t11 * t1779 - 1316700 * t4630 * t1749 + 332640 * t4633 * t1749 + 27720 * t38 * t1779 + 55440 * t4648 * t1779 - 831600 * t4651 * t1779 + 571725 * t4654 * t1779 - 27720 * t4657 * t1779 - 554400 * t18 * t1796 + 190575 * t23 * t1796 - 166320 * t18 * t4569 - 199584 * t4627 * t1808 - 10395 * t50 * t1808 + 33264 * t54 * t1808 - 166320 * t28 * t1823 + 166320 * t33 * t1823 + 51975 * t23 * t4569 + 277200 * t28 * t4556 - 381150 * t33 * t4556
    t4732 = eta * t114
    t4741 = zuv3 * t4513
    t4749 = zuv3 * t31 * eta
    t4758 = t550 * t4504
    t4761 = -14784 * t133 * t1850 - 1155 * t145 * t4585 - 17325 * t16 * t2232 - 14784 * t161 * t4732 + 1155 * t161 * t5 + 55440 * t164 * t4732 + 630 * t4741 * t1640 - 27720 * t167 * t4732 - 17010 * t4758 * t1688 + 1155 * t170 * t4732 + 18480 * t170 * t5 + 36960 * t18 * t1850 + 332640 * t4630 * t1808 - 62370 * t4633 * t1808 - 9240 * t23 * t1850 - 5040 * t224 * t4485 + 5040 * t228 * t4491 - 18480 * t28 * t4585 + 18480 * t33 * t4585 - 10080 * t4749 * t4490
    t4767 = t4514 * t4518
    t4790 = eta * t1979
    t4793 = Q2 * t4
    t4794 = eta * t4793
    t4805 = -17325 * t161 * t1713 + 155925 * t164 * t1713 - 41580 * t173 * t1661 + 69300 * t176 * t1661 + 110880 * t179 * t1661 + 173250 * t167 * t1713 + 15120 * t4749 * t1688 - 9450 * t4741 * t1733 + 47250 * t4767 * t1733 + 9450 * t4749 * t1766 + 28350 * t4758 * t1766 + 34650 * t18 * t4790 - 3465 * t4482 * t1979 - 17325 * t28 * t4794 + 55440 * t4497 * t4793 + 5670 * t4741 * t4510 - 45360 * t4767 * t4510 - 47250 * t4758 * t4529 - 9450 * t4741 * t4543
    t4820 = Q2 * t2018
    t4845 = Q2 * t2078
    t4848 = 173250 * t164 * t1779 + 155925 * t167 * t1779 - 17325 * t170 * t1779 + 69300 * t173 * t1749 + 138600 * t176 * t1749 + 69300 * t179 * t1749 + 47250 * t4767 * t1796 - 55440 * t4497 * t2034 + 103950 * t4501 * t2034 - 103950 * t4505 * t2034 + 55440 * t4482 * t2095 - 346500 * t4515 * t2095 - 173250 * t4520 * t2095 - 20790 * t4482 * t4820 + 173250 * t4501 * t4845 + 103950 * t4515 * t4820 + 277200 * t4520 * t4820 - 28350 * t4767 * t4543 + 9450 * t4749 * t4556 - 47250 * t4758 * t4556
    t4876 = Q2 * t2135
    t4891 = 110880 * t173 * t1808 + 69300 * t176 * t1808 - 41580 * t179 * t1808 - 5670 * t4749 * t1823 + 45360 * t4758 * t1823 + 10080 * t4741 * t1850 + 20790 * t4497 * t2152 - 277200 * t4501 * t2152 - 103950 * t4505 * t2152 - 5040 * t224 * t4578 + 5040 * t228 * t4586 + 55440 * t4482 * t4876 - 55440 * t4497 * t4845 + 346500 * t4505 * t4845 + 103950 * t4515 * t4876 - 103950 * t4520 * t4876 - 15120 * t4741 * t4569 + 17010 * t4767 * t4569 - 630 * t4749 * t4585
    t4900 = eta * t2205
    t4903 = Q2 * t80
    t4904 = eta * t4903
    t4931 = yuv3 * t4513
    t4934 = t31 * t4518
    t4937 = 18480 * t161 * t114 + 1155 * t170 * t114 - 5040 * t224 * t1640 + 630 * t339 * t1640 + 15120 * t239 * t1661 + 5670 * t243 * t1661 - 34020 * t4931 * t1661 + 90720 * t4934 * t1661 - 17010 * t228 * t1688 - 17325 * t21 * t2552 + 15120 * t21 * t4605 - 55440 * t4482 * t2205 - 45360 * t224 * t4510 + 5040 * t228 * t4490 + 17325 * t23 * t4900 - 10080 * t231 * t4490 - 34650 * t33 * t4904 - 15120 * t365 * t4605 + 3465 * t4497 * t4903 + 945 * t8 * t4605
    t4947 = t16 * t4504
    t4958 = eta * t4
    t4979 = -3465 * t133 * t1979 + 55440 * t145 * t4793 - 3465 * t161 * t4958 + 51975 * t164 * t4958 - 55440 * t170 * t4958 - 9450 * t220 * t1713 - 14175 * t4500 * t1713 + 141750 * t4947 * t1713 + 47250 * t224 * t1733 + 9450 * t239 * t1749 - 9450 * t243 * t1749 + 56700 * t4931 * t1749 + 56700 * t4934 * t1749 + 28350 * t228 * t1766 + 34650 * t18 * t1979 + 103950 * t28 * t2034 - 28350 * t224 * t4543 - 47250 * t228 * t4529 - 17325 * t28 * t4793
    t5004 = t21 * eta
    t5021 = 55440 * t11 * t2078 + 9450 * t232 * t1779 + 141750 * t4947 * t1779 - 14175 * t5004 * t1779 + 47250 * t224 * t1796 - 346500 * t18 * t2095 + 103950 * t18 * t4820 + 124740 * t4627 * t2018 - 207900 * t4630 * t2018 - 332640 * t4633 * t2018 - 55440 * t50 * t2018 - 20790 * t54 * t2018 - 103950 * t33 * t2034 - 55440 * t38 * t2078 + 55440 * t4648 * t2078 - 173250 * t23 * t2095 - 47250 * t228 * t4556 + 277200 * t23 * t4820 + 173250 * t28 * t4845 + 346500 * t33 * t4845
    t5062 = 103950 * t18 * t4876 - 5670 * t239 * t1808 - 15120 * t243 * t1808 + 90720 * t4931 * t1808 - 34020 * t4934 * t1808 + 45360 * t228 * t1823 + 10080 * t339 * t1850 - 519750 * t4651 * t2078 - 519750 * t4654 * t2078 + 55440 * t4657 * t2078 - 332640 * t4627 * t2135 - 207900 * t4630 * t2135 + 124740 * t4633 * t2135 + 20790 * t50 * t2135 + 55440 * t54 * t2135 - 277200 * t28 * t2152 - 103950 * t33 * t2152 + 17010 * t224 * t4569 - 103950 * t23 * t4876
    t5083 = eta * t80
    t5090 = zuv3 * t4518
    t5094 = eta * yuv3 * zuv3
    t5107 = -55440 * t133 * t2205 + 3465 * t145 * t4903 + 9450 * t16 * t2235 - 55440 * t161 * t5083 + 175 * t5090 * t1640 + 18900 * t368 * t1661 + 51975 * t167 * t5083 - 1575 * t5094 * t1688 - 3465 * t170 * t5083 - 5040 * t224 * t1850 + 945 * t21 * t4732 - 25200 * t21 * t5 + 17325 * t23 * t2205 + 5040 * t228 * t4585 - 630 * t231 * t4585 - 34650 * t33 * t4903 - 15120 * t365 * t4732 + 1400 * t5094 * t4490 + 1575 * t5090 * t4510 + 15120 * t8 * t4732
    t5147 = 31185 * t16 * t2648 - 3465 * t161 * t4 - 119700 * t371 * t1661 - 15750 * t21 * t1713 - 179550 * t365 * t1713 + 9450 * t8 * t1713 - 1050 * t5090 * t1733 - 100800 * t368 * t1749 - 100800 * t371 * t1749 - 3150 * t5094 * t1766 + 59850 * t4749 * t2034 - 28350 * t4758 * t2034 + 34650 * t21 * t2897 - 9450 * t224 * t4790 + 9450 * t228 * t4794 + 2625 * t5094 * t4529 + 3150 * t5090 * t4543 + 9450 * t4741 * t4820 - 50400 * t4749 * t4793
    t5184 = eta * t2438
    t5189 = 34650 * t161 * t2078 - 124740 * t164 * t2078 - 124740 * t167 * t2078 + 34650 * t170 * t2078 + 83160 * t173 * t2018 + 13860 * t176 * t2018 + 9450 * t21 * t1779 - 179550 * t365 * t1779 - 15750 * t8 * t1779 - 69300 * t179 * t2018 - 2625 * t5090 * t1796 - 13860 * t18 * t5184 - 18900 * t4741 * t2095 + 189000 * t4767 * t2095 - 17325 * t23 * t5184 + 3465 * t4482 * t2438 + 1050 * t5094 * t4556 + 18900 * t4749 * t4845 - 189000 * t4758 * t4845 - 179550 * t4767 * t4820
    t5194 = Q2 * t3
    t5195 = eta * t5194
    t5230 = Q2 * t2503
    t5235 = -69300 * t173 * t2135 + 13860 * t176 * t2135 + 83160 * t179 * t2135 - 119700 * t368 * t1808 + 18900 * t371 * t1808 - 1575 * t5094 * t1823 - 9450 * t4749 * t2152 + 179550 * t4758 * t2152 - 10395 * t4497 * t2520 - 51975 * t4501 * t2520 - 62370 * t4505 * t2520 + 6930 * t28 * t5195 + 34650 * t33 * t5195 + 10395 * t4482 * t5230 + 27720 * t4497 * t5194 + 62370 * t4515 * t5230 + 1575 * t5090 * t4569 - 59850 * t4741 * t4876 + 28350 * t4767 * t4876
    t5251 = t8 * t80
    t5260 = eta * t2593
    t5265 = Q2 * t42
    t5266 = eta * t5265
    t5279 = 51975 * t4520 * t5230 - 1400 * t5090 * t1850 - 175 * t5094 * t4585 - 25200 * t2549 + 9450 * t16 * t2552 + 50400 * t4741 * t2205 - 9450 * t224 * t4900 + 9450 * t228 * t4904 + 34650 * t16 * t5251 + 31185 * t21 * t2833 - 3465 * t170 * t80 - 27720 * t4482 * t2593 - 34650 * t18 * t5260 - 6930 * t23 * t5260 + 17325 * t28 * t5266 + 13860 * t33 * t5266 - 3465 * t4497 * t5265 + 175 * t547 * t1640 + 1400 * t550 * t4490 + 525 * t7 * t4605
    t5287 = yuv3 * t4518
    t5296 = t16 * eta
    t5321 = -4200 * t16 * t4605 - 1575 * t456 * t1661 + 1575 * t460 * t1661 - 9450 * t5287 * t1661 - 1050 * t448 * t1713 - 3150 * t4504 * t1713 + 2625 * t452 * t1713 - 7875 * t5296 * t1713 - 3150 * t456 * t1749 + 3150 * t460 * t1749 - 18900 * t5287 * t1749 - 9450 * t224 * t1979 + 59850 * t239 * t2018 + 9450 * t243 * t2018 - 28350 * t228 * t2034 + 75600 * t21 * t4958 + 9450 * t228 * t4793 - 50400 * t231 * t4793 - 28350 * t365 * t4958
    t5362 = 3465 * t133 * t2438 + 27720 * t145 * t5194 - 2625 * t448 * t1779 - 7875 * t4504 * t1779 + 1050 * t452 * t1779 - 3150 * t5296 * t1779 - 13860 * t18 * t2438 - 56700 * t4931 * t2018 + 359100 * t4934 * t2018 - 18900 * t220 * t2078 + 18900 * t232 * t2078 - 28350 * t4500 * t2078 + 567000 * t4947 * t2078 - 28350 * t5004 * t2078 + 189000 * t224 * t2095 - 179550 * t224 * t4820 - 189000 * t228 * t4845 - 17325 * t23 * t2438 + 6930 * t28 * t5194 + 34650 * t33 * t5194
    t5365 = eta * t3
    t5404 = 3465 * t161 * t5365 - 20790 * t164 * t5365 - 51975 * t167 * t5365 - 27720 * t170 * t5365 + 62370 * t18 * t5230 - 1575 * t456 * t1808 + 1575 * t460 * t1808 - 9450 * t5287 * t1808 - 9450 * t239 * t2135 - 59850 * t243 * t2135 + 359100 * t4931 * t2135 - 56700 * t4934 * t2135 + 179550 * t228 * t2152 + 28350 * t224 * t4876 + 51975 * t23 * t5230 - 10395 * t50 * t2503 + 10395 * t54 * t2503 - 51975 * t28 * t2520 - 62370 * t33 * t2520
    t5441 = eta * t42
    t5446 = -27720 * t133 * t2593 - 3465 * t145 * t5265 + 525 * t16 * t4732 - 27720 * t161 * t5441 - 51975 * t164 * t5441 - 34650 * t18 * t2593 - 1400 * t547 * t1850 - 9450 * t224 * t2205 + 50400 * t339 * t2205 + 9450 * t228 * t4903 - 6930 * t23 * t2593 - 62370 * t4627 * t2503 - 124740 * t4630 * t2503 - 62370 * t4633 * t2503 + 17325 * t28 * t5265 + 13860 * t33 * t5265 - 28350 * t365 * t5083 - 175 * t550 * t4585 - 4200 * t7 * t4732 + 75600 * t8 * t5083
    t5484 = 15225 * t16 * t1713 + 6825 * t16 * t1779 - 35910 * t16 * t2897 + 7875 * t16 * t5 + 16800 * t583 * t1661 - 20790 * t167 * t5441 + 3465 * t170 * t5441 + 6825 * t7 * t1713 + 33600 * t583 * t1749 + 15225 * t7 * t1779 + 525 * t5090 * t1979 - 79380 * t368 * t2018 + 59220 * t371 * t2018 - 8400 * t5094 * t2034 - 3150 * t21 * t4 + 7875 * t5094 * t4793 + 8400 * t5090 * t4820 - 525 * t2235 + 1890 * t2648
    t5525 = -10395 * t16 * t2939 + 3465 * t161 * t3 - 17325 * t170 * t3 + 16800 * t583 * t1808 - 35910 * t21 * t2078 + 136080 * t365 * t2078 - 35910 * t8 * t2078 - 7350 * t5090 * t2095 - 31185 * t21 * t2942 + 59220 * t368 * t2135 - 79380 * t371 * t2135 - 8400 * t5094 * t2152 + 13230 * t224 * t5184 - 13230 * t228 * t5195 - 1890 * t4741 * t2438 + 15120 * t4749 * t2520 + 45360 * t4758 * t2520 - 28350 * t4749 * t5194 + 7350 * t5094 * t4845 + 8400 * t5090 * t4876
    t5565 = -525 * t16 * t114 - 35910 * t16 * t2833 - 31185 * t16 * t3007 - 17325 * t161 * t42 - 41580 * t173 * t2503 - 83160 * t176 * t2503 - 41580 * t179 * t2503 - 10395 * t21 * t3010 + 1890 * t21 * t80 - 7875 * t5090 * t2205 + 13230 * t224 * t5260 - 13230 * t228 * t5266 + 28350 * t4741 * t2593 - 15120 * t4741 * t5230 + 1890 * t4749 * t5265 - 45360 * t4767 * t5230 - 525 * t5094 * t4903 + 7875 * t2552 - 3150 * t5251
    t5568 = Q1 * eta
    t5575 = Q2 * eta
    t5607 = 525 * eta * t1713 + 525 * eta * t1779 - 1155 * t133 * t5568 + 1155 * t145 * t5575 - 23625 * t16 * t4958 + 3465 * t170 * t42 - 2310 * t18 * t5568 + 525 * t547 * t1979 - 8400 * t456 * t2018 + 8400 * t460 * t2018 - 50400 * t5287 * t2018 - 7350 * t448 * t2078 - 22050 * t4504 * t2078 + 7350 * t452 * t2078 - 1155 * t23 * t5568 + 1155 * t28 * t5575 + 2310 * t33 * t5575 + 7875 * t550 * t4793 + 1575 * t7 * t4958 + 175 * t4605
    t5646 = -22050 * t5296 * t2078 + 42525 * t21 * t5365 - 8400 * t456 * t2135 + 8400 * t460 * t2135 - 50400 * t5287 * t2135 - 7875 * t547 * t2205 + 13230 * t224 * t2438 - 45360 * t224 * t5230 + 45360 * t228 * t2520 - 13230 * t228 * t5194 - 28350 * t231 * t5194 + 15120 * t239 * t2503 - 15120 * t243 * t2503 - 1890 * t339 * t2438 + 90720 * t4931 * t2503 + 90720 * t4934 * t2503 + 39690 * t365 * t5365 - 2835 * t8 * t5365 + 175 * t4732
    t5682 = -2310 * t4514 * t15 * Q1 + 1575 * t16 * t5083 + 1155 * t550 * t1663 + 2310 * t231 * t1667 - 2835 * t21 * t5441 + 13230 * t224 * t2593 - 13230 * t228 * t5265 + 1890 * t231 * t5265 + 28350 * t339 * t2593 + 39690 * t365 * t5441 - 1155 * t4519 * t447 - 525 * t550 * t4903 - 23625 * t7 * t5083 + 42525 * t8 * t5441 - 1155 * t11 + 1155 * t38 - 1155 * t4648 - 3465 * t4651 - 3465 * t4654 - 1155 * t4657
    t5717 = 4410 * t16 * t2078 + 22680 * t16 * t2942 + 3045 * t16 * t4 + 3360 * t583 * t2018 + 4410 * t7 * t2078 + 26460 * t21 * t3 + 3360 * t583 * t2135 - 420 * t5090 * t2438 + 60480 * t368 * t2503 + 60480 * t371 * t2503 - 5040 * t5094 * t2520 + 5040 * t5090 * t5230 + 5460 * t5094 * t5194 - 350 * t114 - 1050 * t1713 - 1050 * t1779 + 1365 * t2897 - 3780 * t2939 - 350 * t5
    t5751 = 2100 * eta * t2078 + 22680 * t16 * t3010 + 1365 * t16 * t80 - 3780 * t21 * t42 + 1260 * t224 * t5568 - 1260 * t228 * t5575 - 1260 * t231 * t5575 - 420 * t547 * t2438 - 5460 * t5090 * t2593 + 1260 * t339 * t5568 + 420 * t5094 * t5265 + 5460 * t550 * t5194 - 1260 * t7 * t5365 - 1155 * t161 - 3465 * t164 - 3465 * t167 - 1155 * t170 + 3045 * t2833 + 26460 * t3007 + 1050 * t4958
    t5783 = -10080 * t16 * t3 - 16380 * t16 * t5365 - 1260 * t16 * t5441 - 1260 * t550 * t1667 - 5040 * t456 * t2503 + 5040 * t460 * t2503 - 30240 * t5287 * t2503 - 5460 * t547 * t2593 + 1260 * t4514 * t447 + 420 * t550 * t5265 - 16380 * t7 * t5441 - 420 * t2078 + 1260 * t220 - 1260 * t232 - 210 * t4 + 1890 * t4500 + 3780 * t4947 + 1890 * t5004 + 1050 * t5083
    t5806 = -20160 * t583 * t2503 - 280 * t547 * t5568 + 280 * t550 * t5575 + 80 * eta - 840 * t16 + 1890 * t21 + 480 * t3 - 10080 * t3010 + 3780 * t365 + 480 * t42 - 280 * t448 - 840 * t4504 + 280 * t452 - 840 * t5296 + 840 * t5365 + 840 * t5441 - 840 * t7 + 1890 * t8 - 210 * t80 + 80
    dMLdt = 3/128 * t857 / (1 + eta) * t862 * t864 * t2 * (t4674 + t4583 + t5062 + t4805 + t5484 + t4761 + t5189 + t5021 + t5446 + t5525 + t5646 + t4937 + t5107 + t5235 + t5806 + t5717 + t4626 + t5147 + t4538 + t5279 + t5751 + t5565 + t5362 + t4979 + t5682 + t5783 + t4715 + t5321 + t4848 + t4891 + t5404 + t5607)



    out = mu3*[dsmadt,dP1dt,dP2dt,dQ1dt,dQ2dt,dMLdt]
    return out;
end

############ drag
"""
    dragveragedgauss(k,k1r,skm,smk,m,KK,T1,T2,T3,cdelta,sdelta,cpsih,spsih,sma,P1,P2,eta,Q1,Q2,GG,Idrag,satellite,muPlanet,couplingstrategy)
    Returns the terms of the Gaussian Planetary equations due to the atmospheric drag
    averaged with respect to the orbital mean anomaly 

    The atmospheric drag acceleration model is coherent with the one used in the attitude propagator is the input couplingstrategy is set equal to 1.
    Otherwise, if couplingstrategy=2, the model of acceleration is simplified: the perturbing force is considered directed exactly as the relative incident
    flow. The control parameter couplingstrategy is an internal settings (it cannot be changed by an external user just through the propagator inputs). It must be 
    changed in semianalyticalpropagators_girf.jl, inside the function handleinputs. 

    This function is suitable only for triaxial satellites or axisymmetric satellite with one moment of inertia different from the other two.  

    
    INPUT  
    k,k1r,skm,smk,,m,T1,T2,T3,cdelta,sdelta,cpsih,spsih, = Sadov variables related quantities ***
    sma = semi-major axis
    P1  = ec*sin(om+OM), , ec = eccentricity, om = argument of the pericenter, OM = right ascension of the ascending node
    P2  = ec*cos(om+OM) 
    eta = sqrt(1-P1^2-P2^2)
    Q1  = tan(incl/2)sin(OM), incl = inclination
    Q2  = tan(incl/2)cos(OM)
    GG  = 1+Q1^2+Q2^2
    Idrag = averaged terms over the orbital mean anomaly of the Gauss planetary equations. 
    satellite = Dict with the characteristics of the satellite **
    muPlanet = gravitational parameter of the central body
    couplingstrategy = integer equal to either 1 or 2 to be used to select the drag acceleration model (see comment above)
    
    **
    satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD);
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
    constantstermsdragacc : array of two vectors collecting constants terms depending on the satellite geometry appearing in the Gauss planetary equations
   
    ***
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
    to GV and the xy plane:
    k1r = sqrt(1+k)
    smk = m/(k+m)
    skm = k/(m+k)
    T1 = ((m-1.0)*KK+EE)/m/KK;
    T2 = 2.0*(KK-EE)/KK/m^2.0 + (EE-2.0*KK)/m/KK;
    T3 = EE/KK
    KK,EE = complete elliptic integrals of first and second kind with characteristics m
    cdelta = cos(delta), sdelta = sin(delta), delta = inclination angle between 
    the XY plane and the plane perpendicular to GV
    cpsih = cos(psih), spsih = sin(psih), psih = angle between X axis and N
                
    OUTPUT
    out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt]
    with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

"""
function dragveragedgauss(k,k1r,skm,smk,m,KK,T1,T2,T3,cdelta,sdelta,cpsih,spsih,sma,P1,P2,eta,Q1,Q2,GG,Idrag,satellite,muPlanet,couplingstrategy)
    
    VVfVdrag = satellite["constantstermsdragacc"][1]
    VVfNdrag = satellite["constantstermsdragacc"][2]
    mass = satellite["mass"]

    if couplingstrategy == 1
        sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN = getdragaveragedtermsgausseq_strategy1(VVfVdrag,VVfNdrag,k,k1r,skm,smk,m,KK,T1,T2,T3,cdelta,sdelta,cpsih,spsih,Q1,Q2,Idrag,mass)
    else
        sTLfR,cTLfR,phiLP1fT,phiLM1fT,phiLM1sTLfT,sTLfT,phiLM1cTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN,phiLM1fR = getdragaveragedtermsgausseq_strategy2(VVfVdrag,k,k1r,skm,smk,m,KK,T1,T2,T3,cdelta,sdelta,cpsih,spsih,Q1,Q2,Idrag,mass,satellite,P1,P2,sma)
    end
    out = zeros(6);
    smaOmuS = sqrt(sma/muPlanet);
    out[1] =  2.0/eta*sqrt(sma^3.0/muPlanet)*(-P1*cTLfR + P2*sTLfR + phiLP1fT);
    out[2] =  eta*smaOmuS*(-cTLfR + P1*phiLM1fT + phiLM1sTLfT + sTLfT - P2*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
    out[3] =  eta*smaOmuS*( sTLfR + P2*phiLM1fT + phiLM1cTLfT + cTLfT + P1*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
    out[4] =  0.5*eta*smaOmuS*GG*phiLM1sTLfN;
    out[5] =  0.5*eta*smaOmuS*GG*phiLM1cTLfN;
    out[6] =  -eta*smaOmuS*( (P1*sTLfR+P2*cTLfR +P1*phiLM1cTLfT -P2*phiLM1sTLfT + P1*cTLfT - P2*sTLfT)/(1+eta) + 2*eta*phiLM1fR + Q1*phiLM1cTLfN - Q2*phiLM1sTLfN); 

    # out[abs.(out).<1e-15].=0.0;
    return out;
end

"""
    dragveragedgauss(k,k1r,skm,smk,m,KK,T1,T2,T3,cdelta,sdelta,cpsih,spsih,sma,P1,P2,eta,Q1,Q2,GG,Idrag,satellite,muPlanet,couplingstrategy)
    Returns the terms of the Gaussian Planetary equations due to the atmospheric drag
    averaged with respect to the orbital mean anomaly 

    The atmospheric drag acceleration model is coherent with the one used in the attitude propagator is the input couplingstrategy is set equal to 1.
    Otherwise, if couplingstrategy=2, the model of acceleration is simplified: the perturbing force is considered directed exactly as the relative incident
    flow. The control parameter couplingstrategy is an internal settings (it cannot be changed by an external user just through the propagator inputs). It must be 
    changed in semianalyticalpropagators_girf.jl, inside the function handleinputs. 

    This function is suitable only for axisymmetric satellites with equalt moments of inertia
    
    INPUT  
    k,k1r,skm,smk,,m,T1,T2,T3,cdelta,sdelta,cpsih,spsih, = Sadov variables related quantities ***
    sma = semi-major axis
    P1  = ec*sin(om+OM), , ec = eccentricity, om = argument of the pericenter, OM = right ascension of the ascending node
    P2  = ec*cos(om+OM) 
    eta = sqrt(1-P1^2-P2^2)
    Q1  = tan(incl/2)sin(OM), incl = inclination
    Q2  = tan(incl/2)cos(OM)
    GG  = 1+Q1^2+Q2^2
    Idrag = averaged terms over the orbital mean anomaly of the Gauss planetary equations. 
    satellite = Dict with the characteristics of the satellite **
    muPlanet = gravitational parameter of the central body
    couplingstrategy = integer equal to either 1 or 2 to be used to select the drag acceleration model (see comment above)
    
    **
    satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD);
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
    constantstermsdragacc : array of two vectors collecting constants terms depending on the satellite geometry appearing in the Gauss planetary equations
   
    ***
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
    to GV and the xy plane:
    k1r = sqrt(1+k)
    smk = m/(k+m)
    skm = k/(m+k)
    T1 = ((m-1.0)*KK+EE)/m/KK;
    T2 = 2.0*(KK-EE)/KK/m^2.0 + (EE-2.0*KK)/m/KK;
    T3 = EE/KK
    KK,EE = complete elliptic integrals of first and second kind with characteristics m
    cdelta = cos(delta), sdelta = sin(delta), delta = inclination angle between 
    the XY plane and the plane perpendicular to GV
    cpsih = cos(psih), spsih = sin(psih), psih = angle between X axis and N
                
    OUTPUT
    out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt]
    with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

"""
function dragveragedgauss_cube(csigma,ssigma,clA,slA,cdelta,sdelta,chA,shA,sma,P1,P2,eta,Q1,Q2,GG,Idrag,satellite,muPlanet,couplingstrategy)
    
    VVfVdrag = satellite["constantstermsdragacc"][1]
    VVfNdrag = satellite["constantstermsdragacc"][2]
    mass = satellite["mass"]

    if couplingstrategy == 1
        sTLfR,cTLfR,phiLM1fR,phiLP1fT,phiLM1fT,phiLM1sTLfT,phiLM1cTLfT,sTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN = getdragaveragedtermsgausseq_strategy1_cube(VVfVdrag,VVfNdrag,csigma,ssigma,clA,slA,cdelta,sdelta,chA,shA,Q1,Q2,Idrag,mass)
    else
        sTLfR,cTLfR,phiLP1fT,phiLM1fT,phiLM1sTLfT,sTLfT,phiLM1cTLfT,cTLfT,phiLM1sTLfN,phiLM1cTLfN,phiLM1fR = getdragaveragedtermsgausseq_strategy2_cube(VVfVdrag,csigma,ssigma,clA,slA,cdelta,sdelta,chA,shA,Q1,Q2,Idrag,mass,satellite,P1,P2,sma)
    end
    out = zeros(6);
    smaOmuS = sqrt(sma/muPlanet);
    out[1] =  2.0/eta*sqrt(sma^3.0/muPlanet)*(-P1*cTLfR + P2*sTLfR + phiLP1fT);
    out[2] =  eta*smaOmuS*(-cTLfR + P1*phiLM1fT + phiLM1sTLfT + sTLfT - P2*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
    out[3] =  eta*smaOmuS*( sTLfR + P2*phiLM1fT + phiLM1cTLfT + cTLfT + P1*(Q1*phiLM1cTLfN - Q2*phiLM1sTLfN));
    out[4] =  0.5*eta*smaOmuS*GG*phiLM1sTLfN;
    out[5] =  0.5*eta*smaOmuS*GG*phiLM1cTLfN;
    out[6] =  -eta*smaOmuS*( (P1*sTLfR+P2*cTLfR +P1*phiLM1cTLfT -P2*phiLM1sTLfT + P1*cTLfT - P2*sTLfT)/(1+eta) + 2*eta*phiLM1fR + Q1*phiLM1cTLfN - Q2*phiLM1sTLfN); 

    # out[abs.(out).<1e-15].=0.0;
    return out;
end




############################################### OLD ignore
# # not used
# function dragaveragedgauss_wrapper1(sma,P1,P2,Q1,Q2,ML,satellite,skm,Jg,Jh,psil,psig,psih,muPlanet,IDrag)

#     inverseballisticcoeff = dragballisticcoeff_nonaveraged_sadov(sma,P1,P2,Q1,Q2,ML,satellite,skm,Jg,Jh,psil,psig,psih);

#     # gauss equations
#     eta = sqrt(1-P1^2-P2^2)
#     out = dragaveragedgauss(inverseballisticcoeff,sma,sqrt(muPlanet/sma^3),eta,P1,P2,Q1,Q2,IDrag);

#     return out;
# end

# function dragaveragedgauss_wrapper2(sma,P1,P2,Q1,Q2,ML,satellite,LA,GA,HA,lA,gA,hA,muPlanet,IDrag)
#     inverseballisticcoeff = dragballisticcoeff_nonaveraged_andoyer(sma,P1,P2,Q1,Q2,ML,satellite,LA,GA,HA,lA,gA,hA)
#     # integrals
#     out = dragaveragedgauss(inverseballisticcoeff,sma,sqrt(muPlanet/sma^3),sqrt(1-P1^2-P2^2),P1,P2,Q1,Q2,IDrag);
#     return out;
# end

# # used
# # triaxial and axysimmetric such that (A=B!=C)
# """
#     dragaveragedgauss_wrapper1_AV(eV0mean,sma,P1,P2,Q1,Q2,satellite,k,skm,KK,T1,cdelta,sdelta,cpsih,spsih,muPlanet,IDrag)
#     Returns the terms of the Gaussian Planetary equations due to the atmospheric drag
#     averaged with respect to the orbital mean anomaly

#     INPUT  
#     eV0mean = vectors with averaged terms with respect to the orbital mean anomaly needed to couple dynamics and orbit
#     sma = semi-major axis
#     P1  = ec*sin(om+OM),ec = eccentricity om = argument of the pericenter, OM = right ascension of the ascending node
#     P2  = ec*cos(om+OM) 
#     Q1  = tan(incl/2)sin(OM), incl = inclination
#     Q2  = tan(incl/2)cos(OM)
#     satellite = Dict with the characteristics of the satellite **
#     k,skm,KK,cdelta,sdelta,cpsih,spsih = Sadov variables related quantities ***
#     muPlanet  = gravitational parameter of the central body
#     IDrag   = averaged terms of the gauss equations of motion with respect to the mean anomaly computed numerically
            

#     **
#     satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD);
#     IV = [A,B,C] : principal moment of inertia
#     IM = intrinsic magnetic moment of the satellite [A m^2]
#     n  = number of facets in which the satellite surface is divided
#     facets = [facet1...facetn]
#         where  
#             facetk = [facet_coeff,facet_area,facet_Vinfo,facet_nv]
#             facet_coeff = [cai*psrp,cdi*psrsp,csi*psrp], psrp = solar radiation pressure
#             facet_area  = area of the facet [m^2]
#             facet_Vinfo  = [rhoiv,vv1,vv2,vv3], rhoiv = centre of mass to facet centroid vector [m], vvi = verteces of the facet [m], i=1..3
#             facet_nv    = normal unit vector  
#             with 
#                 cai = (1-rf*sp)
#                 cdi = (2/3*(1-sp)*rf+2/3*(1-rf))
#                 csi = 2*rf*sp*psrp
#                 with rf reflectivity and sp specular coefficient (see Benson&Sheers2021)
#     CD: aerodynamic coefficient

    
#     ***
#     Consider an inertial reference Frame XYZ.
#     Consider a body reference frame xyz in principal axes of inertia. 
#     The reference frame is such that m belongs to in[0,1] and k>0 where
#     k = C/A*(B-A)/(C-B) 
#     m = k*(C-delta)/(delta-A)
#     with
#     A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
#     delta = G^2/2/Phi, with G the angular momentum of the body and Phi the kinetic energy
#     Let GV=diag([A,B,C])omV be angular momentum of the body and G = norm(GV)
#     (omV = [p,q,r]^T is the angular velocity) 
#     Let N be the node between the inertial reference XY plane and the plane
#     perpendicular to GV; let N'' be the node between the plane perpendicular 
#     to GV and the xy plane:
#     skm = k/(m+k)
#     KK  = complete integral of firt kind with characteristic m
#     T1 = ((m-1.0)*KK+EE)/m/KK;
#     cdelta = cos(delta), delta = inclination angle between 
#     the XY plane and the plane perpendicular to GV
#     sdelta = sin(delta)
#     cpsih  = cos(psih), psih = angle between X axis and N
#     spsih  = sin(psih)               

#     OUTPUT
#     out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt]
#     with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

# """
# function dragaveragedgauss_wrapper1_AV(eV0mean,sma,P1,P2,Q1,Q2,satellite,k,skm,KK,T1,cdelta,sdelta,cpsih,spsih,muPlanet,IDrag)
#     inverseballisticcoeff =  dragballisticcoeff_averaged(eV0mean,sma,P1,P2,Q1,Q2,satellite,k,skm,KK,T1,cdelta,sdelta,cpsih,spsih)
#     eta = sqrt(1-P1^2-P2^2)
#     out = dragaveragedgauss(inverseballisticcoeff,sma,sqrt(muPlanet/sma^3),eta,P1,P2,Q1,Q2,IDrag);
#     return out;
# end

# # for axysimmetric (A=B=C)
# """
#     dragaveragedgauss_wrapper2_AV(eV0mean,sma,P1,P2,Q1,Q2,satellite,ssigma,csigma,clA,slA,cdelta,sdelta,cpsih,spsih,muPlanet,IDrag)
#     Returns the terms of the Gaussian Planetary equations due to the atmospheric drag
#     averaged with respect to the orbital mean anomaly

#     INPUT  
#     eV0mean = vectors with averaged terms with respect to the orbital mean anomaly needed to couple dynamics and orbit
#     sma = semi-major axis
#     P1  = ec*sin(om+OM),ec = eccentricity om = argument of the pericenter, OM = right ascension of the ascending node
#     P2  = ec*cos(om+OM) 
#     Q1  = tan(incl/2)sin(OM), incl = inclination
#     Q2  = tan(incl/2)cos(OM)
#     satellite = Dict with the characteristics of the satellite **
#     ssigma,csigma,clA,slA,cdelta,sdelta,cpsih,spsih = Andoyer-Serret variables related quantities ***
#     muPlanet  = gravitational parameter of the central body
#     IDrag   = averaged terms of the gauss equations of motion with respect to the mean anomaly computed numerically
            

#     **
#     satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD);
#     IV = [A,B,C] : principal moment of inertia
#     IM = intrinsic magnetic moment of the satellite [A m^2]
#     n  = number of facets in which the satellite surface is divided
#     facets = [facet1...facetn]
#         where  
#             facetk = [facet_coeff,facet_area,facet_Vinfo,facet_nv]
#             facet_coeff = [cai*psrp,cdi*psrsp,csi*psrp], psrp = solar radiation pressure
#             facet_area  = area of the facet [m^2]
#             facet_Vinfo  = [rhoiv,vv1,vv2,vv3], rhoiv = centre of mass to facet centroid vector [m], vvi = verteces of the facet [m], i=1..3
#             facet_nv    = normal unit vector  
#             with 
#                 cai = (1-rf*sp)
#                 cdi = (2/3*(1-sp)*rf+2/3*(1-rf))
#                 csi = 2*rf*sp*psrp
#                 with rf reflectivity and sp specular coefficient (see Benson&Sheers2021)
#     CD: aerodynamic coefficient

    
#     ***
#     Consider an inertial reference Frame XYZ.
#     Consider a body reference frame xyz in principal axes of inertia. 
#     Let GV=diag([A,B,C])omV be angular momentum of the body 
#     (omV = [p,q,r]^T is the angular velocity) 
#     Let N be the node between the inertial reference XY plane and the plane
#     perpendicular to GV; let N'' be the node between the plane perpendicular
#     to GV and the xy plane.
#     ssigma = sin(sigma), sigma = inclination angle between
#             the plane perpendicular to GV and the xy plane
#     csigma = cos(sigma)
#     clA = cos(lA), lA = angle between N'' and x axis
#     slA = sin(lA)
#     cdelta = cos(delta), delta = inclination angle between 
#             the XY plane and the plane perpendicular to GV.
#     sdelta = sin(delta)
#     cpsih  = cos(hA), hA = angle between X axis and N        
#     spsih  = sin(hA)      

#     OUTPUT
#     out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt]
#     with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

# """
# function dragaveragedgauss_wrapper2_AV(eV0mean,sma,P1,P2,Q1,Q2,satellite,ssigma,csigma,clA,slA,cdelta,sdelta,cpsih,spsih,muPlanet,IDrag)
#     inverseballisticcoeff =  dragballisticcoeff_averaged_axisym(eV0mean,sma,P1,P2,Q1,Q2,satellite,ssigma,csigma,clA,slA,cdelta,sdelta,cpsih,spsih)
#     eta = sqrt(1-P1^2-P2^2)
#     out = dragaveragedgauss(inverseballisticcoeff,sma,sqrt(muPlanet/sma^3),eta,P1,P2,Q1,Q2,IDrag_orbit);
#     return out;
# end

# """
#     dragaveragedgauss(inverseballisticcoeff,sma,n,eta,P1,P2,Q1,Q2,IDrag)
#     Returns the terms of the Gaussian Planetary equations due to the atmospheric drag
#     averaged with respect to the orbital mean anomaly

#     INPUT  
#     sma = semi-major axis 
#     n = meanmotion
#     eta = sqrt(1-P1^2-P2^2)
#     P1  = ec*sin(om+OM), om = argument of the pericenter, OM = right ascension of the ascending node
#     P2  = ec*cos(om+OM) 
#     Q1  = tan(incl/2)sin(OM), incl = inclination
#     Q2  = tan(incl/2)cos(OM)
#     IDrag   = averaged terms of the gauss equations of motion with respect to the mean anomaly

#     OUTPUT
#     out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt], out = out(sma,P1,P2,Q1,Q2; mu, rSunV/norm(rSunV), PsrpCRAOm)
#     with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

# """
# function dragaveragedgauss(inverseballisticcoeff,sma,n,eta,P1,P2,Q1,Q2,IDrag);
#     out = zeros(6)
#     out[1] = 2/eta/n*(P2*IDrag[1]-P1*IDrag[2]+IDrag[3]);
#     out[2] = eta/n/sma*(-IDrag[2]+P1*IDrag[4]+IDrag[5]+IDrag[6]-P2*Q1*IDrag[10]+P2*Q2*IDrag[9])
#     out[3] = eta/n/sma*( IDrag[1]+P2*IDrag[4]+IDrag[7]+IDrag[8]+P1*Q1*IDrag[10]-P1*Q2*IDrag[9])
#     out[4] = eta/2/n/sma*(1+Q1^2+Q2^2)*IDrag[9]
#     out[5] = eta/2/n/sma*(1+Q1^2+Q2^2)*IDrag[10]
#     out[6] =  -eta/2/n*( (P1*IDrag[1]+P2*IDrag[2] +P1*IDrag[7] -P2*IDrag[5] + P1*IDrag[8] - P2*IDrag[6])/(1+eta) + 2*eta*IDrag[11] + Q1*IDrag[10] - Q2*IDrag[9]); 
#     out = out * inverseballisticcoeff/2.0;

#     return out;
# end

# function dragaveragedgaussO(inverseballisticcoeff,sma,n,eta,P1,P2,Q1,Q2,IDrag);
    
#     ec = sqrt(P1^2+P2^2);
#     incl = mod(atan(sqrt(Q1^2.0+Q2^2.0))*2.0,2.0*pi);
#     omPOM = 0.0;    
#     if ec!=0.0 
#         omPOM  = mod(atan(P1,P2),2*pi);
#     elseif ec==0.0 && incl!=0.0
#         omPOM  = mod(atan(Q1,Q2),2*pi);
#     end
   
#     out = zeros(6);
    
#     out[1]  = sma^2 * (ec^2 * IDrag[1] + IDrag[2]);
#     out[2]  = 0.5 * sma * eta^2 * (sin(omPOM) * (ec * IDrag[1] + ec * IDrag[4] + 2 * IDrag[6] + ec * IDrag[7]) + 2 * cos(omPOM) * IDrag[5]);
#     out[3]  = 0.5 * sma * eta^2 * (cos(omPOM) * (ec * IDrag[1] + ec * IDrag[4] + 2 * IDrag[6] + ec * IDrag[7]) - 2 * sin(omPOM) * IDrag[5]);
#     out[6]  = eta^2 *sma * ec *(IDrag[5] - eta*(1+eta) * IDrag[8])/(1+eta);

#     out = out * n * inverseballisticcoeff;

#     return out;

# end

##
# """
#     srpaveragedgauss_wrapper3(sma,eta,P1,P2,Q1,Q2,GG,satellite,LA,GA,HA,lA,gA,hA,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
#     Returns the terms of the Gaussian Planetary equations due to the solar radiation pressure
#     averaged with respect to the orbital mean anomaly

#     INPUT  
#     sma = semi-major axis
#     eta = sqrt(1-ec^2), ec = eccentricity
#     P1  = ec*sin(om+OM), om = argument of the pericenter, OM = right ascension of the ascending node
#     P2  = ec*cos(om+OM) 
#     Q1  = tan(incl/2)sin(OM), incl = inclination
#     Q2  = tan(incl/2)cos(OM)
#     GG  = 1+Q1^2+Q2^2
#     satellite = Dict with the characteristics of the satellite **
#     quat = quaternions
#     mu  = gravitational parameter of the central body
#     rSunV = position vector of the sun
#     includeEclipsesEffects = boolean : if true the eclipses effects are considered, if false the body is considered in sun-light (even if it is not)
#     passageinshadowoccurs = boolean: if true part of the orbit is in shadow; if false the orbit is in sun-light (if includeEclipsesEffects, set passageinshadowoccurs = false)
#     Elinshadow  = value of the eccentric longitude at the entrance of the shadow region (if passageinshadowoccurs, set  Elinshadow = NaN)
#     Eloutshadow = value of the eccentric longitude at the exit from the shadow region (if passageinshadowoccurs, set  Eloutshadow = NaN)

#     **
#     satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>n,"facets"=>facets,"CD"=>CD);
#     IV = [A,B,C] : principal moment of inertia
#     IM = intrinsic magnetic moment of the satellite [A m^2]
#     n  = number of facets in which the satellite surface is divided
#     facets = [facet1...facetn]
#         where  
#             facetk = [facet_coeff,facet_area,facet_Vinfo,facet_nv]
#             facet_coeff = [cai,cdi,csi]
#             facet_area  = area of the facet [m^2]
#             facet_Vinfo  = [rhoiv,vv1,vv2,vv3], rhoiv = centre of mass to facet centroid vector [m], vvi = verteces of the facet [m], i=1..3
#             facet_nv    = normal unit vector  
#             with 
#                 cai = (1-rf*sp)
#                 cdi = (2/3*(1-sp)*rf+2/3*(1-rf))
#                 csi = 2*rf*sp*psrp
#                 with rf reflectivity and sp specular coefficient (see Benson&Sheers2021)
#     CD: aerodynamic coefficient

#     OUTPUT
#     out = [dsma/dt, dP1/dt, dP2/dt, dQ1/dt, dQ2/dt dML/dt]
#     with ML = orbital mean longitude (=M+om+OM, M= orbital mean anomaly)

# """
# function srpaveragedgauss_wrapper3(sma,eta,P1,P2,Q1,Q2,GG,satellite,quat,mu,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow)
    
#     # planet-sun unit vector
#     uVect = rSunV/norm(rSunV);

#     Rb2i = quat2birotmat(quat);
#     Ri2b = Transpose(Rb2i);
#     # srp acc
#     uVect = Ri2b*uVect;
#     numberOfFacets = get(satellite,"numberOfFacets",0);
#     facets = get(satellite,"facets",0.0); 
#     srpforce = [0.0,0.0,0.0];
#     for kk = 1:numberOfFacets
#         niv   = facets[kk][4];
#         cai   = facets[kk][1][1];
#         cdi   = facets[kk][1][2];
#         csi   = facets[kk][1][3];
#         surfi = facets[kk][2];
#         nivuVect = niv[1]*uVect[1]+niv[2]*uVect[2]+niv[3]*uVect[3];
#         gifun =  max(0.0,nivuVect) #1/3/pi + nivuVect/2 + 4*nivuVect^2/3/pi;
#         srpforce = srpforce - (surfi*gifun*cai*uVect + surfi*gifun*cdi*niv + surfi*gifun*csi*nivuVect*niv);
#     end
    
#     srpacc = srpforce/satellite["mass"]/1000.0
#     PsrpCRAOm = -LinearAlgebra.dot(srpacc,uVect);

#     out = srpaveragedgauss(sma,eta,P1,P2,Q1,Q2,GG,mu,PsrpCRAOm,rSunV,includeEclipsesEffects,passageinshadowoccurs,ELinshadow,ELoutofshadow);
#     return out;

# end