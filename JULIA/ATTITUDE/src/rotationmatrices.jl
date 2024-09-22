###################################### attitude ref frame
"""
    euler2ibrotmat(euangles) 
    Returns the rotation matrix from the inertial reference frame
    to the rotating frame attached to the body, given the Euler angles


    Consider an inertial reference frame XYZ.
    Consider a body reference frame xyz.    
    Let N' be the node between the XY and xy planes.
    INPUT: 
    euangles = [phi,theta,psi]
        where
        phi    = angle between X axis and N'
        theta  = inclination between the XY plane and the xy plane
        psi    = angle N' and x axis
    OUTPUT 
        R [3x3] = rotation matrix from inertial to body reference frames. 

"""
function euler2ibrotmat(euangles)     
    phi    = euangles[1];
    theta  = euangles[2];
    psi    = euangles[3];

    Rphi   = [cos(phi) sin(phi) 0.0;
          -sin(phi) cos(phi) 0.0;
          0.0 0.0 1.0];

    Rtheta = [1.0 0.0 0.0; 
          0.0 cos(theta) sin(theta); 
          0.0 -sin(theta) cos(theta)];

    Rpsi   = [cos(psi) sin(psi) 0.0;
          -sin(psi) cos(psi) 0.0; 
          0.0 0.0 1.0];

    R = Rpsi*Rtheta*Rphi;

    return R;

end

"""
    sadov2ibrotmat(k,k1r,skmR,smkR,cn,sn,dn,cg,sg,cpsih,spsih,cdelta,sdelta) 
    Returns the rotation matrix from the inertial reference frame
    to the rotating frame attached to the body, given the Sadov variables

    Consider an inertial reference Frame XYZ.
    Consider a body-fixed reference frame xyz in principal axes of inertia. 
    The body reference frame is such that m belongs to [0,1] and k>0 where
    k = C/A*(B-A)/(C-B) 
    m = k*(C-delta)/(delta-A)
    with
    A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
    delta = G^2/2/Phi, with G the angular momentum of the body and Phi the energy of the associated
    free-torque problem    
    Let N' be the node between the XY and xy planes.

    INPUT: 
    k 
    k1r  = sqrt(1+k)
    skmR = sqrt(k/(m+k))
    smkR = sqrt(m/(m+k))
    cn   = JacobiCos(2/pi*psil*ellF(pi/2.0,m),m)
    sn   = JacobiSin(2/pi*psil*ellF(pi/2.0,m),m)
    dn   = sqrt(1-m*sn^2.0)
    cg   = cos(g), g = psig + sqrt{(k+m)(1+k)/k}(ellP(-k,pi/2.0,m)psil*2.0/pi-ellP(-k,JacobiAM(2/pi*psil*ellF(pi/2.0,m),m),m))
    sg   = sin(g)
    cpsih = cos(psih)
    spsih = sin(psih)
    cdelta = Jh/Jg
    sdelta = sqrt(1-Jh^2.0/Jg^2.0)
    with [Jl,Jg,Jh,psil,psig,psih] = sadov variables

    OUTPUT 
        R [3x3] = rotation matrix from inertial to body reference frames. 

"""
function sadov2ibrotmat(k,k1r,skmR,smkR,cn,sn,dn,cg,sg,cpsih,spsih,cdelta,sdelta)
    Rpsilpsig = zeros(3,3);
    dnk = sqrt(1+k*sn^2);   

    Rpsilpsig[1,1] = -(k1r*sn*cg+skmR*dn*cn*sg)/dnk;
    Rpsilpsig[1,2] = -(k1r*sn*sg-skmR*dn*cn*cg)/dnk;
    Rpsilpsig[1,3] = smkR*cn;
    Rpsilpsig[2,1] = -(cn*cg-skmR*k1r*sn*dn*sg)/dnk;
    Rpsilpsig[2,2] = -(cn*sg+skmR*k1r*sn*dn*cg)/dnk;
    Rpsilpsig[2,3] = -smkR*k1r*sn;
    Rpsilpsig[3,1] = smkR*dnk*sg;
    Rpsilpsig[3,2] = -smkR*dnk*cg;
    Rpsilpsig[3,3] = skmR*dn;

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

    R = Rpsilpsig*Rhdelta;

    return R,Rpsilpsig;   
end

"""
    andoyer2ibrotmat(andoyervar) 
    Returns the rotation matrix from the inertial reference frame
    to the rotating frame attached to the body, given the Andoyer-Serret variables
    
    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    Let A,B,C be the principal moments of inertia with 
    A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm. 
    
    Let GV=diag([A,B,C])omV be angular momentum of the body 
    (omV = [p,q,r]^T is the angular velocity) 
    
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N' be the node between the XY and xy planes; let 
    N'' be the node between the plane perpendicular to GV and the xy plane.
    
    INPUT: 
    andoyervar = [L,G,H,l,g,h]
        where
        L = |GV|cos(sigma) with sigma the inclination angle between
            the plane perpendicular to GV and the xy plane
        G = |GV|
        H = |GV|cos(delta) with delta the inclination angle between 
            the XY plane and the plane perpendicular to GV.
        l = angle between N'' and x axis
        g = angle between N and N''
        h = angle between X axis and N
    
    OUTPUT 
        R [3x3] = rotation matrix from inertial to body reference frames. 

"""
function andoyer2ibrotmat(andoyervar)
    L   = andoyervar[1];
    G   = andoyervar[2];
    H   = andoyervar[3];
    l   = andoyervar[4];
    g   = andoyervar[5];
    h   = andoyervar[6];
    
    sl = sin(l);
    cl = cos(l);
    sg = sin(g);
    cg = cos(g);
    sh = sin(h);
    ch = cos(h);

    csigma = L/G;
    if abs(csigma)>1.0 && abs(abs(csigma)-1.0)<1e-13 
        csigma = csigma/abs(csigma);
    elseif abs(csigma)>1.0
        println(L," ",G)
        Base.error("|cos(sigma)|>1 : wrong inputs");
    end
    ssigma = sqrt(1.0-csigma^2.0);
    cdelta = H/G;   
    if abs(cdelta)>1.0 && abs(abs(cdelta)-1.0)<1e-13 
        cdelta = cdelta/abs(cdelta);
    elseif abs(cdelta)>1.0
        println(H," ",G)
        Base.error("|cos(delta)|>1 : wrong inputs");
    end
    sdelta = sqrt(1.0-cdelta^2.0);

    al1 = ch*(-csigma*sg*sl + cg*cl) + sh*(sdelta*ssigma*sl - cdelta*(cg*csigma*sl + cl*sg));
    al2 = -ch*(cl*csigma*sg + cg*sl) + sh*(sdelta*ssigma*cl + cdelta*(-cg*cl*csigma + sg*sl));
    al3 = ch*sg*ssigma + sh*(cdelta*cg*ssigma + csigma*sdelta);
    bl1 = sh*(-csigma*sg*sl + cg*cl) - ch*(sdelta*ssigma*sl - cdelta*(cg*csigma*sl + cl*sg));
    bl2 = -sh*(cl*csigma*sg + cg*sl) - ch*(sdelta*ssigma*cl + cdelta*(-cg*cl*csigma + sg*sl));
    bl3 = sh*sg*ssigma - ch*(cdelta*cg*ssigma + csigma*sdelta);
    cl1 = sdelta*(cg*csigma*sl + cl*sg) + cdelta*ssigma*sl;
    cl2 = sdelta*(cg*cl*csigma - sg*sl) + cdelta*ssigma*cl;
    cl3 = -cg*sdelta*ssigma + cdelta*csigma;

    R = [al1 bl1 cl1; al2 bl2 cl2; al3 bl3 cl3];
    return R
end

"""
    taitbryan2ibrotmat(taitbryan) 
    Returns the rotation matrix from the inertial reference frame
    to the rotating frame attached to the body, given the Tait-bryan angles


    INPUT: 
        taitbryan     = [yaw,pitch,roll]
        tait-bryan angles
    OUTPUT 
        R [3x3] = rotation matrix from inertial to body reference frames. 

"""
function taitbryan2ibrotmat(taitbryan)
    #  Consider an inertial reference Frame XYZ.
    #  Consider a body reference frame xyz in principal axes of inertia. 
    #  Let A,B.C be the principal moments of inertia with 
    #  A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm. 
     
    #  INPUT: 
    # euangles = [phi,theta,psi]
    #             where
    #             phi    = angle between X axis and N'
    #             theta  = inclination between the XY plane and the xy plane
    #             psi    = angle N' and x axis
    #  OUTPUT: 
    #  R [3x3] = rotation matrix from inertial to body-fixed reference frames. 
     
    yaw    = taitbryan[1];
    pitch  = taitbryan[2];
    roll   = taitbryan[3];

    R1 = [cos(yaw) sin(yaw) 0.0; -sin(yaw) cos(yaw) 0.0; 0.0 0.0 1.0];
    R2 = [cos(pitch) 0.0 -sin(pitch); 0.0 1.0 0.0; sin(pitch) 0.0 cos(pitch)];
    R3 = [1.0 0.0 0.0; 0.0 cos(roll) sin(roll); 0.0 -sin(roll) cos(roll)];
    R = R3*R2*R1;
    return R;
end

"""
    quat2birotmat(q) 
    Returns the rotation matrix from the rotating frame attached
    to the body to the inertial reference frame, given the quaternions


    INPUT: 
        q     = [q1,q2,q3,q4]
        quaternions
    OUTPUT 
        R [3x3] = rotation matrix from inertial to body reference frames. 

"""
function quat2birotmat(q)
    q = q/norm(q);

    m = zeros(3,3);
    m[1] = 2.0*(q[1]^2.0 + q[2]^2) - 1.0;
    m[5] = 2.0*(q[1]^2.0 + q[3]^2) - 1.0;
    m[9] = 2.0*(q[1]^2.0 + q[4]^2) - 1.0;

    m[2] = 2.0*(q[2]*q[3] + q[1]*q[4]);
    m[4] = 2.0*(q[2]*q[3] - q[1]*q[4]);

    m[3] = 2.0*(q[2]*q[4] - q[1]*q[3]);
    m[7] = 2.0*(q[2]*q[4] + q[1]*q[3]);

    m[6] = 2.0*(q[3]*q[4] + q[1]*q[2]);
    m[8] = 2.0*(q[3]*q[4] - q[1]*q[2]);

    return m;
end

"""
    rotmatib2euler(R)
    Returns the Euler angles given the rotation matrix from the inertial reference frame
    to the rotating frame attached to the body

    Consider an inertial reference frame XYZ.
    Consider a body reference frame xyz.    
    Let N' be the node between the XY and xy planes.
    INPUT 
        R [3x3] = rotation matrix from inertial to body reference frames. 
    OUTPUT: 
    euangles = [phi,theta,psi]
        where
        phi    = angle between X axis and N'
        theta  = inclination between the XY plane and the xy plane
        psi    = angle N' and x axis

"""
function rotmatib2euler(R)    
    al1 = R[1,1]; 
    cl1 = R[1,3]; 
    al2 = R[2,1]; 
    cl2 = R[2,3];
    al3 = R[3,1]; 
    bl3 = R[3,2]; 
    cl3 = R[3,3];
    
    ctheta = cl3; 
    if abs(ctheta)>1.0 
        if abs(abs(ctheta)-1.0)<1e-15
            if ctheta>0.0
                ctheta=1.0
            else
                ctheta = -1.0
            end
        else
            println("ctheta = ",ctheta, "|ctheta|-1 = ",abs(abs(ctheta)-1.0) );
            Base.error("|ctheta|>1 : not possible")
        end
    end
    stheta = sqrt(1.0-ctheta^2.0); 
    
    if abs(stheta)<1e-20 
       cphi = 1.0; 
       sphi = 0.0; 
    
       cpsi = al1*cphi-al2*sphi; 
       spsi = -al1*sphi-al2*cphi;
    
    else
        cphi = -bl3/stheta; 
        sphi = al3/stheta; 
    
        cpsi = cl2/stheta; 
        spsi = cl1/stheta;
    end
    
    phi = mod(atan(sphi,cphi),2.0*pi); 
    theta = acos(ctheta);
    psi = mod(atan(spsi,cpsi),2.0*pi);
    
    euangles = [phi,theta,psi];

    return euangles;
end

"""
    rotmatib2taitbryan(R)
    Returns the Tait-Bryan angles given the rotation matrix from the inertial reference frame
    to the rotating frame attached to the body

    INPUT 
        R [3x3] = rotation matrix from inertial to body reference frames. 
    OUTPUT: 
        taitbryan     = [yaw,pitch,roll]
        tait-bryan angles

"""
function rotmatib2taitbryan(R)
    #  Consider an inertial reference Frame XYZ.
    #  Consider a body reference frame xyz in principal axes of inertia. 
    #  Let A,B.C be the principal moments of inertia with 
    #  A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm. 
     
    #  INPUT: 
    #  R [3x3] = rotation matrix from inertial to body-fixed reference frames. 
    #  OUTPUT: 
    # euangles = [phi,theta,psi]
    #             where
    #             phi    = angle between X axis and N'
    #             theta  = inclination between the XY plane and the xy plane
    #             psi    = angle N' and x axis
    
    al1 = R[1,1]; 
    bl1 = R[1,2];
    cl1 = R[1,3]; 
    al2 = R[2,1]; 
    cl2 = R[2,3];
    al3 = R[3,1]; 
    bl3 = R[3,2]; 
    cl3 = R[3,3];
    
    spitch = -cl1; 
    cpitch = sqrt(1.0-spitch^2.0); 
    
    if abs(cpitch)<1e-15 
       cyaw = 1.0; 
       syaw = 0.0; 
    
       croll = -(al2*syaw-bl2*cyaw)
       sroll = (al2*cyaw+bl2*syaw)/spitch;
    
    else
        cyaw = al1/cpitch; 
        syaw = bl1/cpitch; 
    
        sroll = cl2/cpitch; 
        croll = cl3/cpitch;
    end
    
    yaw = atan(syaw,cyaw); 
    pitch = atan(spitch,cpitch);
    roll = atan(sroll,croll); 
    
    taitbryan = [yaw,pitch,roll];

    return taitbryan;
end

"""
    rotmatib2taitbryan(R)
    Returns the quaternions given the rotation matrix from the rotating frame attached to the body
    to the inertial reference frame
    
    INPUT 
        R [3x3] = rotation matrix from inertial to body reference frames. 
    OUTPUT: 
        q     = [q1,q2,q3,q4]
        quaternions

"""
function rotmatbi2quat(m)
    trace = LinearAlgebra.tr(m);
    if trace>=0.0
        s = sqrt(trace+1.0);
        r = 0.5*s;
        s = 0.5/s;
        i = (m[3,2] - m[2,3])*s;
        j = (m[1,3] - m[3,1])*s;
        k = (m[2,1] - m[1,2])*s;
    else
        z = 1;
        if(m[2,2]>m[1,1])
            z = 2;
        end
        if (m[3,3]>m[z,z])
            z = 3;
        end
    
        if z==1
            s = sqrt(m[1,1] - (m[2,2]+m[3,3]) + 1.0);
            i = 0.5*s;
            s = 0.5/s;
            j = (m[1,2] + m[2,1])*s;
            k = (m[3,1] + m[1,3])*s;
            r = (m[3,2] - m[2,3])*s;
        elseif z== 2
            s = sqrt(m[2,2] - (m[1,1]+m[3,3]) + 1.0);
            j = 0.5*s;
            s = 0.5/s;
            k = (m[2,3] + m[3,2])*s;
            i = (m[1,2] + m[2,1])*s;
            r = (m[1,3] - m[3,1])*s;
        else
            s = sqrt(m[3,3] - (m[1,1]+m[2,2]) + 1.0);
            k = 0.5*s;
            s = 0.5/s;
            i = (m[3,1] + m[1,3])*s;
            j = (m[2,3] + m[3,2])*s;
            r = (m[2,1] - m[1,2])*s;
        end
    end
    
    q = [r i j k]';
    q = q/norm(q);
    return q;    
end




###################################### orbital ref frame 
"""
   kep2iprotmat(incl,OM,om)
    Return the rotation matrix from the inertial to the perifocal reference
    frames given the orbital inclination, longitude of the ascending node and
    argument of the pericenter

    INPUT :
        incl = inclination [rad]
        OM   = longitude of the ascending node [rad]
        om   = argument of the pericenter [rad]

    OUTPUT
        R [3x3] = rotation matrix from inertial to perifocal reference frames. 

"""
function kep2iprotmat(incl,OM,om)
    ae1 = cos(om)*cos(OM) - sin(om)*cos(incl)*sin(OM);
    be1 = cos(om)*sin(OM) + sin(om)*cos(incl)*cos(OM);
    ce1 = sin(om)*sin(incl); 
    
    ae2 = -sin(om)*cos(OM) - cos(om)*cos(incl)*sin(OM);
    be2 = -sin(om)*sin(OM) + cos(om)*cos(incl)*cos(OM);
    ce2 = cos(om)*sin(incl);
    
    ae3 = sin(incl)*sin(OM);
    be3 = -sin(incl)*cos(OM);
    ce3 = cos(incl);
    
    R = [ae1 be1 ce1;
         ae2 be2 ce2;
         ae3 be3 ce3];   

    return R;
end

"""
    equi2rotmatirth(Q1,Q2,sTL,cTL)
    Return the rotation matrix from the inertial to the rth reference
    frames given the equinoctial elements Q1, Q2, TL 

    INPUT :                 
        Q1 = tan(i/2)*sin(OM)  i = inclination, OM = longitude of ascending node, om = argument of pericentre  
        Q1 = tan(i/2)*cos(OM)  
        sTL = sin(TL), TL = true longitude (nu+OM+om, nu = true anomaly) [rad]
        cTL = cos(TL)
    OUTPUT
        R [3x3] = rotation matrix from inertial to rth reference frames. 

"""
function equi2rotmatirth(Q1,Q2,sTL,cTL)
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

    R = zeros(3,3); 
    R[1,1] =   cTL*q11 + sTL*q21; 
    R[1,2] =   cTL*q12 + sTL*q22; 
    R[1,3] =   cTL*q13 + sTL*q23; 
    R[2,1] = - sTL*q11 + cTL*q21; 
    R[2,2] = - sTL*q12 + cTL*q22; 
    R[2,3] = - sTL*q13 + cTL*q23; 
    R[3,1] =   q31; 
    R[3,2] =   q32; 
    R[3,3] =   q33;
    
    return R
end