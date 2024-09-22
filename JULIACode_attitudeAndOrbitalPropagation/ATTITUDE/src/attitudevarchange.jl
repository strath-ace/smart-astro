"""
    euler2eulercanonical(euanglesV,omV,IV)
    Transformation from euler angles + the angular velocity to Euler canonical variables

    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    Let A,B,C be the principal moments of inertia with 
    A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm. 
    
    Let omV = [p,q,r]^T be the angular velocity in the body reference
    frame. 
    
    Let N' be the node between the XY and xy planes.
    
    INPUT: 
    euanglesV = [phi,theta,psi]
               where
               phi    = angle between X axis and N'
               theta  = inclination between the XY plane and the xy plane
              psi    = angle N' and x axis
    omV = [p,q,r], angular velocity in the body reference
    frame. 
    IV  = [A,B,C]
    
    OUTPUT: 
    eulervar = [pphi, ptheta, ppsi,phi, theta, psi]
               where
    pphi   = momentum conjugated to phi 
                        (=A*p*sin(psi)*sin(theta)+B*q*cos(psi)*sin(theta)+C*r*cos(theta))
               ptheta = momentum conjugated to theta
                        (=A*p*cos(psi)-B*q*sin(psi)
               ppsi   = momentum conjugated to psi
                        (=C*r)

"""
function euler2eulercanonical(euanglesV,omV,IV)
    theta  = euanglesV[2];
    psi    = euanglesV[3];
    p = omV[1]; q = omV[2]; r = omV[3];
    A = IV[1]; B = IV[2]; C = IV[3];

    stheta = sin(theta);
    ctheta = cos(theta);
    spsi   = sin(psi);
    cpsi   = cos(psi);

    pphi   = (A*p*spsi+B*q*cpsi)*stheta+C*r*ctheta;
    ptheta = A*p*cpsi-B*q*spsi;
    ppsi   = C*r;

    eulervar = [pphi; ptheta; ppsi; euanglesV[1]; euanglesV[2]; euanglesV[3]];
    
    return eulervar;
end

"""
    eulercanonical2euler(eulervar,IV)
    Transformation from Euler canonical variables to Euler angles + the angular velocity

    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    Let A,B,C be the principal moments of inertia with 
    A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm. 
    
    Let omV = [p,q,r]^T be the angular velocity in the body reference
    frame. 
    
    Let N' be the node between the XY and xy planes.
    
    INPUT: 
    eulervar = [pphi, ptheta, ppsi, phi, theta, psi]
               where
               pphi   = momentum conjugated to phi 
                        (=A*p*sin(psi)*sin(theta)+B*q*cos(psi)*sin(theta)+C*r*cos(theta))
               ptheta = momentum conjugated to theta
                        (=A*p*cos(psi)-B*q*sin(psi)
               ppsi   = momentum conjugated to psi
    IV  = [A,B,C]
                        (=C*r)
    OUTPUT: 
    euangles = [phi,theta,psi]
               where
               phi    = angle between X axis and N'
               theta  = inclination between the XY plane and the xy plane
               psi    = angle N' and x axis
    omV = [p,q,r], angular velocity in the body reference frame
    frame. 
"""
function eulercanonical2euler(eulervar,IV)   
    pphi    = eulervar[1]; 
    ptheta  = eulervar[2];
    ppsi    = eulervar[3];
    theta   = eulervar[5];
    psi     = eulervar[6];

    A = IV[1]; B = IV[2]; C = IV[3];

    stheta = sin(theta);
    ctheta = cos(theta);
    spsi   = sin(psi);
    cpsi   = cos(psi);

    p = ((pphi-ppsi*ctheta)*spsi/stheta+ptheta*cpsi)/A;
    q = ((pphi-ppsi*ctheta)*cpsi/stheta-ptheta*spsi)/B;
    r = ppsi/C;

    euanglesV = [eulervar[4]; eulervar[5]; eulervar[6]];
    omV = [p;q;r];

    return euanglesV, omV;
end

"""
        euler2quat(euangles)
        Transformation from Euler angles to quaternions

        Consider an inertial reference Frame XYZ.
        Consider a body reference frame xyz in principal axes of inertia. 
        Let A,B,C be the principal moments of inertia with 
        A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.        
        Let N' be the node between the XY and xy planes.
        
        INPUT: 
        euanglesV = [phi,theta,psi]
                   where
                   phi    = angle between X axis and N'
                   theta  = inclination between the XY plane and the xy plane
                   psi    = angle N' and x axis

        OUTPUT:
        quat      = [q1,q2,q3,q4]
                    quaternions
"""
function euler2quat(euangles)
    Rb2i = transpose(euler2ibrotmat(euangles));
    quat = rotmatbi2quat(Rb2i);
    return quat;
end

"""
        euler2quat(euangles)
        Transformation from Euler angles to quaternions
        
        Consider an inertial reference Frame XYZ.
        Consider a body reference frame xyz in principal axes of inertia. 
        Let A,B,C be the principal moments of inertia with 
        A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.        
        Let N' be the node between the XY and xy planes.
        
        INPUT: 
        quat      = [q1,q2,q3,q4]
        quaternions
        
        OUTPUT:
        euanglesV = [phi,theta,psi]
        where
        phi    = angle between X axis and N'
        theta  = inclination between the XY plane and the xy plane
        psi    = angle N' and x axis
"""
function quat2euler(quat)
    Ri2b = transpose(quat2birotmat(quat));
    euangles = rotmatib2euler(Ri2b);
    return euangles;
end

"""
        euler2taitbryan(euangles)
        Transformation from Euler angles to Tait-Bryan angles 
        
        Consider an inertial reference Frame XYZ.
        Consider a body reference frame xyz in principal axes of inertia. 
        Let A,B,C be the principal moments of inertia with 
        A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.        
        Let N' be the node between the XY and xy planes.
        
        INPUT: 
        euanglesV = [phi,theta,psi]
                   where
                   phi    = angle between X axis and N'
                   theta  = inclination between the XY plane and the xy plane
                   psi    = angle N' and x axis

        OUTPUT:
        tb      = [yaw,pitch,roll]
                    tait-bryan angles
"""
function euler2taitybryan(euangles)
    Ri2b = euler2ibrotmat(euangles);
    tb = rotmatib2taitbryan(Ri2b)
    return tb;
end

"""
        taitbryan2euler(euangles)
        Transformation from Tait-Bryan angles to Euler angles
        
        Consider an inertial reference Frame XYZ.
        Consider a body reference frame xyz in principal axes of inertia. 
        Let A,B,C be the principal moments of inertia with 
        A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.        
        Let N' be the node between the XY and xy planes.
        
        INPUT: 
        tb      = [yaw,pitch,roll]
                  tait-bryan angles
        
        OUTPUT:
        euanglesV = [phi,theta,psi]
        where
        phi    = angle between X axis and N'
        theta  = inclination between the XY plane and the xy plane
        psi    = angle N' and x axis
"""
function taitbryan2euler(tb)
    Ri2b = taitbryan2ibrotmat(tb);
    euangles = rotmatib2euler(Ri2b);
    return euangles;
end

"""
    eulercanonical2andoyer(eulervar)
    Transformation from Euler canonical variables to Andoyer-Serret variables

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
    eulervar = [pphi, ptheta, ppsi,phi, theta, psi]
               where
               pphi   = momentum conjugated to phi 
                        (=A*p*sin(psi)*sin(theta)+B*q*cos(psi)*sin(theta)+C*r*cos(theta))
               ptheta = momentum conjugated to theta
                        (=A*p*cos(psi)-B*q*sin(psi)
               ppsi   = momentum conjugated to psi
                        (=C*r)
               phi    = angle between X axis and N'
               theta  = inclination between the XY plane and the xy plane
               psi    = angle N' and x axis
    
    OUTPUT
    andoyervar = [L,G,H,l,g,h]
                where
                L = |GV|cos(sigma), sigma = inclination angle between
                    the plane perpendicular to GV and the xy plane
                G = |GV|
                H = |GV|cos(delta),  delta = inclination angle between 
                    the XY plane and the plane perpendicular to GV.
                l = angle between N'' and x axis
                g = angle between N and N''
                h = angle between X axis and N
"""
function eulercanonical2andoyer(eulervar)
    phi    = eulervar[4];
    theta  = eulervar[5];
    psi    = eulervar[6];
    pphi   = eulervar[1];
    ptheta = eulervar[2];
    ppsi   = eulervar[3];
    
    stheta = sin(theta);
    ctheta = cos(theta);
    sphi   = sin(phi);
    cphi   = cos(phi);
    spsi   = sin(psi);
    cpsi   = cos(psi);
    
    H = pphi;
    L = ppsi;

    Ap = (pphi-ppsi*ctheta)*spsi/stheta+ptheta*cpsi;
    Bq = (pphi-ppsi*ctheta)*cpsi/stheta-ptheta*spsi;
    Cr = ppsi;
    G = sqrt(Ap^2.0+Bq^2.0+Cr^2.0);
    
    csigma = L/G;
    ssigma = sqrt(1.0-csigma^2.0);
    cdelta = H/G;
    sdelta = sqrt(1.0-cdelta^2.0);
    
    sl = Ap/(G*ssigma);
    cl = Bq/(G*ssigma); 
    l = mod(atan(sl,cl),2.0*pi);
    
    cg = (cdelta*csigma-ctheta)/sdelta/ssigma;
    sg = stheta*(spsi*cl-cpsi*sl)/sdelta;
    g  = mod(atan(sg,cg),2.0*pi);
    
    ch = sphi*sg*ssigma/stheta  + cphi*(csigma*sdelta+ssigma*cg*cdelta)/stheta;
    sh = -cphi*sg*ssigma/stheta + sphi*(csigma*sdelta+ssigma*cg*cdelta)/stheta;
    h  = mod(atan(sh,ch),2.0*pi);
       
    andoyervar = [L;G;H;l;g;h];
    
    return andoyervar;
end

"""
    andoyer2eulercanonical(andoyervar)
    Transformation from Andoyer-Serret variables to Euler canonical variables
    
    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    Let A,B,C be the principal moments of inertia with 
    A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm. 
    
    Let GV=diag([A,B,C])omV be angular momentum of the body 
    (omV = [p,q,r]^T is the angular velocity) 
    
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N' be the node between the XY and xy planes; let 
    N'' be the node between the plane perpendicular to GV and the xy plane.
    
    INPUT
    anoyervar = [L,G,H,l,g,h]
                where
                L = |GV|cos(sigma), sigma = inclination angle between
                    the plane perpendicular to GV and the xy plane
                G = |GV|
                H = |GV|cos(delta), delta = inclination angle between 
                    the XY plane and the plane perpendicular to GV.
                l = angle between N'' and x axis
                g = angle between N and N''
                h = angle between X axis and N
    
    OUTPUT: 
    eulervar = [pphi, ptheta, ppsi,phi, theta, psi]
               where
               pphi   = momentum conjugated to phi 
                        (=A*p*sin(psi)*sin(theta)+B*q*cos(psi)*sin(theta)+C*r*cos(theta))
               ptheta = momentum conjugated to theta
                        (=A*p*cos(psi)-B*q*sin(psi)
               ppsi   = momentum conjugated to psi
                        (=C*r)
               phi    = angle between X axis and N'
               theta  = inclination between the XY plane and the xy plane
               psi    = angle N' and x axis
"""
function andoyer2eulercanonical(andoyervar)       
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
    ssigma = sqrt(1.0-csigma^2.0);
    cdelta = H/G;
    sdelta = sqrt(1.0-cdelta^2.0);
    
    pphi=H;
    ppsi=L;
    
    ctheta = cdelta*csigma-sdelta*ssigma*cg;
    stheta = sqrt(1.0-ctheta^2.0); 
    theta = mod(atan(stheta,ctheta),2.0*pi); 
    
    cphi = -sh*sg*ssigma/stheta + ch*(csigma*sdelta+ssigma*cg*cdelta)/stheta;
    sphi = ch*sg*ssigma/stheta  + sh*(csigma*sdelta+ssigma*cg*cdelta)/stheta;
    phi  = mod(atan(sphi,cphi),2.0*pi);
    
    cpsi = -sl*sg*sdelta/stheta + cl*(cdelta*ssigma+csigma*sdelta*cg)/stheta;
    spsi = cl*sg*sdelta/stheta  + sl*(cdelta*ssigma+csigma*sdelta*cg)/stheta;
    psi  = mod(atan(spsi,cpsi),2.0*pi);
    
    ptheta = G*ssigma*(sl*cpsi-cl*spsi);
    
    eulervar = [pphi;ptheta;ppsi;phi;theta;psi];
    
    return eulervar;
end

"""
    euler2andoyer(euanglesV,omV,IV)
    Transformation from euler angles + the angular velocity to Andoyer-Serret variables

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
    euanglesV = [phi,theta,psi]
        where
        phi    = angle between X axis and N'
        theta  = inclination between the XY plane and the xy plane
        psi    = angle N' and x axis
    omV = [p,q,r], angular velocity in the body reference frame. 
    IV  = [A,B,C]
    
    OUTPUT
    andoyervar = [L,G,H,l,g,h]
                where
                L = |GV|cos(sigma), sigma = inclination angle between
                    the plane perpendicular to GV and the xy plane
                G = |GV|
                H = |GV|cos(delta), delta = inclination angle between 
                    the XY plane and the plane perpendicular to GV.
                l = angle between N'' and x axis
                g = angle between N and N''
                h = angle between X axis and N
"""
function euler2andoyer(euanglesV,omV,IV)
    phi    = euanglesV[1];
    theta  = euanglesV[2];
    psi    = euanglesV[3];

    stheta = sin(theta);
    ctheta = cos(theta);
    spsi   = sin(psi);
    cpsi   = cos(psi);
   
    p = omV[1]; q = omV[2]; r = omV[3];
    A = IV[1]; B = IV[2]; C = IV[3];
    Ap = A*p; Bq = B*q; Cr = C*r;

    L = Cr;
    G = sqrt(Ap^2.0+Bq^2.0+Cr^2.0);
    H = (Ap*spsi+Bq*cpsi)*stheta+Cr*ctheta;

    csigma = L/G;
    if abs(csigma)>1.0 && abs(abs(csigma)-1.0)<1e-15 
        csigma = csigma/abs(csigma);
    elseif abs(csigma)>1.0
        Base.error("|cos(sigma)|>1 : wrong inputs");
    end
    ssigma = sqrt(1.0-csigma^2.0);
    cdelta = H/G;   
    if abs(cdelta)>1.0 && abs(abs(cdelta)-1.0)<1e-15 
        cdelta = cdelta/abs(cdelta);
    elseif abs(cdelta)>1.0
        Base.error("|cos(delta)|>1 : wrong inputs");
    end
    sdelta = sqrt(1.0-cdelta^2.0); 
    
    if ssigma==0.0
        sl = 0.0;
        cl = 1.0;
        l = 0.0;
    else
        sl = Ap/(G*ssigma);
        cl = Bq/(G*ssigma); 
        l = mod(atan(sl,cl),2*pi);
    end

    RlAsigmaT = zeros(3,3);
    RlAsigmaT[1,1] = cl;
    RlAsigmaT[1,2] = -sl;
    RlAsigmaT[2,1] = sl*csigma;
    RlAsigmaT[2,2] = cl*csigma;
    RlAsigmaT[2,3] = -ssigma;
    RlAsigmaT[3,1] = sl*ssigma;
    RlAsigmaT[3,2] = cl*ssigma;
    RlAsigmaT[3,3] = csigma;
    Ri2b = euler2ibrotmat([phi,theta,psi]);
    RT   = RlAsigmaT*Ri2b;

    # println("RlAsigmaT =",RlAsigmaT)
    # println("RT =
    # ")

    if sdelta==0.0
        ch = 1.0;
        sh = 0.0; 
        cg = RT[1,1]*ch + RT[1,2]*sh
        sg = -(RT[2,1]*ch + RT[2,2]*sh)
        g = mod(atan(sg,cg),2*pi);
        h  = mod(atan(sh,ch),2.0*pi);
    else
        sg =  RT[1,3]/sdelta;
        cg =  RT[2,3]/sdelta;
        sh =  RT[3,1]/sdelta;
        ch = -RT[3,2]/sdelta;
        g = mod(atan(sg,cg),2*pi);
        h  = mod(atan(sh,ch),2.0*pi);
    end

      
    andoyervar = [L;G;H;l;g;h];
    
    return andoyervar;
end

"""
    andoyer2euler(andoyervar,IV)
    Transformation from Andoyer-Serret variables to euler angles + the angular velocity 

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
        L = |GV|cos(sigma), sigma = inclination angle between
            the plane perpendicular to GV and the xy plane
        G = |GV|
        H = |GV|cos(delta), delta = inclination angle between 
            the XY plane and the plane perpendicular to GV.
        l = angle between N'' and x axis
        g = angle between N and N''
        h = angle between X axis and N
    IV  = [A,B,C]
    
    OUTPUT:
    euanglesV = [phi,theta,psi]
        where
        phi    = angle between X axis and N'
        theta  = inclination between the XY plane and the xy plane
        psi    = angle N' and x axis
    omV = [p,q,r], angular velocity in the body reference frame. 
"""
function andoyer2euler(andoyervar,IV)
  
    R = andoyer2ibrotmat(andoyervar);
    euangles = ATTITUDE.rotmatib2euler(R);
  
    L   = andoyervar[1];
    G   = andoyervar[2];
    l   = andoyervar[4];
    sl = sin(l);
    cl = cos(l);
 
    csigma = L/G;
    ssigma = sqrt(1.0-csigma^2.0);

    A = IV[1]; B = IV[2]; C = IV[3];

    p = G*sl*ssigma/A;
    q = G*cl*ssigma/B;
    r = L/C;

    omV = [p,q,r];
    return euangles,omV
end

"""
        euler2sadov(euanglesV,omV,IV)
        Transformation from euler angles + the angular velocity to Sadov variables

        Consider an inertial reference frame XYZ.
        Consider a body reference frame xyz in principal axes of inertia. 
        The body-reference frame must be such that m belongs to [0,1] and k>0, where    
        k = C/A*(B-A)/(C-B)
        m = k*(C-delta)/(delta-A)
        with
            A,B,C = the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
            delta = G^2/2/Phi, with G the magnitude of the angular momentum of the body and Phi 
                the kinetic energy

        Let omV = [p,q,r]^T be the angular velocity of the satellite in the body frame
        
        Let N' be the node between the XY and xy planes.
        Let GV be angular momentum of the body (G = norm(GV)). 
        Let N be the node between the inertial reference XY plane and the plane
        perpendicular to GV; let N'' be the node between the plane perpendicular 
        to GV and the xy plane.

        INPUT: 
        euanglesV = [phi,theta,psi]
            where
            phi    = angle between X axis and N'
            theta  = inclination between the XY plane and the xy plane
            psi    = angle N' and x axis
        omV = [p,q,r], angular velocity in the body reference frame. 
        IV  = [A,B,C]
        
        OUTPUT
        sadovvar = [Jl,Jg,Jh,psil,psig,psih]
        with
                Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m))
                Jg   = G
                Jh   = Gcos(delta), delta = inclination angle between the XY plane and the plane perpendicular to GV.
                psil = pi/2 ellK(m)ellF(lambda,m), lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
                psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)), g = angle between N and N''
                psih = angle between X axis and N
        skm = k/(m+k)
"""
function euler2sadov(euanglesV,omV,IV)

    # input
    phi    = euanglesV[1];
    theta  = euanglesV[2];
    psi    = euanglesV[3];
    p = omV[1]; 
    q = omV[2];
    r = omV[3];

    stheta = sin(theta);
    ctheta = cos(theta);
    spsi   = sin(psi);
    cpsi   = cos(psi);
       
    # k
    A = IV[1]; B = IV[2]; C = IV[3];
    k = C/A*(B-A)/(C-B); 
    if (C*(B-A)-A*(C-B)==0.0)
        k=0.0
    end
    
    # Jg
    Ap = A*p; Bq = B*q; Cr = C*r;   
    Jg = sqrt((A*p)^2.0+(B*q)^2.0+(C*r)^2.0);

    # skm
    smk = (Ap/Jg)^2.0+(Bq/Jg/sqrt(1+k))^2.0;
    skm = 1.0-smk;
    mtilde   = smk/skm;
    if k==0.0
        m = 0.0;
    else
        m   = mtilde*k;
    end

    if m<0.0 || m>1.0
        println("m = ",m," mtilde ",mtilde," k ",k)
        Base.error("choice of reference frame wrong ",m)
    end
   
    skmR = sqrt(skm); 
    smkR = sqrt(smk);

    # Jl
    KK = ellF(pi/2.0,m);    #Elliptic.K(m);
    PP = ellP(-k,pi/2.0,m); #Elliptic.Pi(-k,pi/2,m);
    Jl = m2JlJgratioV2(skm,k,m,PP,KK)*Jg;
   
    #Jh
    Jh = (Ap*spsi+Bq*cpsi)*stheta+Cr*ctheta;
    cdelta = Jh/Jg;   
    if abs(cdelta)>1.0 && abs(abs(cdelta)-1.0)<1e-15 
        cdelta = cdelta/abs(cdelta);
    elseif abs(cdelta)>1.0
        Base.error("|cos(delta)|>1 : wrong inputs");
    end
    sdelta = sqrt(1.0-cdelta^2.0);

    # psil
    if Ap==0.0 && Bq==0.0
        cnu = 0.0
        snu = -1.0
    else
        cnu = Ap/(Jg*sqrt(smk))
        snu = -Bq/(Jg*sqrt(smk)*sqrt(1+k))
    end
    lambda = mod(atan(snu,cnu),2.0*pi);
    dnu = sqrt(1.0-m*snu^2.0)
    u = ellF(lambda,m);
    
    if m==1.0
        nlambda = floor(lambda*2.0/pi);
        psil = pi/2.0*nlambda;
    else
        psil = pi*u/2.0/KK;
    end 

    # psig & psih
    dnk = sqrt(1.0+k*snu^2.0);
    RlAsigmaT = zeros(3,3);
    RlAsigmaT[1,1] = -sqrt(1+k)*snu/dnk;
    RlAsigmaT[1,2] = -cnu/dnk;
    RlAsigmaT[2,1] = skmR*cnu*dnu/dnk;
    RlAsigmaT[2,2] = -sqrt(1+k)*skmR*snu*dnu/dnk;
    RlAsigmaT[2,3] = -smkR*dnk;
    RlAsigmaT[3,1] =  smkR*cnu;
    RlAsigmaT[3,2] =  -smkR*sqrt(1+k)*snu;
    RlAsigmaT[3,3] = skmR*dnu;
    Ri2b = euler2ibrotmat([phi,theta,psi]);
    RT   = RlAsigmaT*Ri2b;

    if sdelta==0.0
        ch = 1.0;
        sh = 0.0; 
        cg = RT[1,1]*ch + RT[1,2]*sh
        sg = -(RT[2,1]*ch + RT[2,2]*sh)
        g = mod(atan(sg,cg),2*pi);
        psih  = mod(atan(sh,ch),2.0*pi);
    else
        sg =  RT[1,3]/sdelta;
        cg =  RT[2,3]/sdelta;
        sh =  RT[3,1]/sdelta;
        ch = -RT[3,2]/sdelta;
        g = mod(atan(sg,cg),2*pi);
        psih  = mod(atan(sh,ch),2.0*pi);
    end

    
    if m == 1.0
        psig = mod(g + sqrt(1+k)*atan(sqrt(k)*sin(lambda))/sqrt(mtilde+1)/sqrt(k),2*pi); 
    else
        eP = ellP(-k,lambda,m);
        eF = ellF(lambda,m); #Elliptic.F(lambda,m);
        if k==0.0
            psig = g;
        else
            psig = mod(g + sqrt(1+k)/skmR*(eP-PP*eF/KK),2.0*pi) #mod(g + mod(sqrt(1+k)/skmR*(eP-PP*eF/KK),2.0*pi),2*pi); 
        end
    end

    sadov = [Jl,Jg,Jh,psil,psig,psih];
    
    return sadov,skm
end

"""
        sadov2eulerAndJl(sadovvar,skm,IV)
        Transformation from Sadov variables to euler angles + the angular velocity 

        Consider an inertial reference frame XYZ.
        Consider a body reference frame xyz in principal axes of inertia. 
        The body-reference frame must be such that m belongs to [0,1] and k>0, where    
        k = C/A*(B-A)/(C-B)
        m = k*(C-delta)/(delta-A)
        with
            A,B,C = the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
            delta = G^2/2/Phi, with G the magnitude of the angular momentum of the body and Phi 
                the kinetic energy

        Let omV = [p,q,r]^T be the angular velocity of the satellite in the body frame
        
        Let N' be the node between the XY and xy planes.
        Let GV be angular momentum of the body (G = norm(GV)). 
        Let N be the node between the inertial reference XY plane and the plane
        perpendicular to GV; let N'' be the node between the plane perpendicular 
        to GV and the xy plane.

        INPUT: 
        sadovvar = [Jg,Jh,psil,psig,psih]
        with
            Jg   = G
            Jh   = Gcos(delta), delta = inclination angle between the XY plane and the plane perpendicular to GV.
            psil = pi/2 ellK(m)ellF(lambda,m), lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
            psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)), g = angle between N and N''
            psih = angle between X axis and N
        skm = k/(m+k)
        IV  = [A,B,C]

        OUTPUT:
        euanglesV = [phi,theta,psi]
            where
            phi    = angle between X axis and N'
            theta  = inclination between the XY plane and the xy plane
            psi    = angle N' and x axis
        omV = [p,q,r], angular velocity in the body reference frame. 
        Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m))
        
"""
function sadov2eulerAndJl(sadovvar,skm,IV)
    
    Jg = sadovvar[1]; 
    Jh = sadovvar[2]; 
    psil = sadovvar[3]; 
    psig = sadovvar[4]; 
    psih = sadovvar[5]; 
    
    A = IV[1]; B = IV[2]; C = IV[3];
    k = C/A*(B-A)/(C-B);

    if (C*(B-A)-A*(C-B)==0.0)
        k=0.0
    end

    k1r = sqrt(1+k);
    smk = 1-skm;
    skmR = sqrt(skm);
    smkR = sqrt(smk);

    # m
    if skm == 0.0
        m = 0.0;
    else
        m = smk*k/skm;
        if m>1.0 && abs(m-1)<1e-13
            m = 1.0
        end
    end

    # Jl
    PP = ellP(-k,pi/2,m);  #Elliptic.Pi(-k,pi/2.0,m);
    KK = ellF(pi/2,m);     
    Jl = Jg*m2JlJgratioV2(skm,k,m,PP,KK)
  
    # g
    if m == 1.0
        lambda = psil;
        snu = sin(lambda);
        cnu = cos(lambda);
        dnu = sqrt(1-m*snu^2);
        g = mod(psig - sqrt(1+k)*atan(sqrt(k)*sin(lambda))*skmR/sqrt(k),2*pi); 
    else
        u = 2.0*psil*KK/pi;
        snu = Elliptic.Jacobi.sn(u,m);
        cnu = Elliptic.Jacobi.cn(u,m);
        dnu = Elliptic.Jacobi.dn(u,m);
        lambda = mod(atan(snu,cnu),2.0*pi); 
        eP = ellP(-k,lambda,m); 
        g = mod(psig - sqrt(1+k)/skmR*(eP-PP*u/KK),2.0*pi); #mod(psig - mod(sqrt(1+k)/skmR*(eP-u*PP/KK),2.0*pi),2.0*pi)
    end 

    sg = sin(g); 
    cg = cos(g);

    # rotational matrix
    cpsih = cos(psih);
    spsih = sin(psih);
    cdelta = Jh/Jg;
    sdelta = sqrt(1-cdelta^2);
    Ri2b,Rpsilpsig= sadov2ibrotmat(k,k1r,skmR,smkR,cnu,snu,dnu,cg,sg,cpsih,spsih,cdelta,sdelta);
    

    # euler angles and components of angular velocity
    omV =  [Jg*smkR*cnu/A; -Jg*smkR*k1r*snu/B; Jg*skmR*dnu/C]
    euangles = ATTITUDE.rotmatib2euler(Ri2b);
    return euangles,omV,Jl
end

"""
        sadov2euler(sadovvar,skm,IV)
        Transformation from Sadov variables to euler angles + the angular velocity 

        Consider an inertial reference frame XYZ.
        Consider a body reference frame xyz in principal axes of inertia. 
        The body-reference frame must be such that m belongs to [0,1] and k>0, where    
        k = C/A*(B-A)/(C-B)
        m = k*(C-delta)/(delta-A)
        with
            A,B,C = the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
            delta = G^2/2/Phi, with G the magnitude of the angular momentum of the body and Phi 
                the kinetic energy

        Let omV = [p,q,r]^T be the angular velocity of the satellite in the body frame
        
        Let N' be the node between the XY and xy planes.
        Let GV be angular momentum of the body (G = norm(GV)). 
        Let N be the node between the inertial reference XY plane and the plane
        perpendicular to GV; let N'' be the node between the plane perpendicular 
        to GV and the xy plane.

        INPUT: 
        sadovvar = [Jl,Jg,Jh,psil,psig,psih] or [Jg,Jh,psil,psig,psih]
        with
            Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m))
            Jg   = G
            Jh   = Gcos(delta), delta = inclination angle between the XY plane and the plane perpendicular to GV.
            psil = pi/2 ellK(m)ellF(lambda,m), lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
            psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)), g = angle between N and N''
            psih = angle between X axis and N
        skm = k/(m+k)
        IV  = [A,B,C]

        OUTPUT:
        euanglesV = [phi,theta,psi]
            where
            phi    = angle between X axis and N'
            theta  = inclination between the XY plane and the xy plane
            psi    = angle N' and x axis
        omV = [p,q,r], angular velocity in the body reference frame. 
        
        
"""
function sadov2euler(sadovvar,skm,IV)
    if size(sadovvar)[1] == 5
        euangles,omV,Jl = sadov2eulerAndJl(sadovvar,skm,IV)
    elseif size(sadovvar)[1] == 6
        euangles,omV,Jl = sadov2eulerAndJl(sadovvar[2:end],skm,IV)
    else
        Base.error("not enough input variables");
    end
    return euangles,omV; 
end

"""
    andoyer2sadov(euanglesV,omV,IV)
    Transformation from Andoyer-Serret variables to Sadov variables

    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    The reference frame is such that m belongs to [0,1] and k>0 where
    k = C/A*(B-A)/(C-B) 
    m = k*(C-delta)/(delta-A)
    with
    A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
    delta = G^2/2/Phi, with G the angular momentum of the body and Phi the kinetic energy


    Let GV=diag([A,B,C])omV be angular momentum of the body 
    (omV = [p,q,r]^T is the angular velocity) 
    
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular 
    to GV and the xy plane.
    
    INPUT: 
    # andoyervar = [L,G,H,l,g,h]
    #     where
    #     L = |GV|cos(sigma), sigma = inclination angle between
    #         the plane perpendicular to GV and the xy plane
    #     G = |GV|
    #     H = |GV|cos(delta), delta = inclination angle between 
    #         the XY plane and the plane perpendicular to GV.
    #     l = angle between N'' and x axis
    #     g = angle between N and N''
    #     h = angle between X axis and N
    # IV  = [A,B,C]
    
    OUTPUT
    #  sadovvar = [Jl,Jg,Jh,psil,psig,psih]
    #            Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m))
    #            Jg   = G
    #            Jh   = H
    #            psil = pi/2 ellK(m)ellF(lambda,m)
    #            psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m))
    #            psih = h
    #            where
    #            lambda = atan(cos(l),sqrt{1+k}sin(l)),
    #            ellF(lambda,m)    = int_{0}^{lambda} 1/sqrt{1-msin^2(theta)}dtheta
    #            ellP(-k,lambda,m) = int_{0}^{lambda} 1/sqrt{1-msin^2(theta)}/(1+ksin^2(theta))dtheta
    #            ellK(m) = ellF(pi/2,m), ellP(-k,m)=ellP(-k,pi/2,m)
    #  skm = k/(k+m)
"""
function andoyer2sadov(andoyervar,IV)   

    # input
    L   = andoyervar[1];
    G   = andoyervar[2];
    H   = andoyervar[3];
    l   = andoyervar[4];
    g   = andoyervar[5];
    h   = andoyervar[6];

    sl = sin(l);
    cl = cos(l);

    # k
    A = IV[1]; B = IV[2]; C = IV[3];
    k = C/A*(B-A)/(C-B); 
    if (C*(B-A)-A*(C-B)==0.0)
        k=0.0
    end

    # skm
    smk = (1-L^2/G^2)*(sl^2+cl^2/(1+k));
    skm = 1-smk;
    mtilde = smk/skm;
    if k==0.0
        m = 0.0;
    else
        m   = mtilde*k;
    end

    if m<0.0 || m>1.0
        println("m = ",m," mtilde ",mtilde," k ",k)
        Base.error("choice of reference frame wrong ",m)
    end

    skmR = sqrt(skm); 

    # Jl
    KK = ellF(pi/2.0,m);    #Elliptic.K(m);
    PP = ellP(-k,pi/2.0,m); #Elliptic.Pi(-k,pi/2,m);
    Jl = m2JlJgratioV2(skm,k,m,PP,KK)*G;

    # Jg
    Jg = G;

    # Jh 
    Jh = H;

    # psil
    slambda = -cl/sqrt(1.0+k*sl^2.0);
    clambda = sqrt(1.0+k)*sl/sqrt(1.0+k*sl^2.0);
    lambda  = mod(atan(slambda,clambda),2.0*pi)
    u = ellF(lambda,m);
    if m==1.0
        nlambda = floor(lambda*2.0/pi);
        psil = pi/2.0*nlambda;
    else
        psil = pi*u/2.0/KK;
    end 

    if m == 1.0
        psig = mod(g + sqrt(1+k)*atan(sqrt(k)*sin(lambda))/sqrt(mtilde+1)/sqrt(k),2*pi); 
    else
        eP = ellP(-k,lambda,m);
        eF = ellF(lambda,m); #Elliptic.F(lambda,m);
        if k==0.0
            psig = g;
        else
            psig = mod(g + sqrt(1+k)/skmR*(eP-PP*eF/KK),2.0*pi) #mod(g + mod(sqrt(1+k)/skmR*(eP-PP*eF/KK),2.0*pi),2*pi); 
        end
    end

    # psih 
    psih = h;

    # sadov
    sadovvar = [Jl;Jg;Jh;psil;psig;psih];
    
    return sadovvar,skm;
end

"""
    sadov2andoyerV1(sadovvar,skm,IV)
    Transformation from Sadov variables to Andoyer-Serret variables

    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    The reference frame is such that m belongs to in[0,1] and k>0 where
    k = C/A*(B-A)/(C-B) 
    m = k*(C-delta)/(delta-A)
    with
    A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
    delta = G^2/2/Phi, with G the angular momentum of the body and Phi the kinetic energy

    Let GV=diag([A,B,C])omV be angular momentum of the body 
    (omV = [p,q,r]^T is the angular velocity) 
    
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular 
    to GV and the xy plane.
    
    INPUT: 
    sadovvar = [Jl,Jg,Jh,psil,psig,psih] or [Jg,Jh,psil,psig,psih]
    #            Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m))
    #            Jg   = G
    #            Jh   = H
    #            psil = pi/2 ellK(m)ellF(lambda,m)
    #            psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m))
    #            psih = h
    #            where
    #            lambda = atan(cos(l),sqrt{1+k}sin(l)),
    #            ellF(lambda,m)    = int_{0}^{lambda} 1/sqrt{1-msin^2(theta)}dtheta
    #            ellP(-k,lambda,m) = int_{0}^{lambda} 1/sqrt{1-msin^2(theta)}/(1+ksin^2(theta))dtheta
    #            ellK(m) = ellF(pi/2,m), ellP(-k,m)=ellP(-k,pi/2,m)
    skm = k/(k+m)
    IV  = [A,B,C]
    
    OUTPUT
    andoyervar = [L,G,H,l,g,h]
        where
        L = |GV|cos(sigma), sigma = inclination angle between
            the plane perpendicular to GV and the xy plane
        G = |GV|
        H = |GV|cos(delta), delta = inclination angle between 
            the XY plane and the plane perpendicular to GV.
        l = angle between N'' and x axis
        g = angle between N and N''
        h = angle between X axis and N
"""
function sadov2andoyerV1(sadovvar,skm,IV)

    if (size(sadovvar)[1]==6)
        Jg = sadovvar[2];
        Jh = sadovvar[3];
        psil = sadovvar[4];
        psig = sadovvar[5];
        psih = sadovvar[6];
    elseif (size(sadovvar)[1]==5)
        Jg = sadovvar[1];
        Jh = sadovvar[2];
        psil = sadovvar[3];
        psig = sadovvar[4];
        psih = sadovvar[5];
    else
        Base.error("wrong input");
    end

    A = IV[1]; B = IV[2]; C = IV[3];
    k = C/A*(B-A)/(C-B);

    if (C*(B-A)-A*(C-B)==0.0)
        k=0.0
    end

    k1r = sqrt(1+k);
    smk = 1-skm;
    skmR = sqrt(skm);

    # m
    if skm == 0.0
        m = 0.0;
    else
        m = smk*k/skm;
        if m>1.0 && abs(m-1)<1e-13
            m = 1.0
        end
    end
    KK = ellF(pi/2,m)
    PP = ellP(-k,pi/2,m)

    # lambda and g
    if m == 1.0
        lambda = psil;
        snu = sin(lambda);
        cnu = cos(lambda);
        dnu = sqrt(1-m*snu^2);
        g = mod(psig - sqrt(1+k)*atan(sqrt(k)*sin(lambda))*skmR/sqrt(k),2*pi); 
    else
        u = 2.0*psil*KK/pi;
        snu = Elliptic.Jacobi.sn(u,m);
        cnu = Elliptic.Jacobi.cn(u,m);
        dnu = Elliptic.Jacobi.dn(u,m);
        lambda = mod(atan(snu,cnu),2.0*pi); 
        eP = ellP(-k,lambda,m); 
        g = mod(psig - sqrt(1+k)/skmR*(eP-PP*u/KK),2.0*pi); #mod(psig - mod(sqrt(1+k)/skmR*(eP-u*PP/KK),2.0*pi),2.0*pi)
    end 


    # andoyer
    L = Jg*skmR*dnu; 
    G = Jg; 
    H = Jh; 
    l = mod(atan(cnu,-k1r*snu),2*pi);   
    h = psih;
    
    andoyervar = [L;G;H;l;g;h]; 

    return andoyervar;
end

"""
    sadov2andoyerV2(sadovvar,IV)
    Transformation from Andoyer-Serret variables to Sadov variables
    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    The reference frame is such that m belongs to [0,1] and k>0 where
    k = C/A*(B-A)/(C-B) 
    m = k*(C-delta)/(delta-A)
    with
    A,B,C the moments of inertia, A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm.
    delta = G^2/2/Phi, with G the angular momentum of the body and Phi the kinetic energy

    Let GV=diag([A,B,C])omV be angular momentum of the body 
    (omV = [p,q,r]^T is the angular velocity) 
    
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular 
    to GV and the xy plane.
    
    INPUT: 
    sadovvar = [Jl,Jg,Jh,psil,psig,psih]
    #            Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m))
    #            Jg   = G
    #            Jh   = H
    #            psil = pi/2 ellK(m)ellF(lambda,m)
    #            psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m))
    #            psih = h
    #            where
    #            lambda = atan(cos(l),sqrt{1+k}sin(l)),
    #            ellF(lambda,m)    = int_{0}^{lambda} 1/sqrt{1-msin^2(theta)}dtheta
    #            ellP(-k,lambda,m) = int_{0}^{lambda} 1/sqrt{1-msin^2(theta)}/(1+ksin^2(theta))dtheta
    #            ellK(m) = ellF(pi/2,m), ellP(-k,m)=ellP(-k,pi/2,m)
    IV  = [A,B,C]
    
    OUTPUT
    andoyervar = [L,G,H,l,g,h]
        where
        L = |GV|cos(sigma), sigma = inclination angle between
            the plane perpendicular to GV and the xy plane
        G = |GV|
        H = |GV|cos(delta), delta = inclination angle between 
            the XY plane and the plane perpendicular to GV.
        l = angle between N'' and x axis
        g = angle between N and N''
        h = angle between X axis and N
"""
function sadov2andoyerV2(sadovvar,IV)    

    Jl = sadovvar[1]; 
    Jg = sadovvar[2]; 
    
    A = IV[1]; B = IV[2]; C = IV[3];
    k = C/A*(B-A)/(C-B);
    if (C*(B-A)-A*(C-B)==0.0)
        k=0.0
    end

    # m, skm
    if k == 0 
        m = 0.0;
        mtilde = (Jg/Jl)^2-1.0  
    else
        m=JlJgratio2m_bisection(Jl,Jg,k);    
        mtilde=m/k;
    end

    skm = 1/(1+mtilde);
    andoyervar = sadov2andoyerV1(sadovvar,skm,IV);

    return andoyervar;
end

"""
    sadov2sadovlike(sadovvar,skm)    
    Transformation from Sadov variables to Sadov-like variables 
    (Sadov-like variables = variables regularised with respect to the singularity sin(delta)=0, 
    when the angular momentum is aligned with the inertial Z axis).

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
    
    INPUT: 
    sadov variables = [Jl,Jg,Jh,psil,psig,psih] or [Jg,Jh,psil,psig,psih]
                with
                Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m))
                Jg   = G
                Jh   = Gcos(delta), delta = inclination angle between 
                       the XY plane and the plane perpendicular to GV.
                psil = pi/2 ellK(m)ellF(lambda,m), lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
                psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)), g = angle between N and N''
                psih = angle between X axis and N

    skm = k/(m+k)
    
    OUTPUT: 
    sadov-like variables = [J1,J2,J3,J4,J5,J6,J7]
    with J1 = skm
         J2 = Jg
         J3 = Jh
         J4 = psil
         J5 = psig + Jh/Jg*psih
         J6 = sqrt(Jg^2-Jh^2)*cos(psih)
         J7 = sqrt(Jg^2-Jh^2)*sin(psih)
"""
function sadov2sadovlike(sadovvar,skm)
    if (size(sadovvar)[1]==6)
        Jg = sadovvar[2];
        Jh = sadovvar[3];
        psil = sadovvar[4];
        psig = sadovvar[5];
        psih = sadovvar[6];
    elseif (size(sadovvar)[1]==5)
        Jg = sadovvar[1];
        Jh = sadovvar[2];
        psil = sadovvar[3];
        psig = sadovvar[4];
        psih = sadovvar[5];
    else
        Base.error("wrong input");
    end

    cdelta = Jh/Jg;
    if abs(cdelta)>1.0 && abs(abs(cdelta)-1.0)<1e-15 
        cdelta = cdelta/abs(cdelta);
    elseif abs(cdelta)>1.0
        Base.error("|cos(delta)|>1 : wrong inputs");
    end
    sdelta = sqrt(1.0-cdelta^2.0);

    J1 = skm;
    J2 = Jg;
    J3 = Jh;
    J4 = psil;
    J5 = psig + cdelta*psih; #atan(tan(psih));
    J6 = Jg*sdelta*cos(psih);
    J7 = Jg*sdelta*sin(psih);

    sadovlike = [J1,J2,J3,J4,J5,J6,J7];
    return sadovlike;
end

"""
    sadov2sadovlike(sadovvar,skm)    
    Transformation from Sadov-like variables to Sadov variables 
    (Sadov-like variables = variables regularised with respect to the singularity sin(delta)=0, 
    when the angular momentum is aligned with the inertial Z axis).

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
    
    INPUT: 
    sadov-like variables = [J1,J2,J3,J4,J5,J6,J7]
    with J1 = skm
         J2 = Jg
         J3 = Jh
         J4 = psil
         J5 = psig + Jh/Jg*psih
         J6 = sqrt(Jg^2-Jh^2)*cos(psih)
         J7 = sqrt(Jg^2-Jh^2)*sin(psih)
    where 
    skm = k/(m+k)
    Jg   = G
    Jh   = Gcos(delta), delta = inclination angle between 
           the XY plane and the plane perpendicular to GV.
    psil = pi/2 ellK(m)ellF(lambda,m), lambda = atan(cos(l),sqrt{1+k}sin(l)), l = angle between N'' and x axis
    psig = g - sqrt{(k+m)(1+k)/k}(ellP(-k,m)ellF(lambda,m)/ellK(m)-ellP(-k,lambda,m)),  g = angle between N and N''
    psih = angle between X axis and N

    OUTPUT
    sadov variables = [Jl,Jg,Jh,psil,psig,psih] 
                with
                Jl   = 2G/pi sqrt{k+m} sqrt{(1+k)/k} (ellP(-k,m)-m/(k+m)ellK(m))
    skm
"""
function sadovlike2sadov(sadovlike,IV)
    J1 = sadovlike[1]
    J2 = sadovlike[2]
    J3 = sadovlike[3]
    J4 = sadovlike[4]
    J5 = sadovlike[5]
    J6 = sadovlike[6]
    J7 = sadovlike[7]

    skm = J1;
    Jg  = J2;
    Jh  = J3;
    psil = J4;

    cdelta = Jh/Jg;
    if abs(cdelta)>1.0 && abs(abs(cdelta)-1.0)<1e-15 
        cdelta = cdelta/abs(cdelta);
    elseif abs(cdelta)>1.0
        Base.error("|cos(delta)|>1 : wrong inputs");
    end
    sdelta = sqrt(1.0-cdelta^2.0);
     
    if sdelta==0.0
        psih = 0.0;
        psig = mod(J5,2.0*pi); #mod(J5-2*pi,2*pi);
    else
        psih = mod(atan(J7,J6),2*pi);
        psig = mod(J5 - cdelta*psih,2.0*pi); #atan(tan(psih)),2*pi);
    end

    A = IV[1]; B = IV[2]; C = IV[3];
    k = C/A*(B-A)/(C-B); 
    if (C*(B-A)-A*(C-B)==0.0)
        k=0.0
    end
    if k==0.0 
        m = 0.0
    else
        m = k*(1-skm)/skm;
        if m>1.0 && abs(m-1)<1e-15
            m = 1.0
        end
    end
    PP = ellP(-k,pi/2.0,m);
    KK = ellF(pi/2.0,m);
    Jl = Jg*m2JlJgratioV2(skm,k,m,PP,KK)

    sadov = [Jl,Jg,Jh,psil,psig,psih];
    return sadov,skm;
end

"""
    andoyer2andoyerlike(andoyervar)   
    Transformation from Andoyer-Serret variables to Andoyer-like variables 
    (Andoyer-like variables = variables regularised with respect to the singularity sin(delta)=0, 
    when the angular momentum is aligned with the inertial Z axis)

    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    
    Let GV=diag([A,B,C])omV be angular momentum of the body 
    (omV = [p,q,r]^T is the angular velocity) 
    
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular
    to GV and the xy plane.
    
    INPUT: 
    andoyervar = [L,G,H,l,g,h]
        where
        L = |GV|cos(sigma), sigma = inclination angle between
            the plane perpendicular to GV and the xy plane
        G = |GV|
        H = |GV|cos(delta), delta = inclination angle between 
            the XY plane and the plane perpendicular to GV.
        l = angle between N'' and x axis
        g = angle between N and N''
        h = angle between X axis and N
    
    
    OUTPUT: 
    andoyer-like variables = [J1,J2,J3,J4,J5,J6,J7]
    with J1 = L
         J2 = G
         J3 = H
         J4 = l
         J5 = g + H/G*h
         J6 = sqrt(G^2-H^2)*cos(h)
         J7 = sqrt(G^2-H^2)*sin(h)
"""
function andoyer2andoyerlike(andoyervar)
    LA = andoyervar[1]
    GA = andoyervar[2]
    HA = andoyervar[3]
    lA = andoyervar[4]
    gA = andoyervar[5]
    hA = andoyervar[6]

    cdelta = HA/GA;
    if abs(cdelta)>1.0 && abs(abs(cdelta)-1.0)<1e-15 
        cdelta = cdelta/abs(cdelta);
    elseif abs(cdelta)>1.0
        Base.error("|cos(delta)|>1 : wrong inputs");
    end
    if abs(abs(cdelta)-1.0)<1e-15
        cdelta = cdelta/abs(cdelta)
    end
    sdelta = sqrt(1.0-cdelta^2.0);

    JA1 = LA;
    JA2 = GA;
    JA3 = HA;
    JA4 = lA;
    JA5 = gA + cdelta*hA #atan(tan(psih));
    JA6 = GA*sdelta*cos(hA);
    JA7 = GA*sdelta*sin(hA);

    andoyerlike = [JA1,JA2,JA3,JA4,JA5,JA6,JA7];
    return andoyerlike;
end

"""
    andoyerlike2andoyer(andoyerlike)   
    Transformation from Andoyer-like variables to Andoyer-Serret variables
    (Andoyer-like variables = variables regularised with respect to the singularity sin(delta)=0, 
    when the angular momentum is aligned with the inertial Z axis)

    Consider an inertial reference Frame XYZ.
    Consider a body reference frame xyz in principal axes of inertia. 
    Let A,B,C be the principal moments of inertia with 
    A=int(y^2+z^2)dm, B=int(x^2+z^2)dm and C=int(x^2+y^2)dm. 
    
    Let GV=diag([A,B,C])omV be angular momentum of the body 
    (omV = [p,q,r]^T is the angular velocity) 
    
    Let N be the node between the inertial reference XY plane and the plane
    perpendicular to GV; let N'' be the node between the plane perpendicular
    to GV and the xy plane.
    
    INPUT: 
    andoyer-like variables = [J1,J2,J3,J4,J5,J6,J7]
    with J1 = L
         J2 = G
         J3 = H
         J4 = l
         J5 = g + H/G*h
         J6 = sqrt(G^2-H^2)*cos(h)
         J7 = sqrt(G^2-H^2)*sin(h)
    where
    L = |GV|cos(sigma), sigma = inclination angle between
            the plane perpendicular to GV and the xy plane
    G = |GV|
    H = |GV|cos(delta), delta = inclination angle between 
           the XY plane and the plane perpendicular to GV.
    l = angle between N'' and x axis
    g = angle between N and N''
    h = angle between X axis and N

    OUTPUT: 
    andoyervar = [L,G,H,l,g,h]
"""
function andoyerlike2andoyer(andoyerlike)
    JA1 = andoyerlike[1]
    JA2 = andoyerlike[2]
    JA3 = andoyerlike[3]
    JA4 = andoyerlike[4]
    JA5 = andoyerlike[5]
    JA6 = andoyerlike[6]
    JA7 = andoyerlike[7]

    LA = JA1;
    GA = JA2;
    HA = JA3;
    lA = JA4;

    cdelta = HA/GA;
    if abs(cdelta)>1.0 && abs(abs(cdelta)-1.0)<1e-12 
        cdelta = cdelta/abs(cdelta);
    elseif abs(cdelta)>1.0
        println(cdelta)
        Base.error("|cos(delta)|>1 : wrong inputs");
    end
    sdelta = sqrt(1.0-cdelta^2.0);
     
    if sdelta==0.0
        hA = 0.0;
        gA = mod(JA5,2*pi);
    else
        hA = mod(atan(JA7,JA6),2*pi);
        gA = mod(JA5 - cdelta*hA,2.0*pi); #atan(tan(psih)),2*pi);
    end

    andoyervar = [LA,GA,HA,lA,gA,hA];
    return andoyervar;
end


#  useful functions
function  m2JlJgratioV1(m,k,PP,KK)   
    if m==1.0
        PPmPkMKKm=atan(sqrt(k))*sqrt(k);
    else
        PPmPkMKKm=PP*(m+k)-m*KK;
    end
    return 2.0*sqrt(1.0 + k)*PPmPkMKKm/(pi*sqrt(k)*sqrt(k+m));
end

function  m2JlJgratioV2(skm,k,m,PP,KK)   
    if m==1.0
        ratio = 2.0/pi*atan(sqrt(k));
    else
        if k==0.0
            ratio = sqrt(skm);
        else
            ratio = 2.0/pi*sqrt(1+k)/sqrt(skm)*(PP-(1-skm)*KK);
        end
    end
    return ratio;
end

function JlJgratio2m_bisection(Jl,Jg,k)
   maxiter = 67;
   tol = 1e-20;

   JlJgratio = Jl/Jg;
   
   mA = 1e-30;
   PPA = ellP(-k,pi/2,mA);
   KKA = ellF(pi/2,mA); 
   funA = m2JlJgratioV1(mA,k,PPA,KKA)-JlJgratio;

   if abs(funA) <= tol
    return mA;
   end

   mB = 1.0-1e-30;
   PPB = ellP(-k,pi/2,mB); #Elliptic.Pi(-k,pi/2.0,mB);
   KKB = ellF(pi/2,mB); #Elliptic.K(mB);
   funB = m2JlJgratioV1(mB,k,PPB,KKB)-JlJgratio;

   if abs(funB) <= tol
    return mB;
   end

   m = 0.0;

   for jj = 1:maxiter
    m = (mA+mB)/2.0;
    PP = ellP(-k,pi/2,m); #Elliptic.Pi(-k,pi/2.0,m);
    KK = ellF(pi/2,m); #Elliptic.K(m);
    fun = m2JlJgratioV1(m,k,PP,KK)-JlJgratio;
    if abs(fun) <= tol
        break;
    else
        if fun*funA < 0.0
            mB = m;
        else
            mA = m;
            funA = fun;
        end
    end
   end

   return m;
end

function dJlJgratio_dm(m,k,KK)
    return KK*sqrt(k)*sqrt(1.0+k)/((m+k)^(3.0/2.0))/pi;
end

function JlJgratio2m_newton(Jl,Jg,k,m0)
    maxiter = 1000;
    tol = 1e-20;
    m = m0;
    ier = 0;
    JlJgratio = Jl/Jg;
   
    for jj = 1 : maxiter
        PP = Elliptic.Pi(-k,pi/2.0,m);
        KK = Elliptic.K(m);
        fun = m2JlJgratioV1(m,k,PP,KK)-JlJgratio;
        
        if abs(fun)<=tol
           ier = 1;
            break;
        else
            dfun = dJlJgratio_dm(m,k,KK);
            if dfun!=0.0
                m = max(m-fun/dfun,0.0);
            else
                ier = -1;
                break;
            end
        end
    end

    if ier == 0
        ier = -1
    end
    
    return m,ier;

end


