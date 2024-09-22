"""
        kep2car(kepV,mu)
        Transformation from Cartesian velocity and position vector to Keplerian elements

        INPUT:
        kepV=[a,e,i,o,O,nu], with
                a   semi-major axis                        
                e   orbital eccentricity                   
                i   inclination                            [rad]
                o   argumentum of pericentre               [rad]
                O   right ascension of the ascending node  [rad]
                nu  true anomaly                           [rad]
        mu           planetary constant
        !!!! dimensions of a and mu must be coherent
        OUPUT:
        carV = [v.x]
        x      vector position in the cartesian ref 
        v      vector velocity in the cartesian ref 
"""
function kep2car(kepV,mu)
    # kep2car - calculates position and velocity vectors starting from orbital
    # DESCRIPTION:
    #   This function returns two vectors representing the state of the s/c 
    #   (position and velocity) given the orbital parameters.
    # INPUT:
    #   kepV         =[a,e,i,o,O,nu], with
    #                 a   semi-major axis                        
    #                 e   orbital eccentricity                   
    #                 i   inclination                            [rad]
    #                 o   argumentum of pericentre               [rad]
    #                 O   right ascension of the ascending node  [rad]
    #                 nu  true anomaly                           [rad]
    #   mu           planetary constant
    #   !!!! dimensions of a and mu must be coherent
    # OUPUT:
    #  v[3,1]      vector velocity in the cartesian ref [dimension :
    #               sqrt(dimension of mu)/sqrt(dimension of a)]
    #   x[3,1]      vector position in the cartesian ref [dimension of a]
    #   
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    sma  = kepV[1];
    ec  = kepV[2];
    incl  = kepV[3];
    om  = kepV[4];
    OM  = kepV[5];
    nu = kepV[6];
    
    cnu = cos(nu); snu = sin(nu);
    pval = sma*(1.0-ec^2.0);
    den = 1.0+ec*cnu; 
    xpf = [pval*cnu/den; pval*snu/den; 0.0];
    vpf = sqrt(mu/pval)*[-snu; (ec+cnu); 0];
    ROM = [cos(OM) -sin(OM) 0.0; sin(OM) cos(OM) 0.0; 0.0 0.0 1.0];
    Ri  = [1.0 0.0 0.0; 0.0 cos(incl) -sin(incl); 0.0 sin(incl) cos(incl)];
    Rom = [cos(om) -sin(om) 0.0; sin(om) cos(om) 0.0; 0.0 0.0 1.0];
    RR = ROM*Ri*Rom;
    xv = RR*xpf;
    vv = RR*vpf;
    carV = append!(vv,xv);
    return carV;
end

"""
    car2kep(carVect,mu)
    Transformation from Keplerian elements to Cartesian velocity and position vector 
    
    INPUT:
        carV = [vv,xv]
          vv[3,1]  vector velocity      
          xv[3,1]  vector position                   
        mu       planetary constant                    
    
    OUPUT:
        kepV = [a,e,i,o,O,nu]
                  a        semi-major axis                       [same as xv]
                  e        orbital eccentricity                  
                  i        inclination                           [rad]
                  O        right ascension of the ascending node [rad]
                  o        argumentum of pericentre              [rad]
                  nu       true anomaly                          [rad]
"""
function car2kep(carVect,mu)
    # INPUT:
    #   carV = [vv,xv]
    #      vv[3,1]  vector velocity      
    #      xv[3,1]  vector position                   
    #  mu       planetary constant                    
    #
    # OUPUT:
    #  kepV = [a,e,i,o,O,nu]
    #              a        semi-major axis                       [same as xv]
    #              e        orbital eccentricity                  
    #              i        inclination                           [rad]
    #              O        right ascension of the ascending node [rad]
    #              o        argumentum of pericentre              [rad]
    #              nu       true anomaly                          [rad]
    # ------------------------------------------------------------------------
  
    eps = 1e-12;
   
    # position and velocity
    vv = carVect[1:3];  
    xv = carVect[4:6];

    rn=norm(xv); vn=norm(vv); 
    
    # radial velocity
    vr = dot(xv,vv);

    # angular momentum
    hv= cross(xv,vv); hn=norm(hv);
    
    #energy
    en=0.5*vn^2.0-mu/rn;
    
    # semi-major axis a
    sma = -mu/(2.0*en);
    
    # inclination i
    cincl = hv[3]/hn;
    if abs(cincl)<=1.0 
        incl=acos(cincl);
    elseif cincl>1.0 && cincl-1.0<eps
        incl = 0.0;
    elseif cincl<-1.0 && -cincl-1.0<eps
        incl = pi
    else
        Base.error(["|cos(incl)|>1: something wrong"]);
    end

    # orbital eccentricity (vector ev and modulus e)
    ev=(1.0/mu)*((vn^2.0-mu/rn)*xv - vr*vv);
    ec=norm(ev);
    if ec<eps
        ec = 0.0;
    end

    # line of nodes
    lon = [0.0,0.0,0.0];
    if incl!=0.0 && incl!=pi
        lon[1] = -hv[2]
        lon[2] = hv[1]
    end
    lon = lon;
    lonn = norm(lon);

    # right ascension of the ascending node
    if lonn!=0.0
        cOM = lon[1]/lonn;
        if abs(cOM)<=1.0
            OM = acos(cOM);
        elseif cOM>1.0 && cOM-1.0<eps
            OM = 0.0;
        elseif cOM<-1.0 && -cOM-1.0<eps
            OM = pi;
        else
            Base.error(["|cos(OM)|>1: something wrong"]);
        end
        if lon[2] < 0.0
            OM = 2.0*pi-OM;
        end
    else
        OM = 0.0;
    end

    # argument of pericenter
    if ec == 0.0
        om = 0.0;
    else
        if lonn!=0.0 
            com = dot(ev,lon)/ec/lonn;
            if abs(com)<=1.0
                om = acos(com);
            elseif com>1.0 && com-1.0<eps
                om = 0.0;
            elseif com<-1.0 && -com-1.0<eps
                om = pi;
            else
                Base.error(["|cos(om)|>1: something wrong"]);
            end
            if ev[3]<0
                om = 2.0*pi-om;
            end
        else
            com = ev[1]/ec;
            if abs(com)<=1.0
                om = acos(com);
            elseif com>1.0 && com-1.0<eps
                om = 0.0;
            elseif com<-1.0 && -com-1.0<eps
                om = pi;
            else
                Base.error(["|cos(om)|>1: something wrong"]);
            end
            om = acos(com);
            if ev[2]<0.0
                om = 2.0*pi-om;
            end
        end
    end

    # true anomaly
    if ec!=0
        cnu = dot(xv,ev)/rn/ec;
        if abs(cnu)<=1.0
            nu = acos(cnu);
        elseif cnu>1.0 && cnu-1.0<eps
            nu = 0.0;
        elseif cnu<-1.0 && -cnu-1.0<eps
            nu = pi;
        else
            Base.error(["|cos(nu)|>1: something wrong"]);
        end
        if vr<0
            nu=2.0*pi-nu; 
        end 
    elseif lonn !=0.0
        cnu = dot(xv,lon)/rn/lonn;
        if abs(cnu)<=1.0
            nu = acos(cnu);
        elseif cnu>1.0 && cnu-1.0<eps
            nu = 0.0;
        elseif cnu<-1.0 && -cnu-1.0<eps
            nu = pi;
        else
            Base.error(["|cos(nu)|>1: something wrong"]);
        end
        if dot(xv,cross(hv,lon))<0.0
            nu = 2.0*pi-nu; 
        end
    else
        cnu = xv[1]/rn;
        if abs(cnu)<=1.0
            nu = acos(cnu);
        elseif cnu>1.0 && cnu-1.0<eps
            nu = 0.0;
        elseif cnu<-1.0 && -cnu-1.0<eps
            nu = pi;
        else
            Base.error(["|cos(nu)|>1: something wrong"]);
        end
         if xv[2]<0.0
            nu = 2.0*pi-nu; 
        end
    end

    kepV = [sma,ec,incl,om,OM,nu];
    return kepV;

end   

"""
    kep2equi(kepV,idxLL)
    Transformation from Keplerian elements to equinoctial elements

    INPUT:
    kepV=[a,e,i,o,O,nu], with
            a   semi-major axis                        
            e   orbital eccentricity                   
            i   inclination                            [rad]
            o   argumentum of pericentre               [rad]
            O   right ascension of the ascending node  [rad]
            nu  true anomaly                           [rad]
    idxLL = which longitude in output : 1 - true longitude (nu+o+O)
                                        2 - eccentric longitude (Ean+o+O), Ean eccentric anomaly
                                        3 - mean longitude (M+o+O), M mean anomaly
    OUPUT:
    equiV = [a,P1,P2,Q1,Q2,LL]                      
        P1 = e*sin(OM+om)                   
        P2 = e*cos(OM+om)                         
        Q1 = tan(i/2)*sin(OM)  
        Q2 = tan(i/2)*cos(OM)  
        LL = nu+o+O [rad]  if idxLL = 1
            Ean+o+O [rad] if idxLL = 2
            M+o+O [rad]   if idxLL = 3
"""
function kep2equi(kepV,idxLL)

    sma  = kepV[1];
    ec  = kepV[2];
    incl  = kepV[3];
    om  = kepV[4];
    OM  = kepV[5];
    nu = kepV[6];
    
    P1 = ec*sin(om+OM);
    P2 = ec*cos(om+OM);
    Q1 = tan(incl/2.0)*sin(OM);
    Q2 = tan(incl/2.0)*cos(OM);

    if idxLL == 1 || idxLL == 2
        TL = nu+om+OM;
        LL = TL;
        if idxLL == 2
            EL =  truelong2ecclong(TL,P1,P2);
            LL = EL;
        end
    else
        man = true2mean(nu,ec);
        ML = man + om + OM;
        LL = ML;
    end

    equiV = [sma,P1,P2,Q1,Q2,LL];
    return equiV;

end

"""
    equi2kep(equiV,idxLL)
    Transformation from equinoctial elements to Keplerian elements

    INPUT:
        equiV = [a,P1,P2,Q1,Q2,LL]  
            a   semi-major axis              
            P1 = e*sin(OM+om)                   
            P2 = e*cos(OM+om)                         
            Q1 = tan(i/2)*sin(OM)  
            Q2 = tan(i/2)*cos(OM)  
            LL = nu+o+O [rad]  if idxLL = 1
                 Ean+o+O [rad] if idxLL = 2
                 M+o+O [rad]   if idxLL = 3
            with 
                e   orbital eccentricity                   
                i   inclination                            [rad]
                o   argumentum of pericentre               [rad]
                O   right ascension of the ascending node  [rad]
                nu  true anomaly                           [rad]
                M   mean naomaly                           [rad]
                Ean eccentric anomaly                      [rad]
        idxLL = identifier of the input longitude in the set of equinoctial elements

    

    OUPUT:
    kepV=[a,e,i,o,O,nu]

"""
function equi2kep(equiV,idxLL)

    sma = equiV[1];
    P1  = equiV[2];
    P2  = equiV[3];
    Q1  = equiV[4];
    Q2  = equiV[5];
    if idxLL == 1
        TL  = equiV[6];
    elseif idxLL ==2
        TL = ecclong2truelong(equiV[6],P1,P2);
    else
        TL =  meanlong2truelong(equiV[6],P1,P2)
    end    

    ec = sqrt(P1^2.0+P2^2.0);
    incl = mod(atan(sqrt(Q1^2.0+Q2^2.0))*2.0,2.0*pi);

    if ec!=0.0 && incl!=0.0
        omPOM  = mod(atan(P1,P2),2*pi);
        OM     = mod(atan(Q1,Q2),2*pi);
        om     = mod(omPOM-OM,2*pi);
        nu     = mod(TL-omPOM,2.0*pi);
    elseif ec!=0.0 && incl==0.0
        OM  = 0.0;
        om  = mod(atan(P1,P2),2*pi);
        nu  = mod(TL-om,2.0*pi);
    elseif ec==0.0 && incl!=0.0
        OM  = mod(atan(Q1,Q2),2*pi);
        om  = 0.0;
        nu  = mod(TL-OM,2.0*pi);
    else
        OM  = 0.0;
        om  = 0.0;
        nu  = mod(TL,2.0*pi);
    end

    kepV = [sma,ec,incl,om,OM,nu];
    return kepV;

end

"""
    truelong2ecclong(TL,P1,P2)
    from true to eccentric longitude
    INPUT:
        TL = nu+o+O (true longitude)   [rad]        
        P1 = e*sin(O+o)                   
        P2 = e*cos(O+o)                       
            with 
                e   orbital eccentricity                   
                o   argumentum of pericentre               
                O   right ascension of the ascending node  
                nu  true anomaly                           
    OUTPUT
        EL = E+o+O (eccentric longitue) [rad]
            with E eccentric anomaly
"""
function truelong2ecclong(TL,P1,P2)

    k = 1.0/(1.0+sqrt(1-P1.^2.0-P2.^2.0));
    e2 = P1.^2+P2.^2;
    sTL = sin(TL);
    cTL = cos(TL);
    
    phi = (1.0 - P1.^2.0 - P2.^2.0)./(1.0 + P1.*sTL + P2.*cTL);

    sEL = P1 + phi.*(-sTL + k.*P1.*(cTL.*P2+P1.*sTL))./(-1 + k.*e2);
    cEL = P2 + phi.*(cTL.*(-1+k.*P2.^2)+k.*P1.*P2.*sTL)./(-1+k.*e2);

    EL   = mod(atan(sEL,cEL),2.0*pi);

    return EL;
end

"""
    ecclong2truelong(EL,P1,P2)
    from eccentric to true longitude
    INPUT:
        EL = E+o+O (eccentric longitude)   [rad]        
        P1 = e*sin(O+o)                   
        P2 = e*cos(O+o)                         
            with 
                e   orbital eccentricity                   
                o   argumentum of pericentre              
                O   right ascension of the ascending node 
                E   eccentric anomaly                     
    OUTPUT
        TL = nu+o+O (true longitue) [rad]
            with nu true anomaly
"""
function ecclong2truelong(EL,P1,P2)
    sEL = sin(EL);
    cEL = cos(EL);
    rOsma = 1-P1*sEL-P2*cEL;
    etaP1 = 1.0 + sqrt(1-P1^2-P2^2);

    sL = (-P1*etaP1 + P1*P2*cEL + (etaP1-P2^2)*sEL)/etaP1/rOsma;
    cL = (-P2*etaP1 + P1*P2*sEL + (etaP1-P1^2)*cEL)/etaP1/rOsma;
    TL = mod(atan(sL,cL),2.0*pi);

    return TL;
end

"""
    equi2equationofthecentre(TL,Q1,Q2,P1,P2)
    from the true longitude to the equation of the centre
    
    INPUT:
        TL = nu+om+OM                      
        Q1 = tan(i/2)*sin(O)  
        Q2 = tan(i/2)*cos(O)
        P1 = e*sin(O+o)                   
        P2 = e*cos(O+o)  
        with 
                e   orbital eccentricity                   
                i   inclination                            [rad]
                o   argumentum of pericentre               [rad]
                O   right ascension of the ascending node  [rad]
                nu  true anomaly                           [rad]        
    
    OUTPUT
        equation of the centre = nu-M, M mean anomaly
"""
function equi2equationofthecentre(TL,Q1,Q2,P1,P2)

    ec = sqrt(P1^2.0+P2^2.0);
    incl = mod(atan(sqrt(Q1^2.0+Q2^2.0))*2.0,2.0*pi);

    if ec!=0.0 && incl!=0.0
        omPOM  = mod(atan(P1,P2),2*pi);
        OM     = mod(atan(Q1,Q2),2*pi);
        om     = mod(omPOM-OM,2*pi);
        nu     = mod(TL-omPOM,2.0*pi);
    elseif ec!=0.0 && incl==0.0
        OM  = 0.0;
        om  = mod(atan(P1,P2),2*pi);
        nu  = mod(TL-om,2.0*pi);
    elseif ec==0.0 && incl!=0.0
        OM  = mod(atan(Q1,Q2),2*pi);
        om  = 0.0;
        nu  = mod(TL-OM,2.0*pi);
    else
        OM  = 0.0;
        om  = 0.0;
        nu  = mod(TL,2.0*pi);
    end

    M = true2mean(nu,ec);

    phiTA = nu-M;

    return phiTA;


end

"""
    true2mean(nu,ec)
    from true anomaly to mean anomaly

    INPUT:
        nu = true anomaly [rad]
        ec = eccentricity
    OUTPUT
        M = mean anomaly [rad]

    Reference: Farnocchia, Cioci, Milani, 2012
"""
function true2mean(nu,ec)
    # from Farnocchia, Cioci, Milani 2012
    M = nu;
    if ec==0 
        M = nu;
    elseif ec>0 && ec<1 
        cnu = cos(nu)
        snu = sin(nu)
        dean = 1+ec*cnu;
        if abs(dean)<1e-10
            ean = mod(2*atan(sqrt(1-ec)*tan(nu/2)/sqrt(1+ec)),2*pi);
            sean = sin(ean)
        else
            sean = sqrt(1-ec^2)*snu/dean;
            cean = (ec+cnu)/dean;      
            ean = mod(atan(sean,cean),2*pi)
        end
        M = ean-ec*sean
    else
        Base.error("method suitable only for either circular or elliptic orbit (ec>1)")
    end
    return M;
end

function eccentric2mean(ean,ec)
    # from Farnocchia, Cioci, Milani 2012
    M = ean;
    if ec==0 
        M = ean;
    elseif ec>0 && ec<1 
        sean = sin(ean)
        M = ean-ec*sean
    else
        Base.error("method suitable only for either circular or elliptic orbit (ec>1)")
    end
    return M;
end

"""
    mean2true(M,ec)
    from true anomaly to mean anomaly
    (Farnocchia, Cioci, Milani 2012)
    
    INPUT:
        M = mean anomaly [rad]
        ec = eccentricity
    OUTPUT
        nu = mean anomaly [rad]
        ier =  1 --> conversion succeded
              -1 --> conversion failed


    Authors: Irene Cavallari, Clara Grassi
"""
function mean2trueandeccentric(M,ec)
    # from Farnocchia, Cioci, Milani 2012
    ier = 0
    nu = M
    del = 1e-10;
    if ec>= 1.0
        Base.error("method suitable only for either circular or elliptic orbit (ec>1)")
    end
    if ec==0.0
       # circular =====================================================
       nu = M
       ean = M
       ier = 1
    elseif ec<1
       # elliptic =====================================================
       if(ec>=1-del)    
          ean = acos((1-del)/ec);
          if(M>pi)
             ean = 2*pi-ean;
          end 
          sean = sin(ean)
          chiE = ean-ec*sean
          if(abs(M-chiE)<1e-10) 
             cean = cos(ean)
             dnu = 1-ec*cean
             cnu = (cean-ec)/dnu
             snu = sqrt(1-ec^2)*sean/dnu
             nu = mod(atan(snu,cnu),2*pi);
             ier = 1
          elseif(M>chiE) 
             if(M<=pi) 
                nu,ean,ier = ellipticstrong(M,e)
             else
                nu,ean,ier = ellipticnearparabolic(M,ec);
             end
          else
             if(M<pi) 
                nu,ean,ier = ellipticnearparabolic(M,ec)
             else
                nu,ean,ier = ellipticstrong(M,ec)
             end
          end
       else
        nu,ean,ier =  ellipticstrong(M,ec)
       end
    end
    return nu,ean,ier;
end

function Sfun(ec,D)
    x = (ec-1.0)/(ec+1.0)*D^2.0
    Sval = 0.0
    for kk=0:100
       Sval = Sval + (ec - 1.0/(2.0*kk+3.0))*(x^kk)
    end
    return Sval;
end 

function DSfun(ec,D)
    x = (ec-1.0)/(ec+1.0)*D^2.0
    Sval = 0.0
    for kk=0:100
       Sval = Sval + (ec - 1.0/(2.0*kk+3.0))*(x^kk)*(2.0*kk+3.0)
    end
    return Sval
end


  
function ellipticnearparabolic(M,ec)
    nu = M;
    ier = 0
    if(M>pi) 
       MP = (M-2.0*pi)/ sqrt(2.0*(1.0-ec)^3.0)
    else
       MP = M / sqrt(2.0*(1.0-ec)^3.0)
    end 

    B = 3.0*MP/2.0
    A = (B+sqrt(1.0+B^2.0))^(2.0/3.0)
    D = 2.0*A*B/(1.0+A+A^2.0)

    eccost1 = sqrt(2.0)/sqrt(1.0+ec)
    eccost2 = sqrt(2.0)/sqrt((1+ec)^3.0)
    for jj = 1:20
       S =  Sfun(ec,D);
       fun = eccost1*D + eccost2*D^3.0*S - MP
       if(abs(fun)<=1e-10) 
          ier = 1;
          break;
       end 
       if(my_isnan(fun)) 
          ier = -1
          break
       end
       DS =  DSfun(ec,D)
       dfun = eccost1 + eccost2*D^2.0*DS
       if(dfun==0.0)
          ier = -1
          exit
       end
       D = D - fun/dfun
    end
    if(ier==1)
       f = atan(D)
       if(D<0.0) 
          nu = nu - pi
       end 
       nu = mod(2.0*nu,2.0*pi)
       cnu = cos(nu)
       snu = sin(nu)
       dean = 1.0+ec*cnu
       sean = sqrt(1.0-ec^2.0)*snu/dean
       cean = (ec+cnu)/dean      
       ean = mod(atan(sean,cean),2.0*pi)
    end
    return nu,ean,ier 
end 

function ellipticstrong(M,ec)
    ier = 0
    nu = M
    n_initialguess = 2
    ean = M
    if(M<pi)
        ean0 = [pi, M+ec]
    else
        ean0 = [pi, M-ec]
    end 
    
    for n_initialguess = 1:2
        ean = ean0[n_initialguess]
        for jj = 1:20
            fun = M-ean+ec*sin(ean)
            if (my_isnan(fun)) 
                ier = -1
                exit             
            end 
            if (abs(fun)<1e-10) 
                ier = 1 
                break;
            end 
            dfun = -1+ec*cos(ean)
            if(dfun==0) 
                ier = -1
                break;
            end 
            ean = ean -fun/dfun
        end 
        if(ier==1) 
            break;
        elseif n_initialguess==1
            ier = 0
        end
    end 
    if(ier==1)  
       cean = cos(ean)
       sean = sin(ean)
       dnu = 1.0-ec*cos(ean)
       cnu = (cean-ec)/dnu
       snu = sqrt(1.0-ec^2)*sean/dnu
       nu = mod(atan(snu,cnu),2.0*pi)
    end 
    return nu,ean,ier
end 

"""
    mean2truelongitude(P1,P2,Q1,Q2,M)
    from mean anomaly to true longitude

    INPUT:
        P1 = ec*sin(OM+om), ec = eccentricity, OM = longitude of ascending node, om = argument of pericentre                   
        P2 = ec*cos(OM+om)                         
        Q1 = tan(i/2)*sin(OM)  i = inclination
        Q2 = tan(i/2)*cos(OM)  
        M  = mean anomaly [rad]
    OUTPUT
        TL = true longitude (nu+OM+om, nu = true anomaly) [rad]
        ier =  1 --> conversion succeded
              -1 --> conversion failed

"""
function mean2truelongitude(P1,P2,Q1,Q2,M)
    ec = sqrt(P1^2.0+P2^2.0);

    if ec!=0.0 
        omPOM  = mod(atan(P1,P2),2*pi);
    else
        incl = mod(atan(sqrt(Q1^2.0+Q2^2.0))*2.0,2.0*pi);
        if incl!=0.0
            omPOM  = mod(atan(Q1,Q2),2*pi);
        else
            omPOM  = 0.0;
        end
    end

    nu,ean,ier = mean2trueandeccentric(M,ec);
    if ier == 1
        TL = mod(nu+omPOM,2.0*pi);
    else
        TL = 0.0;
    end

    return TL,ier
end

"""
    meanlong2ecclong(ML,P1,P2)
    returns the value of the eccentric longitude given the
    mean longitude


    INPUT 
    ML = mean longitude (M+om+OM, M = mean anomaly, om = argument of pericenter, OM = right ascension of the ascending node) [rad]
    P1 = ec*sin(om+OM), ec=eccentricity
    P2 = ec*cos(om+OM)

    OUTPUT
    EL = eccentric longitude (E+om+OM, E = eccentric anomaly)
"""
function meanlong2ecclong(ML,P1,P2)
  # uses newtons method to solve keplers equation for equinoctial elements
    error = 1e-8
    EL = ML 

    ratio = 1

    while max(abs(ratio)) > error
        ratio = (EL + P1 .* cos(EL) - P2 .* sin(EL) - ML)./(1.0 - P1 .* sin(EL) - P2 .* cos(EL));
        EL = EL - ratio
    end

    return EL;
end

"""
    meanlong2truelong(ML,P1,P2)
    returns the value of the true longitude given the
    mean longitude


    INPUT 
    ML = mean longitude (M+om+OM, M = mean anomaly, om = argument of pericenter, OM = right ascension of the ascending node) [rad]
    P1 = ec*sin(om+OM), ec=eccentricity
    P2 = ec*cos(om+OM)

    OUTPUT
    TL = true longitude (nu+om+OM, nu = true anomaly)
"""
function meanlong2truelong(ML,P1,P2)

    EL = meanlong2ecclong(ML,P1,P2) # SORT THIS
    TL = ecclong2truelong(EL,P1,P2);  
    return TL;

end


"""
    truelong2meanlong(TL,P1,P2)
    returns the value of the true longitude given the
    mean longitude


    INPUT 
    TL = true longitude (nu+om+OM, nu = true anomaly, om = argument of pericenter, OM = right ascension of the ascending node) [rad]
    P1 = ec*sin(om+OM), ec=eccentricity
    P2 = ec*cos(om+OM)

    OUTPUT
    ML = mean longitude (M+om+OM, M = mean anomaly) [rad]
    
"""
function truelong2meanlong(TL,P1,P2)
    EL =  truelong2ecclong(TL,P1,P2)
    ML = EL + P1*cos(EL) - P2*sin(EL);
    ML = mod(ML,2.0*pi);
    return ML;
end

"""
    eccentriclong2meanlong(EL,P1,P2)
    returns the value of the true longitude given the
    mean longitude


    INPUT 
    EL = eccentric longitude (E+om+OM, E = eccentric anomaly, om = argument of pericenter, OM = right ascension of the ascending node) [rad]
    P1 = ec*sin(om+OM), ec=eccentricity
    P2 = ec*cos(om+OM)

    OUTPUT
    ML = mean longitude (M+om+OM, M = mean anomaly) [rad]
    
"""
function ecclong2meanlong(EL,P1,P2)
    ML = EL + P1*cos(EL) - P2*sin(EL);
    ML = mod(ML,2.0*pi);
    return ML;
end



####### WRONG --- USE SatelliteToolbox instead
"""
    equi2hlonglat_earth(equiEQUATORIAL,idxLL,jd)

    input
    equiEQUATORIAL  = equinoctial elements expressed with respect to the Earth centred inertial frame with x axis 
            towards the vernal equinox
            a   semi-major axis              
            P1 = ec*sin(OM+om)                   
            P2 = ec*cos(OM+om)                         
            Q1 = tan(i/2)*sin(OM)  
            Q2 = tan(i/2)*cos(OM)  
            LL = nu+om+OM [rad]  if idxLL = 1
                 Ean+om+OM [rad] if idxLL = 2
                 M+om+OM [rad]   if idxLL = 3
            with 
                e   orbital eccentricity                   
                i   inclination                            [rad]
                o   argumentum of pericentre               [rad]
                O   right ascension of the ascending node  [rad]
                nu  true anomaly                           [rad]
                M   mean naomaly                           [rad]
                Ean eccentric anomaly                      [rad]
    idxLL       = identifier of the input longitude in the set of equinoctial elements
    mjd200      =  modified julian day 2000
    lattype     = if 1 the output gives the geodetic latitude, otherwise the geocentric latitude

    output 
    hh          = altitude 
    long        = longitude
    lat         = latitude (geocentric or geodetic)
    
    ref 
    Vallado,1997

    Irene Cavallari

"""
function EECIposition2hlonglat(rV,mjd2000,lattype)

    earthparams = astroconstants(3);

    # params
    rEarth = earthparams[2];
    oblateness = earthparams[4];

    # norm
    rN  = norm(rV)
    
    # geocentric latitude
    clat = rV[3]/rN;
    latcentric = acos(clat);

    if latcentric <=pi/2    
        latcentric = pi/2-latcentric
    else
        latcentric = pi-latcentric;
        latcentric = -pi/2+latcentric;
    end
    if lattype == 1 # geodetic
        lat = atan(tan(latcentric)/(1-oblateness)^2.0) 
    else
        lat = atan(tan(latcentric));
    end
    
    # longitude
    rightascension = atan(rV[2]/rN,rV[1]/rN);
    GST = mjd20002greenwichsiderealtime(mjd2000,1);
    long = mod(rightascension - GST,2*pi);
    if long>pi
        long = long - 2*pi;
    end

    # altitude
    hh   = rN-rEarth;

    return hh,lat,long;
end