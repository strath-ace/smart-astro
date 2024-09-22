
"""
    eclipseinout(rsun,sma,P1,P2,q11,q12,q13,q21,q22,q23,rPlanet)
    Returns the true longitude at the entrance and exit from the shadow region
    (cylindrical model)

    INPUTS
        rsun = position vector of the Sun wrt to the central planet [km]
        sma  = sei-major axis [km]
        P1 = ec*sin(OM+om), ec = eccentricity, OM = longitude of ascending node, om = argument of pericentre                   
        P2 = ec*cos(OM+om)       
        q11 = (1.0-Q1^2.0+Q2^2.0)/GG, GG = 1.0+Q1^2.0+Q2^2.0, Q1 = tan(i/2)*sin(OM), Q2 = tan(i/2)*cos(OM),  i = inclination
        q12 = (2.0*Q1*Q2)/GG;
        q13 = (-2.0*Q1)/GG;
        q21 = (2.0*Q1*Q2)/GG;
        q22 = (1.0+Q1^2.0-Q2^2.0)/GG;
        q23 = (2.0*Q2)/GG;
        rPlanet = mean radius of the central planet

    OUTPUTS
        inshadow = true if the orbit is shadowed; false if the orbit is always in sunlight
        TLin     = true longitude at entrance of the shadow region (if inshadow=false, TLin=NaN ) [rad]
        TLout    = true longitude at the exit from the shadow region (if inshadow=false, TLout=NaN ) [rad]
        true longitude = OM + om + nu, nu = true anomaly

    AUTHOR:
            Irene Cavallari
         
"""
function eclipseinout(rsun,sma,P1,P2,q11,q12,q13,q21,q22,q23,rPlanet)

    # init output
    inshadow = false;
    TLin  = NaN;
    TLout = NaN; 

    # direction cosines of the sun (inertial ref frame)
    rsunnorm = norm(rsun);
    cs1 = rsun[1]/rsunnorm;
    cs2 = rsun[2]/rsunnorm;
    cs3 = rsun[3]/rsunnorm;
    
    # (rsun dot r)/norm(rsun)/norm(r) = beta*cos(TL)+xi*sin(TL)
    # computation of beta and xi
    beta = cs1*q11+cs2*q12+cs3*q13;
    xi   = cs1*q21+cs2*q22+cs3*q23;

    # semilatus rectum
    p = sma*(1-P1^2.0-P2^2.0)

    # a0 y^4 + a1 y^3 + a2 y^2 + a3 y + a4 = 0, y = sin(TL)
    a0 = (P1^2.0*rPlanet^2.0 - 2.0*P1*beta*p*rPlanet + P2^2.0*rPlanet^2.0 + 2.0*P2*p*rPlanet*xi + beta^2.0*p^2.0 + p^2.0*xi^2.0)*(P1^2.0*rPlanet^2.0 + 2.0*P1*beta*p*rPlanet + P2^2.0*rPlanet^2.0 - 2.0*P2*p*rPlanet*xi + beta^2.0*p^2.0 + p^2.0*xi^2.0);
    a1 = 4.0*rPlanet^2.0*(P1^3.0*rPlanet^2.0 + P1*P2^2.0*rPlanet^2.0 - P1*beta^2.0*p^2.0 + P1*p^2.0*xi^2.0 + 2.0*P2*beta*p^2.0*xi);
    a2 = -2.0*P1^2.0*P2^2.0*rPlanet^4.0 + 2.0*P1^2.0*beta^2.0*p^2.0*rPlanet^2.0 - 8.0*P1*P2*beta*p^2.0*rPlanet^2.0*xi - 2.0*P2^4.0*rPlanet^4.0 - 4.0*P2^2.0*beta^2.0*p^2.0*rPlanet^2.0 + 2.0*P2^2.0*p^2.0*rPlanet^2.0*xi^2.0 - 2.0*beta^4.0*p^4.0 - 2.0*beta^2.0*p^4.0*xi^2.0 - 2.0*P1^2.0*p^2.0*rPlanet^2.0 + 6.0*P1^2.0*rPlanet^4.0 + 2.0*P2^2.0*p^2.0*rPlanet^2.0 + 2.0*P2^2.0*rPlanet^4.0 + 2.0*beta^2.0*p^4.0 - 2.0*beta^2.0*p^2.0*rPlanet^2.0 - 2.0*p^4.0*xi^2.0 + 2.0*p^2.0*rPlanet^2.0*xi^2.0
    a3 = -4.0*rPlanet^2.0*(P1*P2^2.0*rPlanet^2.0 - P1*beta^2.0*p^2.0 + 2.0*P2*beta*p^2.0*xi + P1*p^2.0 - P1*rPlanet^2.0);
    a4 = (P2^2.0*rPlanet^2.0 + beta^2.0*p^2.0 - 2.0*P2*rPlanet^2.0 - p^2.0 + rPlanet^2.0)*(P2^2.0*rPlanet^2.0 + beta^2.0*p^2.0 + 2.0*P2*rPlanet^2.0 - p^2.0 + rPlanet^2.0)
    
    a0 = a0/rPlanet^2.0;
    a1 = a1/rPlanet^2.0;
    a2 = a2/rPlanet^2.0;
    a3 = a3/rPlanet^2.0;
    a4 = a4/rPlanet^2.0;
    

    nreal,ysol = eq4deg([a0,a1,a2,a3,a4]);
    if nreal>0
        ysolgood = ysol[abs.(ysol).<1.0];
        solposs  = mod.(append!(asin.(ysolgood),pi*ones(length(ysolgood),1)-asin.(ysolgood)),2.0*pi);
        cxiV = beta*cos.(solposs)+xi*sin.(solposs);
        solposs = solposs[cxiV.<0.0];
        cxiV    = cxiV[cxiV.<0.0];
        funV =  p^2.0./((1.0 .+ P1*sin.(solposs) + P2*cos.(solposs)).^2.0).*(1.0 .-cxiV.^2.0) .- rPlanet^2.0;
        sol = solposs[abs.(funV).<1e-4];

        if length(sol)>2
            Base.error("method gives an impossible");
        elseif length(sol)==2
            inshadow = true;
            TL1 = minimum(sol);
            TL2 = maximum(sol);
            TLcheck = (TL1+TL2)/2.0;
            funcheck = p^2.0/((1.0 + P1*sin(TLcheck) + P2*cos(TLcheck))^2.0)*(1.0 - (beta*cos(TLcheck)+xi*sin(TLcheck))^2.0) - rPlanet^2;
                 
           if funcheck<0.0
            TLin = TL1;
            TLout = TL2;
           else
            TLin = TL2;
            TLout = TL1; 
           end
        end
    end
    return inshadow,TLin,TLout
end


"""
    eclipseinoutV2(rsun,sma,P1,P2,Q1,Q2,rPlanet)
    Returns the true longitude at the entrance and exit from the shadow region
    (cylindrical model)

    INPUTS
        rsun = position vector of the Sun wrt to the central planet [km]
        sma  = sei-major axis [km]
        P1 = ec*sin(OM+om), ec = eccentricity, OM = longitude of ascending node, om = argument of pericentre                   
        P2 = ec*cos(OM+om)       
        Q1 = tan(i/2)*sin(OM),  i = inclination
        Q2 = tan(i/2)*cos(OM)
        rPlanet = mean radius of the central planet

    OUTPUTS
        inshadow = true if the orbit is shadowed; false if the orbit is always in sunlight
        TLin     = true longitude at entrance of the shadow region (if inshadow=false, TLin=NaN ) [rad]
        TLout    = true longitude at the exit from the shadow region (if inshadow=false, TLout=NaN ) [rad]
        true longitude = OM + om + nu, nu = true anomaly

    AUTHOR:
            Irene Cavallari
         
"""
function eclipseinoutV2(rsun,sma,P1,P2,Q1,Q2,rPlanet)

    # computation of beta and xi
    GG = 1.0+Q1^2+Q2^2;
    q11 = (1.0-Q1^2+Q2^2)/GG;
    q12 = (2.0*Q1*Q2)/GG;
    q13 = (-2.0*Q1)/GG;
    q21 = (2.0*Q1*Q2)/GG;
    q22 = (1.0+Q1^2-Q2^2)/GG;
    q23 = (2.0*Q2)/GG;

    inshadow,TLin,TLout = eclipseinout(rsun,sma,P1,P2,q11,q12,q13,q21,q22,q23,rPlanet)

    
    return inshadow,TLin,TLout
end


