"""
    astroConstants(in)
    Returns a vector of constants relative to the celestial body
    associated to the id in input

    List of identifiers:
 
    in = 0 : Generic astronomical constants:
            out = [1: Universal gravity constant (G) (from DITAN and Horizon) [km^3/(kg*s^2)]; 
                   2: Astronomical Unit (AU) (from DE405) [km]; 
                   3: Solar radiation pressure at 1 au [kg/(s^2*m)]

    in = 1 : Mercury 
            out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
                   2: Mean radius of the planets [km]]

    in = 2 : Venus 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km]]

    in = 3 : Earth 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km];
            3: Earth-Sun mean distance;
            4: Earth magnetic dipole strenght [kg*km^3/(s^2*A)];
            5: Earth rotation mean motion [rad/s];
            6-9: zonal harmonics J2-J5 ]

    in = 4 : Mars 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km]]

    in = 5 : Jupiter 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km]]

    in = 6 : Saturn 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km]]

    in = 7 : Uranus 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km]]

    in = 8 : Neptune 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km]]

    in = 9 : Pluto 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km]]

    in = 10 : Sun 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km]]

    in = 11 : Moon 
    out = [1: Planetary constants of the planets (mu = mass * G) [km^3/s^2];
            2: Mean radius of the planets [km];
            3: Earth-Moon mean distance]
                
          
    INSPIRED by MATLAB VERSION astroconstants (Version 22/03/2013)
    % AUTHOR:
    %   Matteo Ceriotti, 2006, MATLAB, astroConstants.m
    % LAST REVISION:
    %   22/03/2013, Francesca Letizia

    JULIA VERSION
    16/10/2023 : Irene Cavallari
"""
function astroconstants(nb)
    out = [];

    ############################## General
    if nb == 0   
        # astronomical unit [km]
        au = 149.598*1e+6; 

        # gravitationalconstant [km^3/(kg*s^2)]
        G = 6.67259e-20;

        # solar radiation pressure at 1 au [kg/(s^2*m)]
        psrpAU = 4.56*1e-6;

        out = [G,au,psrpAU]
    
    ############################### SUN
    elseif nb == 10
        # sun gravitational parameter [km^3/(s^2)]
        muSun = 0.19891000000000E+31 * 6.67259e-20;

        # sun mean radius [km]
        rSun  = 700000;
        
        out = [muSun, rSun];  

    ############################## EARTH
    elseif nb == 3
        # Earth gravitational parameter [km^3/(s^2)]
        muEarth   = 398600.4418; 
        
        # Earth magnetic dipole strenght [kg*km^3/(s^2*A)]
        muM  = 1e+8;

        # Earth's mean radius [km]
        rEarth =  0.637813700000000E+04; 
       
        # Earth-Sun mean distance
        dES = 149.598*1e+6; 

        # Earth rotation mean motion [rad/s]
        nE = 7.2921158553*1e-5; # Vallado,1997
      
        # Earth Zonal Harmonics 
        J2E = 0.00108263;
        J3E = -2.5327e-6;
        J4E = -1.6196e-6;
        J5E = -2.2730e-7;

        # oblateness
        oblateness = 0.003352813178;

        out = [muEarth,rEarth,dES,oblateness,muM,nE,J2E,J3E,J4E,J5E]; 

    ########################### MOON
    elseif nb == 11
        # Moon gravitational parameter [km^3/(s^2)]
        muMoon = 0.73476418263373E+23 * 6.67259e-20;

        # Moon mean radius [km]
        rMoon = 0.17380000000000E+04;

        # Earth-Moon mean distance [km]
        dEM = 384401;

        out = [muMoon,rMoon,dEM];
    
    ############################ Solar system planets
    # MERCURY gravitational parameter [km^3/(s^2)] and mean radius [km]
    elseif nb == 1 
        muMercury = 0.33020000000000E+24 * 6.67259e-20; 
        rMercury = 0.24400000000000E+04;

        out = [muMercury,rMercury]; 

    # VENUS gravitational parameter [km^3/(s^2)] and mean radius [km]
    elseif nb == 2
        muVenus = 0.48685000000000E+25 * 6.67259e-20;
        rVenus = 0.60518000000000E+04;

        out = [muVenus,rVenus];
        
    # MARS gravitational parameter [km^3/(s^2)] and mean radius [km]
    elseif nb == 4    
        muMars = 0.64184999247389E+24 * 6.67259e-20;
        rMars = 0.33899200000000E+0;

        out = [muMars,rMars]; 

    # JUPITER gravitational parameter [km^3/(s^2)] and mean radius [km]
    elseif  nb == 5
        muJupiter = 0.18986000000000E+28 * 6.67259e-20; 
        rJupiter = 0.69911000000000E+05;

        out = [muJupiter,rJupiter]; 

    # SATURN gravitational parameter [km^3/(s^2)] and mean radius [km]
    elseif nb == 6
        muSaturn = 0.56846000000000E+27 * 6.67259e-20;
        rSaturn = 0.58232000000000E+05;

        out = [muSaturn,rSaturn]; 

    # URANUS gravitational parameter [km^3/(s^2)] and mean radius [km]
    elseif nb == 7
        muUranus = 0.86832000000000E+26 * 6.67259e-20;
        rUranus = 0.25362000000000E+05;

        out = [muUranus,rUranus];

    # NEPTUNE gravitational parameter [km^3/(s^2)] and mean radius [km]
    elseif nb == 8
        muNeptune = 0.10243000000000E+27 * 6.67259e-20;
        rNeptune = 0.24624000000000E+05;

        out = [muNeptune,rNeptune];

    # PLUTO gravitational parameter [km^3/(s^2)] and mean radius [km]
    elseif  nb == 9
        muPluto = 0.14120000000000E+23 * 6.67259e-20;    
        rPluto = 0.11510000000000E+04;

        out = [muPluto,rPluto];


    end

    return out;

end