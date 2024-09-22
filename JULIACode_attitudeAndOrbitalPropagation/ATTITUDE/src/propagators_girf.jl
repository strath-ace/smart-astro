"""
    attitudeprop_quat_girf(civ,date0,params,tstep,tfinal)

    =========================================================================================================================================================================================================
    INPUT ============================================================================================================================================================
    cvi    = [euangles,omV,carVsun, kepEl], vector of Floats **
    date0  = calendar date at initial time, [year(Float64), month(Float64), day(Float64), hours(Int64), minutes(Int64), seconds(Float64)]
    params = [planetsunparams***;inertialrefframeinfo****;satellite*****,atmosphere******,settings*******] vector of structures of parameters 
    tstep  = timestep [sec]
    tfinal = time of propagation [sec] 
   
    __________________________________________
    ** 
    cvi: 
        euangles = Euler angles [rad], 
        omV      = components of angular velocity [rad/s], 
        carVsun  = velocity vector and position vector of the Sun wrt to Planet [km/s,km]
        kepEl    = [semi-major axis [km], eccentricity, inclination [rad], argument of pericenter [rad], longitude of ascending node [rad], true anomaly [rad]]

    ____________________________________________
    ***
    planetsunparams 
    Dictionary with keys
        	centralBodyIDX : integer number identifying the celestial body
            1:   Mercury
            2:   Venus
            3:   Earth
            4:   Mars
            5:   Jupiter
            6:   Saturn
            7:   Uranus
            8:   Neptune
            9:   Pluto
            10:  Sun
            11:  Moon
        	muPlanet = central Planet's gravitational parameter [km^3/s^2]
        	rPlanet = Planet radius [km]
        	muM = Planet's magnetic dipole strenght [kg km^3/s^2/A]
        	psrp = solar radiation pressure [kg/s^2/m] at a reference distance psrp_refdist from the sun
        	psrp_refdist = reference distance from the sun for the given value of psrp [km]   
        	planetRotation = Planet's rotational angular speed [rad/s]
        	zonalharmonicscoff = [Ji], i>=2 increasing, Ji zonal harmonics
        	perturbingbodiesid = vector containing the ids of the perturbing celestial bodies (third body perturbation)
            1:   Mercury
            2:   Venus
            3:   Earth
            4:   Mars
            5:   Jupiter
            6:   Saturn
            7:   Uranus
            8:   Neptune
            9:   Pluto
            11:  Moon
        	planet_flat = control integer: it is equal to 1 to consider planet's oblateness, equal to 0 otherwise.
        	oblateness = oblateness coefficient 
    
    NB: if the central body is the Earth, it is sufficient to give as an input just a Dictionary with the key centralBodyIDX and the associated value 3 (planetsunparams=Dict(“centralBodyIDX”=>3)). The function will automatically fill the dictionary with the other required fields. In this case, the default value of “planet_flat” is 0. 
    
   
    ____________________________________________
    ****      
    inertialrefframeinfo
    Dictionary with keys
        	ecliptic2inertial =  rotation matrix from ecliptic to girf
        	equatorial2inertial = rotation matrix from equatorial to inertial ref frame
    NB: if the central body is the Earth (planetsunparams=Dict(“centralBodyIDX”=>3)), the user can give as an input an empty dictionary. In this case, the program assumes that the inertial reference frame is the equatorial reference frame and automatically generated the dictionary.



    ____________________________________________
    *****     
    satellite
    Dictionary with keys
        	Moments of Inertia = [A,B,C], principal moment of inertia (A<=B<=C) [Kg m^2]
        	intrinsicMagneticMoment = intrinsic magnetic moment of the satellite [A m^2]
        	mass = mass of the satellite [kg]
        	CD = non-dimensional drag coefficient of the satellite
        	numberOfFacets = number of facets in which the satellite external surface is divided
        	facets = array of arrays, each associated to a facet of the satellite. In particular,
            facets = [facet_1..facet_k..facet_nf]
            with 
            facet_k = [facet_coeff,facet_area,facet_Vinfo_facet_nv]
            where
            facet_coeff = [caj,cdj,csj]
            facet_area  = area of the facet
            facet_Vinfo = [rhojv,vv1,vv2,vv3] or [rhojv,vv1,vv2,vv3,vv4]
                    rhojv = vector from centre of mass to centroid
                    vvi   = vertexes of the facet (only triangular and squared facets accepted)

    NB: if the input dictionary does not contain the key CD this is automatically introduced with the default value 2. 



    ____________________________________________
    ******     
    atmosphere 
    Dictionary with keys
        	typeofmodel = integer number identifying the atmospheric model and consequently the drag atmospheric perturbation model to be exploited
            1 : exponential atmospheric model
            2 : nrlmsise00 
        	AtmM = matrix containing information required by the atmospheric model
    
                If typeodmodel = 1
                AtmM = matrix containing the exponential atmospheric model 
                in each row: [min altitude [km], density [kg/m^3], scale height [km]
                
                if typeofmodel = 2
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
        
    NB: if the dictionary does not contain the key “typeofmodel” this is automatically added with default value 1 (exponential atmospheric model)
    
    NB: if the central body is the Earth (planetsunparams=Dict(“centralBodyIDX”=>3)), and the user wants to use the exponential atmospheric model, it is sufficient to give as an input a dictionary with the only key “typeofmodel” equal to 1. The key AtmM will be automatically generated, using the exponential model by Vallado (1997), (see Section 7). 
    

    ____________________________________________
    *******     
    settings 
    Dictionary with keys 
        	includeGravityTorque = Boolean, true if the perturbation by the gravity gradient torque is to be considered
        	includeMagneticTorque = Boolean, true if the perturbation by the residual magnetic torque is to be considered
        	includeSrpTorque = Boolean, true is the perturbation by the light pressure torque is to be considered
        	includeDragTorque = Boolean, true if the perturbation by the atmospheric drag torque is to be considered
        	includeEclipsesEffectsOnAttitude = Boolean, true if the light pressure torque perturbation must be associated with the eclipses effects
        	includeZonalHarmsAcc = Boolean, true is the zonal harmonics perturbation must be considered
        	maxzonalharmonics = natural number equal to the degree of the maximum zonal harmonics to be considered
        	includeThirdBodyAcc = Boolean, true if the third body perturbation must be considered
        	includeSunGravityAcc = Boolean, true if the perturbation by the Sun gravity must be considered
        	includeSrpAcc = Boolean, true if the perturbation by the light pressure acceleration must be considered
        	includeDragAcc = Boolean, true if the perturbation  by the atmospheric drag acceleration must be considered
        	includeEclipsesEffectsOnOrbit= Boolean, true if the light pressure acceleration perturbation must be associated with the eclipses effects

    

    ___________________________________________
    BE CAREFUL TO UNITS!!!!!!!!!!!!!!!!!!

    =========================================================================================================================================================================================================
    OUTPUTS =================================================================================================================================================================================================
    tV
    res
    with
    tV  = time vector, size (nx1)
    res = matrix (nx40); in each row: (sadov variables, components of angular velocity, euler angles, andoyer variables, equinoctial elements, skm, energy, Jd, savarlike/anvarlike, keplerian elements)


    _______________________________________________
    INFO : propagation performed using quaternions

"""

function attitudeprop_quat_girf(civ,date0,params,tstep,tfinal;print_flag=1)
    
    params =  checks_prop_girf(civ,date0,params,tstep,tfinal);
 
    # inputs
    planetsunparams = deepcopy(params[1]);
    inertialrefframeinfo       = deepcopy(params[2]);
    satellite       = deepcopy(params[3]);
    atmosphere      = deepcopy(params[4]);
    settings        = deepcopy(params[5]);    

    euangles   = civ[1:3];
    omV        = civ[4:6];
    kepEl      = civ[7:12];

    # initial conditions in quaternions
    quat = euler2quat(euangles);

    # orbit 
    equiEl      = kep2equi(kepEl,3); 
        
    # initialisation for srp
    if settings["includeSrpAcc"] || settings["includeSrpTorque"]
        for jj=1:satellite["numberOfFacets"]
            satellite["facets"][jj][1] = satellite["facets"][jj][1]*planetsunparams["psrp"];
        end
    end

    # initialisation for zonal harmonics perturbations
    if settings["includeZonalHarmsAcc"]
        if settings["maxzonalharmonics"][1]<5
            for jj = settings["maxzonalharmonics"][1]:4
                planetsunparams["zonalharmonicscoff"][jj] = 0.0;
            end
        end
    end

    # initialisation for third body perturbation
    if settings["includeThirdBodyAcc"]
        planetsidxslist = planetsunparams["perturbingbodiesid"];
        planetsgplist   = zeros(size(planetsidxslist)[1])
        for jj = 1:size(planetsidxslist)[1]
            planetsgplist[jj] = astroconstants(planetsidxslist[jj])[1]
        end
        planetsunparams = merge(planetsunparams,Dict("perturbingbodiesgravparam"=>planetsgplist));
    end

    if (settings["includeDragTorque"] || settings["includeDragAcc"]) && atmosphere["atmosphericmodel"]==1
        println("An exponential atmospheric model is selected --- the aerodynamic force is in the direction of the incident stream (any potential lift component is neglected) -- the adimensional drag coefficient must be set by the user (satellite[CD]) otherwise it is assumed equal to 2");
    end

    # propagation
    paramsint = [planetsunparams;inertialrefframeinfo;satellite;atmosphere;settings;date2mjd2000(date0);zeros(satellite["numberOfFacets"])];
    ts    = collect(0.0:tstep:tfinal);
    if (settings["includeDragTorque"] || settings["includeDragAcc"])
        u0 = append!(equiEl,quat,omV);
        tspanprop = [0.0 sign(tfinal)*3600];
        tpropagated = 0.0
        resQ = zeros(size(ts)[1],13)
        while abs(tpropagated)<abs(tfinal)
            if (print_flag==1 || print_flag==true)
				println("time propagated [day] : ",tpropagated/24.0/3600.0)
			end
            # cb = VectorContinuousCallback(conditiondrag, affectdrag!,12)
            prob = ODEProblem(orbitAttitudeEv_quat,u0,tspanprop,paramsint);
            sol  = solve(prob,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+8); #,callback=cb);
            tt = sol.t .+ tpropagated;
            idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tt[end]);
            if maximum(idx)
                solsol  = sol(ts[idx].-tpropagated);
                resjj = mapreduce(permutedims,vcat,solsol.u)
                resQ[idx,:] = resjj;
            end
            # update
            paramsint[6] = paramsint[6] +  sol.t[end]/24/3600;
            resG    =  mapreduce(permutedims,vcat,sol.u);
            u0      = resG[end,:];
            tpropagated = tt[end];
        end
        tV = ts;
    else
        u0        = append!(equiEl,quat,omV);
        tspan = [0.0 tfinal];
        prob = ODEProblem(orbitAttitudeEv_quat,u0,tspan,paramsint);
        sol  = solve(prob,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);
        sol  = sol(ts);
        resQ = mapreduce(permutedims,vcat,sol.u);
        tV   = sol.t;    
    end

    res = zeros(size(resQ,1),40);
    IV = satellite["MomentsOfInertia"];
    for jj=1:size(resQ)[1]
        omVjj = resQ[jj,11:13];
        euanglesjj = quat2euler(resQ[jj,7:10]);
        
        euanglesmodjj,omVmodjj,isrefframechanged,RM = selectreferenceframe_prop_girf(IV,euanglesjj[1],euanglesjj[2],euanglesjj[3],omVjj[1],omVjj[2],omVjj[3]);
        
        IVmod = RM*IV;
        IVmod[1] = abs(IVmod[1]);
        IVmod[2] = abs(IVmod[2]);
        IVmod[3] = abs(IVmod[3]); 
        savarjj,skmjj = euler2sadov(euanglesmodjj,omVmodjj,IVmod);
        
        anvar     = euler2andoyer(euanglesjj,omVjj,IV);
        if IVmod[3]*(IVmod[2]-IVmod[1])-IVmod[1]*(IVmod[3]-IVmod[2])==0.0
            kjj = 0.0
        else
            kjj = IVmod[3]/IVmod[1]*(IVmod[2]-IVmod[1])/(IVmod[3]-IVmod[2]);
        end
        # slikeC    = sadov2sadovlikeC(savarjj,skmjj,kjj);

        energy = (IV[1]*omVjj[1]^2+IV[2]*omVjj[2]^2+IV[3]*omVjj[3]^2)/2;
        Jd     = anvar[2]^2/2/energy;

        if abs(IV[2]-IV[1])==0.0 && abs(IV[3]-IV[2])==0.0
            varlike  = andoyer2andoyerlike(anvar);
        else
            varlike = sadov2sadovlike(savarjj,skmjj)
        end

        kepjj  = equi2kep(resQ[jj,1:6],3);         

        res[jj,:] = append!(savarjj,omVjj,euanglesjj,anvar,resQ[jj,1:6],skmjj,energy,Jd,varlike,kepjj); 
    end

    return tV,res;
end

#### deprecated
function attitudeprop_sadov_girf(civ,date0,params,tstep,tfinal)

    params =  checks_prop_girf(civ,date0,params,tstep,tfinal);
 
    # inputs
    planetsunparams = deepcopy(params[1]);
    inertialrefframeinfo  = deepcopy(params[2]);
    satellite       = deepcopy(params[3]);
    atmosphere      = deepcopy(params[4]);
    settings        = deepcopy(params[5]);    

    euangles   = civ[1:3];
    omV        = civ[4:6];
    kepEl      = civ[7:12];

    IV = satellite["MomentsOfInertia"];  
    if IV[1]==IV[2] && IV[2]==IV[3]
        Base.error("not suitable method")
    end

    # choice of body reference frame
    euangles,omV,isrefframechanged,RM = selectreferenceframe_prop_girf(IV,euangles[1],euangles[2],euangles[3],omV[1],omV[2],omV[3])
    if isrefframechanged
        IV = RM*IV;
        IV[1] = abs(IV[1]);
        IV[2] = abs(IV[2]);
        IV[3] = abs(IV[3]); 
        satellite["MomentsOfInertia"] = IV;
        numberOfFacets=satellite["numberOfFacets"];
        for kk = 1:numberOfFacets    
            satellite["facets"][kk][3][1]=RM*satellite["facets"][kk][3][1];
            satellite["facets"][kk][3][2]=RM*satellite["facets"][kk][3][2];
            satellite["facets"][kk][3][3]=RM*satellite["facets"][kk][3][3];
            satellite["facets"][kk][3][4]=RM*satellite["facets"][kk][3][4];
            satellite["facets"][kk][4]=RM*satellite["facets"][kk][4];
        end
        satellite["intrinsicMagneticMoment"] = RM*satellite["intrinsicMagneticMoment"];
    end
    if IV[3]*(IV[2]-IV[1])-IV[1]*(IV[3]-IV[2])==0.0
        kconstant = 0.0
    else
        kconstant = IV[3]/IV[1]*(IV[2]-IV[1])/(IV[3]-IV[2]);
    end
    satellite = merge(satellite,Dict("k_constant"=>kconstant));
    
    # input in sadov variables
    savar,skm = euler2sadov(euangles,omV,IV);    


    if settings["includeSrpAcc"] || settings["includeSrpTorque"]
        for jj=1:satellite["numberOfFacets"]
            satellite["facets"][jj][1] = satellite["facets"][jj][1]*planetsunparams["psrp"];
        end
    end

    # initialisation for zonal harmonics perturbations
    if settings["includeZonalHarmsAcc"]
        if settings["maxzonalharmonics"][1]<5
            for jj = settings["maxzonalharmonics"][1]:4
                planetsunparams["zonalharmonicscoff"][jj] = 0.0;
            end
        end
    end

    # initialisation for third body perturbation
    if settings["includeThirdBodyAcc"]
        planetsidxslist = planetsunparams["perturbingbodiesid"];
        planetsgplist   = zeros(size(planetsidxlist)[1])
        for jj = 1:size(planetsidxlist)[1]
            planetsgplist[jj] = astroconstants(planetsidxslist[jj])[1]
        end
        planetsunparams = merge(planetsunparams,"perturbingbodiesgravparam"=>planetsgplist);
    end
   
    # orbit 
    equiEl      = kep2equi(kepEl,3); 
  
    # attitude propagation
    u0        = append!(equiEl,skm,savar[2:6]);
    paramsint = [planetsunparams;inertialrefframeinfo;satellite;atmosphere;settings;date2mjd2000(date0)];
    tspan = [0.0 tfinal];
    ts = collect(0.0:tstep:tfinal);
    prob = ODEProblem(orbitAttitudeEv_sadov,u0,tspan,paramsint);
    sol = solve(prob,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10);
    sol = sol(ts);
    resS = mapreduce(permutedims,vcat,sol.u);
    tV = sol.t;
 
    # handle sol
    for jj = 1:size(resS)[1]
        resS[jj,10] = mod(resS[jj,10],2*pi)
        resS[jj,11] = mod(resS[jj,11],2*pi)
        resS[jj,12] = mod(resS[jj,12],2*pi)
    end

    # transformation
    if isrefframechanged
        RM = transpose(RM)
    end
    
    res = zeros(size(resS,1),40);
    for jj = 1 : size(resS)[1]
        euanglesjj,omVjj,Jljj = sadov2eulerAndJl(resS[jj,8:12],resS[jj,7],IV);
        savarjj = append!([Jljj],resS[jj,8:12])
        IVmod = IV;
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RM*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RM*omVjj;
            IVmod      = RM*IV;
            IVmod[1]   = abs(IVmod[1]);
            IVmod[2]   = abs(IVmod[2]);
            IVmod[3]   = abs(IVmod[3]);
        end
        anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        energy = (IV[1]*omVjj[1]^2+IV[2]*omVjj[2]^2+IV[3]*omVjj[3]^2)/2;
        Jd     = anvar[2]^2/2/energy;
        kepjj  = equi2kep(resS[jj,1:6],3);    
        varlike = sadov2sadovlike(savarjj,resS[jj,7])
        res[jj,1:40] = append!(savarjj,omVjj,euanglesjj,anvar,resS[jj,1:6],resS[jj,7],energy,Jd,varlike,kepjj); 
    end
    
    return tV, res;
end

##########################################################################
### Functions
##########################################################################

function checks_prop_girf(civ,date0,params,tstep,tfinal)

    ################ check civ

    if size(civ)[1] != 12
        Base.error("error: the length of the input vector must be 18");
    end

    ############### check date0
    if size(date0)[1] != 6
        Base.error("date: [year,month,day,hr,min,sec]");
    end


    ############## check time
    if abs(tstep)>abs(tfinal)
        Base.error("tstep>tfinal");
    end

    
    ############### check params
    if size(params)[1] != 5
        Base.error("error: the length of params must be 5");
    end
    
    p1 = params[1];
    p2 = params[2];
    p3 = params[3];
    p4 = params[4];
    p5 = params[5];

    # check settings
    if haskey(p5,"includeGravityTorque") == false
        p5 = merge(p5,Dict("includeGravityTorque"=>false));
    end

    if haskey(p5,"includeMagneticTorque") == false
        p5 = merge(p5,Dict("includeMagneticTorque"=>false));
    end

    if haskey(p5,"includeSrpTorque") == false
        p5 = merge(p5,Dict("includeSrpTorque"=>false));
    end

    if haskey(p5,"includeDragTorque") == false
        p5 = merge(p5,Dict("includeDragTorque "=>false));
    end
    
    if haskey(p5,"includeZonalHarmsAcc") == false
        p5 = merge(p5,Dict("includeZonalHarmsAcc"=>false));
    end

    if haskey(p5,"maxzonalharmonics") == false &&  p5["includeZonalHarmsAcc"]
        p5 = merge(p5,Dict("maxzonalharmonics"=>[5]));
    else
        p5["maxzonalharmonics"] = [min(p5["maxzonalharmonics"][1],5)];
        if p5["maxzonalharmonics"][1] == 1
            p5["maxzonalharmonics"] = [0.0]
            p5["includeZonalHarmsAcc"] = false;
        end
    end

    if haskey(p5,"includeThirdBodyAcc") == false
        p5 = merge(p5,Dict("includeThirdBodyAcc"=>false));
    end

    if haskey(p5,"includeSunGravityAcc") == false
        p5 = merge(p5,Dict("includeSunGravityAcc"=>false));
    end

    if haskey(p5,"includeSrpAcc") == false
        p5 = merge(p5,Dict("includeSrpAcc"=>false));
    end

    if haskey(p5,"includeDragAcc") == false
        p5 = merge(p5,Dict("includeDragAcc"=>false));
    end

    if haskey(p5,"includeEclipsesEffectsOnOrbit") == false
        if p5["includeSrpAcc"]
            p5 = merge(p5,Dict("includeEclipsesEffectsOnOrbit"=>true));
        else
            p5 = merge(p5,Dict("includeEclipsesEffectsOnOrbit"=>false));
        end
    end

    if haskey(p5,"includeEclipsesEffectsOnAttitude") == false
        if p5["includeSrpTorque"]
            p5 = merge(p5,Dict("includeEclipsesEffectsOnAttitude"=>true));
        else
            p5 = merge(p5,Dict("includeEclipsesEffectsOnAttitude"=>false));
        end
    end
   
    # check satellite characteristics
    if haskey(p3,"MomentsOfInertia")==false
        Base.error("required satellite's moments of inertia");
    end

    if haskey(p3,"mass")==false && (p5["includeDragAcc"] || p5["includeSrpAcc"])
        Base.error("required satellite's mass");
    end

    if haskey(p3,"intrinsicMagneticMoment") == false 
        if p5["includeMagneticTorque"]
            Base.error("required satellite's intrinsic magnetic moment");
        else
            p5 = merge(p5,Dict("intrinsicMagneticMoment"=>[0.0,0.0,0.0]))
        end
    end


    messageerror = "lack of essential info about satellite's geometry and characteristics => the satellite has to be divided in panels and the characteristics of each panel are required :\n 
                        numberOfFacets  = number of facets in which the satellite surface is divided \n facets = [facet_1...facet_n], \n  where  
                        facet_k = [facet_coeff,facet_area,facet_rhov,facet_nv], with \n 
                        facet_coeff = [(1-rf*sp)*psrp,(2/3*(1-sp)*rf+2/3*(1-rf))*psrp,2*rf*sp*psrp],  rf reflectivity and sp specular coefficient (see Benson&Sheers2021),\n 
                        facet_area  = area of the facet [m^2],\n
                        facet_rhov  = centre of mass to facet centroid vector [m],\n 
                        facet_nv = normal unit vector";
    if (haskey(p3,"numberOfFacets")==false || haskey(p3,"facets")==false) 
        if  (p5["includeDragTorque"] || p5["includeSrpTorque"] || p5["includeDragAcc"] || p5["includeSrpAcc"])
            
            Base.error(messageerror);
        else
            p3 = merge(p3,Dict("numberOfFacets"=>0));
        end
    else
        if p5["includeDragTorque"] || p5["includeSrpTorque"] || p5["includeDragAcc"] || p5["includeSrpAcc"]
            if p3["numberOfFacets"] == 0.0
                Base.error(messageerror);
            end
            for jj = 1:p3["numberOfFacets"]
                if size(p3["facets"][jj])[1]!= 4
                    Base.error(messageerror);
                elseif size(p3["facets"][jj][1])[1]!= 3 || (size(p3["facets"][jj][3])[1]!=4 && size(p3["facets"][jj][3])[1]!= 5) || size(p3["facets"][jj][4])[1]!= 3 || typeof(p3["facets"][jj][2])!=Float64 
                    Base.error(messageerror);
                end
            end
        end 
    end

    if haskey(p3,"CD")==false && (p5["includeDragTorque"] || p5["includeDragAcc"])
        p3 = merge(p3,Dict("CD"=>2))
    end


    # check central body parameters 

    if haskey(p1,"centralBodyIDX") == false
        Base.error("missing id of the central body");
    else
        if p1["centralBodyIDX"] != 3 
            if p1<0 || p1>11
                Base.error("unknown central body")
            end

            if haskey(p1,"muPlanet") == false 
                p1 =  merge(p5,Dict("muPlanet"=>astroconstants(p1["centralBodyIDX"])[1]));
            else
                p1["muPlanet"] = astroconstants(p1["centralBodyIDX"])[1]; 
            end

            if haskey(p1,"rPlanet") == false 
                p1 =  merge(p5,Dict("rPlanet"=>astroconstants(p1["centralBodyIDX"])[2]));
            else
                p1["rPlanet"] = astroconstants(p1["centralBodyIDX"])[2]; 
            end
            
            if  (haskey(p1,"muM") == false && p5["includeMagneticTorque"])  
                Base.error("missing info about central body: muM = planet’s magnetic dipole strength");
            end

            if  (haskey(p1,"psrp")==false && p5["includeSrpTorque"]) || (haskey(p1,"psrp_refdist")==false && p5["includeSrpTorque"]) 
                Base.error("missing one of these data about central body: psrp, i.e. solar radiation pressur at the distance psrp_refdist from the Sun; psrp_refdist ");
            end

            if  (haskey(p1,"rPlanet")==false && (p5["includeSrpTorque"] || p5["includeDragTorque"]))
                Base.error("missing info about central body: planet mean radius ");
            end

            if  (haskey(p1,"planetRotation")==false && p5["includeDragTorque"])
                Base.error("missing info about central body: planetRotation: planet's mean rotational rate [s]");
            end

            if haskey(p1,"perturbingbodiesidx") == false && p5["includeThirdBodyAcc"]
                Base.error("missing info about central body: 3bodypert");
            end

            if haskey(p1,"zonalharmonicscoff") == false && p5["includeZonalHarmsAcc"]
                Base.error("required central body zonal harmonics");
            else
                if size(p1["zonalharmonicscoff"])[1]<p5["maxzonalharmonics"]-1 && p5["includeZonalHarmsAcc"]
                    Base.error("missing zonal harmonics coeff");
                end
            end
        else # Earth
            earthconstant = astroconstants(3);
            generalconstant = astroconstants(0);
            p1 = Dict("centralBodyIDX"=>3,"muPlanet"=>earthconstant[1],"rPlanet"=>earthconstant[2],"oblatenesscoeff"=>earthconstant[4],"muM"=>earthconstant[5],"planetRotation"=>earthconstant[6]); 
            p1 = merge(p1,Dict("psrp"=>generalconstant[3],"psrp_refdist"=>generalconstant[2])); 
            p1 = merge(p1,Dict("zonalharmonicscoff"=>earthconstant[7:10]));
            p1 = merge(p1,Dict("perturbingbodiesid"=>[11]));

        end
    end

    # check ref frame info
    if  haskey(p2,"ecliptic2inertial")==false && (p5["includeSrpAcc"] || p5["includeSrpTorque"] || p5["includeSunGravityAcc"] || p5["includeThirdBodyAcc"])
        if p1["centralBodyIDX"] == 3
            incl = 23.43928111*pi/180;
            RM = zeros(3,3);
            RM[1,1] = 1.0;
            RM[2,2] = cos(incl);
            RM[2,3] = -sin(incl);
            RM[3,2] = -RM[2,3];
            RM[3,3] = RM[2,2]
            p2 = merge(p2,Dict("ecliptic2inertial"=>RM)); 
        else
            Base.error("required rotation matrix between ecliptic ref frame 2 inertial ref frame")
        end
    end

    if  haskey(p2,"equatorial2inertial")==false && (p5["includeDragAcc"] || p5["includeDragTorque"] || p5["includeMagneticTorque"] || p5["includeThirdBodyAcc"])
        if p1["centralBodyIDX"] == 3
            p2 = merge(p2,Dict("equatorial2inertial"=>[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])); 
        else
            Base.error("required rotation matrix between ecliptic ref frame 2 inertial ref frame")
        end
    end

    if  haskey(p4,"atmosphericmodel")==false && (p5["includeDragTorque"] || p5["includeDragAcc"])
        p4 = merge(p4,Dict("atmosphericmodel"=>1));
        println("Model of atmosphere not specified -- AUTOMATICALLY IMPOSED EXPONENTIAL MODEL")
    else
        if p4["atmosphericmodel"]<1 || p4["atmosphericmodel"]>2
            Base.error(" only two exponential models supported: exponential (input = 1) or nrlmsise00 (input = 2)");
        end
    end

    if  (p5["includeDragTorque"] || p5["includeDragAcc"]) && p4["atmosphericmodel"]==2 && p1["centralBodyIDX"] != 3 
        Base.error("the nrlmsise00 atmospheric model works only if the central body is the Earth (centralBodyIDX = 3)");
    end

    if  haskey(p4,"AtmM") == false && (p5["includeDragTorque"] || p5["includeDragAcc"])
        if  p4["atmosphericmodel"] == 1
            if p1["centralBodyIDX"] == 3
                p4 = merge(p4,Dict("AtmM"=>getExponentialAtmosphericDensityModelForEarth()));
            else
                Base.error("required exponential model of atmospheric density");
            end
        elseif  p4["atmosphericmodel"] == 2
            Base.errore("required solar flux and geomagnetic data");
        end
    end

    return [p1,p2,p3,p4,p5];
end

function selectreferenceframe_prop_girf(IV,phi,theta,psi,p,q,r)

    RM = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
    isrefframechanged = false;

    if IV[1]==IV[2] && IV[2]==IV[3]
        # condition = p>r || q>r 
        condition = abs(p)>1e-5 || abs(q)>1e-5
    else
        condition = IV[1]*(IV[2]-IV[1])*p^2<IV[3]*(IV[3]-IV[2])*r^2;
    end
    
    if condition
        euangles = [phi;theta;psi];
        omV = [p;q;r];
    else
        Ri2b = euler2ibrotmat([phi;theta;psi]);
        RM = [0.0 0.0 1.0;0.0 -1.0 0.0;1.0 0.0 0.0];
        Ri2b = RM*Ri2b;
        euangles = rotmatib2euler(Ri2b);
        omV      = RM*[p;q;r];
        isrefframechanged = true;
    end

    RM2 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
    if omV[3]<0.0
        Ri2b = euler2ibrotmat(euangles);
        RM2 = [1.0 0.0 0.0;0.0 -1.0 0.0;0.0 0.0 -1.0];
        Ri2b = RM2*Ri2b;
        euangles = rotmatib2euler(Ri2b);
        omV      = RM2*omV ;
        isrefframechanged = true;
    end

    RM = RM2*RM;

    return euangles,omV,isrefframechanged,RM
end




