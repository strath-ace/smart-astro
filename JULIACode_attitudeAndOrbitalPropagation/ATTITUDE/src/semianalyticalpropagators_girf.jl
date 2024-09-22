"""
    attitudeaveragedprop_semianalytical_girf(civ,date0,params,tstep,tfinal)

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
        	planet_flat = control integer: it is equal to 1 to consider planet oblateness, equal to 0 otherwise.
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
    __________________________________________
    ******     
    atmosphere = Dict("atmosphericmodel"=>typeofmodel,"AtmM"=>AtmM);
                typeofmodel = 1 exponential atmospheric model
                              2 nrlmsise00
                              other kind of model still not implemented
                
                             
                AtmM =   
                        if typeofmodel==1
                                AtmM = matrix containing the exponential atmospheric model 
                                in each row: [min altitude [km], density [kg/m^3], scale height [km]
                   
                        
                        if typeofmodel==2
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
            	considerthermalCD = boolean (if false in the drag model the friction effects are neglected)                            
               
    
    NB: if the dictionary does not contain the key “typeofmodel” this is automatically added with default value 1 (exponential atmospheric model)

    NB: if the central body is the Earth (planetsunparams=Dict(“centralBodyIDX”=>3)), and the user wants to use the exponential atmospheric model, it is sufficient to give as an input a dictionary with the only key “typeofmodel” equal to 1. The key AtmM will be automatically generated, using the exponential model by Vallado (1997), (see Section 7). 
    
    NB: if the key considerthermalCD is missing, it is automatically set with the default value “true”.
            

    ____________________________________________
    *******     
    settings = Dictionary with keys 
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
        	checkchaoticity = Boolean, true if the user wants to check during the propagation the proximity to the chaotic region. If close to the chaotic region the semi-analytical propagation is interrupted.
        	Higherordercorrections = Boolena, true if the user wants to include higher order terms in the averaged attitude equations of motion
        	osculatingpropagationinchaoticregion = Boolean: if checkchaoticity = true and osculatingpropagationinchaoticregion= true in the chaotic region the propagation of the dynamics is performed by the osculating propagator. 
        	sunfixed = Boolean, if true the position of the Sun is assumed fixed in the space
        	sunposition = position of the Sun in space if sunfixed=true
    
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
    INFO : 
    semi-analytical propagation performed using modified sadov variables/sadov-like variables for trixial satellite or axisymmetric satellite with one moment of inertia different from the other two.
        
    semi-analytical propagation performed using andoyer variables/andoyer-like variables for axisymmetric satellite with equal moments of inertia 
"""
function attitudeaveragedprop_semianalytical_girf(civ,date0,params,tstep,tfinal)    

    # check inputs
    params=checksandhandleinputs(civ,date0,params,tstep,tfinal);

    # handle inputs
    kconstant,savar,skm,euangles,omV,isrefframechanged,RM,paramsint,equiEl = handleinputs(params,civ,date0,tfinal);
    # savar00 = append!([skm],savar[2:end])
    # println(savar00)
    
    # println("is ref frame changed? ",isrefframechanged)
   
    # propagation
    timevector = collect(0.0:tstep:tfinal);

    
    if kconstant==0.0 && paramsint[3]["MomentsOfInertia"][1]==paramsint[3]["MomentsOfInertia"][2] && paramsint[3]["MomentsOfInertia"][2]==paramsint[3]["MomentsOfInertia"][3]
        #  body with equal principal moments of inertia

        # in this case the action angle variables are the andoyer-serret variables 
        anvar = deepcopy(savar[1:6]);
        anvar[4] = anvar[4]+pi/2;

        # orbitci
        equimean0 = transformationci_keplerian(append!([skm],savar[2:6]),euangles,omV,equiEl,paramsint,1);      

        #attitudeci
        variableskind = 1
        if abs(1.0-abs(anvar[3]/anvar[2]))<1e-13
            #### when in proximity of the singularity sdelta=0.0 in the initial conditions (also case of planar problem)
            anvarmean0 = transformationci_axysimmetric_andoyerlike(anvar,euangles,omV,equimean0,paramsint,1);
            variableskind = 2
        else
            anvarmean0 = transformationci_axysimmetric(anvar,euangles,omV,equimean0,paramsint,1);
        end
        if anvarmean0[1]<0.0
            eumean0,omVmean0 = andoyer2euler(anvarmean0,paramsint[3]["MomentsOfInertia"]);
            eumean0,omVmean0,isrefframechanged,RM2 = selectreferenceframe_girf(IV,eumean0[1],eumean0[2],eumean0[3],omVmean0[1],omVmean0[2],omVmean0[3],satellite);
            IV = satellite["MomentsOfInertia"];
            RM = RM2*RM;
            paramsint[2] = satellite;
            anvarmean0 = euler2andoyer(eumean0,omVmean0,paramsint[3]["MomentsOfInertia"]);   
        end

        u0 = zeros(12,1);
        u0[1:6]  = anvarmean0;
        u0[7:12] = equimean0; 

        # propagation
        if paramsint[5]["includeDragTorque"] || paramsint[5]["includeDragAcc"]
            tV,res = prop_semianalytical_st_girf_axysimm_drag(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
        else
            if variableskind == 1
                # println("QUI")
                try 
                    tV,res = prop_semianalytical_st_girf_axysimm_andoyerlike_V2(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
                catch e
                    # println("se qui")
                    if e==e==ErrorException("singularity")
                        tV,res = prop_semianalytical_st_girf_axysimm_andoyerlike_V2(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
                    else
                        println(e)
                        tV = []; res=[];
                    end
                end
            else
                tV,res = prop_semianalytical_st_girf_axysimm_andoyerlike_V2(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
            end
        end
    else

        # triaxial problem or axisymmetric body with one of the principal moments of inertia different from the other two
        tV  = [];
        res = [];
        mval0 = kconstant*(1-skm)/skm;

        # check proximity to chaoticity region
        if (paramsint[5]["checkchaoticity"] && mval0 > 0.9999) || mval0==1.0
            if paramsint[5]["osculatingpropagationinchaoticregion"]
                println("chaotic region: osculating propagation")
                tV,res = attitudeprop_quat_girf(civ,date0,params,tstep,tfinal);
            else
                println("chaotic region")
            end
            return tV,res;
        end
        
        # orbit ci
        # println("qui ",equiEl)
        equimean0 = transformationci_keplerian(append!([skm],savar[2:6]),euangles,omV,equiEl,paramsint,1);  
        # println("here ",equimean0)
       
               
        # attitude ci
        variableskind = 1;
        if abs(1.0-abs(savar[3]/savar[2]))<1e-13  
            savarmean0 = transformationci_triaxial_sadovlike(kconstant,append!([skm],savar[2:6]),euangles,omV,equimean0,paramsint,paramsint[6],1);
            variableskind = 2
        else
            savarmean0 = transformationci_triaxial(kconstant,append!([skm],savar[2:6]),euangles,omV,equimean0,paramsint,paramsint[6],1);
        end

        # propagation
        if paramsint[5]["includeDragTorque"] || paramsint[5]["includeDragAcc"]
            u0 = append!(savarmean0,equimean0);
            tV,res = prop_semianalytical_st_girf_triax_drag(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
        else
            if variableskind == 1
                u0 = append!(savarmean0,equimean0);
                try
                    tV,res = prop_semianalytical_st_girf_triax(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"])
                catch e
                    if e==ErrorException("singularity")
                        println("savarlike")
                        # savarlike = sadov2sadovlike(savarmean0[2:6],savarmean0[1]);
                        # u0 = zeros(13,1);
                        # u0 = append!(savarlike,equimean0);
                        # tV,res = prop_semianalytical_st_girf_triax_sadovlike(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
                        tV,res = prop_semianalytical_st_girf_triax_sadovlike_V2(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
                    else
                        Base.error(e)
                    end
                end
            else
                # savarlike = sadov2sadovlike(savarmean0[2:6],savarmean0[1]);
                # u0 = zeros(13,1);
                # u0 = append!(savarlike,equimean0);
                # tV,res = prop_semianalytical_st_girf_triax_sadovlike(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
                u0 = append!(savarmean0,equimean0);
                tV,res = prop_semianalytical_st_girf_triax_sadovlike_V2(u0,paramsint,timevector,isrefframechanged,RM,paramsint[3]["MomentsOfInertia"]);
                
            end
        end
        
        # if in chaotic region osculating propagation
        if paramsint[5]["osculatingpropagationinchaoticregion"] && tV[end]<tfinal
            println("chaotic region: osculating propagation")
            
            # mean to osculating attitude
            if abs(1.0-res[end,3]/res[end,2])<1e-13
                savar0new = transformationci_triaxial_savarlike(kconstant,append!([res[end,25]],res[end,2:6]),res[end,10:12],res[end,7:9],res[end,19:24],paramsint,(tV[end]-tV[1])/3600/24+paramsint[6],-1);
            else
                savar0new = transformationci_triaxial(kconstant,append!([res[end,25]],res[end,2:6]),res[end,10:12],res[end,7:9],res[end,19:24],paramsint,(tV[end]-tV[1])/3600/24+paramsint[6],-1);
            end
            eu0new,om0new = sadov2euler(sadvar0new[2:6],savar0new[1],paramsint[3]["MomentsOfInertia"])
            if isrefframechanged
                RMT = Transpose(RM)
                paramsint[3]["MomentsOfInertia"] = RMT*IV;
                numberOfFacets = paramsint[3]["numberOfFacets"]
                if numberOfFacets >0
                    for kk = 1:numberOfFacets    
                        paramsint[3]["facets"][kk][3][1]=RMT*paramsint[3]["facets"][kk][3][1];
                        paramsint[3]["facets"][kk][3][2]=RMT*paramsint[3]["facets"][kk][3][2];
                        paramsint[3]["facets"][kk][3][3]=RMT*paramsint[3]["facets"][kk][3][3];
                        paramsint[3]["facets"][kk][3][4]=RMT*paramsint[3]["facets"][kk][3][4];
                        paramsint[3]["facets"][kk][4]=RMT*paramsint[3]["facets"][kk][4];
                    end
                    paramsint[3]["intrinsicMagneticMoment"] = RMT*paramsint[3]["intrinsicMagneticMoment"];
                end
                Ri2b = RMT*euler2ibrotmat(eu0new);
                eu0new = rotmatib2euler(Ri2b);
                om0new = RMT*om0new;
            end

            # mean to osculating orbital
            kepEl0new = res[end,19:24];

            # osculating propagation
            tVoscu,resosc = attitudeprop_quat_girf(append!(eu0new,om0new,kepEl0new),paramrsint[6]+tV[end]/3600.0/24.0,paramsint,tstep,tfinal-tV[end]);

            tV = append!(tV,tVoscu+tV[end]);
            res = [res;resosc]
        end
        
    end

    # output
   return tV,res;
end

##########################################################################
### Functions
##########################################################################

# checking inputs

function checksandhandleinputs(civ,date0,params,tstep,tfinal)

    ################ check civ

    if size(civ)[1] != 12
        Base.error("error: the length of the input vector must be 18");
    end

    if civ[4] == 0 && civ[5] == 0 && civ[6] == 0 
        Base.error("the semi-analytical propagation can be applied only if the body is rotating: the magnitude of the angular velocity has to be larger than 0");
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

    if haskey(p5,"maxzonalharmonics") == false 
        if p5["includeZonalHarmsAcc"]
            p5 = merge(p5,Dict("maxzonalharmonics"=>[5]));
        end
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

    if haskey(p5,"includeDragAcc") == false
        p5 = merge(p5,Dict("includeDragAcc"=>false));
    end

    # if haskey(p5,"includeEclipsesEffects") == false
    #     if p5["includeSrpTorque"]
    #         p5 = merge(p5,Dict("includeEclipsesEffects"=>true));
    #     else
    #         p5 = merge(p5,Dict("includeEclipsesEffects"=>false));
    #     end
    # end

    if haskey(p5,"higherordercorrections") == false
        p5 = merge(p5,Dict("higherordercorrections"=>false));
    end

    if haskey(p5,"checkchaoticity") == false
        p5 = merge(p5,Dict("checkchaoticity"=>false));
    end
    
    if haskey(p5,"osculatingpropagationinchaoticregion") == false
        p5 = merge(p5,Dict("osculatingpropagationinchaoticregion"=>false));
    end

    if haskey(p5,"sunFixed") == false
        p5 = merge(p5,Dict("sunFixed"=>false));
    end

    if haskey(p5,"sunposition") == false && p5["sunFixed"]
        Base.error("NO")
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
        p3 = merge(p3,Dict("CD"=>2.0));
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
            if haskey(p1,"planet_flat") == false && p5["includeDragAcc"]
                p1 =  merge(p1,Dict("planet_flat"=>0.0));
            end
            if haskey(p1,"oblatenesscoeff") == false && p5["includeDragAcc"]
                if p1["planet_flat"]!=0.0
                    Base.error("missing central body oblateness");
                end
            end

        else # Earth
            earthconstant = astroconstants(3);
            generalconstant = astroconstants(0);
            
            p1 = Dict("centralBodyIDX"=>3,"muPlanet"=>earthconstant[1],"rPlanet"=>earthconstant[2],"oblatenesscoeff"=>earthconstant[4],"muM"=>earthconstant[5],"planetRotation"=>earthconstant[6]); 
            p1 = merge(p1,Dict("psrp"=>generalconstant[3],"psrp_refdist"=>generalconstant[2])); 
            p1 = merge(p1,Dict("zonalharmonicscoff"=>earthconstant[7:10]));
            p1 = merge(p1,Dict("perturbingbodiesid"=>[11]));

            if haskey(p1,"planet_flat") == false && (p5["includeDragAcc"] || p5["includeDragTorque"])
                p1 =  merge(p1,Dict("planet_flat"=>0.0));
            end
        end
    end
    
    # check ref frame info
    if  haskey(p2,"ecliptic2inertial")==false && (p5["includeSrpAcc"] || p5["includeSrpTorque"] || p5["includeSunGravityAcc"]  )
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

    # atmospheric model
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

    if haskey(p4,"considerthermalCD") == false
        p4 = merge(p4,Dict("considerthermalCD"=>true));
    end

    return [p1,p2,p3,p4,p5];
end

# adjusting inputs

function handleinputs(params,civ,date0,tfinal)

    planetsunparams      = deepcopy(params[1]);
    inertialrefframeinfo = deepcopy(params[2]);
    satellite            = deepcopy(params[3]);
    atmosphere           = deepcopy(params[4]);
    settings             = deepcopy(params[5]);

    phi   = deepcopy(civ[1]);
    theta = deepcopy(civ[2]);
    psi   = deepcopy(civ[3]);
    p     = deepcopy(civ[4]);
    q     = deepcopy(civ[5]);
    r     = deepcopy(civ[6]);
    kepEl = deepcopy(civ[7:12]);
    mjd2000_0 = date2mjd2000(date0);

    # check if the reference rotating frame is correct
    IV = satellite["MomentsOfInertia"];
    if ! (IV[1]<=IV[2] && IV[2]<=IV[3])
        Base.error("please for the inputs select a body reference frame in principal axis of inertia such that IV[1]<=IV[2]<=IV[3]: ",IV);
    end
  
    # choice of rotating reference frame
    euangles,omV,isrefframechanged,RM = selectreferenceframe_girf(IV,phi,theta,psi,p,q,r,satellite);
    IV = satellite["MomentsOfInertia"];
    if IV[3]*(IV[2]-IV[1])-IV[1]*(IV[3]-IV[2])==0.0
        kconstant = 0.0
    else
        kconstant = IV[3]/IV[1]*(IV[2]-IV[1])/(IV[3]-IV[2]);
    end
    satellite = merge(satellite,Dict("k_constant"=>kconstant))

    # input in sadov variables
    savar,skm = euler2sadov(euangles,omV,IV);

    # input in equinoctial elements
    equi = kep2equi(kepEl,3);

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

    # initialisation for polar axis
    if settings["includeDragTorque"] || settings["includeMagneticTorque"] || settings["includeDragAcc"]
        # Recef2eci = r_ecef_to_eci(J2000(), PEF(),mjd20002jd(mjd2000_0))
        # paV = (inertialrefframeinfo["equatorial2inertial"]*Recef2eci*[0.0,0.0,1.0]);
        paV = (inertialrefframeinfo["equatorial2inertial"]*[0.0,0.0,1.0]);
        inertialrefframeinfo = merge(inertialrefframeinfo,Dict("cpa1"=>paV[1],"cpa2"=>paV[2],"cpa3"=>paV[3]))
    else
        inertialrefframeinfo = merge(inertialrefframeinfo,Dict("cpa1"=>0.0,"cpa2"=>0.0,"cpa3"=>0.0))
    end

    # multiply optic coeff for reference solar radiation pressure
    if settings["includeSrpAcc"] || settings["includeSrpTorque"]
        for jj=1:satellite["numberOfFacets"]
            satellite["facets"][jj][1] = satellite["facets"][jj][1]*planetsunparams["psrp"];
        end
    end
        
    # precompute constants for srp perturbation
    if settings["includeSrpTorque"]
        numberOfFacets = get(satellite,"numberOfFacets",0);
        facets = get(satellite,"facets",0.0);
        facetsmod = Vector{Vector{Any}}(undef,size(facets,1));
        for kk = 1:numberOfFacets 
            VVT1 = srp_initforaverage_T1(facets[kk][1][1],facets[kk][2],facets[kk][3][1],facets[kk][4]);
            VVT2 = srp_initforaverage_T2(facets[kk][1][2],facets[kk][2],facets[kk][3][1],facets[kk][4]);
            VVT3 = srp_initforaverage_T3(facets[kk][1][3],facets[kk][2],facets[kk][3][1],facets[kk][4]);
            facetmod = [facets[kk][1],facets[kk][2],facets[kk][3],facets[kk][4],VVT1,VVT2,VVT3]
            facetsmod[kk] = facetmod;
        end
        satellite["facets"]=facetsmod;
    end

    # precompute constants for drag perturbation and numerical computation of the mean over M 
    timeUpdate              = tfinal;
    correctiondragJ2        = false;
    couplingstrategydrag    = 2;
    if settings["includeDragTorque"] || settings["includeDragAcc"]
        # if atmosphere["atmosphericmodel"]==2 || settings["includeZonalHarmsAcc"] || settings["includeThirdBodyAcc"] || settings["includeSunGravityAcc"] ||settings["includeSrpAcc"] || settings["includeDragAcc"]
        #     timeUpdate  = pi*sqrt(kepEl[1]^3/planetsunparams["muPlanet"]);
        # end 
        timeUpdate  = pi*sqrt(kepEl[1]^3/planetsunparams["muPlanet"]);
    end
    if settings["includeDragTorque"]
        numberOfFacets = get(satellite,"numberOfFacets",0);
        facets = get(satellite,"facets",0.0);
        facetsmod = Vector{Vector{Any}}(undef,size(facets,1));
        VVT1tot = zeros(30)
        VVT2tot = zeros(60)
        VVT3tot = zeros(105)
        VVT4tot = zeros(63)
 
        for kk = 1:numberOfFacets 
            VVT1 = drag_initforaverage_T1(facets[kk][2],facets[kk][3][1],facets[kk][4]);
            VVT2 = drag_initforaverage_T2(facets[kk][2],facets[kk][3][1],facets[kk][4]);
            VVT3 = drag_initforaverage_T3(facets[kk][2],facets[kk][3][1],facets[kk][4]);
            VVT4 = drag_initforaverage_T4(facets[kk][2],facets[kk][3][1],facets[kk][4]);

            if atmosphere["atmosphericmodel"]==1
                VVT1 = VVT1 * satellite["CD"];
                VVT2 = VVT2 * satellite["CD"];
                VVT3 = VVT3 * satellite["CD"];
                VVT4 = VVT4 * satellite["CD"];
            end
            
            VVT1tot = VVT1tot + VVT1
            VVT2tot = VVT2tot + VVT2
            VVT3tot = VVT3tot + VVT3
            VVT4tot = VVT4tot + VVT4
          
            ll = size(facets[kk],1);
            if ll==7
                facetmod = [facets[kk][1],facets[kk][2],facets[kk][3],facets[kk][4],facets[kk][5],facets[kk][6],facets[kk][7],VVT1,VVT2,VVT3,VVT4];
            else
                facetmod = [facets[kk][1],facets[kk][2],facets[kk][3],facets[kk][4],zeros(30,1),zeros(60,1),zeros(105,1),VVT1,VVT2,VVT3,VVT4];
            end
            facetsmod[kk] = facetmod;
        end
        satellite["facets"]=facetsmod; 
        satellite = merge(satellite, Dict("constantstermsdragtorque"=>[VVT1tot,VVT2tot,VVT3tot,VVT4tot]));  
    end

    if settings["includeDragAcc"]
        numberOfFacets = get(satellite,"numberOfFacets",0);
        facets = get(satellite,"facets",0.0);
        facetsmod = Vector{Vector{Any}}(undef,size(facets,1));
        VVTV = zeros(35)
        VVTN = zeros(56)
        for kk = 1:numberOfFacets 
            VVTV = VVTV + drag_initforaverage_TV(facets[kk][2],facets[kk][4]) 
            VVTN = VVTN + drag_initforaverage_TN(facets[kk][2],facets[kk][4])
        end
        if atmosphere["atmosphericmodel"]==1
            VVTV = VVTV*satellite["CD"];
            VVTN = VVTN*satellite["CD"];
        end
        satellite = merge(satellite,Dict("constantstermsdragacc"=>[VVTV,VVTN]));

    end
    
    if settings["includeDragAcc"]
        if settings["includeZonalHarmsAcc"] && planetsunparams["zonalharmonicscoff"][1] !=0.0
            correctiondragJ2 = true;
        end
    end
    if settings["includeDragAcc"] || settings["includeDragTorque"]
        meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(equi,mjd2000_0,inertialrefframeinfo["cpa1"],inertialrefframeinfo["cpa2"],inertialrefframeinfo["cpa3"],
            planetsunparams["rPlanet"],planetsunparams["muPlanet"],planetsunparams["planetRotation"],planetsunparams["planet_flat"],planetsunparams["oblatenesscoeff"],
            planetsunparams["zonalharmonicscoff"][1],settings["includeDragTorque"],settings["includeDragAcc"],atmosphere,correctiondragJ2,couplingstrategydrag,transpose(inertialrefframeinfo["equatorial2inertial"]));
    else
        meanOverMDrag_A = zeros(91);
        meanOverMDrag_O = zeros(735);
    end

    draginfo = Dict("atmosphericinfo"=>atmosphere,"averageoverM_attitude"=>meanOverMDrag_A,"timeupdate"=>timeUpdate,"correction"=>correctiondragJ2,"averageoverM_orbit"=>meanOverMDrag_O,"couplingstrategy"=>couplingstrategydrag);

    paramsint =[planetsunparams;inertialrefframeinfo;satellite;draginfo;settings;mjd2000_0]; 

    return kconstant,savar,skm,euangles,omV,isrefframechanged,RM,paramsint,equi; 

end
 
# selection of rotating frame
function selectreferenceframe_girf(IV,phi,theta,psi,p,q,r,satellite)

    RM = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
    isrefframechanged = false;
    if IV[1]==IV[2] && IV[2]==IV[3]
        # singularity if p=q=0 and r!=0
        # condition = p>r || q>r 
        condition = abs(p)>1e-5 || abs(q)>1e-5
    else
        # sadov variables well defined only if 0<=m<1
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
        omV      = [r;-q;p];
        isrefframechanged = true;
        A = deepcopy(IV[3]); 
        B = deepcopy(IV[2]); 
        C = deepcopy(IV[1]); 
        IV = [A,B,C];
    end

    # select a reference frame such as omV[3]>0
    RM2 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
    if omV[3]<0.0
        Ri2b = euler2ibrotmat(euangles);
        RM2 = [1.0 0.0 0.0;0.0 -1.0 0.0;0.0 0.0 -1.0];
        Ri2b = RM2*Ri2b;
        euangles = rotmatib2euler(Ri2b);
        omV[2]   = -omV[2];
        omV[3]   = -omV[3];
        isrefframechanged = true;
    end
    
    RM = RM2*RM;

    if isrefframechanged
        satellite["MomentsOfInertia"] = IV;
        numberOfFacets = satellite["numberOfFacets"];
        if numberOfFacets >0
            for kk = 1:numberOfFacets    
                satellite["facets"][kk][3][1]=RM*satellite["facets"][kk][3][1];
                satellite["facets"][kk][3][2]=RM*satellite["facets"][kk][3][2];
                satellite["facets"][kk][3][3]=RM*satellite["facets"][kk][3][3];
                satellite["facets"][kk][3][4]=RM*satellite["facets"][kk][3][4];
                satellite["facets"][kk][4]=RM*satellite["facets"][kk][4];
            end
        end
        satellite["intrinsicMagneticMoment"] = RM*satellite["intrinsicMagneticMoment"];
    end
    # println(isrefframechanged)

    return euangles,omV,isrefframechanged,RM
end

# propagation and variables transformations 

function prop_semianalytical_st_girf_triax_drag(u0,paramsint,timevector,isrefframechanged,RM,IV)
    t0 = timevector[1];
    ts = timevector.-t0;
    params = deepcopy(paramsint);
    tV = deepcopy(ts);
    tpropagated = 0.0;
    u0prop = deepcopy(u0);
    mcheck = true;
    ssign = ts[end]/abs(ts[end]);
    tspanprop = [0.0 ssign*params[4]["timeupdate"]];
    resM = zeros(size(ts)[1],12);
  
    meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[7:12],params[6],params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
            params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
            params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
            params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));
    params[4]["averageoverM_attitude"] = meanOverMDrag_A;
    params[4]["averageoverM_orbit"] =meanOverMDrag_O  ;

    while abs(tpropagated)<abs(ts[end]) && mcheck   
        println(" qui ",tpropagated/24/3600) 
        variableskind = 1;
        # propagation step 1
        solM = [];
        if paramsint[5]["checkchaoticity"]
            try 
                cb = ContinuousCallback(conditionpatch,affectpatch!)
                probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0prop,tspanprop,params);
                solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback=cb);
            catch e
                if e==ErrorException("singularity")
                    cb1 = ContinuousCallback(conditionpatch,affectpatch!)
                    cb2 = ContinuousCallback(conditionjump,affectjump!)
                    cbs = CallbackSet(cb1, cb2)
                    u0propnew = zeros(13,1);
                    savarlike = sadov2sadovlike(u0prop[2:6],u0prop[1]);
                    u0propnew = append!(savarlike,u0prop[7:12]);
                    probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0propnew,tspanprop,params);          
                    solM  = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback = cbs);
                    variableskind = 2;
                else
                    Base.error("there is an error: ",e)
                end
            end
        else
            try
                probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0prop,tspanprop,params);
                solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);
            catch e
                if e==ErrorException("singularity")
                    println("savarlike")
                    u0propnew = zeros(13,1);
                    savarlike = sadov2sadovlike(u0prop[2:6],u0prop[1]);
                    u0propnew = append!(savarlike,u0prop[7:12]);
                    cb    = ContinuousCallback(conditionjump,affectjump!)
                    probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0propnew,tspanprop,params);          
                    solM  = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback = cb);
                    variableskind = 2;
                else
                    Base.error("there is an error: ",e)
                end
            end
        end

        zetaend = solM.u[end][1];
    
        if paramsint[5]["checkchaoticity"] && (zetaend-1.0/(1.0+0.9999/(paramsint[3]["k_constant"])))<1e-12
            mcheck = false;
        end           

        if mcheck
            tM = solM.t .+ tpropagated;
            idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
            if maximum(idx)
                solMsol  = solM(ts[idx].-tpropagated);
                resMjj = mapreduce(permutedims,vcat,solMsol.u)
                if variableskind == 1
                    resM[idx,:] = resMjj;
                else
                    resMtemporary = zeros(size(resMjj)[1],12)
                    for jj = 1 : size(resMjj)[1]
                        savarjj,skmjj = sadovlike2sadov(resMjj[jj,1:7],params[3]["MomentsOfInertia"]);
                        resMtemporary[jj,1] = skmjj
                        resMtemporary[jj,2:6] = savarjj[2:end]
                        resMtemporary[jj,7:12] = resMjj[jj,8:end];
                    end
                    resM[idx,:] = resMtemporary;
                end
            end
            # update
            mjd2000 = params[6] +  solM.t[end]/24/3600;
            resC   =  mapreduce(permutedims,vcat,solM.u);
            if variableskind == 1
                u0prop = resC[end,:];
            else
                savarjj,skmjj = sadovlike2sadov(resC[end,1:7],params[3]["MomentsOfInertia"]);
                u0prop[1] = skmjj
                u0prop[2:6] = savarjj[2:end]
                u0prop[7:end] = resC[end,8:end];
            end
            meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[7:12],mjd2000,params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
                    params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
                    params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
                    params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));

            params[4]["averageoverM_attitude"] = meanOverMDrag_A;
            params[4]["averageoverM_orbit"] =meanOverMDrag_O ;
            params[4]["timeupdate"]   = pi*sqrt(u0prop[7]^3.0/params[1]["muPlanet"]);
            params[6] = mjd2000;
            tpropagated = tM[end];
            tspanprop = [0.0 ssign*params[4]["timeupdate"]];
        else
            tM = solM.t + tpropagated;
            idx = ts.>= tpropagated .&& ts.<=tM[end];
            if maximum(idx)
                solMsol  = solM(ts[idx].-tpropagated);
                resMjj = mapreduce(permutedims,vcat,solMsol.u)
                resM[idx,:] = resMjj;
            end
            idxfinal = length(ts);
            for jj = 1 : size(ts)[1]
                if ts[jj]>tpropagated+tM[end]
                    idxfinal = jj;
                    break;
                end
            end
            resM[idxfinal,:] =mapreduce(permutedims,vcat,solM.u)[end,:];
            tV = append!(ts[1:idxfinal-1],tpropagated+tM[end]);
        end
    end
    tV = tV .+ t0;

    for jj = 1 : size(resM)[1]
        resM[jj,4]=
        mod(resM[jj,4],2.0*pi);
        resM[jj,5]=mod(resM[jj,5],2.0*pi);
        resM[jj,6]=mod(resM[jj,6],2.0*pi);
    end

    # transformation of variables
    res = zeros(size(resM,1),40);
    if isrefframechanged
        RMT = deepcopy(transpose(RM));
    end

    for jj = 1 : size(resM)[1]
        euanglesjj,omVjj,Jljj = sadov2eulerAndJl(resM[jj,2:6],resM[jj,1],IV);
        # savarsol  = transformationci_triaxial(paramsint[3]["k_constant"],resM[jj,1:6],zeros(3),zeros(3),resM[jj,7:12],paramsint,-1);
        # euanglesjj,omVjj,Jljj = sadov2eulerAndJl(savarsol[2:6],savarsol[1],IV);

        IVmod = deepcopy(IV);
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RMT*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RMT*omVjj;
            IVmod      = RMT*IV;
            IVmod[1]   = abs(IVmod[1]);
            IVmod[2]   = abs(IVmod[2]);
            IVmod[3]   = abs(IVmod[3]);
        end

        anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        kepjj = equi2kep(resM[jj,7:12],3);
        savarlikejj = sadov2sadovlike(resM[jj,2:6],resM[jj,1])
        energy = (IV[1]*omVjj[1]^2+IV[2]*omVjj[2]^2+IV[3]*omVjj[3]^2)/2;
        Jd = anvar[2]^2/2/energy;
        res[jj,:] = append!([Jljj],resM[jj,2:6],omVjj,euanglesjj,anvar,resM[jj,7:12],resM[jj,1],energy,Jd,savarlikejj,kepjj);
        # res[jj,:] = append!([Jljj],savarsol[2:6],omVjj,euanglesjj,anvar,resM[jj,7:12],savarsol[1],energy,Jd,savarlikejj,kepjj);
    end
    
    return tV, res;
end

function prop_semianalytical_st_girf_triax(u0,paramsint,timevector,isrefframechanged,RM,IV)
    t0 = timevector[1];
    ts = timevector.-t0;
        
    tspan = [0.0 ts[end]];    
    if paramsint[5]["checkchaoticity"]
        cb = ContinuousCallback(conditionpatch,affectpatch!)
        probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0,tspan,paramsint);
        solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback=cb);
        solMsol = solM(ts[ts.<=solM.t[end]]);
        resM = mapreduce(permutedims,vcat,solMsol.u);
        tV = append!(solMsol.t);
    else
        probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0,tspan,paramsint);
        solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10); #,callback=cb);
        solM = solM(ts);
        resM = mapreduce(permutedims,vcat,solM.u);
        tV = solM.t;
    end



   
    tV = tV .+ t0;
    for jj = 1 : size(resM)[1]
        resM[jj,4]=mod(resM[jj,4],2.0*pi);
        resM[jj,5]=mod(resM[jj,5],2.0*pi);
        resM[jj,6]=mod(resM[jj,6],2.0*pi);
    end

    # transformation of variables
    res = zeros(size(resM,1),40);
    if isrefframechanged
        RMT = deepcopy(transpose(RM));
    end

    for jj = 1 : size(resM)[1]
        euanglesjj,omVjj,Jljj = sadov2eulerAndJl(resM[jj,2:6],resM[jj,1],IV);
        # savarsol  = transformationci_triaxial(paramsint[3]["k_constant"],resM[jj,1:6],zeros(3),zeros(3),resM[jj,7:12],paramsint,-1);
        # euanglesjj,omVjj,Jljj = sadov2eulerAndJl(savarsol[2:6],savarsol[1],IV);

        IVmod = deepcopy(IV);
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RMT*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RMT*omVjj;
            IVmod      = RMT*IV;
            IVmod[1]   = abs(IVmod[1]);
            IVmod[2]   = abs(IVmod[2]);
            IVmod[3]   = abs(IVmod[3]);
        end

        anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        kepjj = equi2kep(resM[jj,7:12],3);
        savarlikejj = sadov2sadovlike(resM[jj,2:6],resM[jj,1])
        energy = (IV[1]*omVjj[1]^2+IV[2]*omVjj[2]^2+IV[3]*omVjj[3]^2)/2;
        Jd = anvar[2]^2/2/energy;
        res[jj,:] = append!([Jljj],resM[jj,2:6],omVjj,euanglesjj,anvar,resM[jj,7:12],resM[jj,1],energy,Jd,savarlikejj,kepjj);
        # res[jj,:] = append!([Jljj],savarsol[2:6],omVjj,euanglesjj,anvar,resM[jj,7:12],savarsol[1],energy,Jd,savarlikejj,kepjj);
    end
    
    return tV, res;
end

function prop_semianalytical_st_girf_triax_sadovlike_V2(u0,paramsint,timevector,isrefframechanged,RM,IV)

    t0 = timevector[1];
    ts = timevector.-t0;
    tspan = [0.0 ts[end]];
    ssign = ts[end]/abs(ts[end]);


    u0prop = deepcopy(u0);
    tpropagated = 0.0;
    tspanprop  = [0.0 ssign*minimum([86400,abs(ts[end])])];
    resM = zeros(size(ts)[1],12);
    mcheck = true;
    tV = deepcopy(ts);
    while abs(tpropagated)<abs(ts[end]) && mcheck
        println("time propagated ", tpropagated/24/3600)
        variableskind = 1;
        solM = [];
        if paramsint[5]["checkchaoticity"]
            try 
                cb = ContinuousCallback(conditionpatch,affectpatch!)
                probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0prop,tspanprop,paramsint);
                solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback=cb);
            catch e
                if e==ErrorException("singularity")
                    cb1 = ContinuousCallback(conditionpatch,affectpatch!)
                    cb2 = ContinuousCallback(conditionjump,affectjump!)
                    cbs = CallbackSet(cb1, cb2)
                    u0propnew = zeros(13,1);
                    savarlike = sadov2sadovlike(u0prop[2:6],u0prop[1]);
                    u0propnew = append!(savarlike,u0prop[7:12]);
                    probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0propnew,tspanprop,paramsint);          
                    solM  = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback = cbs);
                    variableskind = 2;
                else
                    Base.error("there is an error: ",e)
                end
            end
        else
            try 
                probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0prop,tspanprop,paramsint);
                solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);
            catch e
                if e==ErrorException("singularity")
                    u0propnew = zeros(13,1);
                    savarlike = sadov2sadovlike(u0prop[2:6],u0prop[1]);
                    u0propnew = append!(savarlike,u0prop[7:12]);
                    cb    = ContinuousCallback(conditionjump,affectjump!)
                    probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0propnew,tspanprop,paramsint);          
                    solM  = solve(probM,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10,callback = cb);
                    variableskind = 2;
                else
                    Base.error("there is an error: ",e)
                end
            end
        end
        
        zetaend = solM.u[end][1];
    
        if paramsint[5]["checkchaoticity"] && (zetaend-1.0/(1.0+0.9999/(paramsint[3]["k_constant"])))<1e-12
            mcheck = false;
        end  
            
        if mcheck
            tM = solM.t .+ tpropagated;
            # save
            idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
            if maximum(idx)
                solMsol  = solM(ts[idx].-tpropagated);
                resMjj = mapreduce(permutedims,vcat,solMsol.u)
                if variableskind == 1
                    resM[idx,:] = resMjj;
                else
                    resMtemporary = zeros(size(resMjj)[1],12)
                    for jj = 1 : size(resMjj)[1]
                        savarjj,skmjj = sadovlike2sadov(resMjj[jj,1:7],paramsint[3]["MomentsOfInertia"]);
                        resMtemporary[jj,1] = skmjj
                        resMtemporary[jj,2:6] = savarjj[2:end]
                        resMtemporary[jj,7:12] = resMjj[jj,8:end];
                    end
                    resM[idx,:] = resMtemporary;
                end
            end
            # update
            mjd2000 = paramsint[6] +  solM.t[end]/24/3600;
            resC   =  mapreduce(permutedims,vcat,solM.u);
            if variableskind == 1
                u0prop = resC[end,:];
            else
                savarjj,skmjj = sadovlike2sadov(resC[end,1:7],paramsint[3]["MomentsOfInertia"]);
                u0prop[1] = skmjj
                u0prop[2:6] = savarjj[2:end]
                u0prop[7:end] = resC[end,8:end];
            end
            tpropagated = tM[end];
            paramsint[6] = mjd2000;
            tspanprop = [0.0 ssign*minimum([86400,abs(ts[end]-tpropagated)])];
                 
        else
            tM = solM.t + tpropagated;
            idx = ts.>= tpropagated .&& ts.<=tM[end];
            solMsol  = solM(ts[idx]-tpropagated);
            resMjj = mapreduce(permutedims,vcat,solMsol.u)
            resM[idx,:] = resMjj;

            idxfinal = length(ts);
            for jj = 1 : size(ts)[1]
                if ts[jj]>tpropagated+tM[end]
                    idxfinal = jj;
                    break;
                end
            end
            resM[idxfinal,:] =mapreduce(permutedims,vcat,solM.u)[end,:];
            tV = append!(ts[1:idxfinal-1],tpropagated+tM[end]);
        end    
    end
    tV = tV .+t0;

    for jj = 1 : size(resM)[1]
        resM[jj,4]=mod(resM[jj,4],2.0*pi);
    end

    # transformation of variables
    res = zeros(size(resM,1),40);
    if isrefframechanged
        RMT = deepcopy(transpose(RM));
    end

    # savarM = zeros(size(resM,1),7);
    # for jj = 1 : size(resM)[1]
    #     savarjj,skmjj = sadovlike2sadov(resM[jj,1:7],IV);
    #     savarM[jj,:]  = append!(savarjj,skmjj);
    #     resM[jj,5]    = savarM[jj,5]+savarM[jj,3]/savarM[jj,2]*savarM[jj,6];
    # end

    for jj = 1 : size(resM)[1]
        euanglesjj,omVjj,Jljj = sadov2eulerAndJl(resM[jj,2:6],resM[jj,1],IV);
        IVmod = deepcopy(IV);
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RMT*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RMT*omVjj;
            IVmod      = RMT*IV;
            IVmod[1]   = abs(IVmod[1]);
            IVmod[2]   = abs(IVmod[2]);
            IVmod[3]   = abs(IVmod[3]);
        end

        anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        kepjj = equi2kep(resM[jj,7:12],3);
        savarlikejj = sadov2sadovlike(resM[jj,2:6],resM[jj,1])
        energy = (IV[1]*omVjj[1]^2+IV[2]*omVjj[2]^2+IV[3]*omVjj[3]^2)/2;
        Jd = anvar[2]^2/2/energy;
        res[jj,:] = append!([Jljj],resM[jj,2:6],omVjj,euanglesjj,anvar,resM[jj,7:12],resM[jj,1],energy,Jd,savarlikejj,kepjj);
        
    end    
    
    return tV, res;
end

function prop_semianalytical_st_girf_axysimm_drag(u0,paramsint,timevector,isrefframechanged,RM,IV)
    # println("here")
    t0 = timevector[1];
    ts = timevector.-t0;

    params = deepcopy(paramsint);
    tV = deepcopy(ts);
    tpropagated = 0.0;
    u0prop = deepcopy(u0);
    ssign = ts[end]/abs(ts[end]);
    tspanprop = [0.0 ssign*params[4]["timeupdate"]];
    resM = zeros(size(ts)[1],12);

    meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[7:12],params[6],params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
            params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
            params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
            params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));
    params[4]["averageoverM_attitude"] = meanOverMDrag_A;
    params[4]["averageoverM_orbit"] =meanOverMDrag_O  ;

    while abs(tpropagated)<abs(ts[end])  
        println(" time propagated ",tpropagated/24/3600) 
        variableskind = 1;
        solM = [];
        try
            probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm,u0prop,tspanprop,params);
            solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);        
        catch e
            if e==ErrorException("singularity")
                anvarlike = andoyer2andoyerlike(u0prop[1:6]);
                u0propnew = append!(anvarlike,u0prop[7:end]);
                probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm_andoyerlike,u0propnew,tspanprop,params);
                solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);  
                variableskind = 2;
            else
                Base.error("error: ",e);     
            end
        end

        tM = solM.t .+ tpropagated;
        idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
        if maximum(idx)
            solMsol  = solM(ts[idx].-tpropagated);
            resMjj = mapreduce(permutedims,vcat,solMsol.u)
            if variableskind == 1
                resM[idx,:] = resMjj;
            else
                resMtemporary = zeros(size(resMjj)[1],12)
                for jj = 1 : size(resMjj)[1]
                    anvarjj = andoyerlike2andoyer(resMjj[jj,1:7]);
                    resMtemporary[jj,1:6] = anvarjj
                    resMtemporary[jj,7:12] = resMjj[jj,8:end];
                end
                resM[idx,:] = resMtemporary;
            end
        end
        # update
        mjd2000 = params[6] +  solM.t[end]/24/3600;
        resC   =  mapreduce(permutedims,vcat,solM.u);
        if variableskind == 1
            u0prop = resC[end,:];
        else
            anvarjj = sadovlike2sadov(resC[end,1:7]);
            u0prop[1:6] = anvarjj
            u0prop[7:end] = resC[end,8:end];
        end
        meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[7:12],mjd2000,params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
                params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
                params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
                params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));

        params[4]["averageoverM_attitude"] = meanOverMDrag_A;
        params[4]["averageoverM_orbit"] =meanOverMDrag_O ;
        params[4]["timeupdate"]   = pi*sqrt(u0prop[7]^3.0/params[1]["muPlanet"]);
        params[6] = mjd2000;
        tpropagated = tM[end];
        tspanprop = [0.0 ssign*params[4]["timeupdate"]];
        
    end

    tV = tV .+ t0;
    
   
    # handle sol
    for jj = 1 : size(resM)[1]
        resM[jj,4]=mod(resM[jj,4],2.0*pi);
        resM[jj,5]=mod(resM[jj,5],2.0*pi);
        resM[jj,6]=mod(resM[jj,6],2.0*pi);
    end

    # transformation of variables
    res = zeros(size(resM,1),40);
    if isrefframechanged
        RMT = deepcopy(transpose(RM));
    end

    for jj = 1 : size(resM)[1]
        anvar = resM[jj,1:6];
        euanglesjj,omVjj = andoyer2euler(anvar,IV);
        sadovjj = resM[jj,1:6];
        sadovjj[4] = mod(sadovjj[4]-pi/2,2*pi); 
        IVmod = IV;
        anvarlike = andoyer2andoyerlike(anvar);
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RMT*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RMT*omVjj;
            IVmod = RMT*IV;
            IVmod[1] = abs(IVmod[1]);
            IVmod[2] = abs(IVmod[2]);
            IVmod[3] = abs(IVmod[3]);
            anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        end
        energy = (IVmod[1]*omVjj[1]^2+IVmod[2]*omVjj[2]^2+IVmod[3]*omVjj[3]^2)/2;
        Jd = anvar[2]^2/2/energy;
        kepjj = equi2kep(resM[jj,7:12],3);
        res[jj,:] = append!(sadovjj,omVjj,euanglesjj,anvar,resM[jj,7:12],NaN,energy,Jd,anvarlike,kepjj);
    end
    
    return tV, res;
end

function prop_semianalytical_st_girf_axysimm(u0,paramsint,timevector,isrefframechanged,RM,IV)

    
    t0 = timevector[1];
    ts = timevector.-t0;
    tspan = [0.0 ts[end]];
    probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm,u0,tspan,paramsint);
    solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13);
    solM = solM(ts);
    resM = mapreduce(permutedims,vcat,solM.u);
    tV = solM.t;
    tV = tV .+ t0;
    
   
    # handle sol
    for jj = 1 : size(resM)[1]
        resM[jj,4]=mod(resM[jj,4],2.0*pi);
        resM[jj,5]=mod(resM[jj,5],2.0*pi);
        resM[jj,6]=mod(resM[jj,6],2.0*pi);
    end

    # transformation of variables
    res = zeros(size(resM,1),40);
    if isrefframechanged
        RMT = deepcopy(transpose(RM));
    end

    for jj = 1 : size(resM)[1]
        anvar = resM[jj,1:6];
        euanglesjj,omVjj = andoyer2euler(anvar,IV);
        sadovjj = resM[jj,1:6];
        sadovjj[4] = mod(sadovjj[4]-pi/2,2*pi); 
        IVmod = IV;
        anvarlike = andoyer2andoyerlike(anvar);
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RMT*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RMT*omVjj;
            IVmod = RMT*IV;
            IVmod[1] = abs(IVmod[1]);
            IVmod[2] = abs(IVmod[2]);
            IVmod[3] = abs(IVmod[3]);
            anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        end
        energy = (IVmod[1]*omVjj[1]^2+IVmod[2]*omVjj[2]^2+IVmod[3]*omVjj[3]^2)/2;
        Jd = anvar[2]^2/2/energy;
        kepjj = equi2kep(resM[jj,7:12],3);
        res[jj,:] = append!(sadovjj,omVjj,euanglesjj,anvar,resM[jj,7:12],NaN,energy,Jd,anvarlike,kepjj);
    end
    
    return tV, res;
end

function prop_semianalytical_st_girf_axysimm_andoyerlike_V2(u0,paramsint,timevector,isrefframechanged,RM,IV)
   
    t0 = timevector[1];
    ts = timevector.-t0;

    params = deepcopy(paramsint);
    tV = deepcopy(ts);
    tpropagated = 0.0;
    u0prop = deepcopy(u0);
    ssign = ts[end]/abs(ts[end]);
    tspanprop = [0.0 ssign*minimum([86400,abs(ts[end])])];
    resM = zeros(size(ts)[1],12);

    while abs(tpropagated)<abs(ts[end]) 
        println(" time propagated ",tpropagated/24/3600) 
        variableskind = 1;
        solM = [];
        try
            probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm,u0prop,tspanprop,params);
            solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);        
        catch e
            if e==ErrorException("singularity")
                anvarlike = andoyer2andoyerlike(u0prop[1:6]);
                u0propnew = append!(anvarlike,u0prop[7:end]);
                probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm_andoyerlike,u0propnew,tspanprop,params);
                solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);  
                variableskind = 2;
            else
                Base.error("error: ",e);     
            end
        end

        tM = solM.t .+ tpropagated;
        idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
        if maximum(idx)
            solMsol  = solM(ts[idx].-tpropagated);
            resMjj = mapreduce(permutedims,vcat,solMsol.u)
            if variableskind == 1
                resM[idx,:] = resMjj;
            else
                resMtemporary = zeros(size(resMjj)[1],12)
                for jj = 1 : size(resMjj)[1]
                    anvarjj = andoyerlike2andoyer(resMjj[jj,1:7]);
                    resMtemporary[jj,1:6] = anvarjj
                    resMtemporary[jj,7:12] = resMjj[jj,8:end];
                end
                resM[idx,:] = resMtemporary;
            end
        end
        # update
        mjd2000 = params[6] +  solM.t[end]/24/3600;
        resC   =  mapreduce(permutedims,vcat,solM.u);
        if variableskind == 1
            u0prop = resC[end,:];
        else
            anvarjj = andoyerlike2andoyer(resC[end,1:7]);
            u0prop[1:6] = anvarjj
            u0prop[7:end] = resC[end,8:end];
        end
        params[6] = mjd2000;
        tpropagated = tM[end];
        tspanprop = [0.0 ssign*minimum([86400,abs(ts[end]-tpropagated)])];
    end

    tV = tV .+ t0;
    
   
    # handle sol
    for jj = 1 : size(resM)[1]
        resM[jj,4]=mod(resM[jj,4],2.0*pi);
        resM[jj,5]=mod(resM[jj,5],2.0*pi);
        resM[jj,6]=mod(resM[jj,6],2.0*pi);
    end

    # transformation of variables
    res = zeros(size(resM,1),40);
    if isrefframechanged
        RMT = deepcopy(transpose(RM));
    end

    for jj = 1 : size(resM)[1]
        anvar = resM[jj,1:6];
        euanglesjj,omVjj = andoyer2euler(anvar,IV);
        sadovjj = resM[jj,1:6];
        sadovjj[4] = mod(sadovjj[4]-pi/2,2*pi); 
        IVmod = IV;
        anvarlike = andoyer2andoyerlike(anvar);
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RMT*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RMT*omVjj;
            IVmod = RMT*IV;
            IVmod[1] = abs(IVmod[1]);
            IVmod[2] = abs(IVmod[2]);
            IVmod[3] = abs(IVmod[3]);
            anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        end
        energy = (IVmod[1]*omVjj[1]^2+IVmod[2]*omVjj[2]^2+IVmod[3]*omVjj[3]^2)/2;
        Jd = anvar[2]^2/2/energy;
        kepjj = equi2kep(resM[jj,7:12],3);
        res[jj,:] = append!(sadovjj,omVjj,euanglesjj,anvar,resM[jj,7:12],NaN,energy,Jd,anvarlike,kepjj);
    end
    
    return tV, res;
end

# previous version

function prop_semianalytical_st_girf_triax_sadovlike(u0,paramsint,timevector,isrefframechanged,RM,IV)

    t0 = timevector[1];
    ts = timevector.-t0;
    tspan = [0.0 ts[end]];

    if paramsint[5]["includeDragTorque"] || paramsint[5]["includeDragAcc"]
        params = deepcopy(paramsint);
        u0prop = deepcopy(u0);
        tpropagated = 0.0;
        ssign = ts[end]/abs(ts[end]);
        tspanprop  = [0, ssign*params[4]["timeupdate"]];
        resM = zeros(size(ts)[1],13);
        mcheck = true;
        tV = deepcopy(ts);

        meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[8:13],params[6],params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
                params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
                params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
                params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));
        params[4]["averageoverM_attitude"] = meanOverMDrag_A;
        params[4]["averageoverM_orbit"] =meanOverMDrag_O  ;
        
        while abs(tpropagated)<abs(ts[end]) && mcheck
            println("time propagated ", tpropagated/86400.0)
            if paramsint[5]["checkchaoticity"]
                cb    = VectorContinuousCallback(conditionjumppatch,affectjumppatch!,2)
                probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0prop,tspanprop,params);          
                solM  = solve(probM,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10,callback = cb);
            else
                cb    = ContinuousCallback(conditionjump,affectjump!)
                probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0prop,tspanprop,params);          
                solM  = solve(probM,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10,callback = cb);
            end
            
            zetaend = solM[end][1]

            if paramsint[5]["checkchaoticity"] && (zetaend-1.0/(1.0+0.9999/(paramsint[3]["kconstant"])))<1e-12
                mcheck = false;
            end

            if mcheck
                tM = solM.t .+ tpropagated;
                idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
                if maximum(idx)
                    solMsol  = solM(ts[idx].-tpropagated);
                    resMjj = mapreduce(permutedims,vcat,solMsol.u)
                    resM[idx,:] = resMjj;
                end
                # adjustment of J5
                resC   =  mapreduce(permutedims,vcat,solM.u);
                u0prop = resC[end,:];
                mjd2000 = params[6] +  solM.t[end]/24/3600;

                savarlike0updated = u0prop[1:7];

                if abs(savarlike0updated[7])<1e-15 && abs(solM.t[end])<abs(params[4]["timeupdate"])
                    # println("QUI!")
                    cdelta0 = savarlike0updated[3]/savarlike0updated[2];
                    ss = sign(savarlike0updated[7]-resC[end-2,7]);
                    savarlike0updated[5] = mod(savarlike0updated[5]-ss*cdelta0*2.0*pi,2.0*pi);
                    u0prop[1:7] = savarlike0updated;
                end

                meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[8:13],mjd2000,params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
                    params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
                    params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
                    params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));
                params[4]["averageoverM_attitude"] = meanOverMDrag_A;
                params[4]["averageoverM_orbit"] = meanOverMDrag_O ;
                params[4]["timeupdate"]   = pi*sqrt(u0prop[8]^3.0/params[1]["muPlanet"]);
                params[6] = mjd2000;
                tspanprop = [0.0 ssign*params[4]["timeupdate"]];
                tpropagated = tM[end];
            else
                tM = solM.t + tpropagated;
                idx = ts.>= tpropagated .&& ts.<=tM[end];
                solMsol  = solM(ts[idx]-tpropagated);
                resMjj = mapreduce(permutedims,vcat,solMsol.u)
                resM[idx,:] = resMjj;

                idxfinal = length(ts);
                for jj = 1 : size(ts)[1]
                    if ts[jj]>tpropagated+tM[end]
                        idxfinal = jj;
                        break;
                    end
                end

                resM[idxfinal,:] =mapreduce(permutedims,vcat,solM.u)[end,:];
                tV = append!(ts[1:idxfinal-1],tpropagated+tM[end]);
            end    
        end
    else
        u0prop = deepcopy(u0);
        tpropagated = 0.0;
        tspanprop  = deepcopy(tspan);
        resM = zeros(size(ts)[1],13);
        mcheck = true;
        tV = deepcopy(ts);
        while abs(tpropagated)<abs(ts[end]) && mcheck
            println("time to propagate ", tspanprop)
            if paramsint[5]["checkchaoticity"]
                cb    = VectorContinuousCallback(conditionjumppatch,affectjumppatch!,2)
                probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0prop,tspanprop,paramsint);          
                solM  = solve(probM,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10,callback = cb);
            else
                cb    = ContinuousCallback(conditionjump,affectjump!)
                probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0prop,tspanprop,paramsint);          
                solM  = solve(probM,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10,callback = cb);
            end
            
            zetaend = solM[end][1]

            if paramsint[5]["checkchaoticity"] && (zetaend-1.0/(1.0+0.9999/(paramsint[3]["kconstant"])))<1e-12
                mcheck = false;
            end

            
            if mcheck
                tM = solM.t .+ tpropagated;
                idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
                if maximum(idx)
                    solMsol  = solM(ts[idx].-tpropagated);
                    resMjj = mapreduce(permutedims,vcat,solMsol.u)
                    resM[idx,:] = resMjj;
                end
                # adjustment of J5
                resC   =  mapreduce(permutedims,vcat,solM.u);
                u0prop = resC[end,:];
                savarlike0updated = u0prop[1:7];
                cdelta0 = savarlike0updated[3]/savarlike0updated[2];
                ss = sign(savarlike0updated[7]-resC[end-2,7]);
                savarlike0updated[5] = mod(savarlike0updated[5]-ss*cdelta0*2.0*pi,2.0*pi);
                u0prop[1:7] = savarlike0updated;
                tpropagated = tM[end];
                tspanprop = [0.0 tspan[2]-tpropagated];
                
            else
                tM = solM.t + tpropagated;
                idx = ts.>= tpropagated .&& ts.<=tM[end];
                solMsol  = solM(ts[idx]-tpropagated);
                resMjj = mapreduce(permutedims,vcat,solMsol.u)
                resM[idx,:] = resMjj;

                idxfinal = length(ts);
                for jj = 1 : size(ts)[1]
                    if ts[jj]>tpropagated+tM[end]
                        idxfinal = jj;
                        break;
                    end
                end

                resM[idxfinal,:] =mapreduce(permutedims,vcat,solM.u)[end,:];
                tV = append!(ts[1:idxfinal-1],tpropagated+tM[end]);
            end    
        end
    end
    tV = tV .+t0;

    for jj = 1 : size(resM)[1]
        resM[jj,4]=mod(resM[jj,4],2.0*pi);
    end

    # transformation of variables
    res = zeros(size(resM,1),40);
    if isrefframechanged
        RMT = deepcopy(transpose(RM));
    end

  
    savarM = zeros(size(resM,1),7);
    for jj = 1 : size(resM)[1]
        savarjj,skmjj = sadovlike2sadov(resM[jj,1:7],IV);
        savarM[jj,:]  = append!(savarjj,skmjj);
        resM[jj,5]    = savarM[jj,5]+savarM[jj,3]/savarM[jj,2]*savarM[jj,6];
    end

    for jj = 1 : size(resM)[1]
        euanglesjj,omVjj = sadov2euler(savarM[jj,1:6],savarM[jj,7],IV);
        IVmod = IV;
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RMT*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RMT*omVjj;
            IVmod = RMT*IV;
            IVmod[1] = abs(IVmod[1]);
            IVmod[2] = abs(IVmod[2]);
            IVmod[3] = abs(IVmod[3]);
        end
        anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        kepjj = equi2kep(resM[jj,8:13],3);
        energy = (IV[1]*omVjj[1]^2+IV[2]*omVjj[2]^2+IV[3]*omVjj[3]^2)/2;
        Jd = anvar[2]^2/2/energy;
        res[jj,:] = append!(savarM[jj,1:6],omVjj,euanglesjj,anvar,resM[jj,8:13],savarM[jj,7],energy,Jd,resM[jj,1:7],kepjj);
    end

    
    
    return tV, res;
end

function prop_semianalytical_st_girf_axysimm_andoyerlike(u0,paramsint,timevector,isrefframechanged,RM,IV)

    t0 = timevector[1];
    ts = timevector.-t0;
    tspan = [0.0 ts[end]];

    if paramsint[5]["includeDragTorque"] || paramsint[5]["includeDragAcc"]
        params = deepcopy(paramsint);
        u0prop = deepcopy(u0);
        tpropagated = 0.0;
        ssign = ts[end]/abs(ts[end]);
        tspanprop  = [0, ssign*params[4]["timeupdate"]];
        resM = zeros(size(ts)[1],13);
        mcheck = true;
        tV = deepcopy(ts);

        meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[8:13],params[6],params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
            params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
            params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
            params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));
        params[4]["averageoverM_attitude"] = meanOverMDrag_A;
        params[4]["averageoverM_orbit"] =meanOverMDrag_O  ;
        
        while abs(tpropagated)<abs(ts[end]) && mcheck
            println("time to propagate ", tpropagated/86400)
            
            cb = ContinuousCallback(conditionjump,affectjump!)
            probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm_andoyerlike,u0prop,tspanprop,params);          
            solM  = solve(probM,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10,callback = cb);
            
            
            tM = solM.t .+ tpropagated;
            idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
            if maximum(idx)
                solMsol  = solM(ts[idx].-tpropagated);
                resMjj = mapreduce(permutedims,vcat,solMsol.u)
                resM[idx,:] = resMjj;
            end
            
            # adjustment of J5
            u0prop = mapreduce(permutedims,vcat,solM.u)[end,:];
            mjd2000 = params[6] +  solM.t[end]/24/3600;
            anvarlike0updated = u0prop[1:7];
            if abs(anvarlike0updated[7])<1e-15 && abs(solM.t[end])>abs(params[4]["timeupdate"])
                psih0updated = mod(atan(anvarlike0updated[7],anvarlike0updated[6]),2.0*pi);
                cdelta0 = anvarlike0updated[3]/anvarlike0updated[2];
                anvarlike0updated[5] = mod(anvarlike0updated[5]-cdelta0*psih0updated,2.0*pi);
                u0prop[1:7] = anvarlike0updated;
            end

            meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[8:13],mjd2000,params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
                params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
                params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
                params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));
            params[4]["averageoverM_attitude"] = meanOverMDrag_A;
            params[4]["averageoverM_orbit"] = meanOverMDrag_O ;
            params[4]["timeupdate"]   = pi*sqrt(u0prop[8]^3.0/params[1]["muPlanet"]);
            params[6] = mjd2000;
            tspanprop = [0.0 ssign*params[4]["timeupdate"]];
            tpropagated = tM[end];
        end
    else
        tpropagated = 0.0;
        tspanprop  = deepcopy(tspan);
        resM = zeros(size(ts)[1],13);
        u0prop = deepcopy(u0);
        tV = deepcopy(ts);
        while abs(tpropagated)<abs(ts[end])
            cb = ContinuousCallback(conditionjump,affectjump!)
            probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm_andoyerlike,u0prop,tspanprop,paramsint);          
            solM  = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback = cb);
            tM = solM.t .+ tpropagated;
            idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
            solMsol  = solM(ts[idx].-tpropagated);
            resMjj = mapreduce(permutedims,vcat,solMsol.u)
            resM[idx,:] = resMjj;
            
            # adjustment of J5
            u0prop = mapreduce(permutedims,vcat,solM.u)[end,:];
            anvarlike0updated = u0prop[1:7];
            psih0updated = mod(atan(anvarlike0updated[7],anvarlike0updated[6]),2.0*pi);
            cdelta0 = anvarlike0updated[3]/anvarlike0updated[2];
            anvarlike0updated[5] = mod(anvarlike0updated[5]-cdelta0*psih0updated,2.0*pi);
            u0prop[1:7] = anvarlike0updated;
            tpropagated = tM[end];
            tspanprop = [0.0 tspan[2]-tpropagated];
        end
    end
    tV = tV .+ t0;

   
    # handle sol
    anvarM = zeros(size(resM,1),6);
    for jj = 1 : size(resM)[1]
        anvarjj = andoyerlike2andoyer(resM[jj,1:7]);
        anvarM[jj,:] = anvarjj
        anvarM[jj,4] = mod(anvarM[jj,4],2.0*pi)
        anvarM[jj,5] = mod(anvarM[jj,5],2.0*pi)
        anvarM[jj,6] = mod(anvarM[jj,6],2.0*pi)
    end

    # transformation of variables
    res = zeros(size(resM,1),40);
    if isrefframechanged
        RMT = deepcopy(transpose(RM));
    end

    for jj = 1 : size(resM)[1]
        anvar = anvarM[jj,:];
        euanglesjj,omVjj = andoyer2euler(anvar,IV);
        sadovjj = anvarM[jj,1:6];
        sadovjj[4] = mod(sadovjj[4]-pi/2,2*pi); 
        IVmod = IV;
        if isrefframechanged
            Ri2b = euler2ibrotmat(euanglesjj);
            Ri2b = RMT*Ri2b;   
            euanglesjj = rotmatib2euler(Ri2b);
            omVjj      = RMT*omVjj;
            IVmod = RMT*IV;
            IVmod[1] = abs(IVmod[1]);
            IVmod[2] = abs(IVmod[2]);
            IVmod[3] = abs(IVmod[3]);
            anvar = euler2andoyer(euanglesjj,omVjj,IVmod);
        end
        energy = (IVmod[1]*omVjj[1]^2+IVmod[2]*omVjj[2]^2+IVmod[3]*omVjj[3]^2)/2;
        Jd = anvar[2]^2/2/energy;
        kepjj = equi2kep(resM[7:13],3);
        res[jj,:] = append!(sadovjj,omVjj,euanglesjj,anvar,resM[jj,8:13],NaN,energy,Jd,resM[jj,1:7],kepjj);
    end
    
    return tV, res;
end

# osculating2mean (or viceversa)

function transformationci_keplerian(savar0,euangles0,omV0,equi0,params,osculating2mean)    
    paramsmod = deepcopy(params);
    mu     = paramsmod[1]["muPlanet"];
    paramsmod  = deepcopy(params);

    paramsmod[5]["includeGravityTorque"] = false;
    paramsmod[5]["includeMagneticTorque"] = false;
    paramsmod[5]["includeSrpTorque"] = false;
    paramsmod[5]["includeDragTorque"] = false;

    includeZonalHarmsAcc = paramsmod[5]["includeZonalHarmsAcc"];
    includeThirdBodyAcc = paramsmod[5]["includeThirdBodyAcc"];
    includeSunGravityAcc = paramsmod[5]["includeSunGravityAcc"];
    includeSrpAcc = paramsmod[5]["includeSrpAcc"];
    includeDragAcc = paramsmod[5]["includeDragAcc"];

    paramsmod[5]["includeEclipsesEffectsOnOrbit"] = false;

    paramsmodosculating = deepcopy(paramsmod);
    paramsmodosculating[4] = paramsmodosculating[4]["atmosphericinfo"];

    orbitpert = false;
    
    if includeZonalHarmsAcc  || includeSrpAcc || includeThirdBodyAcc || includeSunGravityAcc || includeDragAcc
        orbitpert = true;
    end
    if (orbitpert == false) 
        equi0sol = equi0;
    else 
        # PERIODS COMPUTATION
        if includeDragAcc 
            Tref = paramsmod[4]["timeupdate"]-1e-3;
        else
            nma    = sqrt(mu/equi0[1]^3.0);
            Tma    = 2*pi/nma;
            Tref = Tma
        end
        tstep = Tref/5000.0;
   
        # TIME SETTINGS
        tspan = [0.0 Tref];
        ts  =  collect(0.0:tstep:Tref);

        # oscularing propagation initial conditions
        quat0 = euler2quat(euangles0);
        u0 = deepcopy(equi0);
        u0 = append!(u0,quat0,omV0);
        prob = ODEProblem(orbitAttitudeEv_quat,u0,tspan,paramsmodosculating);
        sol  = solve(prob,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);
        sol  = sol(ts);
        resQ = mapreduce(permutedims,vcat,sol.u);
        tV = sol.t;
        osculatingvar = zeros(size(resQ)[1],7);     
        
        for jj=1:size(tV)[1]
            osculatingvar[jj,:] = append!(resQ[jj,1:5],sin(resQ[jj,6]),cos(resQ[jj,6]));
        end
        
    
        var0sol = deepcopy(equi0);
        timeInt = Inf*ones(7);
        while maximum(abs.(timeInt[1:5]))>1e-5 && maximum(abs.(timeInt[6:7]))>1e-2
            # println("orbitcorrection ", maximum(abs.(timeInt[1:5]))," ",maximum(abs.(timeInt[6:7])))
            # semi-analytical propagation initial conditions
            u0mean    = deepcopy(savar0);
            u0mean    = append!(u0mean,var0sol);
            # println(u0mean," ",tspan)
            probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0mean,tspan,paramsmod);
            solM  = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);
            solM  = solM(ts);
            resM  = mapreduce(permutedims,vcat,solM.u);
            meanvar = zeros(size(resM)[1],7);
            for jj=1:size(tV)[1]
                meanvar[jj,:] = append!(resM[jj,7:11],sin(resM[jj,12]),cos(resM[jj,12]));
            end

            # mean value of the difference of the outputs
            dXV = osculatingvar-meanvar;    
            timeInt = zeros(7,1);
            for jj=1:size(tV)[1]-1
                timeInt = timeInt + (dXV[jj,:]+dXV[jj+1,:])*(tV[jj+1]-tV[jj])/2.0;
            end
            timeInt = timeInt/(tV[end]-tV[1]);     
            timeInt[abs.(timeInt).<1e-15].=0.0; 

            vupdt = append!(var0sol[1:5],sin(var0sol[6]),cos(var0sol[6])) + osculating2mean*timeInt;
            var0sol[1:5] = vupdt[1:5];
            var0sol[6]   = mod(atan(vupdt[6],vupdt[7]),2.0*pi);
        end
        equi0sol = var0sol;
    end

    return equi0sol;
end

function transformationci_triaxial(kconstant,savar0,euangles0,omV0,equi0,params,mjd2000,osculating2mean)
    paramsmod  = deepcopy(params);
    paramsmod[5]["includeZonalHarmsAcc"] = false;
    paramsmod[5]["includeThirdBodyAcc"]  = false;
    paramsmod[5]["includeSunGravityAcc"] = false;
    paramsmod[5]["includeDragAcc"]       = false;
    paramsmod[5]["includeSrpAcc"]        = false;

    includeGravityTorque  = paramsmod[5]["includeGravityTorque"];
    includeMagneticTorque = paramsmod[5]["includeMagneticTorque"];
    includeSrpTorque      = paramsmod[5]["includeSrpTorque"];
    includeDragTorque     = paramsmod[5]["includeDragTorque"];

    # m  = (1-savar0[1])/savar0[1]*kconstant;
    # KK0 = ellF(pi/2.0,m);
    # PP0 = ellP(-kconstant,pi/2.0,m);
    # IV = paramsmod[3]["MomentsOfInertia"]
    # skm0 = savar0[1]
    # npsil  = -pi/2.0/sqrt(1+kconstant)*sqrt(skm0)*savar0[2]/IV[1]/IV[3]*(IV[3]-IV[1])/KK0;
    # npsig  = savar0[2]/IV[1]/IV[3]*(IV[3]-IV[1])/KK0*(PP0-(1-skm0)*KK0) + savar0[2]/IV[1]/IV[3]*(IV[3]*(1.0-skm0)+IV[1]*skm0);

    # println([npsil,npsig,sqrt(paramsmod[1]["muPlanet"]/(equi0[1]^3.0))])
    attitudepert = false;

    if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque
        attitudepert = true;
    end 

    if attitudepert == false
        savar0sol = savar0;
    else
          
        # geneting function gravity and magnetic torque ---> true longitude needed
        TL0 = meanlong2truelong(equi0[6],equi0[2],equi0[3]);
        equi0mod = append!(equi0[1:5],TL0)

        # geneting function srp torque ---> eccentric longitude, sun position and shadow info needed
        if includeSrpTorque 
            mjd2000 = paramsmod[6];
            rSunV   = -paramsmod[2]["ecliptic2inertial"]*celestialbodiesephemeris_position(paramsmod[1]["centralBodyIDX"],mjd2000);
            includeEclipsesEffects = paramsmod[5]["includeEclipsesEffectsOnAttitude"]
            if includeEclipsesEffects
                passageinshadowoccurs,TLin,TLout = eclipseinoutV2(rSunV,equi0mod[1],equi0mod[2],equi0mod[3],equi0mod[4],equi0mod[5],paramsmod[1]["rPlanet"]);
                if passageinshadowoccurs
                    ELin  = truelong2ecclong(TLin,equi0mod[2],equi0mod[3]);
                    ELout = truelong2ecclong(TLout,equi0mod[2],equi0mod[3]);
                else
                    ELin = NaN;
                    ELout = NaN;
                end
            else
                passageinshadowoccurs = false;
                ELin = NaN;
                ELout = NaN;
            end
            EL = truelong2ecclong(equi0mod[6],equi0mod[2],equi0mod[3]);
        end

        if includeDragTorque
            inertial2equatorialRM = transpose(paramsmod[2]["equatorial2inertial"]);
        end

        if osculating2mean == 1           
            # println(savar0)
            W2 = [0.0,0.0,0.0,0.0,0.0,0.0];
            if includeGravityTorque
                WG2 = getWG2(paramsmod[3]["MomentsOfInertia"],kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[1]["muPlanet"]);
                W2 = W2 + WG2;
                # println("gravity WG2",WG2)
            end
            if includeMagneticTorque 
                WM2 = getWM2(paramsmod[3]["MomentsOfInertia"],kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[1]["muM"],paramsmod[3]["intrinsicMagneticMoment"][1],paramsmod[3]["intrinsicMagneticMoment"][2],
                        paramsmod[3]["intrinsicMagneticMoment"][3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"]);
                # println("magnetic WM2",WM2)
                W2  = W2 + WM2;
            end   
            if includeSrpTorque || includeDragTorque 
                m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS = termsWtbWtsForSrpDragW2(savar0[1],kconstant);
            end 
            
            if includeSrpTorque
                Wsrp2 = getWSRP2(kconstant,savar0[2:6],savar0[1],m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSunV,paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout);
                W2  = W2+Wsrp2; 
                # println("srp Wsrp2 ",Wsrp2)
            end
            if includeDragTorque 
                Wdrag2 = getWDRAG2(kconstant,savar0[2:6],savar0[1],m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi0mod,paramsmod[3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000)
                W2  = W2+Wdrag2; 
                # println("drag W2", Wdrag2)
            end
            
            savar1 = savar0 - W2;         
            # println(savar1)     

            W2bis = [0.0,0.0,0.0,0.0,0.0,0.0];
            if includeSrpTorque
                Wsrp2bis = getWSRP2bis(kconstant,savar1[2:6],savar1[1],equi0mod[2],equi0mod[3],EL,rSunV,paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout)
                W2bis  = W2bis + Wsrp2bis;
                # println("srp W2bis ",Wsrp2bis)
            end
            if includeDragTorque # && paramsmod[4]["atmosphericinfo"]["atmosphericmodel"]==1
                Wdrag2bis = getWDRAG2bis(kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[3],paramsmod[4]["averageoverM_attitude"],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],
                 paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000);
                W2bis  = W2bis + Wdrag2bis;
                # println("drag W2bis ",Wdrag2bis)
            end
            savar1 = savar1 - W2bis;
            # println(savar1)

            W3 = [0.0,0.0,0.0,0.0,0.0,0.0];
            if includeGravityTorque
                WG3 = getWG3(paramsmod[3]["MomentsOfInertia"],kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[1]["muPlanet"]);
                W3 = W3 + WG3;
                # println("magnetic WG3",WG3)
            end
            if includeMagneticTorque
                WM3 = getWM3(paramsmod[3]["MomentsOfInertia"],kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[1]["muM"],paramsmod[3]["intrinsicMagneticMoment"][1],
                        paramsmod[3]["intrinsicMagneticMoment"][2],paramsmod[3]["intrinsicMagneticMoment"][3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],paramsmod[1]["muPlanet"]);
                W3  = W3 + WM3;
                # println("magnetic WM3",WM3)
            end
            if includeSrpTorque
                Wsrp3 = getWSRP3(kconstant,savar1[2:6],savar1[1],equi0mod[1],equi0mod[2],equi0mod[3],EL,rSunV,paramsmod[1]["muPlanet"],paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout);
                W3  = W3 + Wsrp3;
                # println("srp Wsrp3 ",Wsrp3)
            end
            if includeDragTorque #&& paramsmod[4]["atmosphericinfo"]["atmosphericmodel"]==1
                Wdrag3 = getWDRAG3(kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[3],paramsmod[4]["averageoverM_attitude"],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000);
                W3  = W3 + Wdrag3;
                # println("drag W3 ",Wdrag3)
            end
            savar0sol = savar1 - W3;
           
        else
            # geneting function gravity and magnetic torque ---> true longitude needed
            W3 = [0.0,0.0,0.0,0.0,0.0,0.0];
            if includeGravityTorque
                WG3 = getWG3(paramsmod[3]["MomentsOfInertia"],kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[1]["muPlanet"]);
                W3 = W3 + WG3;
            end
            if includeMagneticTorque
                WM3 = getWM3(paramsmod[3]["MomentsOfInertia"],kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[1]["muM"],paramsmod[3]["intrinsicMagneticMoment"][1],
                        paramsmod[3]["intrinsicMagneticMoment"][2],paramsmod[3]["intrinsicMagneticMoment"][3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],paramsmod[1]["muPlanet"]);
                W3  = W3 + WM3;
            end
            if includeSrpTorque
                Wsrp3 = getWSRP3(kconstant,savar0[2:6],savar0[1],equi0mod[1],equi0mod[2],equi0mod[3],EL,rSunV,paramsmod[1]["muPlanet"],paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout);
                W3  = W3 + Wsrp3;
            end
            if includeDragTorque
                Wdrag3 = getWDRAG3(kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[3],paramsmod[4]["averageoverM_attitude"],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000);
                W3  = W3 + Wdrag3;
            end
            savar1 = savar0 + W3;

            W2bis = [0.0,0.0,0.0,0.0,0.0,0.0];
            if includeSrpTorque
                Wsrp2bis = getWSRP2bis(kconstant,savar1[2:6],savar1[1],equi0mod[2],equi0mod[3],EL,rSunV,paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout)
                W2bis  = W2bis + Wsrp2bis;
            end
            if includeDragTorque
                Wdrag2bis = getWDRAG2bis(kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[3],paramsmod[4]["averageoverM_attitude"],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],
                paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000);
                W2bis  = W2bis + Wdrag2bis;
            end
            savar1 = savar1 + W2bis;

            W2 = [0.0,0.0,0.0,0.0,0.0,0.0];
            if includeGravityTorque
                WG2 = getWG2(paramsmod[3]["MomentsOfInertia"],kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[1]["muPlanet"]);
                W2 = W2 + WG2;
            end
            if includeMagneticTorque
                WM2 = getWM2(paramsmod[3]["MomentsOfInertia"],kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[1]["muM"],paramsmod[3]["intrinsicMagneticMoment"][1],paramsmod[3]["intrinsicMagneticMoment"][2],
                        paramsmod[3]["intrinsicMagneticMoment"][3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"]);
                W2  = W2 + WM2;
            end    
            if includeSrpTorque || includeDragTorque 
                m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS = termsWtbWtsForSrpDragW2(savar1[1],kconstant);
            end 
            if includeSrpTorque
                Wsrp2 = getWSRP2(kconstant,savar1[2:6],savar1[1],m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSunV,paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout);
                W2  = W2 + Wsrp2; 
            end
            if includeDragTorque
                Wdrag2 = getWDRAG2(kconstant,savar1[2:6],savar1[1],m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi0mod,paramsmod[3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000)
                W2  = W2+Wdrag2; 
            end
            savar0sol = savar1 + W2;
        end
    end       
    return savar0sol;
end

function transformationci_triaxial_numerical(kconstant,savar0,euangles0,omV0,equi0,params,osculating2mean)

    IV     = params[3]["MomentsOfInertia"];
    mu     = params[1]["muPlanet"];
   
    includeGravityTorque  = params[5]["includeGravityTorque"];
    includeMagneticTorque = params[5]["includeMagneticTorque"];
    includeSrpTorque      = params[5]["includeSrpTorque"];
    includeDragTorque     = params[5]["includeDragTorque"];

    paramsmod = deepcopy(params);
    paramsmodosculating = deepcopy(params);
    paramsmodosculating[4] = paramsmodosculating[4]["atmosphericinfo"];

    attitudepert = false;

    if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque 
        attitudepert = true;
    end
    if (attitudepert == false) 
        savar0sol = savar0;
    else 
        # PERIODS COMPUTATION
        skm0   = savar0[1];
        m0     = kconstant*(1-skm0)/skm0;
        KK0    = ellF(pi/2,m0);
        PP0    = ellP(-kconstant,pi/2,m0);
        npsil  = -pi/2.0/sqrt(1+kconstant)*sqrt(skm0)*savar0[2]/IV[1]/IV[3]*(IV[3]-IV[1])/KK0;
        npsig  = savar0[2]/IV[1]/IV[3]*(IV[3]-IV[1])/KK0*(PP0-(1-skm0)*KK0) + savar0[2]/IV[1]/IV[3]*(IV[3]*(1.0-skm0)+IV[1]*skm0);
        nma    = sqrt(mu/equi0[1]^3.0);
      
        Tma    = 2*pi/nma;
        if m0==0.0
            Tref   = 2*pi/(npsil+npsig);
        else
            Tpsil  = 2*pi/abs(npsil);
            Tpsig  = 2*pi/abs(npsig);
            Tref   = max(Tpsil,Tpsig);
        end
        Tref = max(Tref,Tma)
        tstep = Tref/1000.0;
   
        # TIME SETTINGS
        tspan = [0.0 Tref];
        ts  =  collect(0.0:tstep:Tref);

        # oscularing propagation initial conditions
        quat0 = euler2quat(euangles0);
        u0 = deepcopy(equi0);
        u0 = append!(u0,quat0,omV0);
        prob = ODEProblem(orbitAttitudeEv_quat,u0,tspan,paramsmodosculating);
        sol  = solve(prob,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10);
        sol  = sol(ts);
        resQ = mapreduce(permutedims,vcat,sol.u);
        tV = sol.t;
        osculatingvar = zeros(size(resQ)[1],17);
        
        for jj=1:size(tV)[1]
            omVjj = resQ[jj,11:13];
            euanglesjj = quat2euler(resQ[jj,7:10]);
            savarjj,skmjj = euler2sadov(euanglesjj,omVjj,IV);
            osculatingvar[jj,:] = append!([skmjj],savarjj[2:3],cos(savarjj[4]),sin(savarjj[4]),cos(savarjj[5]),sin(savarjj[5]),cos(savarjj[6]),sin(savarjj[6]),omVjj[1:2],cos(euanglesjj[1]),sin(euanglesjj[1]),cos(euanglesjj[2]),sin(euanglesjj[2]),cos(euanglesjj[3]),sin(euanglesjj[3]));
        end
    
        # semi-analytical propagation initial conditions
        savar0var  = deepcopy(savar0)
        m0val = deepcopy(m0);
        maxtimeInt = Inf;
        count = 0;
        # while maxtimeInt>1e-3
            count = count+1;
            # println(count)
            u0mean    = deepcopy(savar0var);
            u0mean    = append!(u0mean,equi0);
            probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0mean,tspan,paramsmod);
            solM  = solve(probM,Feagin14(),reltol=1e-14,abstol=1e-14,maxiters=1e+10);


            if includeDragTorque
                if paramsmod[4]["atmosphericinfo"]["strategy"]==2
                    println("correction")
                    ttvect = collect(range(solM.t[1],solM.t[end],1000));
                    resttvect   =  mapreduce(permutedims,vcat,solM(ttvect).u);
                    resNum = getdragcorrectionattitudeeq(1,ttvect,resttvect,paramsmod[6],
                    paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                    paramsmod[1]["rPlanet"],params[1]["muPlanet"],paramsmod[1]["planetRotation"],paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],
                    paramsmod[1]["zonalharmonicscoff"][1],paramsmod[5]["includeDragTorque"],paramsmod[5]["includeDragAcc"],
                    paramsmod[4]["atmosphericinfo"],paramsmod[4]["correction"],transpose(paramsmod[2]["equatorial2inertial"]),paramsmod[3])
                    paramsmod[4]["numericallyaveragedterms"] = resNum;
                    probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0mean,tspan,paramsmod);
                    solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);
                end
            end
            solM  = solM(ts);
            resM  = mapreduce(permutedims,vcat,solM.u);
            meanvar = zeros(size(resM)[1],17);
            for jj=1:size(tV)[1]
                euangles,omV = sadov2euler(resM[jj,2:6],resM[jj,1],IV)
                meanvar[jj,:] = append!(resM[jj,1:3],cos(resM[jj,4]),sin(resM[jj,4]),cos(resM[jj,5]),sin(resM[jj,5]),cos(resM[jj,6]),sin(resM[jj,6]),omV[1:2],cos(euangles[1]),sin(euangles[1]),cos(euangles[2]),sin(euangles[2]),cos(euangles[3]),sin(euangles[3]));
            end
            # mean value of the difference of the outputs
            dXV = osculatingvar-meanvar;    
            timeInt = zeros(17,1);
            for jj=1:size(tV)[1]-1
                timeInt = timeInt + (dXV[jj,:]+dXV[jj+1,:])*(tV[jj+1]-tV[jj])/2.0;
            end
            timeInt = timeInt/(tV[end]-tV[1]);     
            timeInt[abs.(timeInt).<1e-15].=0.0; 

            # computation of the initial mean values
            var0sol = [savar0[1],savar0[2],savar0[3],cos(savar0[4]),sin(savar0[4]),cos(savar0[5]),sin(savar0[5]),cos(savar0[6]),sin(savar0[6])] + osculating2mean*timeInt[1:9]; 
            savar0var = append!(var0sol[1:3],mod(atan(var0sol[5],var0sol[4]),2.0*pi),mod(atan(var0sol[7],var0sol[6]),2.0*pi),mod(atan(var0sol[9],var0sol[8]),2.0*pi));
            maxtimeInt = maximum(abs.(timeInt[1:9]));
            println(maxtimeInt)
            # if m0 = 0.0 psil and psig are not well defined => used of euler angles and components of angular velocity to determine the initial conditions
            if m0val==0.0
                eucomupdated = [omV0[1],omV0[2],cos(euangles0[1]),sin(euangles0[1]),cos(euangles0[2]),sin(euangles0[3]),cos(euangles0[3]),sin(euangles0[3])] + osculating2mean*timeInt[10:17];
                pupdated     = eucomupdated[1];
                qupdated     = eucomupdated[2];
                phiupdated   = mod(atan(eucomupdated[4],eucomupdated[3]),2.0*pi);
                thetaupdated = mod(atan(eucomupdated[6],eucomupdated[5]),2.0*pi);
                psiupdated   = mod(atan(eucomupdated[8],eucomupdated[7]),2.0*pi);
                # smk = ((IV[1]*pupdated)^2.0+(IV[2]*qupdated)^2.0/(1+kconstant))/(savar0var[2]^2.0);
                # if smk<1e-16
                #     smk = 0.0
                # end
                # skm = 1-smk;
                smk = 0.0
                skm = 1.0

                cdelta = (savar0var[3]/savar0var[2]);
                sdelta = sqrt(1.0-cdelta^2.0);
                chA    = cos(savar0var[6]);
                shA    = sin(savar0var[6]);
                RdhT    = zeros(3,3)
                RdhT[1,1] =  chA
                RdhT[1,2] = -cdelta*shA
                RdhT[1,3] =  sdelta*shA
                RdhT[2,1] =  shA
                RdhT[2,2] =  cdelta*chA
                RdhT[2,3] = -sdelta*chA
                RdhT[3,2] =  sdelta
                RdhT[3,3] =  cdelta
                Ri2b  =  euler2ibrotmat([phiupdated,thetaupdated,psiupdated]);
                RgAslA = Ri2b*RdhT;
                if smk == 0.0
                    lA = 0.0;
                    clA = 1.0;
                    slA = 0.0;
                    cgA = RgAslA[1,1]*clA-RgAslA[2,1]*slA
                    sgA = RgAslA[1,2]*clA-RgAslA[2,2]*slA
                    gA  = mod(atan(sgA,cgA),2.0*pi);
                else
                    ssigma =  sqrt(((IV[1]*pupdated)^2.0+(IV[2]*qupdated)^2.0)/(savar0var[2]^2.0));
                    slA    = RgAslA[1,3]/ssigma
                    clA    = RgAslA[2,3]/ssigma
                    sgA    = RgAslA[3,1]/ssigma
                    cgA    = -RgAslA[3,3]/ssigma
                    gA  = mod(atan(sgA,cgA),2.0*pi);
                    lA  = mod(atan(slA,clA),2.0*pi);
                end 
                
                lambda = mod(atan(-cos(lA),sqrt(1+kconstant)*sin(lA)),2.0*pi);
                mval = smk/skm*kconstant; 
                Fl = ellF(lambda,mval)
                KK = ellF(pi/2.0,mval)
                Pl = ellP(-kconstant,lambda,mval)
                PP = ellP(-kconstant,pi/2.0,mval)
                psil = pi/2.0*Fl/KK;
                psig = gA + sqrt(1+kconstant)/sqrt(skm)*(Pl-2.0/pi*PP*psil);
                savar0var[1] = skm;
                savar0var[4] = psil; 
                savar0var[5] = psig;   
                maxtimeInt = maximum(abs.(append!(timeInt[2:3],timeInt[8:9],timeInt[10:17])));
            end
            m0val = kconstant*(1-savar0var[1])/savar0var[1];
            # display(plot(tV,[osculatingvar[:,2],meanvar[:,2]]));
            # display(plot(tV,[osculatingvar[:,3],meanvar[:,3]]));
            # display(plot(tV,[osculatingvar[:,9],meanvar[:,9]]));
            # display(plot(tV,[osculatingvar[:,8],meanvar[:,8]]));

        end
        savar0sol = savar0var;
    # end

    return savar0sol;
end

function transformationci_triaxial_sadovlike(kconstant,savar0,euangles0,omV0,equi0,params,mjd2000,osculating2mean)
    paramsmod  = deepcopy(params);
    paramsmod[5]["includeZonalHarmsAcc"] = false;
    paramsmod[5]["includeThirdBodyAcc"]  = false;
    paramsmod[5]["includeSunGravityAcc"] = false;
    paramsmod[5]["includeDragAcc"]       = false;
    paramsmod[5]["includeSrpAcc"]        = false;

    includeGravityTorque  = paramsmod[5]["includeGravityTorque"];
    includeMagneticTorque = paramsmod[5]["includeMagneticTorque"];
    includeSrpTorque      = paramsmod[5]["includeSrpTorque"];
    includeDragTorque     = paramsmod[5]["includeDragTorque"];
    inertial2equatorialRM = transpose(paramsmod[2]["equatorial2inertial"]);

    attitudepert = false;

    if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque
        attitudepert = true;
    end 

    if attitudepert == false
        savar0sol = savar0;
    # elseif includeDragTorque 
    #     savar0sol = transformationci_triaxial_sadovlike_numerical(kconstant,savar0,euangles0,omV0,equi0,paramsmod,osculating2mean);
    else

        # geneting function gravity and magnetic torque ---> true longitude needed
        TL0 = meanlong2truelong(equi0[6],equi0[2],equi0[3]);
        equi0mod = append!(equi0[1:5],TL0)

        # geneting function srp torque ---> eccentric longitude, sun position and shadow info needed
        if includeSrpTorque 
            mjd2000 = paramsmod[6];
            rSunV   = -paramsmod[2]["ecliptic2inertial"]*celestialbodiesephemeris_position(paramsmod[1]["centralBodyIDX"],mjd2000);
            includeEclipsesEffects = paramsmod[5]["includeEclipsesEffects"]
            if includeEclipsesEffects 
                passageinshadowoccurs,TLin,TLout = eclipseinoutV2(rSunV,equi0mod[1],equi0mod[2],equi0mod[3],equi0mod[4],equi0mod[5],paramsmod[1]["rPlanet"]);
                if passageinshadowoccurs
                    ELin  = truelong2ecclong(TLin,equi0mod[2],equi0mod[3]);
                    ELout = truelong2ecclong(TLout,equi0mod[2],equi0mod[3]);
                else
                    ELin = NaN;
                    ELout = NaN;
                end
            else
                passageinshadowoccurs = false;
                ELin = NaN;
                ELout = NaN;
            end
            EL = truelong2ecclong(equi0mod[6],equi0mod[2],equi0mod[3]);
        end
        

        if osculating2mean == 1           
            
            sadovlike0 = sadov2sadovlike(savar0[2:6],savar0[1]);
            W2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
            if includeGravityTorque
                WG2 = getWG2sadovlike(paramsmod[3]["MomentsOfInertia"],kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[1]["muPlanet"]);
                W2 = W2 + WG2;
            end
            if includeMagneticTorque
                WM2 = getWM2sadovlike(paramsmod[3]["MomentsOfInertia"],kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[1]["muM"],paramsmod[3]["intrinsicMagneticMoment"][1],paramsmod[3]["intrinsicMagneticMoment"][2],
                        paramsmod[3]["intrinsicMagneticMoment"][3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"]);
                W2  = W2 + WM2;
            end    
            if includeSrpTorque || includeDragTorque 
                m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS = termsWtbWtsForSrpDragW2(savar0[1],kconstant);
            end 
            if includeSrpTorque
                Wsrp2 = getWSRP2sadovlike(kconstant,savar0[2:6],savar0[1],m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSunV,paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout);
                W2  = W2+Wsrp2; 
            end
            if includeDragTorque
                Wdrag2 =  getWDRAG2sadovlike(kconstant,savar0[2:6],savar0[1],m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi0mod,paramsmod[3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                    paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],
                    inertial2equatorialRM,mjd2000);
                W2  = W2+Wdrag2; 
            end
            sadovlike1  = sadovlike0 - W2;
            savar1,skm1 = sadovlike2sadov(sadovlike1,paramsmod[3]["MomentsOfInertia"]);
            savar1[1] = skm1;
                            

            W2bis = [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
            if includeSrpTorque
                Wsrp2bis = getWSRP2bissadovlike(kconstant,savar1[2:6],savar1[1],equi0mod[2],equi0mod[3],EL,rSunV,paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout)
                W2bis  = W2bis + Wsrp2bis;
            end
            if includeDragTorque
                Wdrag2bis = getWDRAG2bissadovlike(kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[3],paramsmod[4]["averageoverM_attitude"],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],
                paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000);
                W2bis  = W2bis + Wdrag2bis;
            end
            sadovlike1 = sadovlike1 - W2bis;
            savar1,skm1 = sadovlike2sadov(sadovlike1,paramsmod[3]["MomentsOfInertia"]);
            savar1[1] = skm1;

            W3 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
            if includeGravityTorque
                WG3 = getWG3sadovlike(paramsmod[3]["MomentsOfInertia"],kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[1]["muPlanet"]);
                W3 = W3 + WG3;
            end
            if includeMagneticTorque
                WM3 = getWM3sadovlike(paramsmod[3]["MomentsOfInertia"],kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[1]["muM"],paramsmod[3]["intrinsicMagneticMoment"][1],
                        paramsmod[3]["intrinsicMagneticMoment"][2],paramsmod[3]["intrinsicMagneticMoment"][3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],paramsmod[1]["muPlanet"]);
                W3  = W3 + WM3;
            end
            if includeSrpTorque
                Wsrp3 = getWSRP3sadovlike(kconstant,savar1[2:6],savar1[1],equi0mod[1],equi0mod[2],equi0mod[3],EL,rSunV,paramsmod[1]["muPlanet"],paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout);
                W3  = W3 + Wsrp3;
            end
            if includeDragTorque
                Wdrag3 = getWDRAG3sadovlike(kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[3],paramsmod[4]["averageoverM_attitude"],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000);
                W3  = W3 + Wdrag3;
            end
            sadovlike0sol = sadovlike1 - W3;
            savar0sol,skmsol = sadovlike2sadov(sadovlike0sol,paramsmod[3]["MomentsOfInertia"]);
            savar0sol[1] = skmsol;

        else
            sadovlike0 = sadov2sadovlike(savar0[2:6],savar0[1]);
            W3 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
            if includeGravityTorque
                WG3 = getWG3sadovlike(paramsmod[3]["MomentsOfInertia"],kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[1]["muPlanet"]);
                W3 = W3 + WG3;
            end
            if includeMagneticTorque
                WM3 = getWM3sadovlike(paramsmod[3]["MomentsOfInertia"],kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[1]["muM"],paramsmod[3]["intrinsicMagneticMoment"][1],
                        paramsmod[3]["intrinsicMagneticMoment"][2],paramsmod[3]["intrinsicMagneticMoment"][3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],paramsmod[1]["muPlanet"]);
                W3  = W3 + WM3;
            end
            if includeSrpTorque
                Wsrp3 = getWSRP3sadovlike(kconstant,savar0[2:6],savar0[1],equi0mod[1],equi0mod[2],equi0mod[3],EL,rSunV,paramsmod[1]["muPlanet"],paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout);
                W3  = W3 + Wsrp3;
            end
            if includeDragTorque
                Wdrag3 = getWDRAG3sadovlike(kconstant,savar0[2:6],savar0[1],equi0mod,paramsmod[3],paramsmod[4]["averageoverM_attitude"],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000);
                W3  = W3 + Wdrag3;
            end
            sadovlike1  = sadovlike0 + W3;
            savar1,skm1 = sadovlike2sadov(sadovlike1,paramsmod[3]["MomentsOfInertia"]);
            savar1[1] = skm1;

            W2bis = [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
            if includeSrpTorque
                Wsrp2bis = getWSRP2bissadovlike(kconstant,savar1[2:6],savar1[1],equi0mod[2],equi0mod[3],EL,rSunV,paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout)
                W2bis  = W2bis + Wsrp2bis;
            end
            if includeDragTorque
                Wdrag2bis = getWDRAG2bissadovlike(kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[3],paramsmod[4]["averageoverM_attitude"],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],
                paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],inertial2equatorialRM,mjd2000);
                W2bis  = W2bis + Wdrag2bis;
            end
            sadovlike1 = sadovlike1 + W2bis;
            savar1,skm1 = sadovlike2sadov(sadovlike1,paramsmod[3]["MomentsOfInertia"]);
            savar1[1] = skm1

            W2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0];
            if includeGravityTorque
                WG2 = getWG2sadovlike(paramsmod[3]["MomentsOfInertia"],kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[1]["muPlanet"]);
                W2 = W2 + WG2;
            end
            if includeMagneticTorque
                WM2 = getWM2sadovlike(paramsmod[3]["MomentsOfInertia"],kconstant,savar1[2:6],savar1[1],equi0mod,paramsmod[1]["muM"],paramsmod[3]["intrinsicMagneticMoment"][1],paramsmod[3]["intrinsicMagneticMoment"][2],
                        paramsmod[3]["intrinsicMagneticMoment"][3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"]);
                W2  = W2 + WM2;
            end    
            if includeSrpTorque || includeDragTorque 
                m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS = termsWtbWtsForSrpDragW2(savar0[1],kconstant);
            end 
            if includeSrpTorque
                Wsrp2 = getWSRP2sadovlike(kconstant,savar1[2:6],savar1[1],m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,EL,rSunV,paramsmod[3],includeEclipsesEffects,passageinshadowoccurs,ELin,ELout);
                W2  = W2+Wsrp2; 
            end
            if includeDragTorque
                Wdrag2 =  getWDRAG2sadovlike(kconstant,savar1[2:6],savar1[1],m,KK,EE,PP,znC1,znC2,znC3,WTb,WTS,equi0mod,paramsmod[3],paramsmod[2]["cpa1"],paramsmod[2]["cpa2"],paramsmod[2]["cpa3"],
                    paramsmod[1]["planet_flat"],paramsmod[1]["oblatenesscoeff"],paramsmod[1]["muPlanet"],paramsmod[4]["atmosphericinfo"],paramsmod[1]["planetRotation"],paramsmod[1]["rPlanet"],
                    inertial2equatorialRM,mjd2000);
                W2  = W2+Wdrag2; 
            end
            sadovlike0sol = sadovlike1 + W2;
            savar0sol,skmsol = sadovlike2sadov(sadovlike0sol,paramsmod[3]["MomentsOfInertia"]);
            savar0sol[1] = skmsol;
        end
        
    end    
    return savar0sol;
end

function transformationci_triaxial_sadovlike_numerical(kconstant,savar0,euangles0,omV0,equi0,params,osculating2mean)

    paramsmod = deepcopy(params);
    IV     = paramsmod[3]["MomentsOfInertia"];
    mu     = paramsmod[1]["muPlanet"];
    paramsmod  = deepcopy(params);
    paramsmod[5]["includeZonalHarmsAcc"] = false;
    paramsmod[5]["includeSrpAcc"] = false;

    includeGravityTorque  = paramsmod[5]["includeGravityTorque"];
    includeMagneticTorque = paramsmod[5]["includeMagneticTorque"];
    includeSrpTorque      = paramsmod[5]["includeSrpTorque"];
    includeDragTorque     = paramsmod[5]["includeDragTorque"];

    paramsmodosculating = deepcopy(paramsmod);
    paramsmodosculating[4] = paramsmodosculating[4]["atmosphericinfo"];
    
    attitudepert = false;
    if includeGravityTorque || includeMagneticTorque || includeSrpTorque || includeDragTorque 
        attitudepert = true;
    end

    if attitudepert == false 
        savar0sol = savar0;
    else 
        savarlike0 =  sadov2sadovlike(savar0[2:6],savar0[1]);
       
        skm0 = savarlike0[1];
        m0             = kconstant*(1-skm0)/skm0;
        KK0            = ATTITUDE.ellF(pi/2,m0);
        PP0            = ATTITUDE.ellP(-kconstant,pi/2,m0);
        npsil          = -pi/2.0/sqrt(1+kconstant)*sqrt(skm0)*savar0[2]/IV[1]/IV[3]*(IV[3]-IV[1])/KK0;
        npsig          = savar0[2]/IV[1]/IV[3]*(IV[3]-IV[1])/KK0*(PP0-(1-skm0)*KK0) + savar0[2]/IV[1]/IV[3]*(IV[3]*(1.0-skm0)+IV[1]*skm0);
        nma            = sqrt(mu/equi0[1]^3.0);
        
        Tma    = 2*pi/nma;
        if m0==0.0
            Tref   = 2*pi/(npsil+npsig);
        else
            Tpsil  = 2*pi/abs(npsil);
            Tpsig  = 2*pi/abs(npsig);
            Tref   = max(Tpsil,Tpsig);
        end
        Tref = max(Tref,Tma)
        tstep = Tref/1000.0;

        # TIME SETTINGS
        tspan = [0.0 Tref];
        ts    = collect(0.0:tstep:Tref);

        # oscularing propagation initial conditions
        quat0 = euler2quat(euangles0);
        u0 = deepcopy(equi0);
        u0 = append!(u0,quat0,omV0);
        prob = ODEProblem(orbitAttitudeEv_quat,u0,tspan,paramsmodosculating);
        sol  = solve(prob,Vern9(),reltol=5e-14,abstol=5e-14,maxiters=1e+10);
        sol  = sol(ts);
        resQ = mapreduce(permutedims,vcat,sol.u);
        tV = sol.t;
        osculatingvar = zeros(size(resQ)[1],16);
        for jj=1:size(tV)[1]
            omVjj = resQ[jj,11:13];
            euanglesjj = quat2euler(resQ[jj,7:10]);
            savarjj,skmjj = euler2sadov(euanglesjj,omVjj,IV);
            savarlikejj = sadov2sadovlike(savarjj,skmjj);
            osculatingvar[jj,:] = append!(savarlikejj[1:3],cos(savarlikejj[4]),sin(savarlikejj[4]),savarjj[5],savarlikejj[6:7],omVjj[1:2],cos(euanglesjj[1]),sin(euanglesjj[1]),cos(euanglesjj[2]),sin(euanglesjj[2]),cos(euanglesjj[3]),sin(euanglesjj[3]));
        end
        
        # J5 var computation
        if m0 != 0.0
            # J5
            corrangle = 0.0;
            psigEV = deepcopy(osculatingvar[:,6]);
            for jj=2:size(tV)[1]
                if sign(npsig)*psigEV[jj]<sign(npsig)*psigEV[jj-1] 
                    corrangle = corrangle+sign(npsig)*2.0*pi;
                end
                osculatingvar[jj,6] = osculatingvar[jj,6] + corrangle;
            end
            deltaJ5O = osculatingvar[:,3]./osculatingvar[:,2].*mod.(atan.(osculatingvar[:,8],osculatingvar[:,7]),2.0*pi)
            osculatingvar[:,6] =   osculatingvar[:,6] + deltaJ5O;
        end

        # averaging prop
        u0mean = deepcopy(savarlike0)
        u0mean    = append!(u0mean,equi0);
        tpropagated = 0.0;
        tspanprop  = deepcopy(tspan);
        resM = zeros(size(ts)[1],13);
        while tpropagated<ts[end]
            cb = ContinuousCallback(conditionjump,affectjump!)
            probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial_sadovlike,u0mean,tspanprop,paramsmod);          
            solM  = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback = cb);
            tM = solM.t .+ tpropagated;
            idx = ts.>= tpropagated .&& ts.<=tM[end];
            solMsol  = solM(ts[idx].-tpropagated);
            resMjj = mapreduce(permutedims,vcat,solMsol.u)
            resM[idx,:] = resMjj;
            
            u0mean = mapreduce(permutedims,vcat,solM.u)[end,:];
            savarlike0updated = u0mean[1:7];
            psih0updated = mod(atan(savarlike0updated[7],savarlike0updated[6]),2.0*pi);
            cdelta0 = savarlike0updated[3]/savarlike0updated[2];
            savarlike0updated[5] = mod(savarlike0updated[5]-cdelta0*psih0updated,2.0*pi);
            u0mean[1:7] = savarlike0updated;
            tpropagated = tM[end];
            tspanprop = [0.0 tspan[2]-tpropagated];
        end

        meanvar = zeros(size(resM)[1],16);
        for jj=1:size(tV)[1]
            savarjj,skmjj = sadovlike2sadov(resM[jj,1:7],IV);
            eujj,omjj = sadov2euler(savarjj,skmjj,IV);
            meanvar[jj,:] = append!(resM[jj,1:3],cos(resM[jj,4]),sin(resM[jj,4]),savarjj[5],resM[jj,6:7],omjj[1:2],cos(eujj[1]),sin(eujj[1]),cos(eujj[2]),sin(eujj[2]),cos(eujj[3]),sin(eujj[3])); 
        end
        
        # J5 var computation
        if m0 != 0.0
            corrangle = 0.0;
            psigEV = deepcopy(meanvar[:,6]);
            for jj=2:size(tV)[1] #-1
                if sign(npsig)*psigEV[jj]<sign(npsig)*psigEV[jj-1]
                    corrangle = corrangle+sign(npsig)*2.0*pi;
                end
                meanvar[jj,6] = meanvar[jj,6] + corrangle;
            end
            deltaJ5 = meanvar[:,3]./meanvar[:,2].*mod.(atan.(meanvar[:,8],meanvar[:,7]),2.0*pi);
            meanvar[:,6] = meanvar[:,6] + deltaJ5;
        end

        # mean of the difference
        dXV = osculatingvar-meanvar;            
        timeInt = zeros(16,1);
        for jj=1:size(tV)[1]-1
            timeInt = timeInt + (dXV[jj,:]+dXV[jj+1,:])/2.0;
        end
        timeInt = timeInt/(length(tV)-1)  
        
        # mean values        
        var0new = [savarlike0[1],savarlike0[2],savarlike0[3],cos(savarlike0[4]),sin(savarlike0[4]),savarlike0[5],savarlike0[6],savarlike0[7]] + osculating2mean*timeInt[1:8]; 
        savarlike0new = append!(var0new[1:3],mod(atan(var0new[5],var0new[4]),2.0*pi),var0new[6:8]);
        savar0sol,skmmean0 = sadovlike2sadov(savarlike0new,IV);    
        savar0sol[1]   = skmmean0;

        # if m0 = 0.0 psil and psig are not well defined => used of euler angles and components of angular velocity to define the initial conditions
        if m0 == 0.0
            eucomupdated = [omV0[1],omV0[2],cos(euangles0[1]),sin(euangles0[1]),cos(euangles0[2]),sin(euangles0[2]),cos(euangles0[3]),sin(euangles0[3])] + osculating2mean*timeInt[9:16];
            pupdated     = eucomupdated[1];
            qupdated     = eucomupdated[2];
            phiupdated   = mod(atan(eucomupdated[4],eucomupdated[3]),2.0*pi);
            thetaupdated = mod(atan(eucomupdated[6],eucomupdated[5]),2.0*pi);
            psiupdated   = mod(atan(eucomupdated[8],eucomupdated[7]),2.0*pi);
            # smk = ((IV[1]*pupdated)^2.0+(IV[2]*qupdated)^2.0/(1+kconstant))/(savar0sol[2]^2.0);
            # if smk<1e-16
            #     smk = 0.0
            # end
            smk = 0.0
            skm = 1-smk;

            cdelta = (savar0sol[3]/savar0sol[2]);
            sdelta = sqrt(1.0-cdelta^2.0);
            chA    = cos(savar0sol[6]);
            shA    = sin(savar0sol[6]);
            RdhT    = zeros(3,3)
            RdhT[1,1] =  chA
            RdhT[1,2] = -cdelta*shA
            RdhT[1,3] =  sdelta*shA
            RdhT[2,1] =  shA
            RdhT[2,2] =  cdelta*chA
            RdhT[2,3] = -sdelta*chA
            RdhT[3,2] =  sdelta
            RdhT[3,3] =  cdelta
            Ri2b  =  euler2ibrotmat([phiupdated,thetaupdated,psiupdated]);
            RgAslA = Ri2b*RdhT;
            if smk == 0.0
                lA = 0.0;
                clA = 1.0;
                slA = 0.0;
                cgA = RgAslA[1,1]*clA-RgAslA[2,1]*slA
                sgA = RgAslA[1,2]*clA-RgAslA[2,2]*slA
                gA  = mod(atan(sgA,cgA),2.0*pi);
            else
                ssigma =  sqrt(((IV[1]*pupdated)^2.0+(IV[2]*qupdated)^2.0)/(savar0sol[2]^2.0));
                slA    = RgAslA[1,3]/ssigma
                clA    = RgAslA[2,3]/ssigma
                sgA    = RgAslA[3,1]/ssigma
                cgA    = -RgAslA[3,3]/ssigma
                gA  = mod(atan(sgA,cgA),2.0*pi);
                lA  = mod(atan(slA,clA),2.0*pi);
            end 
            
            lambda = mod(atan(-cos(lA),sqrt(1+kconstant)*sin(lA)),2.0*pi);
            mval = smk/skm*kconstant; 
            Fl = ellF(lambda,mval)
            KK = ellF(pi/2.0,mval)
            Pl = ellP(-kconstant,lambda,mval)
            PP = ellP(-kconstant,pi/2.0,mval)
            psil = pi/2.0*Fl/KK;
            psig = gA + sqrt(1+kconstant)/sqrt(skm)*(Pl-2.0/pi*PP*psil);
            savar0sol[1] = skm;
            savar0sol[4] = psil; 
            savar0sol[5] = psig;
        end

        # println("here")
        
        # display(plot(tV,[osculatingvar[:,1],meanvar[:,1]],label="skm"))
        # display(plot(tV,[osculatingvar[:,2],meanvar[:,2]],label="Jg"))
        # display(plot(tV,[osculatingvar[:,3],meanvar[:,3]],label="Jh"))
        # display(plot(tV,[osculatingvar[:,4],meanvar[:,4]],label="cpsil"))
        # display(plot(tV,[osculatingvar[:,5],meanvar[:,5]],label="spsil"))
        # display(plot(tV,[osculatingvar[:,6],meanvar[:,6]],label="J5"))
        # display(plot(tV,[osculatingvar[:,7],meanvar[:,7]],label="J6"))
        # display(plot(tV,[osculatingvar[:,8],meanvar[:,8]],label="J7"))
        # display(plot(tV,[osculatingvar[:,9],meanvar[:,9]],label="p"))
        # display(plot(tV,[osculatingvar[:,10],meanvar[:,10]],label="q"))
        # display(plot(tV,[osculatingvar[:,11],meanvar[:,11]],label="cphi"))
        # display(plot(tV,[osculatingvar[:,12],meanvar[:,12]],label="sphi"))
        # display(plot(tV,[osculatingvar[:,13],meanvar[:,13]],label="ctheta"))
        # display(plot(tV,[osculatingvar[:,14],meanvar[:,14]],label="stheta"))
        # display(plot(tV,[osculatingvar[:,15],meanvar[:,15]],label="cpsi"))
        # display(plot(tV,[osculatingvar[:,16],meanvar[:,16]],label="spsi"))
        # display(plot(tV,[dXV[:,1]],label="deltaskm"))
        # display(plot(tV,[dXV[:,2]],label="deltaJg"))
        # display(plot(tV,[dXV[:,3]],label="deltaJh"))
        # display(plot(tV,[dXV[:,4]],label="deltacpsil"))
        # display(plot(tV,[dXV[:,5]],label="deltaspsil"))
        # display(plot(tV,[dXV[:,6]],label="deltaJ5"))
        # display(plot(tV,[dXV[:,7]],label="deltaJ6"))
        # display(plot(tV,[dXV[:,8]],label="deltaJ7"))
        # display(plot(tV,[dXV[:,9]],label="deltap"))
        # display(plot(tV,[dXV[:,10]],label="deltaq"))
        # display(plot(tV,[dXV[:,11]],label="deltacphi"))
        # display(plot(tV,[dXV[:,12]],label="deltasphi"))
        # display(plot(tV,[dXV[:,13]],label="deltactheta"))
        # display(plot(tV,[dXV[:,14]],label="deltastheta"))
        # display(plot(tV,[dXV[:,15]],label="deltacpsi"))
        # display(plot(tV,[dXV[:,16]],label="deltaspsi"))
       
    end

    # println(savar0sol)

    return savar0sol;
end

function transformationci_axysimmetric(anvar0,euangles0,omV0,equi0,params,osculating2mean)
    paramsmod  = deepcopy(params);
    paramsmod[5]["includeZonalHarmsAcc"] = false;
    paramsmod[5]["includeThirdBodyAcc"]  = false;
    paramsmod[5]["includeSunGravityAcc"] = false;
    paramsmod[5]["includeDragAcc"]       = false;
    paramsmod[5]["includeSrpAcc"]        = false;
    IV     = paramsmod[3]["MomentsOfInertia"];
    mu     = paramsmod[1]["muPlanet"];

    includeMagneticTorque = paramsmod[5]["includeMagneticTorque"];
    includeSrpTorque      = paramsmod[5]["includeSrpTorque"];
    includeDragTorque     = paramsmod[5]["includeDragTorque"];

    paramsmodosculating = deepcopy(paramsmod);
    paramsmodosculating[4] = paramsmodosculating[4]["atmosphericinfo"];

    
    attitudepert = false;   
    if includeMagneticTorque || includeSrpTorque || includeDragTorque 
        attitudepert = true;
    end
    if (attitudepert == false && orbitpert == false) 
        anvar0sol = anvar0;
    else 
        # PERIODS COMPUTATION
        ngA  = anvar0[2]/IV[1]
        nma    = sqrt(mu/equi0[1]^3.0);
        TgA  = 2*pi/abs(ngA);
        Tma    = 2*pi/nma;
        Tref   = max(TgA,Tma);
        tstep = Tref/1000.0;
        
        # TIME SETTINGS
        tspan = [0.0 Tref];
        ts    = collect(0.0:tstep:Tref);

        # oscularing propagation initial conditions
        quat0 = euler2quat(euangles0);
        u0 = deepcopy(equi0);
        u0 = append!(u0,quat0,omV0);
        prob = ODEProblem(orbitAttitudeEv_quat,u0,tspan,paramsmodosculating);
        sol  = solve(prob,Feagin14(),reltol=1e-12,abstol=1e-12,maxiters=1e+10);
        sol  = sol(ts);
        resQ = mapreduce(permutedims,vcat,sol.u);
        tV = sol.t;
        osculatingvar = zeros(size(resQ)[1],9);
        for jj=1:size(tV)[1]
            omVjj = resQ[jj,11:13];
            euanglesjj = quat2euler(resQ[jj,7:10]);
            anvar = euler2andoyer(euanglesjj,omVjj,IV);
            osculatingvar[jj,:] = append!(anvar[1:3],cos.(anvar[4:6]),sin.(anvar[4:6]));
        end


        # semi-analytical propagation
        u0mean    = deepcopy(anvar0);
        u0mean    = append!(u0mean,equi0);
        probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm,u0mean,tspan,paramsmod);
        solM  = solve(probM,Feagin14(),reltol=1e-12,abstol=1e-12,maxiters=1e+10);
        solM  = solM(ts);
        resM  = mapreduce(permutedims,vcat,solM.u);
        meanvar = zeros(size(resM)[1],9);
        for jj=1:size(tV)[1]
            meanvar[jj,:] = append!(resM[jj,1:3],cos.(resM[jj,4:6]),sin.(resM[jj,4:6]));
        end
           
        # mean value of the difference
        dXV = osculatingvar-meanvar;    
        timeInt = zeros(9,1);
        for jj=1:size(tV)[1]-1
            timeInt = timeInt + (dXV[jj,:]+dXV[jj+1,:])*(tV[jj+1]-tV[jj])/2.0;
        end
        timeInt = timeInt/(tV[end]-tV[1]);
            
        var0sol = append!(anvar0[1:3],cos.(anvar0[4:6]),sin.(anvar0[4:6]))  + osculating2mean*timeInt[1:9]; ## + [0.0, 0.0, 0.0, psil0, psig0, 0.0];
        anvar0sol = append!(var0sol[1:3],mod(atan(var0sol[7],var0sol[4]),2.0*pi),mod(atan(var0sol[8],var0sol[5]),2.0*pi),mod(atan(var0sol[9],var0sol[6]),2.0*pi));
        

        # display(plot(tV/3600.0/24.0,[osculatingvar[:,1],meanvar[:,1]]))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,2],meanvar[:,2]]))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,3],meanvar[:,3]]))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,4],meanvar[:,4]]))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,5],meanvar[:,5]]))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,6],meanvar[:,6]]))          
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,7],meanvar[:,7]]))  
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,8],meanvar[:,8]]))  
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,9],meanvar[:,9]]))  
    end
    return anvar0sol;
end

function transformationci_axysimmetric_andoyerlike(anvar0,euangles0,omV0,equi0,params,osculating2mean)

    paramsmod  = deepcopy(params);
    paramsmod[5]["includeZonalHarmsAcc"] = false;
    paramsmod[5]["includeThirdBodyAcc"]  = false;
    paramsmod[5]["includeSunGravityAcc"] = false;
    paramsmod[5]["includeDragAcc"]       = false;
    paramsmod[5]["includeSrpAcc"]        = false;
    IV     = paramsmod[3]["MomentsOfInertia"];
    mu     = paramsmod[1]["muPlanet"];

    includeMagneticTorque = paramsmod[5]["includeMagneticTorque"];
    includeSrpTorque      = paramsmod[5]["includeSrpTorque"];
    includeDragTorque     = paramsmod[5]["includeDragTorque"];

    paramsmodosculating = deepcopy(paramsmod);
    paramsmodosculating[4] = paramsmodosculating[4]["atmosphericinfo"];

    attitudepert = false;

    if includeMagneticTorque || includeSrpTorque || includeDragTorque 
        attitudepert = true;
    end

    if attitudepert == false 
        anvar0sol = anvar0;
    else 
        anvarlike0 = andoyer2andoyerlike(anvar0)

        # PERIODS COMPUTATION
        ngA  = anvarlike0[2]/IV[1]
        nma    = sqrt(mu/equi0[1]^3.0);
        TgA  = 2*pi/abs(ngA);
        Tma    = 2*pi/nma;
        Tref   = max(TgA,Tma);
        tstep = Tref/1000.0;
        
        # TIME SETTINGS
        tspan = [0.0 Tref];
        ts    = collect(0.0:tstep:Tref);

        # oscularing propagation initial conditions
        quat0 = euler2quat(euangles0);
        u0 = deepcopy(equi0);
        u0 = append!(u0,quat0,omV0);
        prob = ODEProblem(orbitAttitudeEv_quat,u0,tspan,paramsmodosculating);
        sol  = solve(prob,Feagin14(),reltol=1e-12,abstol=1e-12,maxiters=1e+10);
        sol  = sol(ts);
        resQ = mapreduce(permutedims,vcat,sol.u);
        tV = sol.t;
        osculatingvar = zeros(size(resQ)[1],9);
        for jj=1:size(tV)[1]
            omVjj = resQ[jj,11:13];
            euanglesjj = quat2euler(resQ[jj,7:10]);
            anvarjj = euler2andoyer(euanglesjj,omVjj,IV);
            anvarlikejj = andoyer2andoyerlike(anvarjj);
            osculatingvar[jj,:] = append!(anvarlikejj[1:3],cos(anvarlikejj[4]),sin(anvarlikejj[4]),cos(anvarjj[5]),sin(anvarjj[5]),anvarlikejj[6:7]);
        end

        # averaging perturbations
        u0mean = deepcopy(anvarlike0)
        u0mean    = append!(u0mean,equi0);
        tpropagated = 0.0;
        tspanprop  = deepcopy(tspan);
        resM = zeros(size(ts)[1],13);
        while tpropagated<ts[end]
            cb = ContinuousCallback(conditionjump,affectjump!)
            probM = ODEProblem(attitudeMeanEv_semianalytical_axisymm_andoyerlike,u0mean,tspanprop,paramsmod);          
            solM  = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback = cb);
            tM = solM.t .+ tpropagated;
            idx = ts.>= tpropagated .&& ts.<=tM[end];
            solMsol  = solM(ts[idx].-tpropagated);
            resMjj = mapreduce(permutedims,vcat,solMsol.u)
            resM[idx,:] = resMjj;
            
            u0mean = mapreduce(permutedims,vcat,solM.u)[end,:];
            savarlike0updated = u0mean[1:7];
            psih0updated = mod(atan(savarlike0updated[7],savarlike0updated[6]),2.0*pi);
            cdelta0 = savarlike0updated[3]/savarlike0updated[2];
            savarlike0updated[5] = mod(savarlike0updated[5]-cdelta0*psih0updated,2.0*pi);
            u0mean[1:7] = savarlike0updated;
            tpropagated = tM[end];
            tspanprop = [0.0 tspan[2]-tpropagated];
        end
        
        meanvar = zeros(size(resM)[1],9);
        for jj=1:size(tV)[1]
            anvarjj = andoyerlike2andoyer(resM[jj,1:7])
            meanvar[jj,:] = append!(resM[jj,1:3],cos(resM[jj,4]),sin(resM[jj,4]),cos(anvarjj[5]),sin(anvarjj[5]),resM[jj,6:7]);
        end        

        # meanvalue mean
        dXV = osculatingvar-meanvar;    
        timeInt = zeros(9,1);
        for jj=1:size(tV)[1]-1
            timeInt = timeInt + (dXV[jj,:]+dXV[jj+1,:])/2.0;
        end
        timeInt = timeInt/(length(tV)-1.0)
           
        var0new = append!(anvarlike0[1:3],cos(anvarlike0[4]),sin(anvarlike0[4]),cos(anvar0[5]),sin(anvar0[5]),anvarlike0[6:7]) + osculating2mean*timeInt[1:9];
        anvarlike0new = append!(var0new[1:3],mod(atan(var0new[5],var0new[4]),2.0*pi),mod(atan(var0new[7],var0new[6]),2.0*pi)+var0new[3]/var0new[2]*mod(atan(var0new[9],var0new[8]),2.0*pi),var0new[8:9]);
        anvar0sol = andoyerlike2andoyer(anvarlike0new);

        # display(plot(tV/3600.0/24.0,[osculatingvar[:,1],meanvar[:,1]],label="LA"))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,2],meanvar[:,2]],label="GA"))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,3],meanvar[:,3]],label="HA"))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,4],meanvar[:,4]],label="clA"))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,5],meanvar[:,5]],label="slA"))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,6],meanvar[:,6]],label="cgA"))          
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,7],meanvar[:,7]],label="sgA"))  
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,8],meanvar[:,8]],label="JA6"))
        # display(plot(tV/3600.0/24.0,[osculatingvar[:,9],meanvar[:,9]],label="JA7"))  

    end


    return anvar0sol;
end




# ##########################################################################
# # VERSION 2 - different handling of eclipses for srp torque --- IGNORE

# function attitudeaveragedprop_semianalytical_girf_version2_onlytriaxialO(civ,date0,params,tstep,tfinal)    

#     includeSrpTorque = get(params[5],"includeSrpTorque",false);
#     includeEclipsesEffectsOnAttitude = get(params[5],"includeEclipsesEffectsOnAttitude",false);

#     if ! (includeSrpTorque && includeEclipsesEffectsOnAttitude)
#         tV,res = attitudeaveragedprop_semianalytical_girf(civ,date0,params,tstep,tfinal);
#         return tV,res;
#     end

#     params=checksandhandleinputs(civ,date0,params,tstep,tfinal);
#     kconstant,savar,skm,euangles,omV,isrefframechanged,RM,paramsint,equiEl = handleinputs(params,civ,date0,tfinal);
#     paramsint[5] = merge(paramsint[5],Dict("srpTorqueAnalytical"=>true));
#     paramsint[4] = merge(paramsint[4],Dict("srpmeanv"=>zeros(1980*12)));

#     updatedragmean = false;
#     if paramsint[5]["includeDragTorque"] 
#         if paramsint[4]["timeupdate"]<tfinal
#             paramsint[4]["timeupdate"] = 2.0*tfinal;
#             updatedragmean = true;
#         end
#     end
    
#     # propagation #####################################################################################
#     timeVector = collect(0.0:tstep:tfinal);

#     # check chaoticity
#     mval0 = kconstant*(1-skm)/skm;
#     if (paramsint[5]["checkchaoticity"] && mval0 > 0.9999) || mval0==1.0
#         if paramsint[5]["osculatingpropagationinchaoticregion"]
#             println("chaotic region: osculating propagation")
#             tV,res = attitudeprop_quat_girf(civ,date0,params,tstep,tfinal);
#         else
#             tV = [];
#             res = [];
#             println("chaotic region")
#         end
#         return tV,res;
#     end

#     # init output
#     tV = [0.0]
#     res = zeros(1,40);

#     # init settings
#     propagationdirection = sign(tstep);
#     paramsint[5]["includeEclipsesEffectsOnAttitude"]=false;  #### the shadow is managed step by step and not directly included in the equations of motion
#     instableregion = true; 
#     tStart = 0.0;
#     mjd2000Start = deepcopy(paramsint[6]); 
#     adjustmentrequired = false;

#     # conversion equinoctial elements
#     equimean0 = transformationci_keplerian(append!([skm],savar[2:6]),euangles,omV,equiEl,paramsint,1);

#     # conversion initial sadov variables
#     satinlight = true;
#     rSunV = -paramsint[2]["ecliptic2inertial"]*celestialbodiesephemeris_position(paramsint[1]["centralBodyIDX"],paramsint[6])
#     inshadow,TLin,TLout = eclipseinoutV2(rSunV,equimean0[1],equimean0[2],equimean0[3],equimean0[4],equimean0[5],paramsint[1]["rPlanet"])
#     if inshadow
#         MLin  = truelong2meanlong(TLin,equimean0[2],equimean0[3]);
#         MLout = truelong2meanlong(TLout,equimean0[2],equimean0[3]);
#         MLcurrent = mod(equimean0[6],2.0*pi);
#         if MLin<MLout && MLcurrent>=MLin && MLcurrent<MLout 
#             satinlight = false;
#         elseif MLin>MLout && (MLcurrent>=MLin || MLcurrent<MLout)
#             satinlight = false;
#         end
#     end
#     paramsint[5]["includeSrpTorque"]=satinlight;
#     if abs(1.0-abs(savar[3]/savar[2]))<1e-13 
#         savarmean0 = transformationci_triaxial_sadovlike(kconstant,append!([skm],savar[2:6]),euangles,omV,equimean0,paramsint,1);
#     else
#         savarmean0 = transformationci_triaxial(kconstant,append!([skm],savar[2:6]),euangles,omV,equimean0,paramsint,1);
#     end
   
#     # propagation
#     while abs(tStart)<abs(timeVector[end]) && instableregion
        
#         # define time interval
#         meanmotion  = sqrt(paramsint[1]["muPlanet"]/equimean0[1]^3)
#         currentOPeriod = 2.0*pi/meanmotion;
#         ts = deepcopy(timeVector);
#         ts = ts[abs.(ts).<=abs(tStart+propagationdirection*currentOPeriod) .&& abs.(ts).>=abs(tStart)];
#         ts = append!([tStart],ts,tStart+propagationdirection*currentOPeriod);

#         # eclipse
#         rSunV = -paramsint[2]["ecliptic2inertial"]*celestialbodiesephemeris_position(paramsint[1]["centralBodyIDX"],paramsint[6])
#         inshadow,TLin,TLout = eclipseinoutV2(rSunV,equimean0[1],equimean0[2],equimean0[3],equimean0[4],equimean0[5],paramsint[1]["rPlanet"])
        
#         if inshadow
#             # define light-shadow sequence
#             MLin  = truelong2meanlong(TLin,equimean0[2],equimean0[3]);
#             MLout = truelong2meanlong(TLout,equimean0[2],equimean0[3]);
#             MLcurrent = mod(equimean0[6],2.0*pi);
#             ts1,ts2,ts3,sequence = defineshadowsequence(ts,MLin,MLout,MLcurrent,meanmotion,propagationdirection);
            
#             ## propagation part 1
#             paramsint[5]["srpTorqueAnalytical"] = true;
#             paramsint[5]["includeSrpTorque"]    = sequence[1]
#             paramsint[6] = mjd2000Start;
#             tV1,res1 = launchpropagationO(ts1,savarmean0,equimean0,paramsint,isrefframechanged,RM);
            
#             # update ci
#             savarmean0    = deepcopy(append!([res1[end,25]],vec(res1[end,2:6])));
#             equimean0     = deepcopy(vec(res1[end,19:24]));
#             mjd2000Start  = mjd2000Start + (tV1[end]-tV1[1])/3600/24;

#             # check caoticity
#             if abs(tV1[end])<abs(ts1[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                 instableregion = false;
#             end
            
#             # determine whether adjustment is necessary
#             npsil,npsig = getmeanmotionpsilpsig(kconstant,savarmean0[1],savarmean0[2],paramsint[3]["MomentsOfInertia"][1],paramsint[3]["MomentsOfInertia"][3]);

#             if abs(npsil)>abs(meanmotion) || abs(npsig)>abs(meanmotion)
#                 adjustmentrequired = true;
#             end

#             # adjustment
#             if instableregion && adjustmentrequired
#                 println("adjustment")
#                 if abs(npsil)>abs(npsig)
#                     angleidx = 4;
#                     otherangleidx = 5;
#                     nangle = npsil;
#                 else
#                     angleidx = 5;
#                     otherangleidx = 4;
#                     nangle = npsig;
#                 end
                
#                 # define the time interval of the adjustment  
#                 tcurrent     = tV1[end]
#                 anglecurrent = mod(res1[end,angleidx]-res1[1,angleidx],2.0*pi);
#                 otheranglecurrent = res1[end,otherangleidx]
#                 t0   = -anglecurrent/nangle+tcurrent;
#                 t2pi = (2*pi-anglecurrent)/nangle+tcurrent;
#                 tback = propagationdirection*min(propagationdirection*t0,propagationdirection*t2pi);
#                 tforw = propagationdirection*max(propagationdirection*t0,propagationdirection*t2pi);
#                 # println("anglecurrent ", anglecurrent*180/pi," n: ",nangle, "anglelimits: ",[res1[end,angleidx],res1[1,angleidx]]*180/pi);
#                 # println("tback: ",tback,", tcurrent: ",tcurrent,", tforw: ",tforw)

#                 # backward propagation
#                 paramsint[5]["srpTorqueAnalytical"] = true;
#                 paramsint[5]["includeSrpTorque"]    = sequence[1]
#                 paramsint[6] = mjd2000Start;
#                 tVb,resb = launchpropagationO([tcurrent, tback],savarmean0,equimean0,paramsint,isrefframechanged,RM);
 
#                 # update ci
#                 savarmean0    = deepcopy(append!([resb[end,25]],vec(resb[end,2:6])));
#                 equimean0     = deepcopy(vec(resb[end,19:24]));
#                 mjd2000Start  = mjd2000Start + (tVb[end]-tVb[1])/3600/24;

#                 # correct tV1
#                 idx = size(tV1)[end];
#                 for kk = 1:size(tV1)[1]
#                     if abs(tV1[kk])>=abs(tVb[end])
#                         break;
#                     end
#                     idx = kk;
#                 end
              
#                 tV1  = append!(tV1[1:idx],tVb[end]);
#                 res1[idx+1,:] = resb[end,:];
#                 res1 = res1[1:idx+1,:];
             
#                 # define successive time step
#                 tend = propagationdirection*min(propagationdirection*tforw,propagationdirection*ts3[end]);
#                 tsA = deepcopy(append!([tVb[end]],ts[propagationdirection*ts.>=propagationdirection*tVb[end] .&& propagationdirection*ts.<=propagationdirection*tend],[tend]));                  

#                 # compute numerical mean
#                 if sequence[1]>0
#                     if tend >= abs(ts2[end])
#                         anglelimit = anglecurrent + (ts2[end]-tcurrent)*nangle;
#                     else
#                         anglelimit = anglecurrent + (tend-tcurrent)*nangle;
#                     end
#                 else
#                     anglelimit = anglecurrent + (tback-tcurrent)*nangle;
#                 end
#                 if abs(anglelimit)<1e-15
#                     anglelimit =0.0;
#                 end
#                 if anglecurrent<anglelimit
#                     anglein  = anglecurrent;
#                     angleout = anglelimit;
#                 else
#                     anglein  = anglelimit
#                     angleout = anglecurrent;
#                 end
#                 srpmeanvector = getMeanOverpsilpsig_srp(paramsint[3]["MomentsOfInertia"][1],paramsint[3]["MomentsOfInertia"][3],kconstant,resb[1,25],resb[1,2],true,angleidx,anglein,angleout);
#                 # srpmeanvector = getMeanOverpsilpsig_srpV2(paramsint[3]["MomentsOfInertia"][1],paramsint[3]["MomentsOfInertia"][3],kconstant,resb[1,25],resb[1,2],true,angleidx,anglein,angleout,mod(otheranglecurrent-npsil/npsig*anglecurrent,2.0*pi),npsil,npsig);

#                 # propagation
#                 paramsint[4]["srpmeanv"] = srpmeanvector;
#                 paramsint[5]["srpTorqueAnalytical"] = false;
#                 paramsint[5]["includeSrpTorque"]    = true;
#                 paramsint[6] = mjd2000Start;
#                 tVA,resA = launchpropagationO(tsA,savarmean0,equimean0,paramsint,isrefframechanged,RM);
                
#                 # update ci
#                 savarmean0    = deepcopy(append!([resA[end,25]],vec(resA[end,2:6])));
#                 equimean0     = deepcopy(vec(resA[end,19:24]));
#                 mjd2000Start  = mjd2000Start + (tVA[end]-tVA[1])/3600/24;
                
#                 # check caoticity
#                 if abs(tV1[end])<abs(tsA[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                     instableregion = false;
#                 end

#                 # adjust result first part
#                 tV1 = append!(tV1[1:end-1],tVA[2:end])
#                 res1 = [res1[1:end-1,:];resA[2:end,:]];
                
#                 # adjust successive propagation time intervals
#                 ts2 = ts2[propagationdirection*ts2.>propagationdirection*tend];
#                 if !isempty(ts2)
#                     ts2 = append!([tforw],ts2);
#                 else
#                     ts3 = ts3[propagationdirection*ts3.>propagationdirection*tend];
#                     if !isempty(ts3)
#                         ts3 = append!([tend],ts3);
#                     end
#                 end
#                 # update
#                 adjustmentrequired == false;
#             end

#             tVcurrent = tV1;
#             rescurrent = res1;

#             ## propagation part 2
#             if instableregion && !isempty(ts2)
#                 paramsint[5]["srpTorqueAnalytical"] = true;
#                 paramsint[5]["includeSrpTorque"]    = sequence[2]
#                 paramsint[6] = mjd2000Start;
#                 tV2,res2 = launchpropagationO(ts2,savarmean0,equimean0,paramsint,isrefframechanged,RM);
                
#                 # ci update
#                 savarmean0    = deepcopy(append!([res2[end,25]],vec(res2[end,2:6])));
#                 equimean0     = deepcopy(vec(res2[end,19:24]));
#                 mjd2000Start  = mjd2000Start + (tV2[end]-tV2[1])/3600/24;

#                 if abs(tV2[end])<abs(ts2[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                     instableregion = false;
#                 end

#                 # determine whether adjustment is necessary
#                 npsil,npsig = getmeanmotionpsilpsig(kconstant,savarmean0[1],savarmean0[2],paramsint[3]["MomentsOfInertia"][1],paramsint[3]["MomentsOfInertia"][3]);
#                 if abs(npsil)>meanmotion || abs(npsig)>meanmotion
#                     adjustmentrequired = true;
#                 end

#                 # adjustment
#                 if instableregion && adjustmentrequired
#                     println("adjustment part 2")
#                     # definition of the fastest angle
#                     if abs(npsil)>abs(npsig)
#                         angleidx = 4;
#                         nangle = npsil;
#                     else
#                         angleidx = 5;
#                         nangle = npsig;
#                     end
                    
#                     # define the time interval of the adjustment  
#                     tcurrent     = tV2[end]
#                     anglecurrent = mod(res2[end,angleidx]-res2[1,angleidx],2.0*pi);
#                     otheranglecurrent = res2[end,otherangleidx];
#                     t0   = -anglecurrent/nangle+tcurrent;
#                     t2pi = (2*pi-anglecurrent)/nangle+tcurrent;
#                     tback = propagationdirection*min(propagationdirection*t0,propagationdirection*t2pi);
#                     tforw = propagationdirection*max(propagationdirection*t0,propagationdirection*t2pi);
    
#                     # backward propagation
#                     paramsint[5]["srpTorqueAnalytical"] = true;
#                     paramsint[5]["includeSrpTorque"]    = sequence[2]
#                     paramsint[6] = mjd2000Start;
#                     tVb,resb = launchpropagationO([tcurrent, tback],savarmean0,equimean0,paramsint,isrefframechanged,RM);
                
#                     # update ci
#                     savarmean0    = deepcopy(append!([resb[end,25]],vec(resb[end,2:6])));
#                     equimean0     = deepcopy(vec(resb[end,19:24]));
#                     mjd2000Start  = mjd2000Start + (tVb[end]-tVb[1])/3600/24;
    
#                     # correct tV1
#                     idx = size(tV2)[end];
#                     for kk = 1:size(tV2)[1]
#                         if abs(tV2[kk])>=abs(tVb[end])
#                             break;
#                         end
#                         idx = kk;
#                     end
#                     tV2  = append!(tV2[1:idx],tVb[end]);
#                     res2[idx+1,:] = resb[end,:];
#                     res2 = res2[1:idx+1,:];
    
#                     # define successive time step
#                     tend = propagationdirection*min(propagationdirection*tforw,propagationdirection*ts3[end]);
#                     tsA = deepcopy(append!([tVb[end]],ts[propagationdirection*ts.>=propagationdirection*tVb[end] .&& propagationdirection*ts.<=propagationdirection*tend],[tend]));                  

#                     # compute numerical mean
#                     if sequence[1]>0
#                         anglelimit = anglecurrent + (tend-tcurrent)*nangle;
#                     else
#                         anglelimit = anglecurrent + (tback-tcurrent)*nangle;
#                     end
#                     if abs(anglelimit)<1e-15
#                         anglelimit = 0.0;
#                     end
#                     if anglecurrent<anglelimit
#                         anglein  = anglecurrent;
#                         angleout = anglelimit;
#                     else
#                         anglein  = anglelimit
#                         angleout = anglecurrent;
#                     end
#                     srpmeanvector = getMeanOverpsilpsig_srp(paramsint[3]["MomentsOfInertia"][1],paramsint[3]["MomentsOfInertia"][3],kconstant,resb[1,25],resb[1,2],true,angleidx,anglein,angleout);
#                     # srpmeanvector = getMeanOverpsilpsig_srpV2(paramsint[3]["MomentsOfInertia"][1],paramsint[3]["MomentsOfInertia"][3],kconstant,resb[1,25],resb[1,2],true,angleidx,anglein,angleout,mod(otheranglecurrent-npsil/npsig*anglecurrent,2.0*pi),npsil,npsig);

#                     # propagation
#                     paramsint[4]["srpmeanv"] = srpmeanvector;
#                     paramsint[5]["srpTorqueAnalytical"] = false;
#                     paramsint[5]["includeSrpTorque"]    = true;
#                     paramsint[6] = mjd2000Start;
#                     tVA,resA = launchpropagationO(tsA,savarmean0,equimean0,paramsint,isrefframechanged,RM);
                    
#                     # update ci
#                     savarmean0    = deepcopy(append!([resA[end,25]],vec(resA[end,2:6])));
#                     equimean0     = deepcopy(vec(resA[end,19:24]));
#                     mjd2000Start  = mjd2000Start + (tVA[end]-tVA[1])/3600/24;
                    
#                     # check caoticity
#                     if abs(tV2[end])<abs(tsA[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                         instableregion = false;
#                     end
    
#                     # adjust result first part
#                     tV2 = append!(tV2[1:end-1],tVA[2:end])
#                     res2 = [res2[1:end-1,:];resA[2:end,:]];
                    
#                     # adjust successive propagation time intervals
#                     ts3 = ts3[propagationdirection*ts3.>propagationdirection*tend];
#                     if !isempty(ts3)
#                         ts3 = append!([tforw],ts3);
#                     end

#                     # update
#                     adjustmentrequired == false;
#                 end

#                 tVcurrent = append!(tVcurrent[1:end-1],tV2[2:end]); 
#                 rescurrent = [rescurrent[1:end-1,:];res2[2:end,:]];
#             end

#             ## propagation part 3            
#             if instableregion && !isempty(ts3)
#                 paramsint[5]["srpTorqueAnalytical"] = true;
#                 paramsint[5]["includeSrpTorque"]    = sequence[3]
#                 paramsint[6] = mjd2000Start;
#                 tV3,res3 = launchpropagationO(ts3,savarmean0,equimean0,paramsint,isrefframechanged,RM);
                
#                 # update ci
#                 savarmean0    = deepcopy(append!([res3[end,25]],vec(res3[end,2:6])));
#                 equimean0     = deepcopy(vec(res3[end,19:24]));
#                 mjd2000Start  = mjd2000Start + (tV3[end]-tV3[1])/3600/24;

#                 # check caoticity
#                 if abs(tV3[end])<abs(ts3[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                     instableregion = false;
#                 end

#                 tVcurrent = append!(tVcurrent[1:end-1],tV3[2:end]); 
#                 rescurrent = [rescurrent[1:end-1,:];res3[2:end,:]];

#             end
#         else
#             # propagation
#             paramsint[5]["srpTorqueAnalytical"] = true;
#             paramsint[5]["includeSrpTorque"]    = true
#             paramsint[6] = mjd2000Start;
#             tVcurrent,rescurrent = launchpropagationO(ts,savarmean0,equimean0,paramsint,isrefframechanged,RM);

#             # update ci
#             savarmean0    = deepcopy(append!([rescurrent[end,25]],vec(rescurrent[end,2:6])));
#             equimean0     = deepcopy(vec(rescurrent[end,19:24]));
#             mjd2000Start  = mjd2000Start + (tVcurrent[end]-tVcurrent[1])/3600/24;
#             if abs(tVcurrent[end])<abs(ts[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                 instableregion = false;
#             end
#         end

#         if updatedragmean
#             P1 = equimean0[2];
#             P2 = equimean0[3];
#             Q1 = equimean0[4];
#             Q2 = equimean0[5];
#             GG = 1.0+Q1^2.0+Q2^2.0;
#             q11 = (1.0-Q1^2.0+Q2^2.0)/GG;
#             q12 = (2.0*Q1*Q2)/GG;
#             q13 = (-2.0*Q1)/GG;
#             q21 = (2.0*Q1*Q2)/GG;
#             q22 = (1.0+Q1^2.0-Q2^2.0)/GG;
#             q23 = (2.0*Q2)/GG;
#             q31 = (2.0*Q1)/GG;
#             q32 = (-2.0*Q2)/GG;
#             q33 = (1.0-Q1^2.0-Q2^2.0)/GG;  
#             ec  = sqrt(equimean0[2]^2+equimean0[3]^2);
#             if ec!=0.0
#                 Ri2p = [P2/ec P1/ec 0.0; -P1/ec P2/ec 0.0; 0.0 0.0 1.0]*RMq;
#             else
#                 ti2 = sqrt(Q1^2+Q2^2);
#                 if ti2!= 0.0
#                     Ri2p = [Q2/ti2 Q1/ti2 0.0; -Q1/ti2 Q2/ti2 0.0; 0.0 0.0 1.0]*RMq;
#                 else
#                     Ri2p = RMq;
#                 end
#             end
#             Zeqp  = Ri2p*[inertialrefframeinfo["cpa1"],inertialrefframeinfo["cpa2"],inertialrefframeinfo["cpa3"]];
#             meanOverMDrag_p = getMeanOverM_drag(sma,eta,u[14],planetsunparams["planetRotation"],Zeqp[1],Zeqp[2],Zeqp[3],mu,draginfo["AtmM"],planetsunparams["rPlanet"])/2.0/pi;
#             meanOverMDrag   = getmeaninertial(meanOverMDrag_p,transpose(Ri2p))*1000.0;
#             paramsint[4]["averageoverM_attitude"] = meanOverMDrag;
#         end

#         if tStart == 0.0
#            tV = vec(tVcurrent[2:end-1,:]);
#            res = rescurrent[2:end-1,:];
#         else
#            if size(tVcurrent)[1]>2
#                 tV = append!(tV, vec(tVcurrent[2:end-1,:]));
#                 res = [res; rescurrent[2:end-1,:]];
#            end
#         end
        
#         tStart = tStart + currentOPeriod;
  
#     end

#     # osculating propagation
#     if paramsint[5]["osculatingpropagationinchaoticregion"] && tV[end]<tfinal
#         println("chaotic region: osculating propagation")
             
#         paramsint[6] = mjd2000Start;
#         paramsint[5]["srpTorqueAnalytical"] = true;
        
#         # mean to osculating attitude
#         if abs(1.0-res[end,3]/res[end,2])<1e-13
#             savar0new = transformationci_triaxial_savarlike(kconstant,append!([res[end,25]],res[end,2:6]),res[end,10:12],res[end,7:9],res[end,19:24],paramsint,-1);
#         else
#             savar0new = transformationci_triaxial(kconstant,append!([res[end,25]],res[end,2:6]),res[end,10:12],res[end,7:9],res[end,19:24],paramsint,-1);
#         end

#         eu0new,om0new = sadov2euler(sadvar0new[2:6],savar0new[1],IV)
#         if isrefframechanged
#             RMT = Transpose(RM)
#             Ri2b = RMT*euler2ibrotmat(eu0new);
#             eu0new = rotmatib2euler(Ri2b);
#             om0new = RMT*om0new;
#             paramsint[3] = params[3];            
#         end

#         # mean to osculating orbital
#         kepEl0new = res[end,19:24];

#         paramsint["includeEclipsesEffectsOnAttitude"] = true;
#         paramsint[5]["includeSrpTorque"] = true;

#         # osculating propagation
#         tVoscu,resosc = attitudeprop_quat_girf(append!(eu0new,om0new,kepEl0new),paramrsint[6]+tV[end]/3600.0/24.0,paramsint,tstep,tfinal-tV[end]);

#         tV = append!(tV,tVoscu+tV[end]);
#         res = [res;resosc]
#     end


#     # output
#    return tV,res;
# end

# function getmeanmotionpsilpsig(kconstant,skm,Jg,A,C)

#     m   = (1-skm)*kconstant/skm;
#     KK  = ellF(pi/2.0,m)
#     PP  = ellP(-kconstant,pi/2.0,m);
#     npsil   = -pi/sqrt(1+kconstant)*sqrt(skm)*Jg/2/A/C*(C-A)/KK;
#     npsig   = Jg/A/C*(C-A)*(PP-(1-skm)*KK)/KK + Jg/A/C*(C*(1-skm)+A*skm);

#     return npsil,npsig;
# end

# function launchpropagationO(ts,savarmean0,equimean0,params,isrefframechanged,RM)
#     if abs(1.0-abs(savarmean0[3]/savarmean0[2]))<1e-10
#         savarlike = sadov2sadovlike(savarmean0[2:6],savarmean0[1]);
#         u0 = append!(savarlike,equimean0);
#         tV,res = prop_semianalytical_st_girf_triax_sadovlike(append!(savarlike,equimean0),params,ts,isrefframechanged,RM,params[3]["MomentsOfInertia"]);
#     else
#         try
#             println(" triaxial normal")
#             tV,res = prop_semianalytical_st_girf_triax(append!(savarmean0,equimean0),params,ts,isrefframechanged,RM,params[3]["MomentsOfInertia"])
#         catch e
#             if e==ErrorException("singularity")
#                 savarlike = sadov2sadovlike(savarmean0[2:6],savarmean0[1]);
#                 u0 = zeros(13,1);
#                 u0 = append!(savarlike,equimean0);
#                 tV,res = prop_semianalytical_st_girf_triax_sadovlike(append!(savarlike,equimean0),params,ts,isrefframechanged,RM,params[3]["MomentsOfInertia"]);
#             else
#                 Base.error(e)
#             end
#         end
#     end
    
#     return tV,res;
# end

# function attitudeaveragedprop_semianalytical_girf_version2_onlytriaxial(civ,date0,params,tstep,tfinal)    

#     includeSrpTorque = get(params[5],"includeSrpTorque",false);
#     includeEclipsesEffectsOnAttitude = get(params[5],"includeEclipsesEffectsOnAttitude",false);

#     if ! (includeSrpTorque && includeEclipsesEffectsOnAttitude)
#         tV,res = attitudeaveragedprop_semianalytical_girf(civ,date0,params,tstep,tfinal);
#         return tV,res;
#     end

#     params=checksandhandleinputs(civ,date0,params,tstep,tfinal);
#     kconstant,savar,skm,euangles,omV,isrefframechanged,RM,paramsint,equiEl = handleinputs(params,civ,date0,tfinal);
#     paramsint[5] = merge(paramsint[5],Dict("srpTorqueAnalytical"=>true));
#     paramsint[4] = merge(paramsint[4],Dict("srpmeanv"=>zeros(19800*12)));

#     updatedragmean = false;
#     if paramsint[5]["includeDragTorque"] 
#         if paramsint[4]["timeupdate"]<tfinal
#             paramsint[4]["timeupdate"] = 2.0*tfinal;
#             updatedragmean = true;
#         end
#     end
    
#     # propagation #####################################################################################
#     timeVector = collect(0.0:tstep:tfinal);

#     # check chaoticity
#     mval0 = kconstant*(1-skm)/skm;
#     if (paramsint[5]["checkchaoticity"] && mval0 > 0.9999) || mval0==1.0
#         if paramsint[5]["osculatingpropagationinchaoticregion"]
#             println("chaotic region: osculating propagation")
#             tV,res = attitudeprop_quat_girf(civ,date0,params,tstep,tfinal);
#         else
#             tV = [];
#             res = [];
#             println("chaotic region")
#         end
#         return tV,res;
#     end

#     # init
#     paramsint[5]["includeEclipsesEffectsOnAttitude"]=false;
#     instableregion = true; 
#     tV = [0.0]
#     res = zeros(1,40);
#     tStart = 0.0;
#     propagationdirection = sign(tstep);
#     mjd2000Start     = deepcopy(paramsint[6]); 

#     # conversion equinoctial elements
#     equimean0 = transformationci_keplerian(append!([skm],savar[2:6]),euangles,omV,equiEl,paramsint,1);

#     # conversion initial sadov variables
#     # rSunV = -paramsint[2]["ecliptic2inertial"]*celestialbodiesephemeris_position(paramsint[1]["centralBodyIDX"],paramsint[6])
#     # inshadow,TLin,TLout = eclipseinoutV2(rSunV,equimean0[1],equimean0[2],equimean0[3],equimean0[4],equimean0[5],paramsint[1]["rPlanet"])
#     # MLin  = truelong2meanlong(TLin,equimean0[2],equimean0[3]);
#     # MLout = truelong2meanlong(TLout,equimean0[2],equimean0[3]);
#     # MLcurrent = mod(equimean0[6],2.0*pi);
#     # satinlight = true;
#     # if MLin<MLout && MLcurrent>=MLin && MLcurrent<MLout 
#     #     satinlight = false;
#     # elseif MLin>MLout && (MLcurrent>=MLin || MLcurrent<MLout)
#     #     satinlight = false;
#     # end
#     satinlight = true;
#     paramsint[5]["includeSrpTorque"]=satinlight;
#     savar0 = deepcopy(savar);
#     savar0[1] = deepcopy(skm);
#     if abs(1.0-abs(savar0[3]/savar0[2]))<1e-13 
#         savarmean0 = transformationci_triaxial_sadovlike(kconstant,append!([skm],savar[2:6]),euangles,omV,equimean0,paramsint,1);
#     else
#         savarmean0 = transformationci_triaxial(kconstant,append!([skm],savar[2:6]),euangles,omV,equimean0,paramsint,1);
#     end
#     euanglesmean0,omVmean0 = sadov2euler(savarmean0[2:6],savarmean0[1],paramsint[3]["MomentsOfInertia"]);
   
#     # propagation
#     while abs(tStart)<abs(timeVector[end]) && instableregion
        
#         # define time interval
#         meanmotion  = sqrt(paramsint[1]["muPlanet"]/equimean0[1]^3)
#         currentOPeriod = 2.0*pi/meanmotion;
#         ts = deepcopy(timeVector);
#         ts = ts[abs.(ts).<=abs(tStart+propagationdirection*currentOPeriod) .&& abs.(ts).>=abs(tStart)];
#         ts = append!([tStart],ts,tStart+propagationdirection*currentOPeriod);

#         # eclipse
#         rSunV = -paramsint[2]["ecliptic2inertial"]*celestialbodiesephemeris_position(paramsint[1]["centralBodyIDX"],paramsint[6])
#         inshadow,TLin,TLout = eclipseinoutV2(rSunV,equimean0[1],equimean0[2],equimean0[3],equimean0[4],equimean0[5],paramsint[1]["rPlanet"])
#         println()
#         if inshadow
#             MLin  = truelong2meanlong(TLin,equimean0[2],equimean0[3]);
#             MLout = truelong2meanlong(TLout,equimean0[2],equimean0[3]);
#             MLcurrent = mod(equimean0[6],2.0*pi);
#             ts1,ts2,ts3,sequence = defineshadowsequence(ts,MLin,MLout,MLcurrent,meanmotion,propagationdirection);
#             println("seq ",sequence[1]," ",sequence[2]," ",sequence[3])
#             println(ts1[1],"-",ts1[end])
#             println(ts2[1],"-",ts2[end])
#             println(ts3[1],"-",ts3[end])
#             ## propagation part 1
#             # println("PARTE1")
#             tV1,res1 = launchpropagation(ts1,sequence[1],savarmean0,euanglesmean0,omVmean0,satinlight,equimean0,paramsint,kconstant,true,mjd2000Start,isrefframechanged,RM)
            
#             # update ci
#             savarmean0    = deepcopy(append!([res1[end,25]],vec(res1[end,2:6])));
#             euanglesmean0 = deepcopy(vec(res1[end,10:12]));
#             omVmean0      = deepcopy(vec(res1[end,7:9]));
#             equimean0     = deepcopy(vec(res1[end,19:24]));
#             mjd2000Start  = mjd2000Start + (tV1[end]-tV1[1])/3600/24;
#             satinlight    = sequence[1]

#             # check 
#             if abs(tV1[end])<abs(ts1[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                 tVcurrent  = tV1[1:end]; 
#                 rescurrent = res1[1:end,:];
#                 instableregion = false;
#             end

#             ## propagation part 2
#             if instableregion
#                 # println("PARTE2",sequence[2],satinlight)
#                 tV2,res2 = launchpropagation(ts2,sequence[2],savarmean0,euanglesmean0,omVmean0,satinlight,equimean0,paramsint,kconstant,true,mjd2000Start,isrefframechanged,RM)
                
#                 # update
#                 savarmean0    = deepcopy(append!([res2[end,25]],vec(res2[end,2:6])));
#                 euanglesmean0 = deepcopy(vec(res2[end,10:12]));
#                 omVmean0      = deepcopy(vec(res2[end,7:9]));
#                 equimean0     = deepcopy(vec(res2[end,19:24]));
#                 mjd2000Start  = mjd2000Start + (tV2[end]-tV2[1])/3600/24;
#                 satinlight    = sequence[2]

#                 if abs(tV2[end])<abs(ts2[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                     tVcurrent = append!(tV1[1:end-1],tV2[2:end]); 
#                     rescurrent = [res1[1:end-1,:];res2[2:end,:]];
#                     instableregion = false;
#                 end
#             end
            
#             if instableregion
#                 ## propagation part 3
#                 # println("PARTE3")
#                 tV3,res3 = launchpropagation(ts3,sequence[3],savarmean0,euanglesmean0,omVmean0,satinlight,equimean0,paramsint,kconstant,true,mjd2000Start,isrefframechanged,RM)
                
#                 # update
#                 savarmean0    = deepcopy(append!([res3[end,25]],vec(res3[end,2:6])));
#                 euanglesmean0 = deepcopy(vec(res3[end,10:12]));
#                 omVmean0      = deepcopy(vec(res3[end,7:9]));
#                 equimean0     = deepcopy(vec(res3[end,19:24]));
#                 mjd2000Start  = mjd2000Start + (tV3[end]-tV3[1])/3600/24;
#                 satinlight= sequence[3] 

#                 if abs(tV3[end])<abs(ts3[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                     instableregion = false;
#                 end
#                 tVcurrent = append!(tV1[1:end-1],tV2[2:end-1],tV3[2:end]); 
#                 rescurrent = [res1[1:end-1,:];res2[2:end-1,:];res3[2:end,:]];

#             end
#             println("STEP IN SHADOW ENDED")
#             println()
#         else
#             tVcurrent,rescurrent = launchpropagation(ts,true,savarmean0,euanglesmean0,omVmean0,satinlight,equimean0,paramsint,kconstant,true,mjd2000Start,isrefframechanged,RM)
#             savarmean0    = deepcopy(append!([rescurrent[end,25]],vec(rescurrent[end,2:6])));
#             euanglesmean0 = deepcopy(vec(rescurrent[end,10:12]));
#             omVmean0      = deepcopy(vec(rescurrent[end,7:9]));
#             equimean0     = deepcopy(vec(rescurrent[end,19:24]));
#             mjd2000Start  = mjd2000Start + (tVcurrent[end]-tVcurrent[1])/3600/24;
#             satinlight= true
#             if abs(tVcurrent[end])<abs(ts[end]) && paramsint[5]["checkchaoticity"] && kconstant*(1-savarmean0[1])/savarmean0[1]>0.9999
#                 instableregion = false;
#             end
#         end

#         if updatedragmean
#             # sma = equimean0[1];
#             # P1  = equimean0[2];
#             # P2  = equimean0[3];
#             # Q1  = equimean0[4];
#             # Q2  = equimean0[5];
#             meanOverMDrag,eV0mean,Idrag = dragaverageOverMUpdate(equimean0,paramsint[2]["cpa1"],paramsint[2]["cpa2"],paramsint[2]["cpa3"],
#             paramsint[1]["rPlanet"],paramsint[1]["muPlanet"],paramsint[1]["planetRotation"],paramsint[1]["planet_flat"],paramsint[1]["oblatenesscoeff"],
#             paramsint[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],paramsint[4]["atmosphericinfo"],paramsint[4]["couplingsimplification"],paramsint[4]["correction"]);
#             paramsint[4]["averageoverM_attitude"] = meanOverMDrag;
#             paramsint[4]["averageoverM_orbit"] = Idrag;
#             paramsint[4]["averageoverM_coupling"] = eV0mean;
#         end

#         if tStart == 0.0
#            tV = vec(tVcurrent[2:end-1,:]);
#            res = rescurrent[2:end-1,:];
#         else
#            if size(tVcurrent)[1]>2
#                 tV = append!(tV, vec(tVcurrent[2:end-1,:]));
#                 res = [res; rescurrent[2:end-1,:]];
#            end
#         end
        
#         tStart = tStart + currentOPeriod;
  
#     end

#     # propagation
#     if paramsint[5]["osculatingpropagationinchaoticregion"] && tV[end]<tfinal
#         println("chaotic region: osculating propagation")
        
#         if satinlight
#             paramsint[5]["includeSrpTorque"] = true
#         else
#             paramsint[5]["includeSrpTorque"] = false
#         end
        
#         # mean to osculating attitude
#         if abs(1.0-res[end,3]/res[end,2])<1e-13
#             savar0new = transformationci_triaxial_savarlike(kconstant,append!([res[end,25]],res[end,2:6]),res[end,10:12],res[end,7:9],res[end,19:24],paramsint,-1);
#         else
#             savar0new = transformationci_triaxial(kconstant,append!([res[end,25]],res[end,2:6]),res[end,10:12],res[end,7:9],res[end,19:24],paramsint,-1);
#         end

#         eu0new,om0new = sadov2euler(sadvar0new[2:6],savar0new[1],IV)
#         if isrefframechanged
#             RMT = Transpose(RM)
#             paramsint[3]["MomentsOfInertia"] = RMT*IV;
#             numberOfFacets = paramsint[3]["numberOfFacets"]
#             if numberOfFacets >0
#                 for kk = 1:numberOfFacets    
#                     paramsint[3]["facets"][kk][3]=RMT*paramsint[3]["facets"][kk][3];
#                     paramsint[3]["facets"][kk][4]=RMT*paramsint[3]["facets"][kk][4];
#                 end
#                 paramsint[3]["intrinsicMagneticMoment"] = RMT*paramsint[3]["intrinsicMagneticMoment"];
#             end
#             Ri2b = RMT*euler2ibrotmat(eu0new);
#             eu0new = rotmatib2euler(Ri2b);
#             om0new = RMT*om0new;
#         end

#         # mean to osculating orbital
#         kepEl0new = res[end,19:24];

#         paramsint["includeEclipsesEffectsOnAttitude"] = true;
#         paramsint[5]["includeSrpTorque"] = true;

#         # osculating propagation
#         tVoscu,resosc = attitudeprop_quat_girf(append!(eu0new,om0new,kepEl0new),paramrsint[6]+tV[end]/3600.0/24.0,paramsint,tstep,tfinal-tV[end]);

#         tV = append!(tV,tVoscu+tV[end]);
#         res = [res;resosc]
#     end


#     # output
#    return tV,res;
# end

# function launchpropagation(ts,inlight,savarmean0,euanglesmean0,omVmean0,satinlight,equimean0,paramsint,kconstant,includeSrpTorque,mjd2000_0,isrefframechanged,RM)
#     params = deepcopy(paramsint);
#     params[6] = mjd2000_0;

    
#     if satinlight==inlight
#         change = false;
#     else
#         change = true;
#     end

#     if change
#         if abs(1.0-abs(savarmean0[3]/savarmean0[2]))<1e-13 
#             params[5]["includeSrpTorque"] = includeSrpTorque*(!inlight);
#             savarosc = transformationci_triaxial_sadovlike(kconstant,savarmean0,euanglesmean0,omVmean0,equimean0,params,-1);
#             euanglesosc,omVosc = sadov2euler(savarosc[2:6],savarosc[1],params[3]["MomentsOfInertia"]);
#             params[5]["includeSrpTorque"] = includeSrpTorque*inlight;
#             savarmean0 = transformationci_triaxial_sadovlike(kconstant,savarosc,euanglesosc,omVosc,equimean0,params,1);
#         else
#             println("change ",includeSrpTorque*(!inlight)," ",includeSrpTorque*inlight)
#             params[5]["includeSrpTorque"] = includeSrpTorque*(!inlight);
#             savarosc = transformationci_triaxial(kconstant,savarmean0,euanglesmean0,omVmean0,equimean0,params,-1);
#             euanglesosc,omVosc = sadov2euler(savarosc[2:6],savarosc[1],params[3]["MomentsOfInertia"]);
#             println(savarosc,euanglesosc,omVosc)
#             params[5]["includeSrpTorque"] = includeSrpTorque*inlight;
#             savarmean0 = transformationci_triaxial(kconstant,savarosc,euanglesosc,omVosc,equimean0,params,1);
#         end
#         # if inlight
#         #     if abs(1.0-abs(savarmean0[3]/savarmean0[2]))<1e-13 
#         #         params[5]["includeSrpTorque"] = includeSrpTorque*false;
#         #         savarosc = transformationci_triaxial_sadovlike(kconstant,savarmean0,euanglesmean0,omVmean0,equimean0,params,-1);
#         #         euanglesosc,omVosc = sadov2euler(savarosc[2:6],savarosc[1],params[3]["MomentsOfInertia"]);
#         #         params[5]["includeSrpTorque"] = includeSrpTorque*true;
#         #         savarmean0 = transformationci_triaxial_sadovlike(kconstant,savarosc,euanglesosc,omVosc,equimean0,params,1);
#         #     else
#         #         params[5]["includeSrpTorque"] = includeSrpTorque*false;
#         #         savarosc = transformationci_triaxial(kconstant,savarmean0,euanglesmean0,omVmean0,equimean0,params,-1);
#         #         euanglesosc,omVosc = sadov2euler(savarosc[2:6],savarosc[1],params[3]["MomentsOfInertia"]);
#         #         params[5]["includeSrpTorque"] = includeSrpTorque*true;
#         #         savarmean0 = transformationci_triaxial(kconstant,savarosc,euanglesosc,omVosc,equimean0,params,1);
#         #     end
#         # else
#         #     # println("change!!!!!!!!!")
#         #     #     println("savar0 ",savarmean0)
#         #     if abs(1.0-abs(savarmean0[3]/savarmean0[2]))<1e-13 
#         #         params[5]["includeSrpTorque"] = includeSrpTorque*true;
#         #         savarosc = transformationci_triaxial_sadovlike(kconstant,savarmean0,euanglesmean0,omVmean0,equimean0,params,-1);
#         #         euanglesosc,omVosc = sadov2euler(savarosc[2:6],savarosc[1],params[3]["MomentsOfInertia"]);
#         #         params[5]["includeSrpTorque"] = includeSrpTorque*false;
#         #         savarmean0 = transformationci_triaxial_sadovlike(kconstant,savarosc,euanglesosc,omVosc,equimean0,params,1);
#         #     else
                
#         #         params[5]["includeSrpTorque"] = includeSrpTorque*true;
#         #         savarosc = transformationci_triaxial(kconstant,savarmean0,euanglesmean0,omVmean0,equimean0,params,-1);
#         #         euanglesosc,omVosc = sadov2euler(savarosc[2:6],savarosc[1],params[3]["MomentsOfInertia"]);
#         #         params[5]["includeSrpTorque"] = includeSrpTorque*false;
#         #         savarmean0 = transformationci_triaxial(kconstant,savarosc,euanglesosc,omVosc,equimean0,params,1);
#         #     end
#         #     # println("osc ",savarosc)
#         #     # println("savarmean0new ",savarmean0)
#         # end
#     else
#         params[5]["includeSrpTorque"] = includeSrpTorque*inlight;
#     end

#     if abs(1.0-abs(savarmean0[3]/savarmean0[2]))<1e-10
#         savarlike = sadov2sadovlike(savarmean0[2:6],savarmean0[1]);
#         u0 = append!(savarlike,equimean0);
#         tV,res = prop_semianalytical_st_girf_triax_sadovlike(append!(savarlike,equimean0),params,ts,isrefframechanged,RM,params[3]["MomentsOfInertia"]);
#     else
#         try
#             tV,res = prop_semianalytical_st_girf_triax(append!(savarmean0,equimean0),params,ts,isrefframechanged,RM,params[3]["MomentsOfInertia"])
#         catch e
#             if e==ErrorException("singularity")
#                 savarlike = sadov2sadovlike(savarmean0[2:6],savarmean0[1]);
#                 u0 = zeros(13,1);
#                 u0 = append!(savarlike,equimean0);
#                 tV,res = prop_semianalytical_st_girf_triax_sadovlike(append!(savarlike,equimean0),params,ts,isrefframechanged,RM,params[3]["MomentsOfInertia"]);
#             else
#                 Base.error(e)
#             end
#         end
#     end
#     return tV,res;
# end

# function defineshadowsequence(ts,MLin,MLout,MLcurrent,meanmotion,propagationdirection)
#     ts1 = deepcopy(ts);
#     ts2 = deepcopy(ts);
#     ts3 = deepcopy(ts);

#     if MLin<MLout
#         if MLcurrent<MLin
#             if propagationdirection>0.0
#                 tinshadow = (MLin-MLcurrent)/meanmotion + ts[1];
#                 toutofshadow = (MLout-MLin)/meanmotion  + tinshadow;
#             else
#                 tinshadow = -(2*pi-MLout+MLcurrent)/meanmotion + ts[1];
#                 toutofshadow = -(MLout-MLin)/meanmotion  + tinshadow;
#             end
#             ts1 = append!(ts1[ts1.<=tinshadow],tinshadow);
#             ts2 = append!([tinshadow],ts2[ts2.>tinshadow .&& ts2.<=toutofshadow],toutofshadow);
#             ts3 = append!([toutofshadow],ts3[ts3.>toutofshadow]);
#             sequence = [true,false,true];
#         elseif MLcurrent>=MLin && MLcurrent<MLout
#             if propagationdirection >0
#                 toutofshadow =  (MLout-MLcurrent)/meanmotion + ts[1];
#                 tinshadow    =  (2*pi-MLout+MLin)/meanmotion + toutofshadow;
                
#             else
#                 toutofshadow =  -(MLcurrent-MLin)/meanmotion + ts[1];
#                 tinshadow    =  -(2*pi-MLout+MLin)/meanmotion + toutofshadow;
#             end
#             ts1 = append!(ts1[ts1.<=toutofshadow],toutofshadow);
#             ts2 = append!([toutofshadow],ts2[ts2.>toutofshadow .&& ts2.<=tinshadow],tinshadow);
#             ts3 = append!([tinshadow],ts3[ts3.>tinshadow]);
#             sequence = [false,true,false];
#         else # case MLcurrent>=MLout
#             if propagationdirection>0
#                 tinshadow    = (2*pi-MLcurrent+MLin)/meanmotion + ts[in];
#                 toutofshadow = (MLout-MLin)/meanmotion  + tinshadow;
#             else
#                 tinshadow    = -(MLcurrent-MLout)/meanmotion + ts[in];
#                 toutofshadow = -(MLout-MLin)/meanmotion  + tinshadow;
#             end
#             ts1 = append!(ts1[ts1.<=tinshadow],tinshadow);
#             ts2 = append!([tinshadow],ts2[ts2.>tinshadow .&& ts2.<=toutofshadow],toutofshadow);
#             ts3 = append!([toutofshadow],ts3[ts3.>toutofshadow]);
#             sequence = [true,false,true];
#         end
#     else
#         if MLcurrent>=MLout && MLcurrent< MLin
#             if propagationdirection>0.0
#                 tinshadow = (MLin-MLcurrent)/meanmotion + ts[1];
#                 toutofshadow = (2*pi+MLout-MLin)/meanmotion  + tinshadow;
#             else
#                 tinshadow = -(MLcurrent-MLout)/meanmotion + ts[1];
#                 toutofshadow = -(2*pi+MLout-MLin)/meanmotion  + tinshadow;
#             end
#             ts1 = append!(ts1[ts1.<=tinshadow],tinshadow);
#             ts2 = append!([tinshadow],ts2[ts2.>tinshadow .&& ts2.<=toutofshadow],toutofshadow);
#             ts3 = append!([toutofshadow],ts3[ts3.>toutofshadow]);
#             sequence = [true,false,true];
#         elseif MLcurrent<MLout
#             if propagationdirection >0.0
#                 toutofshadow =  (MLout-MLcurrent)/meanmotion + ts[1];
#                 tinshadow    =  (MLin-MLout)/meanmotion + toutofshadow;
#             else
#                 toutofshadow =  -(2*pi-MLin+MLcurrent)/meanmotion + ts[1];
#                 tinshadow    =  -(MLin-MLout)/meanmotion + toutofshadow;
#             end
#             ts1 = append!(ts1[ts1.<=toutofshadow],toutofshadow);
#             ts2 = append!([toutofshadow],ts2[ts2.>toutfshadow .&& ts2.<=tinshadow],tinshadow);
#             ts3 = append!([tinshadow],ts3[ts3.>tinshadow]);
#             sequence = [false,true,false];
#         else # case MLcurrent>=MLin
#             if propagationdirection >0.0
#                 toutofshadow =  (2*pi+MLout-MLcurrent)/meanmotion + ts[1];
#                 tinshadow    =  (MLin-MLout)/meanmotion + toutofshadow;
#             else
#                 toutofshadow =  -(MLcurrent-MLin)/meanmotion + ts[1];
#                 tinshadow    =  -(MLin-MLout)/meanmotion + toutofshadow;
#             end
#             ts1 = append!(ts1[ts1.<=tinshadow],tinshadow);
#             ts2 = append!([tinshadow],ts2[ts2.>tinshadow .&& ts2.<=toutofshadow],toutofshadow);
#             ts3 = append!([toutofshadow],ts3[ts3.>toutofshadow]);
#             sequence = [true,false,true];
#         end
#     end
#     # println("out: ",toutofshadow,", in: ",tinshadow)
#     return ts1,ts2,ts3,sequence;


# end








################## olde pieces of code
 # if paramsint[5]["includeDragTorque"] || paramsint[5]["includeDragAcc"]

    #     params = deepcopy(paramsint);
    #     tV = deepcopy(ts);
    #     tpropagated = 0.0;
    #     u0prop = deepcopy(u0);
    #     mcheck = true;
    #     ssign = ts[end]/abs(ts[end]);
    #     tspanprop = [0.0 ssign*params[4]["timeupdate"]];
    #     resM = zeros(size(ts)[1],12);

    #     meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[7:12],params[6],params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
    #             params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
    #             params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
    #             params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));
    #     params[4]["averageoverM_attitude"] = meanOverMDrag_A;
    #     params[4]["averageoverM_orbit"] =meanOverMDrag_O  ;

    #     while abs(tpropagated)<abs(ts[end]) && mcheck   
    #         println(" qui ",tpropagated/24/3600) 

    #         # propagation step 1
    #         if paramsint[5]["checkchaoticity"]
    #             cb = ContinuousCallback(conditionpatch,affectpatch!)
    #             probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0prop,tspanprop,params);
    #             solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback=cb);
    #         else
    #             probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0prop,tspanprop,params);
    #             solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10);
    #         end

    #         zetaend = solM.u[end][1];
    
    #         if paramsint[5]["checkchaoticity"] && (zetaend-1.0/(1.0+0.9999/(paramsint[3]["k_constant"])))<1e-12
    #             mcheck = false;
    #         end           

    #         if mcheck
    #             tM = solM.t .+ tpropagated;
    #             idx = abs.(ts).>= abs(tpropagated) .&& abs.(ts).<=abs(tM[end]);
    #             if maximum(idx)
    #                 solMsol  = solM(ts[idx].-tpropagated);
    #                 resMjj = mapreduce(permutedims,vcat,solMsol.u)
    #                 resM[idx,:] = resMjj;
    #             end
    #             # update
    #             mjd2000 = params[6] +  solM.t[end]/24/3600;
    #             resC   =  mapreduce(permutedims,vcat,solM.u);
    #             u0prop = resC[end,:];
    #             meanOverMDrag_A,meanOverMDrag_O = dragnumericallyaveragedterms(u0prop[7:12],mjd2000,params[2]["cpa1"],params[2]["cpa2"],params[2]["cpa3"],
    #                     params[1]["rPlanet"],params[1]["muPlanet"],params[1]["planetRotation"],params[1]["planet_flat"],params[1]["oblatenesscoeff"],
    #                     params[1]["zonalharmonicscoff"][1],paramsint[5]["includeDragTorque"],paramsint[5]["includeDragAcc"],
    #                     params[4]["atmosphericinfo"],params[4]["correction"],params[4]["couplingstrategy"],transpose(params[2]["equatorial2inertial"]));

    #             params[4]["averageoverM_attitude"] = meanOverMDrag_A;
    #             params[4]["averageoverM_orbit"] =meanOverMDrag_O ;
    #             params[4]["timeupdate"]   = pi*sqrt(u0prop[7]^3.0/params[1]["muPlanet"]);
    #             params[6] = mjd2000;
    #             tpropagated = tM[end];
    #             tspanprop = [0.0 ssign*params[4]["timeupdate"]];
    #         else
    #             tM = solM.t + tpropagated;
    #             idx = ts.>= tpropagated .&& ts.<=tM[end];
    #             if maximum(idx)
    #                 solMsol  = solM(ts[idx].-tpropagated);
    #                 resMjj = mapreduce(permutedims,vcat,solMsol.u)
    #                 resM[idx,:] = resMjj;
    #             end
    #             idxfinal = length(ts);
    #             for jj = 1 : size(ts)[1]
    #                 if ts[jj]>tpropagated+tM[end]
    #                     idxfinal = jj;
    #                     break;
    #                 end
    #             end
    #             resM[idxfinal,:] =mapreduce(permutedims,vcat,solM.u)[end,:];
    #             tV = append!(ts[1:idxfinal-1],tpropagated+tM[end]);
    #         end
    #     end
    # else
    #     tspan = [0.0 ts[end]];    
    #     if paramsint[5]["checkchaoticity"]
    #         cb = ContinuousCallback(conditionpatch,affectpatch!)
    #         probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0,tspan,paramsint);
    #         solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10,callback=cb);
    #         solMsol = solM(ts[ts.<=solM.t[end]]);
    #         resM = [mapreduce(permutedims,vcat,solMsol.u), mapreduce(permutedims,vcat,solM.u)[end,:]];
    #         tV = append!(solMsol.t,solM.t[end]);
    #     else

    #         probM = ODEProblem(attitudeMeanEv_semianalytical_triaxial,u0,tspan,paramsint);
    #         solM = solve(probM,Feagin14(),reltol=1e-13,abstol=1e-13,maxiters=1e+10); #,callback=cb);
    #         solM = solM(ts);
    #         resM = mapreduce(permutedims,vcat,solM.u);
    #         tV = solM.t;
    #     end
    # end