# needed packages
import Pkg

Pkg.add("Revise")
using Revise

push!(LOAD_PATH,"../")
import ATTITUDE
using ATTITUDE

Pkg.add("LinearAlgebra")
using LinearAlgebra

Pkg.add("Plots")
using Plots

Pkg.add("JLD2");
using JLD2;

Pkg.add("TimerOutputs")
using TimerOutputs

Pkg.add("LaTeXStrings")
using LaTeXStrings

Pkg.add("Printf")
using Printf

myrad2deg = 180.0/pi;
mydeg2rad = pi/180.0;

######## assumptions
# central body = earth
# inertial frame = equatorial

######## IMPORTANT
##  The osculating propagation works only
# considering a rotating propagator comoving with the satellite
# in principal axes

########################################################################################################
#### USER INPUTS
#######################################################################################################

######### initial state
# angular velocity [deg/s]
angVel = [0.0,0.1,4.0]; 
# euler angles 313 [deg]
euAng  = [0.0,20.0,80.0];
# semi-major axis [km]
sma = 20000 
# eccentricity
ec = 0.1
# inclination [deg]
incl = 20.0
# argument of perigee [deg]
om = 60.0
# longitude of the node [deg]
OM = 40.0
# true anomaly [deg]
nu = 0
# initial date [year-month-day-hour-min-sec]
date0 = [2023,01,01,0,0,0]

######### satellite
satname    = "satellite1"
# specify type of facets: triangular or square
typeoffacets = "triangular"

######### settings
# attitude pert
includeGravityTorque   = true;
includeMagneticTorque  = true;
includeSrpTorque       = true;
includeDragTorque      = false;
includeEclipsesEffectsOnAttitude = true;

# orbit pert
includeZonalHarmsAcc   = true;
maxzonalharmonics      = 5;
includeThirdBodyAcc    = true;
includeSunGravityAcc   = true;
includeSrpAcc          = true;
includeDragAcc         = false;
includeEclipsesEffectsOnOrbit= true;

# other settings
higherordercorrections = true
osculatingpropagationinchaoticregion = false;
checkchaoticity = false;

# atmospheric model
## atmospheric model = 1 --> exponential [RECOMMENDED]
## atmospheric model = 2 --> nrlmsise0
atmosphericmodel = 1
solarfluxandgeomagneticindexMAT = load_object("DATA/solarfluxandgeomagneticindexMAT.jld2");
AtmM = solarfluxandgeomagneticindexMAT;

######### time of propagation [sec]
timePropagation = 24.0*3600*1;
tstep = 30; 

######### outputs
folderouput = "folderoutput_example";
fileoutputname = "example2"

##############################################################################################################
# inputs managment
###############################################################################################################

############## settings ######################################################################
settings = Dict("includeMagneticTorque"=>includeMagneticTorque,"includeGravityTorque"=>includeGravityTorque,"includeSrpTorque"=>includeSrpTorque,"includeDragTorque"=>includeDragTorque);
settings = merge(settings,Dict("includeZonalHarmsAcc"=>includeZonalHarmsAcc,"includeThirdBodyAcc"=>includeThirdBodyAcc,"includeSunGravityAcc"=>includeSunGravityAcc,"includeSrpAcc"=>includeSrpAcc,"includeDragAcc"=>includeDragAcc));
settings = merge(settings,Dict("maxzonalharmonics"=>[maxzonalharmonics]));
settings = merge(settings,Dict("includeEclipsesEffectsOnAttitude"=>includeEclipsesEffectsOnAttitude,"includeEclipsesEffectsOnOrbit"=>includeEclipsesEffectsOnOrbit));
settings = merge(settings,Dict("higherordercorrections"=>higherordercorrections,"osculatingpropagationinchaoticregion"=>osculatingpropagationinchaoticregion,"checkchaoticity"=>checkchaoticity))
settings = merge(settings,Dict("sunposition"=>[0.5;0.3;0.4]))

############## satellite ###################################################################
fidsat = open(string("SATELLITES\\",satname,".txt"),"r");
fidsatl = readlines(fidsat);
close(fidsat)

if typeoffacets == "triangular"
    mass = parse(Float64,split(fidsatl[1],"= ")[2]);
    moiL = split(split(fidsatl[2],"= ")[2]," ");
    A    = parse(Float64,moiL[1]); 
    B    = parse(Float64,moiL[2]); 
    C    = parse(Float64,moiL[3]);
    IV = [A,B,C];
    numberOfFacets = 0;
    facets = Array{Vector{Any}}(undef,numberOfFacets,1)
    if size(fidsatl)[1]>3
        global numberOfFacets = parse(Int64,split(fidsatl[3],"= ")[2]);
        global facets = Array{Vector{Any}}(undef,numberOfFacets,1);
        for ll = 1:numberOfFacets
            coeffl = split(fidsatl[3+7*(ll-1)+1],"= "); coeffl = split(coeffl[2]," "); coeffl = coeffl[2:4];
            coeff  = [parse(Float64,coeffl[1]),parse(Float64,coeffl[2]),parse(Float64,coeffl[3])];
            vv1l   = split(fidsatl[3+7*(ll-1)+2],"= ");  vv1l = split(vv1l[2]," "); vv1l = vv1l[2:4];
            vv1    = [parse(Float64,vv1l[1]),parse(Float64,vv1l[2]),parse(Float64,vv1l[3])];
            vv2l   = split(fidsatl[3+7*(ll-1)+3],"= "); vv2l = split(vv2l[2]," "); vv2l = vv2l[2:4];
            vv2    = [parse(Float64,vv2l[1]),parse(Float64,vv2l[2]),parse(Float64,vv2l[3])];
            vv3l   = split(fidsatl[3+7*(ll-1)+4],"= "); vv3l = split(vv3l[2]," "); vv3l = vv3l[2:4];
            vv3    = [parse(Float64,vv3l[1]),parse(Float64,vv3l[2]),parse(Float64,vv3l[3])];
            nivl   = split(fidsatl[3+7*(ll-1)+5],"= ");  nivl = split(nivl[2]," "); nivl = nivl[2:4];
            niv    = [parse(Float64,nivl[1]),parse(Float64,nivl[2]),parse(Float64,nivl[3])];
            rhoivl = split(fidsatl[3+7*(ll-1)+6],"= ");  rhoivl = split(rhoivl[2]," "); rhoivl = rhoivl[2:4];
            rhoiv  = [parse(Float64,rhoivl[1]),parse(Float64,rhoivl[2]),parse(Float64,rhoivl[3])];
            area   = parse(Float64,split(fidsatl[3+7*(ll-1)+7],"= ")[2]);
            global facets[ll] = [coeff,area,[rhoiv,vv1,vv2,vv3],niv];
        end
    end
    if length(fidsatl)>3+7*(numberOfFacets-1)+7 
        cdval = parse(Float64,split(fidsatl[3+7*(numberOfFacets-1)+8],"= ")[2]);
    else
        cdval = 0.0;
    end
    if length(fidsatl)>3+7*(numberOfFacets-1)+8
        IML = split(split(fidsatl[3+7*(numberOfFacets-1)+9],"= ")[2]," ");
        IM  = [parse(Float64,IML[1]), parse(Float64,IML[2]),parse(Float64,IML[3])];
    else
        IM = [0.0,0.0,0.0];
    end
    satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>numberOfFacets,"facets"=>facets,"CD"=>cdval,"mass"=>mass);
elseif typeoffacets == "square"
    mass = parse(Float64,split(fidsatl[1],"= ")[2]);
    moiL = split(split(fidsatl[2],"= ")[2]," ");
    A    = parse(Float64,moiL[1]); 
    B    = parse(Float64,moiL[2]); 
    C    = parse(Float64,moiL[3]);
    IV = [A,B,C];
    numberOfFacets = 0;
    facets = Array{Vector{Any}}(undef,numberOfFacets,1)
    if size(fidsatl)[1]>3
        global numberOfFacets = parse(Int64,split(fidsatl[3],"= ")[2]);
        global facets = Array{Vector{Any}}(undef,numberOfFacets,1);
        for ll = 1:numberOfFacets
            coeffl = split(fidsatl[3+8*(ll-1)+1],"= "); coeffl = split(coeffl[2]," "); coeffl = coeffl[2:4];
            coeff  = [parse(Float64,coeffl[1]),parse(Float64,coeffl[2]),parse(Float64,coeffl[3])];
            vv1l   = split(fidsatl[3+8*(ll-1)+2],"= ");  vv1l = split(vv1l[2]," "); vv1l = vv1l[2:4];
            vv1    = [parse(Float64,vv1l[1]),parse(Float64,vv1l[2]),parse(Float64,vv1l[3])];
            vv2l   = split(fidsatl[3+8*(ll-1)+3],"= "); vv2l = split(vv2l[2]," "); vv2l = vv2l[2:4];
            vv2    = [parse(Float64,vv2l[1]),parse(Float64,vv2l[2]),parse(Float64,vv2l[3])];
            vv3l   = split(fidsatl[3+8*(ll-1)+4],"= "); vv3l = split(vv3l[2]," "); vv3l = vv3l[2:4];
            vv3    = [parse(Float64,vv3l[1]),parse(Float64,vv3l[2]),parse(Float64,vv3l[3])];
            vv4l   = split(fidsatl[3+8*(ll-1)+5],"= "); vv4l = split(vv4l[2]," "); vv4l = vv4l[2:4];
            vv4    = [parse(Float64,vv4l[1]),parse(Float64,vv4l[2]),parse(Float64,vv4l[3])];
            nivl   = split(fidsatl[3+8*(ll-1)+6],"= ");  nivl = split(nivl[2]," "); nivl = nivl[2:4];
            niv    = [parse(Float64,nivl[1]),parse(Float64,nivl[2]),parse(Float64,nivl[3])];
            rhoivl = split(fidsatl[3+8*(ll-1)+7],"= ");  rhoivl = split(rhoivl[2]," "); rhoivl = rhoivl[2:4];
            rhoiv  = [parse(Float64,rhoivl[1]),parse(Float64,rhoivl[2]),parse(Float64,rhoivl[3])];
            area   = parse(Float64,split(fidsatl[3+8*(ll-1)+8],"= ")[2]);
            global facets[ll] = [coeff,area,[rhoiv,vv1,vv2,vv3,vv4],niv];
        end
    end
    if length(fidsatl)>3+8*(numberOfFacets-1)+8 
        cdval = parse(Float64,split(fidsatl[3+8*(numberOfFacets-1)+9],"= ")[2]);
    else
        cdval = 0.0;
    end
    if length(fidsatl)>3+8*(numberOfFacets-1)+9
        IML = split(split(fidsatl[3+8*(numberOfFacets-1)+10],"= ")[2]," ");
        IM  = [parse(Float64,IML[1]), parse(Float64,IML[2]),parse(Float64,IML[3])];
    else
        IM = [0.0,0.0,0.0];
    end

    satellite=Dict("MomentsOfInertia"=>IV,"intrinsicMagneticMoment"=>IM,"numberOfFacets"=>numberOfFacets,"facets"=>facets,"CD"=>cdval,"mass"=>mass);
else 
    Base.error("only triangular or square facets")
end



#################### initial state######################################################################

incl    = incl*mydeg2rad;
om      = om*mydeg2rad;
OM      = OM*mydeg2rad;
nu      = nu*mydeg2rad;
euangles0  = euAng *mydeg2rad;
angVel     = angVel* mydeg2rad;

kepEl     = [sma,ec,incl,om,OM,nu];

########### central body #####################################################################################
planetparams = Dict("centralBodyIDX"=>3);

########## inertial ref frame ################################################################################
inertialrefframe = Dict();

######### atmosphere info #####################################
if atmosphericmodel == 1
    atmosphere = Dict("atmosphericmodel"=>atmosphericmodel);
else
    atmosphere = Dict("atmosphericmodel"=>atmosphericmodel,"AtmM"=>AtmM);
end

##############################################################################################################
# initialisation of time of propagation
tfinal = timePropagation;

################################################################################################################
###################### propagation #####################################################
civ    = (append!(deepcopy(euangles0),angVel,kepEl));
params = [planetparams;inertialrefframe;satellite;atmosphere;settings];

println("semi-analytical solution under computation")
to = TimerOutput();
enable_timer!(to)
tV3, resPre = ATTITUDE.attitudeaveragedprop_semianalytical_girf(civ,date0,params,tstep,tfinal)
res = zeros(size(tV3)[1],32)
res[:,1] = tV3;
res[:,2:7] = resPre[:,7:12];
res[:,8:12] = resPre[:,19:23];
for jj=1:size(tV3)[1]
    Ri2b = ATTITUDE.euler2ibrotmat(resPre[jj,10:12]);
    res[jj,13:15] = ATTITUDE.rotmatib2taitbryan(Ri2b);
end
res[:,16:21] = resPre[:,1:6];
res[:,22] = resPre[:,25];
res[:,23] = resPre[:,24];
res[:,24] = sqrt.(res[:,2].^2.0+res[:,3].^2.0+res[:,4].^2.0); 
res[:,25] = acos.(res[:,4]./res[:,24]); 
res[:,26] = mod.(atan.(res[:,3],res[:,2]),2.0*pi);
res[:,27:32] = resPre[:,35:40];

disable_timer!(to)
println(to)

ftest = open(string(folderouput,"\\semianalyticaljulia_",fileoutputname,".txt"),"w");
for kk=1:size(tV3)[1]
    @printf(ftest,"%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e\n", res[kk,1], res[kk,2], res[kk,3], res[kk,4], res[kk,5], res[kk,6], res[kk,7], res[kk,8], res[kk,9], res[kk,10], res[kk,11], res[kk,12], res[kk,13],res[kk,14],res[kk,15],res[kk,16],res[kk,17],res[kk,18],res[kk,19],res[kk,20],res[kk,21],res[kk,22],res[kk,23],res[kk,24],res[kk,25],res[kk,26],res[kk,27],res[kk,28],res[kk,29],res[kk,30],res[kk,31],res[kk,32]);
end
close(ftest);
fper = open(string(folderouput,"\\semianalyticaljulia_",fileoutputname,"_performance.txt"),"w");
@printf(fper,"%s\n",to) 
close(fper);

fid = open(string(folderouput,"\\myinput_semianalyticaljulia.txt"),"w");
@printf(fid,"initial date = %d %d %d %d %d %.4e \n",date0[1],date0[2],date0[3],date0[4],date0[5],date0[6]);
@printf(fid,"a = %20.16e \n", sma)
@printf(fid,"e = %20.16e \n", ec)
@printf(fid,"incl = %20.16e \n", incl*myrad2deg)
@printf(fid,"om = %20.16e \n", om*myrad2deg)
@printf(fid,"OM = %20.16e \n", OM*myrad2deg)
@printf(fid,"nu = %20.16e\n", nu*myrad2deg)
@printf(fid,"angularvelocity = %20.16e %20.16e %20.16e \n",angVel[1]*myrad2deg,angVel[2]*myrad2deg,angVel[3]*myrad2deg)
@printf(fid,"anglestype = eu313\n")
@printf(fid,"angles = %20.16e %20.16e %20.16e\n",euangles0[1]*myrad2deg, euangles0[2]*myrad2deg, euangles0[3]*myrad2deg);
close(fid)



##################################### !!!!! IN THE OUTPUT file
## in each row:

#######
# time
# x component angular velocity
# y component angular velocity
# z component angular velocity
# euler angles 313
# semi-major axis
# equinoctial element P1 = ec*sin(OM+om) (ec=eccentricity, OM=longitude of the node, om=argumentof perigee)
# equinoctial element P2 = ec*cos(OM+om)
# equinoctial element Q1 = tan(incl/2)*sin(OM) (incl=inclination)
# equinoctial element Q2 = tan(incl/2)*cos(OM)
# Tait-Bryan angles [yaw pitch roll]
# Sadov variables
# Sadov parameter skm
# orbital mean longitude
# magnitude angular velocity
# angle between angular velocity and z axis of rotating frame
# angle between the projection of the angular velocity over the xy plane of the rotating frame and the x axis
# semi-major axis
# eccentricity
# inclination
# argument of the perigee
# longitude of the node
# true anomaly



