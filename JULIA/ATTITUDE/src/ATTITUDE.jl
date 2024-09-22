module ATTITUDE

using LinearAlgebra
using Elliptic
using DifferentialEquations
using Plots
using DataInterpolations
using SatelliteToolboxAtmosphericModels
using SatelliteToolboxTransformations
using SpecialFunctions
using Polynomials
using HypergeometricFunctions
using FastGaussQuadrature

include("functionslibrary.jl") 
include("astroconstants.jl") 
include("myElliptic.jl") 
include("polynomialequations.jl") 
include("atmosphere.jl"); 
include("orbitvarchange.jl"); 
include("attitudevarchange.jl"); 
include("rotationmatrices.jl"); 
include("timevarchange.jl") 
include("ephemeris.jl") 

### for averaged eq of motion - torque
include("srp_average.jl"); 
include("drag_average.jl");
include("eclipse.jl") 
include("averagedorbitpert.jl") 

# generatingfunction
include("generatingfunction_srp_drag_termsP1.jl") 
include("generatingfunction_srp_drag_termsP2.jl") 
include("generatingfunction_srp_drag_termsP3.jl") 

include("generatingfunctionsrpdrag_p1.jl"); 
include("generatingfunctionsrp_p2.jl");
include("generatingfunctiondrag_p2.jl");
include("generatingfunctiongravity.jl");
include("generatingfunctionmagnetic.jl");
include("generatingfunctionsterms_casesmk0.jl") 


#### averaged eq of moyion
include("sam_higherorderterms.jl");
include("eqMotion_averagedmodel_girf.jl");

#### eq of motion
include("eqMotion.jl");

#### propagators
include("semianalyticalpropagators_girf.jl")
include("propagators_girf.jl");






#OLD
#### for averaged eq of motion
# include("gravitymagnetic_average.jl");
# include("gravitymagnetic_average_girf.jl");
#### averaged eq of moyion
# include("eqMotion_averagedmodel.jl");
# include("eqMotion_averagedmodel_girf.jl");
#### eq motion
# include("eqMotion_withorbitalcoupling.jl");
#### propagators
# include("propagators.jl");
# include("semianalyticalpropagators.jl")

end # module ATTITUDE



