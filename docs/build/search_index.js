var documenterSearchIndex = {"docs":
[{"location":"metrics.html#Metrics","page":"Metrics","title":"Metrics","text":"","category":"section"},{"location":"metrics.html#Kerr-Schild-metric","page":"Metrics","title":"Kerr-Schild metric","text":"","category":"section"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"The Kerr-Schild metric takes the form:","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"    g_mu nu = eta_mu nu + f  k_mu  k_nu ","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"where ","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"  k_mu = left( fracr  x + a  yr^2+a^2\n                 fracr  y - a  xr^2+a^2 fraczrright)","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"  f = frac2  G  M  r^3r^4+a^2  z^2","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"and r is implicitly defined by:","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"fracx^2+y^2r^2+a^2 + fracz^2r^2 = 1 ","category":"page"},{"location":"metrics.html#Weak-field-metric","page":"Metrics","title":"Weak field metric","text":"","category":"section"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"For solar system and terrestrial positioning, the weak field metric suffices. The weak field metric has the form:","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"g_mu nu = eta_mu nu - 2  V  delta_mu nu","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"where V is the gravitational potential.","category":"page"},{"location":"metrics.html#Gordon-metric","page":"Metrics","title":"Gordon metric","text":"","category":"section"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"To incorporate atmospheric and ionospheric effects, one uses the analogue Gordon metric, which takes the form (with g_mu nu being the gravitational metric):","category":"page"},{"location":"metrics.html","page":"Metrics","title":"Metrics","text":"    barg_mu nu = g_mu nu + \n    left(1-frac1n^2right) u_mu u_nu ","category":"page"},{"location":"geodesics.html#Geodesic-solver","page":"Geodesic solver","title":"Geodesic solver","text":"","category":"section"},{"location":"geodesics.html#Null-Geodesics:-Hamilton's-equations","page":"Geodesic solver","title":"Null Geodesics: Hamilton's equations","text":"","category":"section"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"Null geodesics in a spacetime geometry described by a metric g_μν may be described by the Hamiltonian:","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"H = frac12 g^μν p_μ p_ν","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"and the associated Hamilton equations:","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"fracdx^μdλ = fracHp_μ","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"fracdp_μdλ = - fracHx^μ","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"The conjugate momenta p_μ are defined as:","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"p_μ = g_μν fracdx^νdλ","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"and for null geodesics, the initial data satisfies:","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"leftδ^i_μ fracdx^μdλright_λ=0 = bf v^i","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"left g_μν(x_0) \nfracdx^μdλfracdx^νdλ right_λ=0 = 0","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"The solution to Hamilton's equations is formally given by x^μ=x^μ(λx_0bf v) (with x_0 denoting the initial position of the geodesic), and since λ (being an affine parameter) can be redefined linearly, it is appropriate to set up the problem so that λ01, with λ=0 being the initial point and λ=1 is the final point.","category":"page"},{"location":"geodesics.html#Implementation","page":"Geodesic solver","title":"Implementation","text":"","category":"section"},{"location":"geodesics.html#Hamiltonian","page":"Geodesic solver","title":"Hamiltonian","text":"","category":"section"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"The Hamiltonian function takes the form:","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"squirrel.HamGeo","category":"page"},{"location":"geodesics.html#squirrel.HamGeo","page":"Geodesic solver","title":"squirrel.HamGeo","text":"HamGeo( Z::RealVec , gfunc::Function )\n\nThe HamGeo function computes the geodesic Hamiltonian  H = frac12 g^μν p_μ p_ν from the metric g_μν provided in the function gfunc and the  spacetime position x and momentum p encoded in Z=(x,p).\n\n\n\n\n\n","category":"function"},{"location":"geodesics.html#Hamilton's-equations","page":"Geodesic solver","title":"Hamilton's equations","text":"","category":"section"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"Hamilton's equations may be written in terms of the phase space coordinate z^α, where z=(xp).","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"fracdz^αdλ = J^αβ fracHz^β","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"where J^αβ is the symplectic matrix, which has the block matrix form:","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"J =\nleft\n  beginarraycccc\n     O    I  \n     -I    O  \n  endarray\nright","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"where O is a 44 matrix of zeros and I is the identity matrix. The symplectic matrix is implemented as an operator acting on a vector Hz^β:","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"squirrel.Jsympl","category":"page"},{"location":"geodesics.html#squirrel.Jsympl","page":"Geodesic solver","title":"squirrel.Jsympl","text":"Jsympl( Zarg::RealVec )\n\nThe Jsympl function effectively applies the symplectic matrix J^αβ to the vector Zarg (the vector Hz^β).\n\n\n\n\n\n","category":"function"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"The quantity J^αβ fracHz^β is evaluated in the following function: ","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"squirrel.ZdotGeo","category":"page"},{"location":"geodesics.html#squirrel.ZdotGeo","page":"Geodesic solver","title":"squirrel.ZdotGeo","text":"ZdotGeo( Z::RealVec , gfunc::Function )\n\nThe ZdotGeo function computes the gradient of the Hamiltonian and applies the symplectic operator to the gradient.\n\n\n\n\n\n","category":"function"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"The symplectic operator ","category":"page"},{"location":"geodesics.html#Geodesic-solver-function","page":"Geodesic solver","title":"Geodesic solver function","text":"","category":"section"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"Geodesics are solved with the following function, which outputs the endpoint (λ=1) of the solution to Hamilton's equations:","category":"page"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"squirrel.solveZ","category":"page"},{"location":"geodesics.html#squirrel.solveZ","page":"Geodesic solver","title":"squirrel.solveZ","text":"solveZ( Z0::RealVec , gfunc::Function , δ1::Real , δ2::Real ,\n        integrator=AutoVern7(Rodas5()) , δt=0 )\n\nThe function solveZ performs the integration of geodesics given the initial data Z0, and the metric function gfunc. The integration is performed using the integrator specified by the variable integrator using the relative tolerance parameter δ1 and the absolute tolerance  parameter δ2. The parameter δt determines the timestep parameter  dt in ODEProblem. \n\nThe function solveZ outputs only the endpoint of the solution for the geodesic equation.\n\n\n\n\n\n","category":"function"},{"location":"geodesics.html","page":"Geodesic solver","title":"Geodesic solver","text":"Following the recommendations in the ODE Solver documentation, the integrators AutoVern7(Rodas5()) and AutoVern9(Rodas5()) are used in squirrel.jl.","category":"page"},{"location":"index.html#Home","page":"Home","title":"Home","text":"","category":"section"},{"location":"index.html#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Relativistic positioning refers to the concept of establishing spacetime positions from proper time broadcasts emitted by a system of satellites. Central to relativistic positioning is the relativistic location location problem, which is the problem of finding the intersection of future pointing light cones from a collection of at least four emission points. Algorithms for relativistic location in flat spacetime are provided in the cereal.jl package. This package, squirrel.jl, contains a collection of functions for the relativistic location problem in slightly curved spacetime geometries.","category":"page"},{"location":"index.html#Short-tutorial","page":"Home","title":"Short tutorial","text":"","category":"section"},{"location":"index.html#Setup","page":"Home","title":"Setup","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The squirrel.jl code was written for and tested in Julia 1.6; we recommend Julia 1.6 or newer. To add the code, run the following command in the package manager for the Julia REPL:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"pkg> add https://github.com/justincfeng/squirrel.jl/","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Once added, one may access the cereal module with the following  command:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> using squirrel","category":"page"},{"location":"index.html#Positioning-examples","page":"Home","title":"Positioning examples","text":"","category":"section"},{"location":"index.html#Vacuum-case","page":"Home","title":"Vacuum case","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"One begins by defining a metric. The Kerr-Schild metric for a rotating object with the mass and angular momentum of the Earth is given by the following function:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> gks = squirrel.metric.ge","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The resulting function gks takes a four-component vector xc=[t,x,y,z] (representing spacetime coordinates) as an argument, and returns a 44 matrix of metric components. The units are normalized to Earth mass (1 length unit = 0.4435 cm), and chosen so that the speed of light is 1. One can randomly generate a target point Xtar and a set of 5 emission points X on the past light cone of Xtar with the following:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> (X,Xtar) = squirrel.seval.pgen(6e9,gks,1e-14,5) ;","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The first argument in the function on the left hand side specifies (roughly) the spatial radius of the emission points from the origin, and the third argument specified the tolerance for the geodesic integrator. The target point Xtar is placed on the WGS84 reference ellipsoid (defined with respect to the Cartesian Kerr-Schild coordinates).","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The intersection of future light cones from the emission points is computed with the squirrel.locator function:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> Xs = squirrel.locator(X,gks,1e-10)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"The third argument is the tolerance for the geodesic solvers; the tolerance is looser here to minimize computation time. The accuracy of the result may be estimated by comparing Xs and Xtar:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> ΔX = Xs-Xtar","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Upon multiplying by the conversion factor 0.4435 to obtain the result in units of centimeters (0.4435*ΔX), one typically obtains a result well in the submillimeter range. The relative error may be obtained by running:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> squirrel.norm(ΔX)/squirrel.norm(Xtar)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"and one typically obtains errors on the order of sim 10^-12 - 10^-13.","category":"page"},{"location":"index.html#Including-atmospheric-and-ionospheric-effects","page":"Home","title":"Including atmospheric and ionospheric effects","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"One can incorporate the effects of the atmosphere and ionosphere with the following metric:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> g = squirrel.metric.g","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"This metric is the Gordon metric for light propagatation through media; here, a simple model for the atmosphere and ionosphere is implemented. One may repeat the steps of the vacuum case to obtain the errors:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"julia> (X,Xtar) = squirrel.seval.pgen(6e9,g,1e-14,5) ;\n\njulia> Xs = squirrel.locator(X,g,1e-10)\n\njulia> ΔX = Xs-Xtar\n\njulia> 0.4435*ΔX\n\njulia> squirrel.norm(ΔX)/squirrel.norm(Xtar)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"In this case, one typically finds that the errors are larger than the vacuum case by an order of magnitude, though still in the submillimeter range.","category":"page"},{"location":"index.html#Evaluation","page":"Home","title":"Evaluation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Scripts are provided for a more comprehensive evaluation of the accuracy and performance of the functions in squirrel.jl.","category":"page"},{"location":"index.html#References","page":"Home","title":"References","text":"","category":"section"},{"location":"evaluation.html#Evaluation-functions","page":"Evaluation","title":"Evaluation functions","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"The functions described here are contained in the seval submodule. These functions are encapsulated so that they provide a somewhat independent evaluation of the locator functions in the squirrel.jl code. At present, the parameters of the evaluation functions have values corresponding to terrestrial positioning.","category":"page"},{"location":"evaluation.html#Basic-strategy","page":"Evaluation","title":"Basic strategy","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"The functions in seval, particularly seval.gen and seval.main, may be used to evaluate the accuracy of the locator functions implemented in squirrel.jl according to the following strategy:","category":"page"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"Generate (stochastically) a target point Xtar.\nGenerate ne emission points X on the past light cone of Xtar  by solving the geodesic equation for some metric g.\nFeed the emission points X and metric g into the  squirrel.locator or squirrel.locator4 function and compare with  target point Xtar.","category":"page"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"The function seval.gen is implements steps 1. and 2., and the function seval.main implements step 3.","category":"page"},{"location":"evaluation.html#Structs-and-format-conversion-functions","page":"Evaluation","title":"Structs and format conversion functions","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"The evaluation submodule seval defines two datatypes, one for the generation of test cases (TestCases), and one for storing the results of the evaluation (TestData).","category":"page"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.TestCases","category":"page"},{"location":"evaluation.html#squirrel.seval.TestCases","page":"Evaluation","title":"squirrel.seval.TestCases","text":"The TestCases datatype may be populated by the associated function of the form:\n\nTestCases( par::NTuple , N::Int , ne::Int , X::Array{RealMtx,1} , Xtar::Array{RealVec,1} )\n\nThe variables are defined in the following way:\n\npar     ::  NTuple              # Parameter tuple\nN       ::  Int                 # Number of generated test cases\nne      ::  Int                 # Number of emission points\nX       ::  Array{RealMtx,1}    # Emission points\nXtar    ::  Array{RealVec,1}    # Target point\n\n\n\n\n\n","category":"type"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.TestData","category":"page"},{"location":"evaluation.html#squirrel.seval.TestData","page":"Evaluation","title":"squirrel.seval.TestData","text":"The TestData datatype may be populated by the associated function of the form:\n\nTestData(   par::NTuple , N::Int    , \n            X::Array{RealMtx,1}     ,   Xtar::Array{RealVec,1}     ,\n            Xc::Array{RealVec,1}    ,   Xsc::Array{RealVec,1}      ,\n            erh::RealVec    ,   erv::RealVec    ,   err::RealVec   ,\n            erhC::RealVec   ,   ervC::RealVec   ,   errC::RealVec  ,\n            Xc2::Array{RealVec,1}   ,   Xsc2::Array{RealVec,1}     ,\n            erh2::RealVec   ,   erv2::RealVec   ,   err2::RealVec  ,\n            erhC2::RealVec  ,   ervC2::RealVec  ,   errC2::RealVec )\n\nThe quantities are defined in the following way:\n\npar     ::  NTuple              # Parameter tuple\nN       ::  Int                 # Number of cases\nX       ::  Array{RealMtx,1}    # Emission points\nXtar    ::  Array{RealVec,1}    # Target point\nXc      ::  Array{RealVec,1}    # Flat spacetime position\nXsc     ::  Array{RealVec,1}    # Squirrel position\nerh     ::  RealVec             # Horizontal error relative to Xsc\nerv     ::  RealVec             # Vertical error relative to Xsc\nerr     ::  RealVec             # Total error relative to Xsc\nerhC    ::  RealVec             # Horizontal error relative to Xc\nervC    ::  RealVec             # Vertical error relative to Xc\nerrC    ::  RealVec             # Total error relative to Xc\nXc2     ::  Array{RealVec,1}    # Flat spacetime position     (Aux.)\nXsc2    ::  Array{RealVec,1}    # Squirrel position      (Auxiliary)\nerh2    ::  RealVec             # Horizontal errors      (Auxiliary)\nerv2    ::  RealVec             # Vertical error         (Auxiliary)\nerr2    ::  RealVec             # Total error            (Auxiliary)\nerhC2   ::  RealVec             # Horizontal errors      (Auxiliary)\nervC2   ::  RealVec             # Vertical error         (Auxiliary)\nerrC2   ::  RealVec             # Total error            (Auxiliary)\n\nThe quantities suffixed with 2 are auxiliary quantities which are needed in the four-point case, in which the special relativistic location algorithms generally suffer from the bifurcation problem.\n\n\n\n\n\n","category":"type"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"The functions seval.gen and seval.main can generate large amounts of data. The data may be saved to a file using Julia's built in serializer, but it is recommended that the data be converted to a tuple. The following functions may be used to convert the between tuples and the  TestCases/TestData datatypes:","category":"page"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.tc2tup","category":"page"},{"location":"evaluation.html#squirrel.seval.tc2tup","page":"Evaluation","title":"squirrel.seval.tc2tup","text":"tc2tup( tc::TestCases )\n\nThe tc2tup function changes tc from an object of type TestCases to  a tuple.\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.td2tup","category":"page"},{"location":"evaluation.html#squirrel.seval.td2tup","page":"Evaluation","title":"squirrel.seval.td2tup","text":"td2tup( td::TestData )\n\nThe td2tup function changes td from an object of type TestData to  a tuple.\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.tup2tc","category":"page"},{"location":"evaluation.html#squirrel.seval.tup2tc","page":"Evaluation","title":"squirrel.seval.tup2tc","text":"tup2tc( tctup::Tuple )\n\nThe tc2tup function changes tctup from a Tuple to an object of type  TestCases.\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.tup2td","category":"page"},{"location":"evaluation.html#squirrel.seval.tup2td","page":"Evaluation","title":"squirrel.seval.tup2td","text":"tup2td( tdtup::Tuple )\n\nThe td2tup function changes tdtuple from a Tuple to an object of  type TestData.\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"In some cases, one may wish to change to a higher precision floating point type for the generated samples. For this purpose, the following function is provided:","category":"page"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.tcfl","category":"page"},{"location":"evaluation.html#squirrel.seval.tcfl","page":"Evaluation","title":"squirrel.seval.tcfl","text":"tcfl( tc::TestCases , tpfl::DataType=Float64 )\n\nThe tcfl function changes the floating point type of tc to tpfl.\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html#Test-case-generation","page":"Evaluation","title":"Test case generation","text":"","category":"section"},{"location":"evaluation.html#Point-generation","page":"Evaluation","title":"Point generation","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"Emission points are generated on the past light cone of a target point. The idea is to stochastically generate a target point Xtar, and then generate initial data for N past-directed null geodesics, which are then integrated to obtain the emission points.","category":"page"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.pgen","category":"page"},{"location":"evaluation.html#squirrel.seval.pgen","page":"Evaluation","title":"squirrel.seval.pgen","text":"pgen( Rs::Real , gfunc::Function , tol::Real , ne::Int=4 , \n      Δψ::Real=Δψ0 )\n\nThe pgen function generates a target point Xtar and ne emission points in a 4×ne matrix X such that the points in X lie on the past light cone of Xtar with respect to the metric gfunc, and have spatial radii values of ~Rs (defined with respect to spatial origin). The parameter tol is the tolerance parameter for the integration. This function returns the tuple (X,Xtar).\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html#Multiple-case-generator","page":"Evaluation","title":"Multiple case generator","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"The following function generates N sets of target points Xtar and emission points X. It returns a quantity of datatype TestCases.","category":"page"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.gen","category":"page"},{"location":"evaluation.html#squirrel.seval.gen","page":"Evaluation","title":"squirrel.seval.gen","text":"gen( N::Int , g::Function , ne::Int=4 , Δψ::Real=Δψ0 \n          , tpfl::DataType=Float64 \n          , REs ::  Real=1.4365276211950395e9       # Earth radius\n          , RR  ::  Real=1.4365277e9    # Just above Earth surface\n          , Rs  ::  Real=6e9            # Satellite radius\n          , tolh::  Real=1e-14          # Tolerance for ODE solvers\n          , tol ::  Real=1e-10          # Tolerance for ODE solvers\n        )\n\nThe gen function generates N sets of target points Xtar and  emission points X with respect to the metric gfunc. This function returns a quantity of the TestCases datatype.\n\n\n\n\n\n","category":"function"},{"location":"evaluation.html#Evaluation-function","page":"Evaluation","title":"Evaluation function","text":"","category":"section"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"The following is the main evaluation function. It returns a quantity of datatype TestData.","category":"page"},{"location":"evaluation.html","page":"Evaluation","title":"Evaluation","text":"squirrel.seval.main","category":"page"},{"location":"evaluation.html#squirrel.seval.main","page":"Evaluation","title":"squirrel.seval.main","text":"main(  tc::TestCases , sloc::Function , g::Function , Nx::Int=-1 , \n       tpfl::DataType=Double64 , tol::Real=1e-10 , ξ::Real=2e1 , \n       nb::Int=24 , ne::Real=6 )\n\nThis is the main evaluation function. It takes as input a collection of test cases contained in a quantity of datatype TestCases, and evaluates the squirrel locator function sloc on each test case relative to a given metric function g. The function main returns results in a quantity of type TestData.\n\nThe variable Nx is a limiter for the number of test cases to run. The variable tpfl specifies the floating point precision for the calculations. variables tol, ξ, nb,  and ne are inputs for the sloc function, and correspond respectively to the integration tolerance, outlier detection threshold, number of steps for the Broyden algorithm, and number of emission points.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html#Squirrel-algorithm","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"","category":"section"},{"location":"squirrel.html#Description","page":"Squirrel algorithm","title":"Description","text":"","category":"section"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"The squirrel algorithm takes the coordinates of ne ≥ 4 emission points, arranged as column vectors in a 4×ne matrix X,  and computes the intersection of future pointing null geodesics defined with respect to a (slightly curved) spacetime metric g. ","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"Null geodesics may be described as solutions to the geodesic equation (implemented here in the form of Hamilton's equations), which may be formally written as the functions:","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"    x^μ_I(λX_Ibf v_I)","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"where the indices I12n_e distinguish the emission points X_I and their associated null geodesics, and bf v_I denotes the spatial components of the initial four-velocity vector for the null geodesic (the time component is determined by the requirement that the four-velocity is null). ","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"Given a collection of four geodesic functions x_1x_2x_3x_4 (each having the form x^μ_I=x^μ_I(λx_0bf v)), the condition that they intersect is the vanishing of the following vector-valued function:","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"f = left( x_1 - x_2  x_1 - x_3  x_1 - x_4 right)","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"which may be formally written as f=f(X_1X_2X_3X_4bf v_1bf v_2bf v_3bf v_4). The squirrel algorithm seeks to find the roots of the function f.","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"The squirrel algorithm is then summarized:","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"Apply the flat spacetime algorithm to the emission points X to  obtain an initial guess for the zeros of the function f.\nApply a root finding algorithm to the function f.","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"A quasi-Newton Broyden algorithm (which will be described in detail below) is employed to do the root-finding; in such a method, the Jacobian for f is computed once in the first iteration of the root finding algorithm, and is updated in the subsequent iterations. The function f is computed by way of numerical integration of geodesics; if the numerical integration is performed using native Julia libraries, one can compute the Jacobian by way of automatic differentiation.","category":"page"},{"location":"squirrel.html#Geodesic-endpoint-function","page":"Squirrel algorithm","title":"Geodesic endpoint function","text":"","category":"section"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"To compute the function f, the geodesic endpoint function x^μ_I=x^μ_I(λX_Ibf v)_λ=0 is implemented as:","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.gsolve","category":"page"},{"location":"squirrel.html#squirrel.gsolve","page":"Squirrel algorithm","title":"squirrel.gsolve","text":"gsolve( Xi::RealVec , Vi::RealVec , g::Function , tol::Real , integrator=AutoVern7(Rodas5()) )\n\nThe gsolve function takes an initial point Xi and four velocity Vi and computes the endpoint of a future pointing null geodesic in the metric func g. The variable integrator specifies the integration scheme, and tol is the tolerance parameter.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"With the geodesic endpoint functions x^μ_I=x^μ_I(λX_Ibf v)_λ=0","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.zF","category":"page"},{"location":"squirrel.html#squirrel.zF","page":"Squirrel algorithm","title":"squirrel.zF","text":"zF( Vid::RealVec , Zi::RealMtx , g::Function , tol::Real )\n\nThe function zF returns differences between the endpoints of four geodesics for the initial data encoded in Vid and Zi, and the metric function g. The variable tol is the tolerance parameter for the integration. This function vanishes when the four geodesics intersect.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"To find the roots of the function f, the squirrel algorithm first computes the Jacobian of f. ","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.gejac","category":"page"},{"location":"squirrel.html#squirrel.gejac","page":"Squirrel algorithm","title":"squirrel.gejac","text":"gejac( Xi::RealVec , Vi::RealVec , g::Function , δ::Real )\n\nThe function gejac computes the endpoint of a geodesic and its Jacobian. The variables have the same meaning as those in gsolve.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.geocJ","category":"page"},{"location":"squirrel.html#squirrel.geocJ","page":"Squirrel algorithm","title":"squirrel.geocJ","text":"geocJ( Zi::RealMtx , g::Function , δ::Real )\n\nThe function geocJ computes the Jacobian of the function zF from the Jacobian endpoints computed in gejac.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.idf","category":"page"},{"location":"squirrel.html#squirrel.idf","page":"Squirrel algorithm","title":"squirrel.idf","text":"idf( Zi::RealMtx , gfunc::Function , tol::Real , nb::Int )\n\nThe function idf applies the Broyden algorithm to the function zF, using the Jacobian initially computed with geocJ.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html#Root-finding-algorithm","page":"Squirrel algorithm","title":"Root finding algorithm","text":"","category":"section"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"Given some function f(x), the standard approach is the Newton method, which finds the roots according to the prescription:","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"    x_i+1 = x_i + bf J^-1_i f(x_i)","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"where bf J^-1_i is the inverse of the Jacobian matrix bf J for f(x) evaluated at x_i. However, in situations which the computation of the Jacobian matrix bf J becomes expensive, one may instead employ the Broyden method, which is a quasi-Newton root finding method. The Broyden algorithm involves first computing the Jacobian matrix bf J and its inverse for f(x). However, instead of computing the Jacobian matrix at each iteration (as is done in the Newton method), the Broyden algorithm updates the (inverse) Jacobian matrix according to the Sherman-Morrison update formula:","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"    bf J^-1_i+1 \n    = \n        bf J^-1_i\n        +\n        fracΔv^T_i - bf J^-1_i Δf_i \n        Δv^T_i bf J^-1_i Δf_i \n        Δv^T_i bf J^-1_i","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"which is implemented in the following function:","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.JiSMU","category":"page"},{"location":"squirrel.html#squirrel.JiSMU","page":"Squirrel algorithm","title":"squirrel.JiSMU","text":"JiSMU( ΔF::RealVec , Δx::RealVec , Ji::RealMtx )\n\nThe function JiSMU implements the Sherman-Morrison update formula.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"The Broyden update formula is implemented in the following function:","category":"page"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.bsolve","category":"page"},{"location":"squirrel.html#squirrel.bsolve","page":"Squirrel algorithm","title":"squirrel.bsolve","text":"bsolve( F::Function , J::RealMtx , f0::RealVec , x0::RealVec ,\n        nb::Int=24 )\n\nThe function bsolve implements the Broyden algorithm; in particular, it finds the roots of the function F(x), given the Jacobian matrix J and the initial guesses f0 and x0. The parameter nb specifies the maximum number of iterations.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html#Locator","page":"Squirrel algorithm","title":"Locator","text":"","category":"section"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.locator4","category":"page"},{"location":"squirrel.html#squirrel.locator4","page":"Squirrel algorithm","title":"squirrel.locator4","text":"locator4( X::RealMtx , Xc::RealVec , gfunc::Function , tol::Real , nb::Int=24 , idv::Bool=false , V::RealMtx=zeros(Float64,4,4) )\n\nThe function locator4 computes the intersection point from a set of four emission points X using the guess Xc. The intersection point is computed by first finding the initial data using the function idf, then integrating geodesics to obtain the result. If an improved guess for the four-velocities is available, one can set ìdv=true and specify the four-velocities as column vectors in the matrix V.\n\n\n\n\n\n","category":"function"},{"location":"squirrel.html","page":"Squirrel algorithm","title":"Squirrel algorithm","text":"squirrel.locator","category":"page"},{"location":"squirrel.html#squirrel.locator","page":"Squirrel algorithm","title":"squirrel.locator","text":"locator( X::RealMtx , gfunc::Function , δ::Real , nb::Int=24 , tpflc::DataType=Double64 , outthresh::Real=1e1 , ne::Int=5 )\n\nThe function locator computes the intersection point from a set of ne>4 emission points X by applying locator4 to all combinations of 4 points out of ne in X. A basic outlier detection algorithm (implemented in the function odetc) is applied to reduce errors.\n\n\n\n\n\n","category":"function"}]
}
