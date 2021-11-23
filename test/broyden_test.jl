#-----------------------------------------------------------------------
#   FULL TEST OF BROYDEN CODE
#-----------------------------------------------------------------------

d = 12

x0 = 2*rand(d)/1000
x1 = rand(d)
x2 = -π*rand(d)/10

f0 = x->2*π*(x-x0)
Jac0 = ForwardDiff.jacobian(f0,zeros(d));

@test (bsolve(f0,Jac0,f0(zeros(d)),zeros(d),2*d)-x0)./x0 ≈ zeros(d)

f1 = x->(x-x0).*(x-x1)
Jac1 = ForwardDiff.jacobian(f1,zeros(d));

@test (bsolve(f1,Jac1,f1(zeros(d)),zeros(d),2*d)-x0)./x0 ≈ zeros(d)

f2 = x->(x-x0).*(x-x1).*(x-x2)
Jac2 = ForwardDiff.jacobian(f2,zeros(d));

@test (bsolve(f2,Jac2,f2(zeros(d)),zeros(d),2*d)-x0)./x0 ≈ zeros(d)
