
function result = thermal_integrand_sphere_InGaAs(k,p,T,MatParams)
hbar = 1.05457173e-34;
kb = 1.3806488e-23;
prefactor = 1/(2*pi)^3;
c = MatParams.vs;
%c = [8000,4000,4000];
cp = c(p);
omega = cp*k;

vg2 = cp^2;
tau = 1./Get_TauInv_kappa_InGaAs(k,p,T,MatParams); %10e-9/cp;%./sin(phi);
x = hbar.*omega/(kb*T);
f = 1./(exp(x)-1);

dfdT = f.*(f+1).*x/T;
hbaromega_dfdT = hbar.*omega.*dfdT;
zero_logic = (x==0); %if omega = 0
hbaromega_dfdT(zero_logic) = kb;

dV = (k.^2);
result = prefactor.*hbaromega_dfdT.*vg2.*tau.*dV;
