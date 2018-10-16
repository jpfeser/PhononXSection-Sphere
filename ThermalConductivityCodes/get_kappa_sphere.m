function [ktot,kp] = get_kappa_sphere(T,MatParams)
if nargin==0
    PropertiesForSiGe
    T = 300;
    MatParams.eta_NP =0;
end

kp = zeros(3,1);

for p = 1:3  %all polarizations
    integrand = @(k,theta,phi) (k.^2.*abs(sin(phi)));
    kmax = MatParams.kmax(p);
    anglefactor = 2*pi*integral(@(phi) abs(sin(phi)).*cos(phi).^2,0,pi);
    kp(p,1)=anglefactor*integral(@(k) thermal_integrand_sphere_InGaAs(k,p,T,MatParams),0,kmax,'RelTol',1,'AbsTol',1e-8);
end

ktot = sum(kp);
%         sum(k_p)                              %change back if needed
end
