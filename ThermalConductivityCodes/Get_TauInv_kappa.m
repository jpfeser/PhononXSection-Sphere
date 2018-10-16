function TauInv_tot = Get_TauInv_kappa(k,p,T,MatParams)

%%
if nargin<3
    T = 300;
end
if nargin<4 %if MatParams not provided in function call
    %then we need to input something...lets use properties for Si
    MatParams.vs = [3970 2770 2770];
    MatParams.A=876e-021; %Anharmonic scattering terms
    MatParams.B=85.6; %Anharmonic scattering terms
    MatParams.F = 2e-30; 
    MatParams.Lb = 300e-6;
    MatParams.xalloy = 0;
end

%% Unpackage the necessary materials parameters (stored as an object):
vs=MatParams.vs; %Sound velocity vector [vL vT1 vT2]

%phonon-phonon scattering parameters:  tau^(1) = A*T*omega^2*exp(-B/T)
Si_A=MatParams.A_Si; 
Si_B=MatParams.B_Si; 
Ge_A=MatParams.A_Ge;
Ge_B=MatParams.B_Ge;

% phonon-impurity scattering tau^(1) = F*(1-x)*x*omega^4/vs^3; 
F=MatParams.F;
xalloy = MatParams.xalloy;

%phonon-boundary scattering tau^(1) = vs*Lb;
Lb=MatParams.Lb; %Boundary Scattering Length

%phonon-np scattering:->  complicated function of (k, theta, phi, p)

%% Calculate the scattering times
omega = vs(p)*k;
TauInv_phph_Si = Si_A*T.*(omega.^2)*exp(-Si_B/T);
TauInv_phph_Ge = Ge_A*T.*(omega.^2)*exp(-Ge_B/T);
TauInv_phph = xalloy.*TauInv_phph_Ge+(1-xalloy).*TauInv_phph_Si;
% TauInv_phph = Ge_A*T.*omega.^3*exp(-Ge_B/T);
TauInv_alloy = F*xalloy*(1-xalloy)*(omega.^4)/(vs(p)^3);
TauInv_boundary = vs(p)/Lb;
TauInv_np = Get_TauInv_sphere(k,p,MatParams);

TauInv_tot = TauInv_phph + TauInv_alloy + TauInv_boundary + TauInv_np;

end

