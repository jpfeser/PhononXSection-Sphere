function TauInv_tot = Get_TauInv_kappa_InGaAs(k,p,T,MatParams)

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
Lcutoff = MatParams.Lcutoff; %min MFP to consider (set tau_inv = Inf for stuff above that)
Wavecutoff = MatParams.Lambda_cuttoff;

%phonon-phonon scattering parameters:  tau^(1) = A*T*omega^2*exp(-B/T)
InAs_AU=MatParams.AU_InAs; 
InAs_BU=MatParams.BU_InAs; 
GaAs_AU=MatParams.AU_GaAs;
GaAs_BU=MatParams.BU_GaAs;
%InAs_AN=MatParams.AN_InAs;
%GaAs_AN=MatParams.AN_GaAs;

% phonon-impurity scattering tau^(1) = F*(1-x)*x*omega^4/vs^3; 
F=MatParams.F;
xalloy = MatParams.xalloy;

%phonon-boundary scattering tau^(1) = vs*Lb;
Lb=MatParams.Lb; %Boundary Scattering Length

%phonon-np scattering:->  complicated function of (k, theta, phi, p)

%% Calculate the scattering times
omega = vs(p)*k;
TauInv_U_InAs = InAs_AU*T*omega.^2*exp(-InAs_BU/T);
TauInv_U_GaAs = GaAs_AU*T*omega.^2*exp(-GaAs_BU/T);
% TauInv_U_InAs = InAs_AU*(T*(1-exp(-3*T/MatParams.BU_InAs)))*(omega.^3); %old InAs_A*T*omega.^2*exp(-InAs_B/T)
% TauInv_U_GaAs = GaAs_AU*(T*(1-exp(-3*T/MatParams.BU_GaAs)))*(omega.^3);
TauInv_U = xalloy.*TauInv_U_GaAs+(1-xalloy).*TauInv_U_InAs;

% TauInv_N_InAs = InAs_AN*(T*(1-exp(-3*T/MatParams.BU_InAs)))*(omega.^2); %old InAs_A*T*omega.^2*exp(-InAs_B/T)
% TauInv_N_GaAs = GaAs_AN*(T*(1-exp(-3*T/MatParams.BU_GaAs)))*(omega.^2);
% TauInv_N = xalloy.*TauInv_N_GaAs+(1-xalloy).*TauInv_N_InAs;

% TauInv_phph = GaAs_A*T.*omega.^3*exp(-GaAs_B/T);
TauInv_alloy = F*xalloy*(1-xalloy)*(omega.^4)/(vs(p)^3);
TauInv_boundary = vs(p)/(Lb);
TauInv_np = Get_TauInv_sphere(k,p,MatParams);

TauInv_tot = TauInv_U + TauInv_alloy + TauInv_boundary + TauInv_np; %TauInv_N?

TauInv_min = vs(p)/Lcutoff;
cutoff_logic = (TauInv_tot<TauInv_min); %if scat rate is higher than this, set scat rate to Inf;
TauInv_tot(cutoff_logic) = Inf;


cutoff_logic = (k<(2*pi./Wavecutoff)); %if scat rate is higher than this, set scat rate to Inf;
TauInv_tot(cutoff_logic) = Inf;


end

