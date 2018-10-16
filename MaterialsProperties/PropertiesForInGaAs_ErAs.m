function MatParams =  PropertiesForInGaAs_ErAs(x_GaAs,T)

if nargin<1
    x_GaAs = 0.5;
elseif nargin<2
    T = 300
end

MatParams.xalloy = x_GaAs; %alloy percentage ,x: perc GaAs in alloy
MatParams.F = 4e-30; %alloy scattering
MatParams.Lb = 2e-6;%200e-6; %boundary scattering / film thickness
MatParams.Lcutoff = Inf; %The longest MFP considered (helps figure out accumulation function)
MatParams.Lambda_cuttoff = Inf; %The longest wavelength considered (helps figure out accumulation function)

% Build Materials Properties of a Matrix / NP system
% Matrix Properties
% InAs
MatParams.InAs_vs = [3830 2640 2640];
MatParams.InAs_omega_max = [31 10.7 10.7]*1e12;%from ioffe
MatParams.AU_InAs = 5.06e-19; %  <---still need to fit (Good for A*omega^2*exp(-B/T), A=5.06e-19, B = 50;
MatParams.BU_InAs = 50; %  <---still need to fit
%MatParams.AN_InAs = 50/300^3;
InAs_F = 1; %  <---still need to fit
MatParams.InAs_rho = 5680; %  <---still need to fit
MatParams.InAs_alat = 6.08e-10;

%GaAs
MatParams.GaAs_vs = [4266 2984 2984]%[5171 3180  2780];
MatParams.GaAs_omega_max = 2*pi*[6.44 2.02 2.02]*1e12;%rad/s from ioffe/Weber (1977)
MatParams.AU_GaAs = 0.9*3.67e-19; %<---(Good for A*omega^2*exp(-B/T), A=3.67e-19, B = 65;)
MatParams.BU_GaAs = 65;
% MatParams.AU_GaAs = 1.45*3.35953483679161e-030; 
% MatParams.BU_GaAs = 360;%<---still need fit
% MatParams.AU_GaAs = MatParams.AU_GaAs/(300*(1-exp(-3*300/MatParams.BU_GaAs)));
% MatParams.AN_GaAs = 1.45*21.1085799254870e-018/(300*(1-exp(-3*300/MatParams.BU_GaAs)));
GaAs_F = 1;%<---still need fit
MatParams.GaAs_rho = 5317; %kg/m3
MatParams.GaAs_alat = 5.74e-10;

%% Alloy Code
%Alloy Properties                                                              
MatParams.vs = MatParams.xalloy.*MatParams.GaAs_vs+(1-MatParams.xalloy).*MatParams.InAs_vs; %sound speeds
MatParams.omega_max = MatParams.xalloy.*MatParams.GaAs_omega_max+(1-MatParams.xalloy).*MatParams.InAs_omega_max; %max frequency
MatParams.alat =  MatParams.xalloy.*MatParams.GaAs_alat+(1-MatParams.xalloy).*MatParams.InAs_alat;
%MatParams.kmax = pi/MatParams.alat;
MatParams.kmax = MatParams.omega_max./MatParams.vs;


MatParams.rho = MatParams.xalloy*MatParams.GaAs_rho+(1-MatParams.xalloy)*MatParams.InAs_rho; %density                                                                             

%% Nanoparticle Properties ErAs
MatParams.a_NP = 1.5e-9; %nanoparticle radius
C_ErAs = [230.5 40.5 40.5]*1e9;
MatParams.rho_NP_Material = 242.18*4*1.6726e-27/(5.74e-10)^3; %nanoparticle density (8567 kg/m3)
MatParams.vs_NP_Material = sqrt(C_ErAs/MatParams.rho_NP_Material);
MatParams.VolFrac_NP = 0.0; %volume fraction of nanocylinders                                             %change back to 0.05
MatParams.eta_NP = MatParams.VolFrac_NP/(4/3*pi*MatParams.a_NP^3); %number density (#/m3) of nanocylinders.
