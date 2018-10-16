function [sigma,scat_eff] = GetSigmaSphere(k,p,MatParams)

if p==1
    sigma = GetSigma_SphereComp(k,MatParams);
else
    sigma = GetSigma_SphereTrans(k,MatParams);
end
%    % geometric scattering only!!!  Uncomment for Debug only!
%    sigma = 2*(pi*MatParams.a_NP^2)*ones(size(k));

%    % density based pertubation theory
% Crude (incorrect) model, good for comparing. Uncomment for Debug only!
%drho = (MatParams.rho - MatParams.rho_NP_Material)/MatParams.rho
%sigma_LW = pi*MatParams.a_NP^2*(4/9)*drho^2*(k*MatParams.a_NP).^4;
%sigma_SW = 2*pi*MatParams.a_NP^2;
%sigma = 1./(1./sigma_LW+1./sigma_SW);

    scat_eff = sigma/(pi*MatParams.a_NP^2);
end