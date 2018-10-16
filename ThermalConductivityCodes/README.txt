get_kappa_sphere:  is the top-level function for calculating thermal conductivity of stuff with embedded nanoparticles.
thermal_integrand_sphere_InGaAs: is a function that evaluates the thermal conductivity integrand.
Get_TauInv_kappa_InGaAs: function, calculates all the scattering rates (including the nanoparticle scattering rate).

as part of "Get_TauInv_kappa_InGaAs" the scattering cross sections of the nanoparticles are calculated using "GetSigma_Sphere" which can be found in another folder.

