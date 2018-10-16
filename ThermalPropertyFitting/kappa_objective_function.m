function k = kappa_objective_function(X,T,MatParams)

MatParams.a_NP = X;
MatParams.eta_NP = MatParams.VolFrac_NP/(4/3*pi*MatParams.a_NP^3);

k = get_kappa_sphere(T,MatParams);

[X,k];

end