function sigma = GetSigma_SphereComp(k,MatParams)

if nargin==0
    a = 2e-9;
    omega_vect = 2*pi*logspace(9,14,200)';
    
    
    cL_1 = 8500;
    cT_1 = 5000;
    rho_1 = 2200;
    k = omega_vect/cT_1;
    
    cL_2 = 0.7*8500;
    cT_2 = 0.7*5000;
    rho_2 = 2200;
    
else
    %unpack some things
    a = MatParams.a_NP;
    vs=MatParams.vs;
    rho_NP = MatParams.rho_NP_Material;
    vs_NP = MatParams.vs_NP_Material;
    cL_NP = vs_NP(1);
    cT_NP = vs_NP(2);
    
    rho_MATRIX = MatParams.rho;
    cL_MATRIX = vs(1);
    cT_MATRIX = vs(2);
    
    omega_mat = vs(1)*k;
    omega_vect = omega_mat(:);
    omega_zero_logic = (omega_vect==0);
end
C11_NP = cL_NP^2*rho_NP;
C44_NP = cT_NP^2*rho_NP;

C11_MATRIX =cL_MATRIX^2*rho_MATRIX;
C44_MATRIX = cT_MATRIX^2*rho_MATRIX;

[GammaN] = TruellXSection_JPFCorrections(omega_vect,a,C11_MATRIX,C44_MATRIX,rho_MATRIX,C11_NP,C44_NP,rho_NP);

    sigma = GammaN*pi*a^2;
    [n,m]=size(k);
    sigma = reshape(sigma,n,m);   
end