function sigma = GetSigma_SphereTrans(k,MatParams)

%calculation is based on:
%Iwashimazu, Journal of Sound and Vibration, 40(2), 267-271, 1975
%
% NOTE: subscript 2 is the Nanoparticle in this code!!!
%
if nargin==0
    a = 2e-9;
    omega_vect = 2*pi*logspace(9,14,200)';
    
    
    cL_1 = 8500;
    cT_1 = 5000;
    rho_1 = 2200;
    k = omega_vect/cT_1
    
    cL_2 = 0.7*8500;
    cT_2 = 0.7*5000;
    rho_2 = 2200;
    
else
    %unpack some things
    a = MatParams.a_NP;
    vs=MatParams.vs;
    rho_2 = MatParams.rho_NP_Material;
    vs_NP = MatParams.vs_NP_Material;
    cL_2 = vs_NP(1);
    cT_2 = vs_NP(2);
    
    rho_1 = MatParams.rho;
    cL_1 = vs(1);
    cT_1 = vs(2);
    
    omega_mat = vs(2)*k;
    omega_vect = omega_mat(:);
    omega_zero_logic = (omega_vect==0);
end




C11_1 = cL_1^2*rho_1;
C44_1 = cT_1^2*rho_1;

C11_2 =cL_2^2*rho_2;
C44_2 = cT_2^2*rho_2;


mu_1 = C44_1;
lambda_1 = C11_1 - 2*mu_1;
mu_2 = C44_2;
lambda_2 = C11_2 - 2*mu_2;

w = warning('off','MATLAB:nearlySingularMatrix');
w = warning('off','MATLAB:illConditionedMatrix');
w = warning('off','MATLAB:SingularMatrix');
%QT1n = zeros(200,1);
%QT2n =zeros(200,1);
%QLn =zeros(200,1);
for nomega = 1:length(omega_vect)
    nomega;
    Ka_L_1 = omega_vect(nomega)/cL_1*a;
    Ka_T_1 = omega_vect(nomega)/cT_1*a;

    Ka_L_2 = omega_vect(nomega)/cL_2*a;
    Ka_T_2 = omega_vect(nomega)/cT_2*a;
    
    p = mu_2/mu_1;
    kappa = Ka_L_1/Ka_T_1;
    
    QTot_loop = 0;
    relerr= Inf;
    n=0;
    while relerr > 1e-13;
        n = n+1;
  
        J_ka_1 = Ka_T_1*jsphere(n+1,Ka_T_1)/jsphere(n,Ka_T_1);
        H_ka_1 = Ka_T_1*hsphere(n-1,Ka_T_1)/hsphere(n,Ka_T_1);
        J_Ka_1 = Ka_L_1*jsphere(n+1,Ka_L_1)/jsphere(n,Ka_L_1);
        H_Ka_1 = Ka_L_1*hsphere(n-1,Ka_L_1)/hsphere(n,Ka_L_1);
        J_ka_2 = Ka_T_2*jsphere(n+1,Ka_T_2)/jsphere(n,Ka_T_2);
        H_ka_2 = Ka_T_2*hsphere(n-1,Ka_T_2)/hsphere(n,Ka_T_2);
        J_Ka_2 = Ka_L_2*jsphere(n+1,Ka_L_2)/jsphere(n,Ka_L_2);
        H_Ka_2 = Ka_L_2*hsphere(n-1,Ka_L_2)/hsphere(n,Ka_L_2);
        
        
        A1 = [-(n+1) + H_Ka_1;
            1
            2*(n+1)*(n+2)-((Ka_T_1)^2+4*H_Ka_1);
            -2*(n+2)+2*H_Ka_1];
        
        A2 = [-(n+1);
            1-H_ka_1/n;
            2*(n+1)*(n+2)-2*(n+1)*H_ka_1;
            -2*(n+2)+((Ka_T_1)^2+2*H_ka_1)/n];
        
        A3 = [n-J_Ka_1;
            1;
            2*n*(n-1)*p + (-Ka_T_2^2 + 4*J_Ka_2)*p;
            2*(n-1)*p - 2*J_Ka_2*p];
        
        A4 = [n;
            1-J_ka_2/(n+1);
            2*n*(n-1)*p-2*n*J_ka_2*p;
            2*(n-1)*p+(-Ka_T_2^2 + 2*J_ka_2)*p/(n+1)];
        
        B = [1;
            1/n - J_ka_1/(n*(n+1));
            2*(n-1)-2*J_ka_1;
            2*(n-1)/n-((Ka_T_1)^2-2*J_ka_1)/(n*(n+1))];
        
        C1 = [1;
            -(n+2) + H_ka_1];
        
        C2 = [1;
            (n-1)*p - p*J_ka_2];
        
        D = [1;
            (n-1)-J_ka_1];
        
        detA = det([A1,A2,A3,A4]);
        detC = det([C1,C2]);
        phi_1 = 1i^(n+1)*(2*n+1)*(jsphere(n,Ka_T_1)/(Ka_T_1*hsphere(n,Ka_L_1)))*det([B,A2,A3,A4])/detA;
        psi_1 = 1i^(n+1)*(2*n+1)/n*(jsphere(n,Ka_T_1)/(Ka_T_1*hsphere(n,Ka_T_1)))*det([A1,B,A3,A4])/detA;
        omega_1 = 1i^(n)*(2*n+1)/(n*(n+1))*(jsphere(n,Ka_T_1)/(hsphere(n,Ka_T_1)))*det([D,C2])/detC;
        
        phi_2 = -1i^(n+1)*(2*n+1)*(jsphere(n,Ka_T_1)/(Ka_T_1*jsphere(n,Ka_L_2)))*det([A1,A2,B,A4])/detA;
        psi_2 = 1i^(n+1)*((2*n+1)/n+1)*(jsphere(n,Ka_T_1)/(Ka_T_1*jsphere(n,Ka_T_2)))*det([A1,A2,A3,B])/detA;
        omega_2 = -1i^(n)*(2*n+1)/(n*(n+1))*(jsphere(n,Ka_T_1)/(jsphere(n,Ka_T_2)))*det([C1,D])/detC;
        
        
        pref = 2*n*(n+1)/(2*n+1);
        QLn(n,1) = pref*(Ka_T_1/Ka_L_1)*phi_1*conj(phi_1); 
        QT1n(n,1) = pref*n*(n+1)*psi_1*conj(psi_1);
        QT2n(n,1) = pref*n*(n+1)*omega_1*conj(omega_1)/(Ka_T_1)^2;
        
        QTotprev = QTot_loop;
        QTot_loop = QTot_loop + (QT1n(n,1) + QT2n(n,1) + QLn(n,1));
        relerr = (QTot_loop-QTotprev)/QTotprev;
    end
    QT1(nomega,1)=sum(real(QT1n));
    QT2(nomega,1)=sum(real(QT2n));
    QL(nomega,1)=sum(real(QLn));
    
        
    
end
QTot = QT1 + QT2 + QL;

[n,m]=size(k);
sigma = reshape(QTot*pi*a^2,n,m);   

if nargin==0
    figure(1);
    loglog(omega_vect/cT_1*a,QL,omega_vect/cT_1*a,QT1,omega_vect/cT_1*a,QT2,omega_vect/cT_1*a,QTot)
    legend('Longitudinal','Transverse 1','Transverse 2','Total')
    
    figure(2)
    plot(cumsum(real(QT2n)))
    
    figure(3);
    plot(omega_vect/cT_1*a,QL,omega_vect/cT_1*a,QT1,omega_vect/cT_1*a,QT2,omega_vect/cT_1*a,QTot)
    legend('Longitudinal','Transverse 1','Transverse 2','Total')
    axis([0 pi*a/0.27e-10 0 4])
end

end