function [GammaN] = TruellXSection(omega,a,C11_1,C44_1,rho_1,C11_2,C44_2,rho_2)

% NOTE the nanoparticle is material index #2!!!!

%clear all
% % ------ GaSb ----------
% C11_1 = 88e9; %lambda + 2*mu
% C44_1 = 43.5e9; %mu
% rho_1 =  5.61e3;
% 
% % ------ ErSb -----------
% C11_2 = 159e9;
% C44_2 = 25.8e9;
% rho_2 = 5.61e3;%8.49e3;

Mmax = 200;

%vL = (C11_1/rho_1)^(0.5);
%kcrit = 1/a; %changeover from Rayleigh to Geometric
%omegacrit = kcrit*vL;

Nomega = length(omega);
GammaN = zeros(size(omega));
parfor n=1:Nomega
        w = warning('off','MATLAB:nearlySingularMatrix');
w = warning('off','MATLAB:illConditionedMatrix');
w = warning('off','MATLAB:SingularMatrix');
    %calculate wavevector in two materials
    k1 = omega(n)*(C11_1/rho_1)^(-0.5);
    kap1 = omega(n)*(C44_1/rho_1)^(-0.5);
    k2 = omega(n)*(C11_2/rho_2)^(-0.5);
    kap2 = omega(n)*(C44_2/rho_2)^(-0.5);
    
    k1a = k1*a;
    kap1a = kap1*a;
    k2a = k2*a;
    kap2a = kap2*a;
    
    kapokap = kap1/kap2;
    ror = rho_1/rho_2;
    mom = C44_1/C44_2;
    
    %First calculate the m=0 coeficients
    Eterm = 1/ror*((kap1a)^2/(1-k2a*cot(k2a))-4*(kapokap)^2);
    FE= 4+Eterm;
    num = ((kap1a)^2-FE)*sin(k1a)+FE*k1a*cos(k1a);
    denom = ((kap1a^2-FE)^2+FE^2*k1a^2)^0.5;
    multiplier = exp(1i*(k1a-atan((kap1a^2-FE)/(FE*k1a))));

    A=zeros(Mmax,1);
    B=zeros(Mmax,1);
    C=zeros(Mmax,1);
    D=zeros(Mmax,1);
    A(1) = 1/(k1a)*num/denom*multiplier; % analytic result from Johnson/Truell
    
    %m==0 case
    %create the matrix components
    m=0;
        AA = zeros(8);
        delta = zeros(8,1);
   
        AA(1,1) = k1a*jsphere(m+1,k1a);
        AA(1,2) = k1a*ysphere(m+1,k1a);
        AA(1,3) = 0;
        AA(1,4) = 0;
        AA(1,5) = -k2a*jsphere(m+1,k2a);
        AA(1,6) = 0;
        AA(1,7) = 0;
        AA(1,8) = 0;
        
        AA(3,1) = jsphere(m,k1a);
        AA(3,2) = ysphere(m,k1a);
        AA(3,3) = -((m+1)*jsphere(m,kap1a)-kap1a*jsphere(m+1,kap1a));
        AA(3,4) = -((m+1)*ysphere(m,kap1a)-kap1a*ysphere(m+1,kap1a));
        AA(3,5) = -jsphere(m,k2a);
        AA(3,6) = 0;
        AA(3,7) = (m+1)*jsphere(m,kap2a)-kap2a*jsphere(m+1,kap2a);
        AA(3,8) = 0;
        
        AA(5,1) = kap1a^2*jsphere(m,k1a)-2*(m+2)*k1a*jsphere(m+1,k1a);
        AA(5,2) = kap1a^2*ysphere(m,k1a)-2*(m+2)*k1a*ysphere(m+1,k1a);
        AA(5,3) = 0;
        AA(5,4) = 0;
        AA(5,5) = -1/mom*(kap2a^2*jsphere(m,k2a)-2*(m+2)*k2a*jsphere(m+1,k2a));
        AA(5,6) = 0;
        AA(5,7) = 0;
        AA(5,8) = 0;
        
        AA(7,1) = (m-1)*jsphere(m,k1a)-k1a*jsphere(m+1,k1a);
        AA(7,2) = (m-1)*ysphere(m,k1a)-k1a*ysphere(m+1,k1a);
        AA(7,3) = -((m^2-1-kap1a^2/2)*jsphere(m,kap1a)+kap1a*jsphere(m+1,kap1a));
        AA(7,4) = -((m^2-1-kap1a^2/2)*ysphere(m,kap1a)+kap1a*ysphere(m+1,kap1a));
        AA(7,5) = -1/mom*((m-1)*jsphere(m,k2a)-k2a*jsphere(m+1,k2a));
        AA(7,6) = 0;
        AA(7,7) = 1/mom*((m^2-1-kap2a^2/2)*jsphere(m,kap2a)+kap2a*jsphere(m+1,kap2a));
        AA(7,8) = 0;
        
        for kk = 1:4
            for LL=1:8
                if mod(LL,2)==1 %odd number
                    AA(2*kk,LL) = -AA(2*kk-1,LL+1);
                else  
                    AA(2*kk,LL) = AA(2*kk-1,LL-1);
                end
            end
        end
        for LL = 1:8
            if mod(LL,2)==1 %odd number
                delta(LL) = (-1)^m*(1/k1a)*AA(LL,1);
            else  %even numbers
                delta(LL) = 0;
            end
        end
        
        %Now solve
        X = AA\delta;
        A(m+1) = X(1) + 1i*X(2);
        B(m+1) = X(3) + 1i*X(4);
        C(m+1) = X(5) + 1i*X(6);
        D(m+1) = X(7) + 1i*X(8);
        
        summand = zeros(Mmax,1);
        summand(1) = 4*(abs(A(1)))^2;
    
    %Now do m~=0
    relerr = 1;
    m =0;
    while (relerr > 1e-4)
        m = m+1;
        %create the matrix components
        AA = zeros(8);
        delta = zeros(8,1);
   
        AA(1,1) = k1a*jsphere(m+1,k1a)-m*jsphere(m,k1a);
        AA(1,2) = k1a*ysphere(m+1,k1a)-m*ysphere(m,k1a);
        AA(1,3) = m*(m+1)*jsphere(m,kap1a);
        AA(1,4) = m*(m+1)*ysphere(m,kap1a);
        AA(1,5) = -k2a*jsphere(m+1,k2a)+m*jsphere(m,k2a);
        AA(1,6) = 0;
        AA(1,7) = -m*(m+1)*jsphere(m,kap2a);
        AA(1,8) = 0;
        
        AA(3,1) = jsphere(m,k1a);
        AA(3,2) = ysphere(m,k1a);
        AA(3,3) = -((m+1)*jsphere(m,kap1a)-kap1a*jsphere(m+1,kap1a));
        AA(3,4) = -((m+1)*ysphere(m,kap1a)-kap1a*ysphere(m+1,kap1a));
        AA(3,5) = -jsphere(m,k2a);
        AA(3,6) = 0;
        AA(3,7) = (m+1)*jsphere(m,kap2a)-kap2a*jsphere(m+1,kap2a);
        AA(3,8) = 0;
        
        AA(5,1) = kap1a^2*jsphere(m,k1a)-2*(m+2)*k1a*jsphere(m+1,k1a);
        AA(5,2) = kap1a^2*ysphere(m,k1a)-2*(m+2)*k1a*ysphere(m+1,k1a);
        AA(5,3) = m*(kap1a^2*jsphere(m,kap1a)-2*(m+2)*kap1a*jsphere(m+1,kap1a));
        AA(5,4) = m*(kap1a^2*ysphere(m,kap1a)-2*(m+2)*kap1a*ysphere(m+1,kap1a));
        AA(5,5) = -1/mom*(kap2a^2*jsphere(m,k2a)-2*(m+2)*k2a*jsphere(m+1,k2a));
        AA(5,6) = 0;
        AA(5,7) = -1/mom*m*(kap2a^2*jsphere(m,kap2a)-2*(m+2)*kap2a*jsphere(m+1,kap2a));
        AA(5,8) = 0;
        
        AA(7,1) = (m-1)*jsphere(m,k1a)-k1a*jsphere(m+1,k1a);
        AA(7,2) = (m-1)*ysphere(m,k1a)-k1a*ysphere(m+1,k1a);
        AA(7,3) = -((m^2-1-kap1a^2/2)*jsphere(m,kap1a)+kap1a*jsphere(m+1,kap1a));
        AA(7,4) = -((m^2-1-kap1a^2/2)*ysphere(m,kap1a)+kap1a*ysphere(m+1,kap1a));
        AA(7,5) = -1/mom*((m-1)*jsphere(m,k2a)-k2a*jsphere(m+1,k2a));
        AA(7,6) = 0;
        AA(7,7) = 1/mom*((m^2-1-kap2a^2/2)*jsphere(m,kap2a)+kap2a*jsphere(m+1,kap2a));
        AA(7,8) = 0;
        
        for kk = 1:4
            for LL=1:8
                if mod(LL,2)==1 %odd number
                    AA(2*kk,LL) = -AA(2*kk-1,LL+1);
                else  
                    AA(2*kk,LL) = AA(2*kk-1,LL-1);
                end
            end
        end
        for LL = 1:8
            if mod(LL,2)==1 %odd number
                delta(LL) = (-1)^m*(1/k1a)*AA(LL,1);
            else  %even numbers
                delta(LL) = 0;
            end
        end
        
        %Now solve
        X = AA\delta;
        A(m+1) = X(1) + 1i*X(2);
        B(m+1) = X(3) + 1i*X(4);
        C(m+1) = X(5) + 1i*X(6);
        D(m+1) = X(7) + 1i*X(8);
        summand(m+1) = 4*(2*m+1)*((abs(A(m+1)))^2+m*(m+1)*(k1/kap1)*(abs(B(m+1)))^2);
        relerr = abs(summand(m+1)-summand(m))/summand(m);
    end
    
    %now evaluate for the scattering cross section (normalized by pi*a^2)
    GammaN(n) = sum(summand);
end
        
end
        
        
        
    
    

