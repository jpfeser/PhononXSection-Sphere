clear all
close all
x_GaAs = 0.47;

% T = 300;
% MatParams = PropertiesForInGaAs_ErAs(x_GaAs,T);
% MPFcutoff = logspace(-10,-2,21);
% figure(1)
% for i=1:length(MPFcutoff)
%     MatParams.Lcutoff = MPFcutoff(i);
%     k(i) = get_kappa_sphere(T,MatParams)
%     semilogx(MPFcutoff(1:i),k(1:i),'k-')
%     figure(gcf)
%     pause(1)
% end
% semilogx(MPFcutoff,k,'k-')
% figure(gcf)

figure(2)
T = 300;


MatParams = PropertiesForInGaAs_ErAs(x_GaAs,T);
MatParams.Lb = 500e-6; %boundary scattering / film thickness
% properties of the metal (choose one)

% % Ge
C_NP = [120.6 67.7 67.7]*1e9; % from ioffe
MatParams.rho_NP_Material = 5323; % from ioffe

% % ErAs
% C_NP = [230.5 40.5 40.5]*1e9;
% MatParams.rho_NP_Material = 242.18*4*1.6726e-27/(5.74e-10)^3; %nanoparticle density (8567 kg/m3)

% % Cu
% C_NP = [171 75.6 75.6]*1e9;
% MatParams.rho_NP_Material = 8960; %nanoparticle density (8567 kg/m3)

% % GaIr
%C_NP = [320 62 62]*1e9;
%MatParams.rho_NP_Material = 15.17e3; %nanoparticle density (8567 kg/m3)

 % NiGa
 %C_NP = [173 81.4 81.4]*1e9;
 %MatParams.rho_NP_Material = 8.62e3;

% % NiAl
% C_NP = [207 116 116]*1e9;
% MatParams.rho_NP_Material = 5.92e3;

% % Fe
%C_NP = [231 116 116]*1e9;
%MatParams.rho_NP_Material = 7.87e3;

% % InAs
%C_NP = [83.4 43.2 43.2]*1e9;
%MatParams.rho_NP_Material = 5.68e3;

% GaAs
%C_NP = [96.76 47.34 47.34]*1e9;
%MatParams.rho_NP_Material = 5317;

% % FeAl
% C_NP = [248 137 137]*1e9;
% MatParams.rho_NP_Material = 5.79e3;

% % IrAl
% C_NP = [366 1.05*125 1.05*125]*1e9;
% MatParams.rho_NP_Material = 13.24e3;

% % CoAl
%C_NP = [301 139 139]*1e9;
%MatParams.rho_NP_Material = 6.14e3;

% % Al
% C_NP = [107 28.3 28.3]*1e9;
% MatParams.rho_NP_Material = 2.7e3;

% % W
% C_NP = [523 161 161]*1e9;
% MatParams.rho_NP_Material = 19.27e3;

% % Au
% C_NP = [192 42 42]*1e9;
% MatParams.rho_NP_Material = 19.3e3;

temp = [C_NP(1)/1e9 C_NP(2)/1e9 sqrt(C_NP(1)/MatParams.rho_NP_Material) sqrt(C_NP(2)/MatParams.rho_NP_Material)]
fprintf('& %f & %f & %f & %f\n',temp)
% 
%MatParams.a_NP = 1.5e-9; %nanoparticle radius
MatParams.VolFrac_NP = 0.01; %volume fraction of nanocylinders 

%
MatParams.vs_NP_Material = sqrt(C_NP/MatParams.rho_NP_Material);                                          %change back to 0.05
MatParams.eta_NP = MatParams.VolFrac_NP/(4/3*pi*MatParams.a_NP^3); %number density (#/m3) of nanocylinders.

avect = logspace(-10,-7,41);
%avect = 0.85e-9;%1.05e-9;

deltaC11oC11 = linspace(-0.99, 5, 41);
deltarhoorho = linspace(-0.99, 3, 41);
leg_ent = string(deltaC11oC11);

for i=1:length(deltaC11oC11)
    % readout C11 and v for NP
    C11matrix = MatParams.vs(1)^2*MatParams.rho;
    C11NP = (1+deltaC11oC11(i))*C11matrix;
    C44NP = C11NP/2;
    
    for j = 1:length(deltarhoorho)
        [i,j]
        MatParams_copy = MatParams;
        MatParams_copy.rho_NP_Material = (1+deltarhoorho(j))*MatParams_copy.rho;
        MatParams_copy.vs_NP_Material = sqrt([C11NP,C44NP,C44NP]/MatParams_copy.rho_NP_Material);  
%        options = optimset('TolFun',1e-4,'TolX',1e-4);
%        [amin(i,j),kmin(i,j),exitflag(i,j)]=fminsearch(@(X) kappa_objective_function(X,T,MatParams_copy),3e-9,options)
        
        MatParams_copy.a_NP = 3e-9; %amin(i,j);
        MatParams_copy.eta_NP = MatParams_copy.VolFrac_NP/(4/3*pi*MatParams_copy.a_NP^3);
        
        kmax = MatParams_copy.kmax(1);
        ka = 10^(-3);
        kvect = ka/MatParams_copy.a_NP;
        
        [sigma_L(i,j),scat_eff_L(i,j)] = GetSigmaSphere(kvect,1,MatParams_copy);
        [sigma_T(i,j),scat_eff_T(i,j)] = GetSigmaSphere(kvect,2,MatParams_copy);
        
    end 

    
%     if i>2
%         figure(3)
%         [XX,YY]=meshgrid(deltaC11oC11(1:i),deltarhoorho(:));
%             save('ContrastMap')
%         [c,h]=contour(XX,YY,kmin(1:i,:)',[0.9:0.1:2,2.5,3:4],'k-')
%         axis('equal')
%         clabel(c,h)
%         saveas(gcf,'kmin_contrast_map','epsc')
%         
%         figure(4)
%         [c4,h4]=contour(XX,YY,amin(1:i,:)',[1,2,4,6,8,12,16]*1e-9,'k-')
%         axis('equal')
%         clabel(c4,h4)
%         saveas(gcf,'amin_contrast_map','epsc')
%     end
end

[XX,YY] = meshgrid(deltaC11oC11,deltarhoorho)
 figure(1)
 v = 0.5*2.^(0:4)
[c1,h1] = contour(XX,YY,scat_eff_L'/(ka)^4,v,'k','Linewidth',2)
% loglog(kvect*MatParams_copy.a_NP,scat_eff,'LineWidth',3)
% axis([1e-1 1e2 1e-2 10])
% legend([leg_ent],'Location','SouthEast')
% set(gca,'FontSize',16)

 xlabel('\Delta C_{11}/C_{11}')
 ylabel('\Delta \rho/\rho')
 set(gca,'FontSize',16)
 clabel(c1,h1,'fontsize',16)
% ylabel('\gamma/\pi R^2')
% title('Longitudinal')
% %hold on
% loglog([1e-4 1e4],[2 2],'-.k')
% legend([leg_ent],'Location','SouthEast')


 figure(2)
  v2 = 0.2*2.^(0:4)
[c2,h2] = contour(XX,YY,scat_eff_T'/(ka)^4,v2,'k','Linewidth',2)
  xlabel('{\Delta}C_{44}/C_{44}','fontsize',16)
 ylabel('{\Delta}\rho/\rho','fontsize',16)
 title('Transverse','fontsize',16)
  set(gca,'FontSize',16)
 clabel(c2,h2,'fontsize',16)
% loglog(kvect_T*MatParams_copy.a_NP,scat_eff_T,'LineWidth',3)
% legend([leg_ent],'Location','SouthEast')
% axis([1e-1 1e2 1e-2 10])
% set(gca,'FontSize',16)
% xlabel('ka')
% ylabel('\gamma/\pi R^2')
% title('Transverse')
% %hold on
% loglog([1e-4 1e4],[2 2],'-.k')
% legend([leg_ent],'Location','SouthEast')
