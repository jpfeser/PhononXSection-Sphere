clear all
x_GaAs = 0.47;
T = 300;
MatParams = PropertiesForInGaAs_ErAs(x_GaAs,T);
volvect = logspace(-4,-1,20);
for i=1:length(volvect)
    MatParams.a_NP = 0.5e-9;
    MatParams.VolFrac_NP = volvect(i)
    MatParams.eta_NP = MatParams.VolFrac_NP/(4/3*pi*MatParams.a_NP^3);
    k1(i) = get_kappa_sphere(T,MatParams)
    semilogx(volvect(1:i),k1(1:i),'k-')
    figure(gcf)
    pause(1)
end
hold on
for i=1:length(volvect)
    MatParams.a_NP = 2.5e-9;
    MatParams.VolFrac_NP = volvect(i)
    MatParams.eta_NP = MatParams.VolFrac_NP/(4/3*pi*MatParams.a_NP^3);
    k2(i) = get_kappa_sphere(T,MatParams)
    semilogx(volvect(1:i),k2(1:i),'k-')
    figure(gcf)
    pause(1)
end
%%
for i=1:length(volvect)
    MatParams.a_NP = 1e-9;
    MatParams.VolFrac_NP = volvect(i)
    MatParams.eta_NP = MatParams.VolFrac_NP/(4/3*pi*MatParams.a_NP^3);
    k3(i) = get_kappa_sphere(T,MatParams)
    semilogx(volvect(1:i),k3(1:i),'k-')
    figure(gcf)
    pause(1)
end
semilogx(volvect,k,'k-')
figure(gcf)