clear all
x_GaAs = 0;
T = 300;
MatParams = PropertiesForInGaAs_ErAs(x_GaAs,T)
MatParams.VolFrac_NP = 0;
Tvect = linspace(300,800,200);
for i=1:length(Tvect)
    k(i) = get_kappa_sphere(Tvect(i),MatParams);
end
Tdata = [306.88
329.459
355.426
375.032
391.504
417.434
447.257
476.25
505.461
538.436
561.174
604.331
639.47
688.13
725.891];

kdata =[26.728
25.748
24.189
22.637
21.996
20.423
18.422
17.29
16.15
15.191
14.291
12.993
12.136
11.199
10.612];

loglog(Tvect,k,'k-',Tdata,kdata,'ro')
figure(gcf)