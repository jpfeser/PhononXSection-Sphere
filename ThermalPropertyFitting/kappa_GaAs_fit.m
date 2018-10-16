clear all
x_GaAs = 1;
T = 300;
MatParams = PropertiesForInGaAs_ErAs(x_GaAs,T)
MatParams.VolFrac_NP = 0;
Tvect = linspace(200,600,400);
for i=1:length(Tvect)
    k(i) = get_kappa_sphere(Tvect(i),MatParams);
end
Tdata = [196.309
238.79
266.431
290.75
311.08
343.897
352.963
443.653
461.813
496.2
523.36
594.294];

kdata =[82.144
63.011
55.236
52.383
47.525
38.785
36.754
29.301
27.986
25.485
23.749
20.884];

loglog(Tvect,k,'k-',Tdata,kdata,'ro')
figure(gcf)