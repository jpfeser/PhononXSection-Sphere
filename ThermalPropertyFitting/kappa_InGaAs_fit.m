clear all
x_GaAs = 0;
T = 300;

xvect = linspace(0,1,100);
for i=1:length(xvect)
    MatParams = PropertiesForInGaAs_ErAs(xvect(i),T);
    MatParams.VolFrac_NP = 0;
    k(i) = get_kappa_sphere(T,MatParams);
end
xdata = [0.16
0.239
0.311
0.362
0.56
0.659
0.776
0.88];
xdata = 1-xdata;

kdata = [8.99847026
7.079646018
5.734602592
5.352459455
4.784688995
5.44751321
6.561249262
7.881462799]; %Dat from Abeles

semilogy(xvect,k,'k-',xdata,kdata,'ro')
figure(gcf)