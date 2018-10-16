function [results] = hsphere(m,xi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nu = m+1/2;
results = sqrt(pi/2)*besselh(nu,1,xi)./sqrt(xi); 

end

