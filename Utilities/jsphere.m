function [results] = jsphere(m,xi)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
nu = m +1/2;

results = sqrt(pi/2)*besselj(nu,xi)./sqrt(xi);

end

