
function [Q]=sig(V)
Q=340./(1+exp(-(V-3.4*3.8e-3)./3.8e-3));

