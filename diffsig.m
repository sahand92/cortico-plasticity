function [sigprime]=diffsig(Q)
sigprime=(Q/3.8e-3).*(1-(Q/340));

