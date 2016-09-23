function [V]=invsig(phi)
V=0.0204-0.006*log((340-phi)./phi);
