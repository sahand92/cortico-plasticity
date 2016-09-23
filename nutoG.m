function [G]=nutoG(nu)
S=load('/suphys/sahanda/cortico plasticity/data/pdb_all.mat');

nuee=S.nus_final(1,1);
nuei=S.nus_final(1,2);
nues=S.nus_final(1,3);
nuse=S.nus_final(1,4);
nusr=S.nus_final(1,5);
nusn=S.nus_final(1,6);
nure=S.nus_final(1,7);
nurs=S.nus_final(1,8);

%Expression for rho:
Qee=0.6832;   %this is the zero of the steadyphie function
phiee=Qee;
phir=sig(nure*phiee + (nurs/nues).*(invsig(phiee)-(nuee+nuei).*phiee));
%phis=sig(nuse*phie+nusr*phir+nusn*phin);
phis=(1/nues).*(invsig(phiee)-(nuee+nuei)*phiee);
Qr=phir;
Qs=phis;
rhoe=(Qee./3.8e-3).*(1-(Qee/340));
rhor=(Qr./3.8e-3).*(1-(Qr/340));
rhos=(Qs./3.8e-3).*(1-(Qs/340));
G=[rhoe, rhoe, rhoe, rhor, rhor, rhor, rhos, rhos].*nu;
