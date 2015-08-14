clear
load('/suphys/sahanda/neurofield/corticothalamic-model/example_parameters.mat');
run('/suphys/sahanda/steadyparams.m');

%plasticity window
tp=0.01; %plasticity timescale

A_plus=1;
%A_minus=1; %CDP
A_minus=-1;%STDP
H0=(A_plus + A_minus)*tp; 
H1=(A_plus - A_minus)*tp;

fmax=45;
nw=1000;
f = linspace(0,fmax,nw);
w=2*pi*f;
Hw=(H0+1i*w.*tp*H1)./(1+w*tp).^2;
dw=w(2)-w(1);
%G_ab's 5 states, 8 parameters
G=zeros(12,8);
nu=zeros(12,8);
% %wake
% for j=1:8
%     G(j)=wgab(j); %wake 
% end
%N1 sleep
wgab = [2.074,-4.110,0.772,7.768,-3.301,8.097,0.656,0.196];

   
for j=1:8
     G(1,j)=n1.gab(j); %N1 sleep 
     nu(1,j)=n1.nus(j);
    
% end
% Hw=(H0+1i*w.*tp.*H1)./(1+(w.*tp).^2);
% %N2 sleep
% for j=1:8
     G(2,j)=n2.gab(j); %n2
     nu(2,j)=n2.nus(j);
     
% end
% %N3 sleep
% for j=1:8
     G(3,j)=n3.gab(j); %n3
     nu(3,j)=n3.nus(j);
% end
%N2S sleep

%for j=1:8
    G(4,j)=n2s.gab(j); %n2s
    nu(4,j)=n2s.nus(j);
    G(5,j)=wgab(j); %wake
    nu(5,j)=ec.nus(j);
 %REM
  G(6,j)=rem.gab(j);
  nu(6,j)=rem.nus(j);
  %EO
  G(7,j)=eo.gab(j);
  nu(7,j)=eo.gab(j);
   
end
%steadyparams
G(8,:)=Gs(1,:);
G(9,:)=Gs(2,:);
G(10,:)=Gs(3,:);
G(11,:)=Gs(4,:);
G(12,:)=Gs(5,:);
nu(8,:)=nus(1,:);
nu(9,:)=nus(2,:);
nu(10,:)=nus(3,:);
nu(11,:)=nus(4,:);
nu(12,:)=nus(5,:);
 %XYZs
 


for c=1:11


%paramters wake
G_ee=G(c,1);
G_ei=G(c,2);
G_es=0;
G_re=0;
G_rs=0;
G_sr=0;
G_sn=0;
G_se=0;

%nus:
nu_ee=nu(c,1);
nu_ei=nu(c,2);
nu_es=0;
nu_re=0;
nu_rs=0;
nu_sr=0;
nu_sn=0;
nu_se=0;



% %paramters wake
% G_ee=G(5,1);
% G_ei=G(5,2);
% G_es=G(5,3);
% G_re=G(5,7);
% G_rs=G(5,8);
% G_sr=G(5,5);
% G_sn=G(5,6);
% G_se=G(5,4);
% end

%Random connectivity assumption
G_ie=G_ee;
G_ii=G_ei;
G_is=G_es;

nu_ie=nu_ee;
nu_ii=nu_ei;
nu_is=nu_es;


G_ese=G_es*G_se;
%G_ese=5.9943;
G_erse=G_es*G_sr*G_re;
%G_erse=-1.6712;
G_srs=G_sr*G_rs;
%G_srs=-0.6474;
G_esn=G_es*G_sn;





alpha=83.3;
beta=770;
gamma_ee=116;
L=((1-1i*w./alpha).*(1-1i*w./beta)).^-1;

t_0=0.085;
tau_re=t_0/2;
tau_se=t_0/2;
tau_es=t_0/2;
%XYZs
X(c)=G_ee/(1-G_ei);
Y(c)=(G_ese+G_erse)/((1-G_srs)*(1-G_ei));
Z(c)=-G_srs*(alpha*beta)/(alpha+beta)^2;


r_ee=0.086;
Gamma_e=(1-1i*w/gamma_ee).^2+1;
J_ee=G_ee*L.*Gamma_e;
J_ei=G_ei*L;
J_es=G_es*L.*exp(1i*w.*tau_es);
J_re=G_re*L.*exp(1i*w.*tau_re);
J_se=G_se*L.*exp(1i*w.*tau_se);
J_rs=G_rs*L;
J_sr=G_sr*L;
J_sn=G_sn*L;

phi_n=0.001;
Ns=phi_n*G_sn*L;
A=(1-J_ee-J_ei).*(1-J_sr.*J_rs)-J_es.*(J_se+J_sr.*J_re);


Qe=(Ns./A).*(J_ei.*J_es+J_es-J_ee.*J_ei);
Qi=(Ns./A).*(J_es.*J_es+J_es-J_es.*J_ee);
Qr=(Ns./A).*(J_ei.*J_ee+J_es.*J_re-J_es.*J_re.*J_ei+J_rs-J_rs.*J_ei-J_rs.*J_ee+J_rs.*J_ee.*J_ei+J_ei.*J_es.*J_re);
Qs=(Ns./A).*(1+J_ee.*J_ei-J_ee-J_ei-J_ei.*J_ee);
Qn=phi_n*ones(1,length(w));
Q=[Qe;Qi;Qr;Qs];
phi_e=Qe./Gamma_e;
for i=2:length(w)
    %_see
    dsee(1)=(1/2*pi)*conj(Hw(1))*Qe(1)*conj(Gamma_e(1))*conj(Qe(1));
    dsee(i)=(1/2*pi)*conj(Hw(i))*Qe(i)*conj(Gamma_e(i))*conj(Qe(i));
    Iee(1)=dsee(1)*dw;
    Iee(i)=Iee(i-1)+dsee(i)*dw;
    %s_ei
    dsei(1)=(1/2*pi)*conj(Hw(1))*Qe(1)*conj(Qi(1));
    dsei(i)=(1/2*pi)*conj(Hw(i))*Qe(i)*conj(Qi(i));
    Iei(1)=dsei(1)*dw;
    Iei(i)=Iei(i-1)+dsei(i)*dw;
 
    
end

%G_ab after adding terms to make it the same dimension as Imax/sab's
last=length(w);

Gfinal(c,1:2)=[G_ee,G_ei];
nufinal(c,1:2)=[nu_ee,nu_ei,nu_es,nu_ie,nu_ii,nu_is,nu_re,nu_rs,nu_se,nu_sr,nu_sn];
rhos=Gfinal./nufinal;
Ilast(c,1:2)=[Iee(last),Iei(last);
dGdt==rhos.*Ilast; % 5x11 matrix of dG/dts

dG_ee=dGdt(c,1);
dG_ei=dGdt(c,2);
dG_es=dGdt(c,3);
dG_se=dGdt(c,9);
dG_re=dGdt(c,7);
dG_sr=dGdt(c,10);
dG_rs=dGdt(c,8);
dG_ese=dG_es*G_se + G_es*dG_se;
dG_erse=G_es*G_sr*dG_re + G_re*(G_es*dG_sr + dG_es*G_sr);
dG_srs=dG_sr*G_rs + dG_rs*G_sr;


%dXdt(c)=((1-G_ei).*dGdt(c,1)+G_ee.*dGdt(c,2))./(1-G_ei).^2;
dXdt(c)=((1-G_ei)*dG_ee + G_ee*dG_ei)/(1-G_ei)^2;
%dYdt(c)=((1-G_srs).*(1-G_ei).*(dGdt(c,3)*G_se + G_es*dGdt(c,9)+G_es*G_sr*dGdt(c,7)+G_re*(G_es*dGdt(c,10)+dGdt(c,3)*G_sr))-(G_ese+G_erse).*(-(dGdt(c,10)*G_rs+G_sr*dGdt(c,8))-dGdt(c,2)+(dGdt(c,10)*G_rs+G_sr*dGdt(c,8))*G_ei+G_srs*dGdt(c,2)))./((1-G_srs).*(1-G_ei))^.2;
dYdt(c)=((1-G_srs)*(1-G_ei)*(dG_ese+dG_erse)-(G_ese+G_erse)*(-dG_srs-dG_ei+dG_srs*G_ei+G_srs*dG_ei))/((1-G_srs)*(1-G_ei))^2;
dZdt(c)=-(dGdt(c,10)*G_rs+G_sr*dGdt(c,8))*(alpha*beta)/(alpha+beta)^2;

end


% quiver(X,Y,dXdt,dYdt,0.5)
% hold on
% x=linspace(0.25,1,10);
% y=1-x;
% plot(x,y,'--','color','green')
% hold on
% scatter(X,Y)
%quiver(Gfinal(:,1),Gfinal(:,2),dGdt(:,1),dGdt(:,2))
quiver(Gfinal(:,1),Gfinal(:,2),dGdt(:,1),dGdt(:,2))

% sigma=0.006;
% theta=0.0204;
% Qmax=340;
% 
% rho_a=Q/sigma .*(1-Q/Qmax);
 
% %vab=G(1,:)./rho_a;
% Q_0=linspace(0,340,1000);
% %Q_x=-sigma*log((Qmax./Q_0) -1)+theta-Q_0;
% %plot(Q_x,Q_0)
% plot(Q_0,-sigma*log((Qmax./Q_0)-1)+theta-Q_0)
%XYZ space
