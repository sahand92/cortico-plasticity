clear
load('/suphys/sahanda/neurofield/corticothalamic-model/example_parameters.mat');
run('/suphys/sahanda/cortico plasticity/steadyparams.m');

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

 [nuee, nuei]=meshgrid(linspace(0.0003,0.06,10),linspace(0,-0.6,10));
for j=1:8
     G(1,j)=n1.gab(j); %N1 sleep 
   
% end
% Hw=(H0+1i*w.*tp.*H1)./(1+(w.*tp).^2);
% %N2 sleep
% for j=1:8
     G(2,j)=n2.gab(j); %n2
     
     
% end
% %N3 sleep
% for j=1:8
     G(3,j)=n3.gab(j); %n3
     
%N2S sleep

%for j=1:8
    G(4,j)=n2s.gab(j); %n2s
    
    G(5,j)=wgab(j); %wake
    
 %REM
  G(6,j)=rem.gab(j);
  nu(6,j)=rem.nus(j);
  %EO
  G(7,j)=eo.gab(j);
  
   
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
 
 %Gfinal=zeros(12,11);
Ilast=zeros(10,10);
%nufinal=zeros(12,11);   


for c=1:10
 for m=1:10

%paramters wake
G_ee=G(5,1);
G_ei=G(5,2);
G_es=G(5,3);
G_re=G(5,7);
G_rs=G(5,8);
G_sr=G(5,5);
G_sn=G(5,6);
G_se=G(5,4);

%nus:
nu_ee=nuee(c,m);
nu_ei=nuei(c,m);
%nu_ei=nuei(c,2);
% nu_es=nu(c,3);
% nu_re=nu(c,7);
% nu_rs=nu(c,8);
% nu_sr=nu(c,5);
% nu_sn=nu(c,6);
% nu_se=nu(c,4);



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

% nu_ie=nu_ee;
% nu_ii=nu_ei;
% nu_is=nu_es;


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
%     %s_es
%     dses(1)=(1/2*pi)*conj(Hw(1))*Qe(1)*conj(Qs(1));
%     dses(i)=(1/2*pi)*conj(Hw(i))*Qe(i)*conj(Qs(i));
%     Ies(1)=dses(1)*dw;
%     Ies(i)=Ies(i-1)+dses(i)*dw;
%     %s_ie
%     dsie(1)=(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qe(1));
%     dsie(i)=(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qe(i));
%     Iie(1)=dsie(1)*dw;
%     Iie(i)=Iie(i-1)+dsie(i)*dw;
%     %s_ii
%     dsii(1)=(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qi(1));
%     dsii(i)=(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qi(i));
%     Iii(1)=dsii(1)*dw;
%     Iii(i)=Iii(i-1)+dsii(i)*dw;
%     %s_is
%     dsis(1)=(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qs(1));
%     dsis(i)=(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qs(i));
%     Iis(1)=dsis(1)*dw;
%     Iis(i)=Iis(i-1)+dsis(i)*dw;
%     %s_re
%     dsre(1)=(1/2*pi)*conj(Hw(1))*Qr(1)*conj(Qe(1));
%     dsre(i)=(1/2*pi)*conj(Hw(i))*Qr(i)*conj(Qe(i));
%     Ire(1)=dsre(1)*dw;
%     Ire(i)=Ire(i-1)+dsre(i)*dw;
%     %s_rs
%     dsrs(1)=(1/2*pi)*conj(Hw(1))*Qr(1)*conj(Qs(1));
%     dsrs(i)=(1/2*pi)*conj(Hw(i))*Qr(i)*conj(Qs(i));
%     Irs(1)=dsrs(1)*dw;
%     Irs(i)=Irs(i-1)+dsrs(i)*dw;
%     %s_se
%     dsse(1)=(1/2*pi)*conj(Hw(1))*Qs(1)*conj(Qe(1));
%     dsse(i)=(1/2*pi)*conj(Hw(i))*Qs(i)*conj(Qe(i));
%     Ise(1)=dsse(1)*dw;
%     Ise(i)=Ise(i-1)+dsse(i)*dw;
%     %s_sr
%     dssr(1)=(1/2*pi)*conj(Hw(1))*Qs(1)*conj(Qr(1));
%     dssr(i)=(1/2*pi)*conj(Hw(i))*Qs(i)*conj(Qr(i));
%     Isr(1)=dssr(1)*dw;
%     Isr(i)=Isr(i-1)+dssr(i)*dw;
%     %s_sn
%     dssn(1)=(1/2*pi)*conj(Hw(1))*Qs(1)*conj(Qn(1));
%     dssn(i)=(1/2*pi)*conj(Hw(i))*Qs(i)*conj(Qn(i));
%     Isn(1)=dssn(1)*dw;
%     Isn(i)=Isn(i-1)+dssn(i)*dw;
    
end

%G_ab after adding terms to make it the same dimension as Imax/sab's
last=length(w);

Gfinal(c,1:11)=[G_ee,G_ei];

rhos=Gfinal./nufinal;
Ilast(c,m)=[Iee(last),Iei(last)];
dGdt=rhos.*Ilast; % 5x11 matrix of dG/dts

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
end
 radius=sqrt(real(dXdt).^2 + real(dYdt).^2);
 dUdt=(real(dXdt))./real(radius);
 dVdt=(real(dYdt))./real(radius);
 
 quiver(real(X),real(Y),dUdt,dVdt,0.25)
 hold on
 x=linspace(0.25,1,10);
 y=1-x;
% 
 plot(x,y,'--','color','green')
 hold on
scatter(X,Y)

%quiver(Gfinal(:,1),Gfinal(:,2),real(dGdt(:,1))./sqrt(real(dGdt(:,1)).^2+real(dGdt(:,2).^2)),real(dGdt(:,2))./sqrt(real(dGdt(:,1)).^2+real(dGdt(:,2)).^2),0.25)


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

