clear
S=load('/suphys/sahanda/cortico plasticity/data/pdb_sleep.mat');
X=S.xyz_final(:,1);
Y=S.xyz_final(:,2);
Z=S.xyz_final(:,3);

% U=transpose(S.u);
% V=transpose(S.v);
% W=transpose(S.w);
%scatter(X,Y,'.')
% newX=X(X<0.2&X>0);
% newY=Y(Y>0.4&Y<0.8);
% M=[X Y];
% 
% newX=X(M(:,1)<0.2 & M(:,2)>0.4)
% newY=Y(M(:,2)>0.4 & M(:,1)<0.2) 
% 
% 
% 
% %blocks (squares) are M
% M=zeros(size(M))
% M=[newX newY];
n=0.03;  %square size is n x n
%index=zeros(400,2);

for i=0:ceil(0.9/n)
    for j=0:ceil(1.8/n)
   
l=find(X>0+n*i&X<0+n*(i+1) & Y>-1+n*j&Y<-1+n*(j+1),1);
if isempty(l)==1;
     l=0;  
 end 
 %disp(l)
 
%  if l~=[];
%      x=l;
%  else
%      x=0;
%  end
index(i+1,j+1)=l;

% M2=[X(l) Y(l)];
% scatter(X(l),Y(l),'.')
% hold on
  end  

end
a=numel(index);
index=reshape(index,1,a);
index=index(index~=0);
clear i
clear j

tp=0.01; %plasticity timescale

A_plus=1;
%A_minus=1; %CDP
A_minus=-0.4;%STDP

H0=-(A_plus + A_minus)*tp; 
H1=(A_plus - A_minus)*tp;

 fmax=45;
 nw=1000;
 f = linspace(0,fmax,nw);
 w=2*pi*f;
 Hw=(H0+1i*w.*tp*H1)./(1+(w*tp).^2);


 %Window function for RS-SR population------------------------
 
 A_plusrs=1;
%A_minus=1; %CDP
A_minusrs=-0.2;%STDP

H0rs=-(A_plusrs + A_minusrs)*tp; 
H1rs=(A_plusrs - A_minusrs)*tp;

 Hwrs=(H0rs+1i*w.*tp*H1rs)./(1+(w*tp).^2);
 
 
 
 %triphasic H(w) ------------------------------------------------


% a=0.00002;
% b=0.0002;
% Aminus=-0.1;
% Apos=0.25;
% alphaH=0.005;
% 
% Hw=sqrt(pi).*exp(-1i*w.*alphaH).*(Apos.*sqrt(a).*exp(-0.25.*w.^2 .*a) + Aminus.*sqrt(b).*exp(-0.25.*w.^2 .*b));


%H(w) STDP + CDP ------------------------------------------------
% a=1;
% b=-0.5;
% tstpd=0.01;
% tcdp=0.01;
% Hw=a.*(2.*1i.*w.*tstpd.^2)./(1+(w.*tstpd).^2)+ (b.*(2*tcdp)./(1+(w.*tcdp).^2));

dw=w(2)-w(1);
%G_ab's 5 states, 8 parameters
% G=zeros(30,8);
% nu=zeros(30,8);


   

     G=S.gab_final(index,:);  
     nu=S.nus_final(index,:);
   


 %XYZs
 
 Gfinal=zeros(length(G(:,1)),11);
Ilast=zeros(length(G(:,1)),11);
nufinal=zeros(length(G(:,1)),11);   
N=zeros(length(G(:,1)),11);

for c=1:length(G(:,1))


%paramters wake
G_ee=G(c,1);
G_ei=G(c,2);
G_es=G(c,3);
G_re=G(c,7);
G_rs=G(c,8);
G_sr=G(c,5);
G_sn=G(c,6);
G_se=G(c,4);

%nus:
nu_ee=nu(c,1);
nu_ei=nu(c,2);
nu_es=nu(c,3);
nu_re=nu(c,7);
nu_rs=nu(c,8);
nu_sr=nu(c,5);
nu_sn=nu(c,6);
nu_se=nu(c,4);



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
% X(c)=G_ee/(1-G_ei);
% Y(c)=(G_ese+G_erse)/((1-G_srs)*(1-G_ei));
% Z(c)=-G_srs*(alpha*beta)/(alpha+beta)^2;


r_ee=0.086;
Gamma_e=(1-1i*w/gamma_ee).^-2;
% J_ee=G_ee*L.*Gamma_e;
% J_ei=G_ei*L;
% J_es=G_es*L.*exp(1i*w.*tau_es);
% J_re=G_re*L.*exp(1i*w.*tau_re);
% J_se=G_se*L.*exp(1i*w.*tau_se);
% J_rs=G_rs*L;
% J_sr=G_sr*L;
% J_sn=G_sn*L;
% 
 phi_n=0.001;
% Ns=phi_n*G_sn*L;
% A=(1-J_ee-J_ei).*(1-J_sr.*J_rs)-J_es.*(J_se+J_sr.*J_re);
% 
% 
% Qe=(Ns./A).*(J_ei.*J_es+J_es-J_ee.*J_ei);
% 
% Qi=(Ns./A).*(J_es.*J_es+J_es-J_es.*J_ee);
% Qr=(Ns./A).*(J_ei.*J_ee+J_es.*J_re-J_es.*J_re.*J_ei+J_rs-J_rs.*J_ei-J_rs.*J_ee+J_rs.*J_ee.*J_ei+J_ei.*J_es.*J_re);
% Qs=(Ns./A).*(1+J_ee.*J_ei-J_ee-J_ei-J_ei.*J_ee);
% Qn=phi_n*ones(1,length(w));
% Q=[Qe;Qi;Qr;Qs];
% phi_e=Qe./Gamma_e;

%expression for phi's:
phi_e=((L.^2 *G_esn.*exp(1i.*w.*tau_es).*phi_n)./(1-L.^2 .*G_srs)) ./ (Gamma_e.*(1-G_ei.*L)-L.*G_ee -(L.^2 .*G_ese.*exp(1i*w*(tau_es+tau_se)) +L.^3.*G_erse.*exp(1i*w*(tau_es+tau_re)))./(1-L.^2.*G_srs));
phi_s=((L.*G_se.*phi_e.*exp(1i*w*tau_se)+L.^2.*G_sr.*G_re.*phi_e.*exp(1i*w*tau_re)) +L.*G_sn.*phi_n) ./ (1-L.^2.*G_srs);
phi_i=(L.*G_ie.*phi_e + L.*G_is.*phi_s.*exp(1i*w*tau_es)) ./ (1-G_ii*L);
phi_r=L.*(G_re.*phi_e.*exp(1i*w*tau_re)+G_rs.*phi_s);



%expression for Q's:
Qe=L.*(G_ee*phi_e+G_ei*phi_i+G_es.*exp(1i.*w.*tau_es).*phi_s);
Qi=L.*(G_ii.*phi_i + G_ie.*phi_e + G_is.*exp(1i.*w.*tau_es).*phi_s);
Qr=L.*(G_re*phi_e.*exp(1i.*w.*tau_re)+G_rs.*phi_s);
Qs=L.*(G_se*phi_e.*exp(1i.*w.*tau_se)+G_sr*phi_r+G_sn*phi_n);
Qn=phi_n.*ones(1,length(w));


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
    %s_es
    dses(1)=(1/2*pi)*conj(Hw(1))*Qe(1)*conj(Qs(1));
    dses(i)=(1/2*pi)*conj(Hw(i))*Qe(i)*conj(Qs(i));
    Ies(1)=dses(1)*dw;
    Ies(i)=Ies(i-1)+dses(i)*dw;
    %s_ie
    dsie(1)=(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qe(1));
    dsie(i)=(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qe(i));
    Iie(1)=dsie(1)*dw;
    Iie(i)=Iie(i-1)+dsie(i)*dw;
    %s_ii
    dsii(1)=(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qi(1));
    dsii(i)=(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qi(i));
    Iii(1)=dsii(1)*dw;
    Iii(i)=Iii(i-1)+dsii(i)*dw;
    %s_is
    dsis(1)=(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qs(1));
    dsis(i)=(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qs(i));
    Iis(1)=dsis(1)*dw;
    Iis(i)=Iis(i-1)+dsis(i)*dw;
    %s_re
    dsre(1)=(1/2*pi)*conj(Hw(1))*Qr(1)*conj(Qe(1));
    dsre(i)=(1/2*pi)*conj(Hw(i))*Qr(i)*conj(Qe(i));
    Ire(1)=dsre(1)*dw;
    Ire(i)=Ire(i-1)+dsre(i)*dw;
    %s_rs
    dsrs(1)=(1/2*pi)*conj(Hw(1))*Qr(1)*conj(Qs(1));
    dsrs(i)=(1/2*pi)*conj(Hw(i))*Qr(i)*conj(Qs(i));
    Irs(1)=dsrs(1)*dw;
    Irs(i)=Irs(i-1)+dsrs(i)*dw;
    %s_se
    dsse(1)=(1/2*pi)*conj(Hw(1))*Qs(1)*conj(Qe(1));
    dsse(i)=(1/2*pi)*conj(Hw(i))*Qs(i)*conj(Qe(i));
    Ise(1)=dsse(1)*dw;
    Ise(i)=Ise(i-1)+dsse(i)*dw;
    %s_sr
    dssr(1)=(1/2*pi)*conj(Hw(1))*Qs(1)*conj(Qr(1));
    dssr(i)=(1/2*pi)*conj(Hw(i))*Qs(i)*conj(Qr(i));
    Isr(1)=dssr(1)*dw;
    Isr(i)=Isr(i-1)+dssr(i)*dw;
    %s_sn
    dssn(1)=(1/2*pi)*conj(Hw(1))*Qs(1)*conj(Qn(1));
    dssn(i)=(1/2*pi)*conj(Hw(i))*Qs(i)*conj(Qn(i));
    Isn(1)=dssn(1)*dw;
    Isn(i)=Isn(i-1)+dssn(i)*dw;
    
end

%G_ab after adding terms to make it the same dimension as Imax/sab's
last=length(w);

Gfinal(c,1:11)=[G_ee,G_ei,G_es,G_ie,G_ii,G_is,G_re,G_rs,G_se,G_sr,G_sn];
nufinal(c,1:11)=[nu_ee,nu_ei,nu_es,nu_ie,nu_ii,nu_is,nu_re,nu_rs,nu_se,nu_sr,nu_sn];
rhos=Gfinal./nufinal;
Ilast(c,1:11)=[Iee(last),Iei(last),Ies(last),Iie(last),Iii(last),Iis(last),Ire(last),Irs(last),Ise(last),Isr(last),Isn(last)];
%number of synapses
N(c,1:11)=[10^4,1600,800,10^4,1600,700,2700,700,1100,550,450];
dGdt=N.*rhos.*Ilast; % nx11 matrix of dG/dts


dG_ee=dGdt(c,1);
dS_ee(c,1)=Ilast(c,1);

dG_ei=dGdt(c,2);
dS_ei(c,1)=Ilast(c,2);


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

dGeedt(c,1)=real(dG_ee);
dGeidt(c,1)=real(dG_ei);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
%   dXdts=smoothn(real(dXdt));
%   dYdts=smoothn(real(dYdt));
     radius=sqrt(real(dXdt).^2 + real(dYdt).^2);
     dUdt=(real(dXdt))./real(radius);
     dVdt=(real(dYdt))./real(radius);
     dXdts=smoothn(real(dUdt),1);
     dYdts=smoothn(real(dVdt),1);
      %quiver(X(index),Y(index),transpose(dUdt),transpose(dVdt),0.3,'blue')
     rad=sqrt(real(dXdts).^2 + real(dYdts).^2);
     dXdt_sm=dXdts./rad;
     dYdt_sm=dYdts./rad;
     quiver(X(index),Y(index),transpose(dXdt_sm),transpose(dYdt_sm),0.3,'blue')
     
     
     
     
     
     
     
     
     xlabel('X')
     ylabel('Y')
%        hold on 
%      tent.compute
     hold on
       tent.draw_blobs({'ec','n1'},0.01)
       hold on
%      
%      hold on
      %quiver3(real(X),real(Y),real(Z),dUdt,dVdt,dWdt,0.3,'black')
      grid off
      xlabel('\bf{X}')
      ylabel('\bf{Y}')
      hold on
       patch([1.2 1.2 0],[-0.2 1 1],[0.9 0.9 0.9],'EdgeAlpha',0.2)
%      legend
%   
      %hold on
      % tent.surface
     view(0,90)
%      
     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     