clear
load('/suphys/sahanda/phd/corticothalamic-model/example_parameters.mat');


load('/suphys/sahanda/cortico plasticity/indexdatabig.mat');
xyz_data=S.xyz;
gab_data=S.gab;
nus_data=S.nus;


%plasticity window



nw=45000;
fmax=50;
f = linspace(0,fmax,nw);
w=2*pi*f;
%w=-500:0.1:500;



% Pure STDP ---------------------------------------------------
tp=0.01; %plasticity timescale
A_plus=1;
%A_minus=1; %CDP
A_minus=-0.4;%STDP
%A_minus=-0.4;

H0=(A_plus + A_minus)*tp; 
H1=(A_plus - A_minus)*tp;

Hw=(H0+1i*w.*tp*H1)./(1+(w*tp).^2);


% Pure CDP ---------------------------------------------------
% tp=0.1; %plasticity timescale
% A_plus=1;
% A_minus=1; %CDP
% %A_minus=-1;%STDP
% H0=(A_plus + A_minus)*tp; 
% H1=(A_plus - A_minus)*tp;
% 
% Hw=(H0+1i*w.*tp*H1)./(1+w*tp).^2;


%triphasic H(w) ------------------------------------------------


% a=0.00002;
% b=0.0002;
% Aminus=-0.1;
% Apos=0.25;
% alphaH=0.005;
% 
% Hw=sqrt(pi).*exp(-1i*w.*alphaH).*(Apos.*sqrt(a).*exp(-0.25.*w.^2 .*a) + Aminus.*sqrt(b).*exp(-0.25.*w.^2 .*b));



% STDP + CDP ---------------------------------------------------
% a=1;
% b=-0.03937;
% tstpd=0.01;
% tcdp=0.01;
% Hw=a.*(2.*1i.*w.*tstpd.^2)./(1+(w.*tstpd).^2)+ (b.*(2*tcdp)./(1+(w.*tcdp).^2));





dw=w(2)-w(1);
%tau=linspace(-1,1,4500);
nu(1,:)=ec.nus;
nu(2,:)=n2.nus;
nu(3,:)=n3.nus;
nu(4,:)=n2s.nus;
nu(5,:)=n1.nus;
nu(6,:)=eo.nus;



G=zeros(8,8);
% %wake
% for j=1:8
%     G(j)=wgab(j); %wake 
% end
%N1 sleep
for j=1:8
     G(5,j)=n1.gab(j); %N1 sleep 
% end
% Hw=(H0+1i*w.*tp.*H1)./(1+(w.*tp).^2);
% %N2 sleep
% for j=1:8
     G(2,j)=n2.gab(j); %n2
% end
% %N3 sleep
% for j=1:8
     G(3,j)=n3.gab(j); %n3
% end
%N2S sleep

%for j=1:8
    G(4,j)=n2s.gab(j); %n2s
    G(1,j)=ec.gab(j); %w
    %G(6,1:8)=[10.593825332523812,-15.428458435110269,9.069340086179279,2.768960650292482,-0.345586010850280,2.966146825565251,12.029707107870626,2.660317515202559];
    %min G(6,1:8)=[1.460054689899175,-6.433444266744750,3.416884723262748 ,2.163056508698735, -0.145951540377034, 0.006282944866352,4.162976108515869,3.462423553355372];
     G(6,1:8)=eo.gab;%[9.174490544063291,-14.512659326324547,4.754422470509468,5.332531962375183,-1.161263288195806,7.313578776683109,3.252438507388800,0.181085848746451];
end

nu(7,:)=nus_data(1185,:);
G(7,:)=gab_data(1185,:);

nu(8,:)=nus_data(1187,:);
G(8,:)=gab_data(1187,:);


for c=6:6


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
% Qi=(Ns./A).*(J_es.*J_es+J_es-J_es.*J_ee);
% Qr=(Ns./A).*(J_ei.*J_ee+J_es.*J_re-J_es.*J_re.*J_ei+J_rs-J_rs.*J_ei-J_rs.*J_ee+J_rs.*J_ee.*J_ei+J_ei.*J_es.*J_re);
% Qs=(Ns./A).*(1+J_ee.*J_ei-J_ee-J_ei-J_ei.*J_ee);
% Qn=phi_n*ones(1,length(w));
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

X(c)=G_ee/(1-G_ei);
Y(c)=(G_ese+G_erse)/((1-G_srs)*(1-G_ei));
Z(c)=-G_srs*(alpha*beta)/(alpha+beta)^2;

Q=[Qe;Qi;Qr;Qs];
phi_e=Qe./Gamma_e;
I=zeros(1,length(w));
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
Integrand=Iee(1,length(w));
I=zeros(10,length(w));
% I(1,:)=Iee;
% I(2,:)=Iei;
% I(3,:)=Ies;
% I(4,:)=Iie;
% I(5,:)=Iii;
% I(6,:)=Iis;
% I(7,:)=Ire;
% I(8,:)=Irs;
% I(4,:)=Ise;
% I(5,:)=Isr;
%Isn

k0=10;
Lx=0.5;
Ly=0.5;

%frequency range and resolution
%nw=1000;
f1 = linspace(0,45,nw);
w1=2*pi*f1;
dw1=w1(2)-w1(1);
kmax=6;
dk=2*pi/Lx;
m_rows=-kmax:kmax;
n_cols=-kmax:kmax;
[kxa,kya]=meshgrid(dk*m_rows,dk*n_cols);
k2=kxa.^2+kya.^2;
k2u=unique(k2(:));
k2u=[k2u histc(k2(:),k2u)];

P=zeros(size(w1));

%for m=1:length(m_rows)
%k=1  
   % for n=1:length(n_cols)
     %   k=k2(m,n);
        
   for j=1:size(k2u,1)
       k=k2u(j,1);
       
 

       
         %k=sqrt((k_x)^2 + (k_y)^2);
%L(w)
L=((1-1i*w1./alpha).*(1-1i*w1./beta)).^-1;
%q^2*r_ee^2
q_re=(1-1i*w1./gamma_ee).^2 - (1./(1-G_ei*L)).*(L.*G_ee + ((L.^2.*G_ese + L.^3 .*G_erse).*exp(1i*w1.*t_0))./(1-L.^2.*G_srs));
%Transfer function for e:
T=(L.^2).*G_esn.*exp(1i*w1.*t_0*0.5)./((1-(L.^2)*G_srs).*(1-G_ei*L).*(((k^2) .* (r_ee.^2))+(q_re)));
%Power spectrum
        P=P+k2u(j,2)*abs(T).^2 .* abs(phi_n).^2.*exp(-k/k0^2);
        
        
        
   end
    
P = P.*dk.^2; % Multiply by dk then sum to get the integral over k
    P = P(:).'*2*pi; % Convert to P(f)

alpha_freq = 1/(t_0 + 1/alpha + 1/beta);
freq=w/(2*pi);
%plot(w,Iee,'red')
%hold on 
%plot(w,Iei,'blue')
%hold on
%plot(w,Ies,'green')
subplot(1,2,1)
plot(freq,(dsee),freq,dsei,freq,dses,freq,dsre,freq,dsrs,'--',freq,dsse,'--',freq,dssr,'--',freq,dssn,':','LineWidth',2)
%plot(freq,log(dsee),freq,log(dsei),freq,dses,freq,dsie,freq,dsii,freq,dsis,freq,dsre,freq,dsrs,'--',freq,dsse,'--',freq,dssr,'--',freq,dssn,':')
hold on

xlabel('f(Hz)')
ylabel('Integrand')

%legend('ds_{ee}/dt','ds_{ei}/dt','ds_{es}/dt','ds_{re}/dt','ds_{rs}/dt','ds_{se}/dt','ds_{sr}/dt','ds_{sn}/dt')
box off
handaxes2 = axes('Position', [0.327 0.795 0.124 0.107]);

  loglog(f,P)
 xlabel('f(Hz)')
 ylabel('P(s^{-1})')
 hold on
 plot(alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
 hold on
 plot(2*alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
 hold on
 plot(3*alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
 text(alpha_freq+1,10^-8*1.1,'f_\alpha')
 text(2*alpha_freq+1,10^-8*1.1,'2f_\alpha')
  text(3*alpha_freq+1,10^-8*1.1,'3f_\alpha')
 set(gca,'XLim',[1 50],'YLim',[min(P) max(P)]);
 
 
box off
subplot(1,2,2)

plot(freq,Iee,freq,Iei,freq,Ies,freq,Ire,freq,Irs,'--',freq,Ise,'--',freq,Isr,'--',freq,Isn,':','LineWidth',2)
legend('I_{ee}','I_{ei}','I_{es}','I_{re}','I_{rs}','I_{se}','I_{sr}','I_{sn}')
xlabel('f_{max}(Hz)')
ylabel('Integral')
box off
hold on
 last=length(w);

Gfinal(c,1:11)=[G_ee,G_ei,G_es,G_ie,G_ii,G_is,G_re,G_rs,G_se,G_sr,G_sn];
nufinal(c,1:11)=[nu_ee,nu_ei,nu_es,nu_ie,nu_ii,nu_is,nu_re,nu_rs,nu_se,nu_sr,nu_sn];
rhos=Gfinal./nufinal;
%rhos=ones(size(nufinal));
Ilast(c,1:11)=[Iee(last),Iei(last),Ies(last),Iie(last),Iii(last),Iis(last),Ire(last),Irs(last),Ise(last),Isr(last),Isn(last)];
%number of synapses
N(c,1:11)=[10^4,1600,800,10^4,1600,700,2700,700,1100,550,450];
%N=1;
dGdt=Ilast.*N.*rhos; % nx11 matrix of dG/dts


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

dGeedt(c,1)=dG_ee;
dGeidt(c,1)=dG_ei;
dGrsdt(c,1)=dG_rs; 




end

radius=sqrt(real(dXdt).^2 + real(dYdt).^2);
     dUdt=(real(dXdt))./real(radius);
     dVdt=(real(dYdt))./real(radius);

%END****************************




