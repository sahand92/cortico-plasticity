clear
load('/suphys/sahanda/neurofield/corticothalamic-model/example_parameters.mat');
%plasticity window
tp=0.01; %plasticity timescale

A_plus=1;
%A_minus=1; %CDP
A_minus=-1;%STDP
H0=(A_plus + A_minus)*tp; 
H1=(A_plus - A_minus)*tp;

nw=450000;
fmax=45;
f = linspace(0,fmax,nw);
w=2*pi*f;
Hw=(H0+1i*w.*tp*H1)./(1+w*tp).^2;
dw=w(2)-w(1);
%tau=linspace(-1,1,4500);
wgab = [2.074,-4.110,0.772,7.768,-3.301,8.097,0.656,0.196];

G=zeros(6,8);
% %wake
% for j=1:8
%     G(j)=wgab(j); %wake 
% end
%N1 sleep
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
% end
%N2S sleep

%for j=1:8
    G(4,j)=n2s.gab(j); %n2s
    G(5,j)=wgab(j); %w
    
    
end



for c=5:5


%paramters wake
G_ee=G(c,1);
G_ei=G(c,2);
G_es=G(c,3);
G_re=G(c,7);
G_rs=G(c,8);
G_sr=G(c,5);
G_sn=G(c,6);
G_se=G(c,4);
%Random connectivity assumption
G_ie=G_ee;
G_ii=G_ei;
G_is=G_es;

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
Gamma_e=(1-1i*w/gamma_ee).^2;
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
f = linspace(0,45,nw);
w=2*pi*f;
dw=w(2)-w(1);
kmax=6;
dk=2*pi/Lx;
m_rows=-kmax:kmax;
n_cols=-kmax:kmax;
[kxa,kya]=meshgrid(dk*m_rows,dk*n_cols);
k2=kxa.^2+kya.^2;
k2u=unique(k2(:));
k2u=[k2u histc(k2(:),k2u)];

P=zeros(size(w));

%for m=1:length(m_rows)
%k=1  
   % for n=1:length(n_cols)
     %   k=k2(m,n);
        
   for j=1:size(k2u,1)
       k=k2u(j,1);
       
 

       
         %k=sqrt((k_x)^2 + (k_y)^2);
%L(w)
L=((1-1i*w./alpha).*(1-1i*w./beta)).^-1;
%q^2*r_ee^2
q_re=(1-1i*w./gamma_ee).^2 - (1./(1-G_ei*L)).*(L.*G_ee + ((L.^2.*G_ese + L.^3 .*G_erse).*exp(1i*w.*t_0))./(1-L.^2.*G_srs));
%Transfer function for e:
T=(L.^2).*G_esn.*exp(1i*w.*t_0*0.5)./((1-(L.^2)*G_srs).*(1-G_ei*L).*(((k^2) .* (r_ee.^2))+(q_re)));
%Power spectrum
        P=P+k2u(j,2)*abs(T).^2 .* abs(phi_n).^2.*exp(-k/k0^2);
        
        
        
   end
    
P = P.*dk.^2; % Multiply by dk then sum to get the integral over k
    P = P(:).'*2*pi; % Convert to P(f)

alpha_freq = 1/(t_0 + 1/alpha + 1/beta);
%plot(w,Iee,'red')
%hold on 
%plot(w,Iei,'blue')
%hold on
%plot(w,Ies,'green')
subplot(1,2,1)
freq=w/(2*pi);
plot(freq,Iee,freq,Iei,freq,Ies,freq,Iie,freq,Iii,freq,Iis,freq,Ire,freq,Irs,'--',freq,Ise,'--',freq,Isr,'--',freq,Isn,':')
legend('I_{ee}','I_{ei}','I_{es}','I_{ie}','I_{ii}','I_{is}','I_{re}','I_{rs}','I_{se}','I_{sr}','I_{sn}')
% plot(freq,Iee);
% hold on 
% plot(w/(2*pi),Iei)
% hold on
% plot(w/(2*pi),Ies)
% hold on
% plot(w/(2*pi),Iie)
% hold on
% plot(w/(2*pi),Iii)
% hold on
% plot(w/(2*pi),Iis)
% hold on
% plot(w/(2*pi),Ire)
% hold on
% plot(w/(2*pi),Irs)
% hold on
% plot(w/(2*pi),Ise)
% hold on
% plot(w/(2*pi),Isr)
% hold on
  xlabel('f(Hz)')
  ylabel('Integral')
 subplot(1,2,2)
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
 set(gca,'XLim',[1 50],'YLim',[10^-8 10^-2]);
 hold on
end

%END****************************




