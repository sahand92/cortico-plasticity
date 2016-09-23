clear 
load('/suphys/sahanda/phd/corticothalamic-model/example_parameters.mat');
%load('/suphys/sahanda/cortico plasticity/indexdatabig.mat');
S=load('/suphys/sahanda/cortico plasticity/data/pdb_all.mat');
xyz_data=S.xyz_final(1:100:3000000,:);
gab_data=S.gab_final(1:100:3000000,:);
nus_data=S.nus_final(1:100:3000000,:);
%n=1; %number of points is 3000000/n
% xyz_data=S.xyz;
% gab_data=S.gab;
% nus_data=S.nus;
%plasticity window



nw=450;
fmax=45;
f = linspace(0,fmax,nw);
w=2*pi*f;

% Pure CDP ---------------------------------------------------
% tp=0.01; %plasticity timescale
% A_plus=1;
% %A_minus=1; %CDP
% A_minus=-1;%STDP
% H0=(A_plus + A_minus)*tp; 
% H1=(A_plus - A_minus)*tp;
% 
% Hw=(H0+1i*w.*tp*H1)./(1+w*tp).^2;


% % Pure STDP ---------------------------------------------------
tp=0.01; %plasticity timescale
A_plus=1;
%A_minus=1; %CDP
A_minus=-0.9;%STDP
H0=-(A_plus + A_minus)*tp; 
H1=(A_plus - A_minus)*tp;

Hw=(H0+1i*w.*tp*H1)./(1+(w*tp).^2);

%triphasic H(w) ------------------------------------------------


% a=0.00002;
% b=0.0002;
% Aminus=-0.1;
% Apos=0.25;
% alphaH=0.001;
% 
% Hw=sqrt(pi).*exp(-1i*w.*alphaH).*(Apos.*sqrt(a).*exp(-0.25.*w.^2 .*a) + Aminus.*sqrt(b).*exp(-0.25.*w.^2 .*b));



% STDP + CDP ---------------------------------------------------
% a=1;
% b=-0.3;
% tstpd=0.01;
% tcdp=0.01;
% Hw=a.*(2.*1i.*w.*tstpd.^2)./(1+(w.*tstpd).^2)+ (b.*(2*tcdp)./(1+(w.*tcdp).^2));



dw=w(2)-w(1);
%tau=linspace(-1,1,4500);
wgab = [2.074,-4.110,0.772,7.768,-3.301,8.097,0.656,0.196];



    

% Hw=(H0+1i*w.*tp.*H1)./(1+(w.*tp).^2);


        Gint(3,:)=n1.gab; %N1 sleep 
%      Gint(2,:)=n2.gab; %n2
      Gint(2,:)=n3.gab; %n3
      Gint(4,:)=n2s.gab; %n2s
      Gint(5,:)=wgab; %w
      Gint(1,:)=ec.gab;
%     
%     
% 
%  
   Gint(6,:)=Gint(5,:)+rand(size(Gint(5,:)))/5;
   Gint(7,:)=Gint(5,:)-rand(size(Gint(5,:)))/5;
   Gint(8,:)=gab_data(1059,:);
%      Gint(1,:)=gab_data(934,:);
%      nus=nus_data(934,:);
%    Gint(14,:)=gab_data(934,:)+rand(size(Gint(5,:)))/10;
%    Gint(15,:)=gab_data(934,:)-rand(size(Gint(5,:)))/10;
%    Gint(11,:)=gab_data(1113,:);
%    Gint(12,:)=gab_data(1113,:)+rand(size(Gint(5,:)))/5;
%    Gint(13,:)=gab_data(1113,:)-rand(size(Gint(5,:)))/5;

%-------------------------------------------------
    %Gint(1,:)=gab_data(934,:);
%     Gint(2,:)=gab_data(934,:)+round(randn(1,8))/10;
%     Gint(3,:)=gab_data(934,:)+round(randn(1,8))/10;
%     Gint(4,:)=gab_data(934,:)+round(randn(1,8))/10;
%     Gint(5,:)=gab_data(934,:)+round(randn(1,8))/10;
%     Gint(6,:)=gab_data(934,:)+round(randn(1,8))/10;
%     Gint(7,:)=gab_data(934,:)+round(randn(1,8))/10;
%     Gint(8,:)=gab_data(934,:)+round(randn(1,8))/10;
%     Gint(9,:)=gab_data(934,:)+round(randn(1,8))/10;
    
%-------------------------------------------------

 %colorRed: [1 0.5 0.5]
%Blue: [0.5 0.5 1]

%Gint=gab_data;

C=zeros(1,11);
C(3)=1;
C(1)=1;
C(5)=0.5;
C(2)=1;
C(4)=0.5;
C(6)=1;
C(7)=0.5;
C(8)=0.1;
C(9)=0.8;
C(10)=0.2;
C(11)=0.7;
for p=1:8
%initial G
 G=Gint(p,:);

for c=1:20


%paramters wake
G_ee=G(1);
G_ei=G(2);
G_es=G(3);
G_re=G(7);
G_rs=G(8);
G_sr=G(5);
G_sn=G(6);
G_se=G(4);
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

X=G_ee/(1-G_ei);
Y=(G_ese+G_erse)/((1-G_srs)*(1-G_ei));
Z=-G_srs*(alpha*beta)/(alpha+beta)^2;

% J_ee=G_ee*L.*Gamma_e;
% J_ei=G_ei*L;
% J_es=G_es*L.*exp(1i*w.*tau_es);
% J_re=G_re*L.*exp(1i*w.*tau_re);
% J_se=G_se*L.*exp(1i*w.*tau_se);
% J_rs=G_rs*L;
% J_sr=G_sr*L;
% J_sn=G_sn*L;
% 
% % default phi_n=0.001;
 phi_n=0.001;%+0.0002*c*normpdf(f,9,1)+c*0.0001*normpdf(f,18,1);
% Ns=phi_n.*G_sn.*L;
% A=(1-J_ee-J_ei).*(1-J_sr.*J_rs)-J_es.*(J_se+J_sr.*J_re);
% 
% 
% Qe=(Ns./A).*(J_ei.*J_es+J_es-J_ee.*J_ei);
% Qi=(Ns./A).*(J_es.*J_es+J_es-J_es.*J_ee);
% Qr=(Ns./A).*(J_ei.*J_ee+J_es.*J_re-J_es.*J_re.*J_ei+J_rs-J_rs.*J_ei-J_rs.*J_ee+J_rs.*J_ee.*J_ei+J_ei.*J_es.*J_re);
% Qs=(Ns./A).*(1+J_ee.*J_ei-J_ee-J_ei-J_ei.*J_ee);
% Qn=phi_n.*ones(1,length(w));
% Q=[Qe;Qi;Qr;Qs];

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


Integrand(c,1:8)=real([Iee(length(w)),Iei(length(w)),Ies(length(w)),Ise(length(w)),Isr(length(w)),Isn(length(w)),Ire(length(w)),Irs(length(w))]);

%normalize 

%calculate distance to nearest nus
% nusm=zeros(2,8);
% nusm(1,:)=nus;
% nusm(2,:)=nusm(1,:)+Integrand(1,1:8);
% dist=zeros(1,length(nus_data));
%     for l=1:length(nus_data(:,1))
%         delta=abs(nus_data(l,:)-nusm(2,:));
%         dist(l)=sqrt(delta(1)^2+delta(2)^2+delta(3)^2+delta(4)^2+delta(5)^2+delta(6)^2+delta(7)^2+delta(8)^2);
%         ind=find(dist==min(dist));
%         inde=sort(dist);
%         
%         
%         
%         
%        
%     end
%     
%           if find((dist==min(dist))==find(nus(:,1)==nus_data(:,1)));
%           index=inde(2);
%           else index=inde(1);
%           end

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
freq=w/(2*pi);
%plot(w,Iee,'red')
%hold on 
%plot(w,Iei,'blue')
%hold on
%plot(w,Ies,'green')
%subplot(1,2,1)
%plot(freq,(dsee),'color',[C(c),0.5,C(c+1)],'LineWidth',2)%,freq,dsei,freq,dses,freq,dsre,freq,dsrs,'--',freq,dsse,'--',freq,dssr,'--',freq,dssn,':')
%plot(freq,log(dsee),freq,log(dsei),freq,dses,freq,dsie,freq,dsii,freq,dsis,freq,dsre,freq,dsrs,'--',freq,dsse,'--',freq,dssr,'--',freq,dssn,':')

%hold on


% xlabel('f(Hz)','FontSize',14)
% ylabel('Integrand','FontSize',14)

%legend('ds_{ee}/dt','ds_{ei}/dt','ds_{es}/dt','ds_{re}/dt','ds_{rs}/dt','ds_{se}/dt','ds_{sr}/dt','ds_{sn}/dt')
%box off
%handaxes2 = axes('Position', [0.327 0.795 0.124 0.107]);


% 
%    loglog(f,P,'black' ,'color',[C(c),0.5,C(c+1)])
%   xlabel('f(Hz)','FontSize',14)
%   ylabel('P(s^{-1})','FontSize',14)
%  
 
%scatter(X,Y,'.','black')
%save XYZ point

Xdata(p,c)=X;
Ydata(p,c)=Y;
Zdata(p,c)=Z;

%  hold on
%  plot(alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
%  hold on
%  plot(2*alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
%  hold on
%  plot(3*alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
%  text(alpha_freq+1,10^-8*1.1,'f_\alpha')
%  text(2*alpha_freq+1,10^-8*1.1,'2f_\alpha')
%   text(3*alpha_freq+1,10^-8*1.1,'3f_\alpha')
 %set(gca,'XLim',[1 50],'YLim',[2*min(P) 2*max(P)]);
 
 
box off
%subplot(1,2,2)

%plot(freq,Iee,'color',[C(c),0.5,C(c+1)],'Linewidth',2)%,freq,Iei,freq,Ies,freq,Ire,freq,Irs,'--',freq,Ise,'--',freq,Isr,'--',freq,Isn,':')

% legend('I_{ee}-Sleep','I_{ee}-Wake')
% xlabel('f_{max}(Hz)','FontSize',14)
% ylabel('Integral','FontSize',14)
% box off

hold on
rhos=1;
%number of synapses
N=[10^4,1600,800,700,2700,1100,550,450];
%G=N.*G;
%time interval
tc=1;

if G(1)>0 && G(1)<20
if real(Iee(end))>0;
    G(1)=G(1)+tc.*(Iee(end)).*N(1);%*sqrt(20-abs(G(1)));
else 
    G(1)=G(1)-tc.*(Iee(end)).*N(1);%*sqrt(abs(G(1)));
end
end

if G(2)<0 && G(2)>-20
if real(Iei(end))>0;
    G(2)=G(2)+tc.*Iei(end).*N(2);%*sqrt(abs(G(2)));
else
    G(2)=G(2)-tc.*Iei(end).*N(2);%*sqrt(20-abs(G(2)));
end
end

if G(3)>0 && G(3)<20
if real(Ies(end))>0;
    G(3)=G(3)+tc.*Ies(end).*N(3);%*sqrt(20-abs(G(3)));
else
    G(3)=G(3)-tc.*Ies(end).*N(3);%*sqrt(abs(G(3)));
end
end

if G(7)>0 && G(7)<20
if real(Ire(end))>0;
    G(7)=G(7)+tc.*Ire(end).*N(7);%*sqrt(20-abs(G(7)));
else
    G(7)=G(7)-tc.*Ire(end).*N(7);%*sqrt(abs(G(7)));
end
end

if G(8)>0 && G(8)<20
if real(Irs(end))>0;
    G(8)=G(8)+tc.*Irs(end).*N(8);%*sqrt(20-abs(G(8)));
else
    G(8)=G(8)-tc.*Irs(end).*N(8);%*sqrt(abs(G(8)));
end
end

if G(5)<0 && G(5)>-20
if real(Isr(end))>0;
    G(5)=G(5)+tc.*Isr(end).*N(5);%*sqrt(abs(G(5)));
else
    G(5)=G(5)-tc.*Isr(end).*N(5);%*sqrt(20-abs(G(5)));
end
end

if G(6)>0 && G(6)<20
if real(Isn(end))>0;
    G(6)=G(6)+tc.*Isn(end).*N(6);%*sqrt(20-abs(G(6)));
else
    G(6)=G(6)-tc.*Isn(end).*N(6);%*sqrt(abs(G(6)));
end
end

if G(4)>0 && G(4)<20
if real(Ise(end))>0;
    G(4)=G(4)+tc.*Ise(end).*N(4);%*sqrt(20-abs(G(4)));
else
    G(4)=G(4)-tc.*Ise(end).*N(4);%*sqrt(abs(G(4)));
end



end
% scatter(c,G(1),'.','blue')
% hold on
% scatter(c,G(2),'.','green')
% hold on
% scatter(c,G(3),'.','red')
% hold on
% scatter(c,G(4),'.','black')
% hold on
% scatter(c,G(5),'.','magenta')
% hold on
% scatter(c,G(6),'.','cyan')
% hold on
% scatter(c,G(7),'o','blue')
% hold on
% scatter(c,G(8),'o','red')
end
hold on
end
% hold on
%  patch([1.2 1.2 0],[-0.2 1 1],[0.9 0.9 0.9],'EdgeAlpha',0.2)
%   hold on
%   plot(Xdata(1,:),Ydata(1,:))%,Zdata(1,:))
% hold on
% plot3(Xdata(2,:),Ydata(2,:),Zdata(2,:),'red')
% hold on
% %wgab XYZ
% scatter(0.405870841487280,0.513912441451432);
% %Z:   0.056994615874845
% hold on
% 
% %N1 sleep XYZ
% scatter(0.800697032322903,-0.015804458936975);
% 
% % Z= 0.157290544441891
% hold on
% 
% %N3
% scatter(0.941341620064136,-0.039443508995960);
% % Z:0.026003435853375
% hold on
% 
% %N2
% scatter(0.890587116453429,-0.059580983971334)
% 
% %Z=0.103411319825169)
% hold on
% 
% scatter(0.928147769223030,-0.011785971675854)

%Z=0.393130442657386
%H(w) plot -------------------------------------------------
% figure
% 
% 
% plot(freq,Hw,freq,imag(Hw),'red','LineWidth',2)
% hold on
% plot(freq,(real(Hw)+imag(Hw)),'--','color','black')
% xlabel('f(Hz)','FontSize',14)
% ylabel('H','FontSize',14)
% legend('H_{CDP}(\omega)','H_{STDP}(\omega)','\Sigma H(\omega)')
%  xL = xlim;
% yL = ylim;
% line(xL, [0 0],'color','k','linewidth',1) %x-axis
% line([0 0], yL,'color','k','linewidth',1) %y-axis
% 
 ws=-1000:0.01:1000;
 freqs=ws/(2*pi);
% %STDP+CDP---------------------------------------------------------------------
%  Hws=a.*(2.*1i.*ws.*tstpd.^2)./(1+(ws.*tstpd).^2)+ (b.*(2*tcdp)./(1+(ws.*tcdp).^2));
% %------------------------------------------------------------------------------

%Pure STDP------------------------
tp=0.01; %plasticity timescale

H0=0*tp; 
H1=2*tp;

Hws=(1i*ws*tp*H1)./(1+(ws*tp).^2);
%----------------------------------------


% Pure CDP ---------------------------------------------------
% tp=0.01; %plasticity timescale
% H0=2*tp; 
% H1=0;
% 
% Hws=(H0+1i*ws.*tp*H1)./(1+ws*tp).^2;
%-------------------------------------------------


% %Tri-phasic-----------------------------------------------------------------
% a=0.00002;
% b=0.0002;
% Aminus=-0.1;
% Apos=0.25;
% 
% Hws=sqrt(pi).*exp(-1i*ws.*alphaH).*(Apos.*sqrt(a).*exp(-0.25.*ws.^2 .*a) + Aminus.*sqrt(b).*exp(-0.25.*ws.^2 .*b));
% %-------------------------------------------------------------------------------
% 

figure

plot(freq,Hw,freq,imag(Hw),'red','LineWidth',2)
hold on
%plot(freqs,sqrt((real(Hws).^2 +imag(Hws).^2)),'--','color','black')
xlabel('f(Hz)','FontSize',14)
ylabel('H','FontSize',14)
legend('Real[H(\omega)]','Imag[H(\omega)]',' |H(\omega)|')
 xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',1) %x-axis
line([0 0], yL,'color','k','linewidth',1) %y-axis
hold on
line([45 45], yL,'color','g') %y-axis


figure
for ii=1:p
plot(Xdata(ii,:),Ydata(ii,:),'black')
xlabel('X','FontSize',14)
ylabel('Y','FontSize',14)
hold on
scatter(Xdata(ii,1),Ydata(ii,1),500,'green','.')
hold on
scatter(Xdata(ii,end),Ydata(ii,end),500,'red','.')
hold on
end
hold on
patch([1.2 1.2 0],[-0.2 1 1],[0.9 0.9 0.9],'EdgeAlpha',0.2)
hold on

hold on
tent.draw_blobs({'ec','n1'},0.1)
%END****************************
