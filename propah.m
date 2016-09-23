
%trajectory of nus via plasticity.

%run evol first?
%------------------------------
A_minusrs=-0.4;

tp=0.01; %plasticity timescale

A_plus=1;
%A_minus=1; %CDP
A_minus=-0.8;%STDP

H0=(A_plus + A_minus)*tp; 
H1=(A_plus - A_minus)*tp;

 fmax=45;
 nw=1000;
 f = linspace(0,fmax,nw);
 w=2*pi*f;
 Hw=(H0+1i*w.*tp*H1)./(1+(w*tp).^2);
 
 
 
 %Window function for RS-SR population------------------------
 
 A_plusrs=1;
%A_minus=1; %CDP
%A_minusrs=0.5;%STDP

H0rs=-(A_plusrs + A_minusrs)*tp; 
H1rs=(A_plusrs - A_minusrs)*tp;

 Hwrs=(H0rs+1i*w.*tp*H1rs)./(1+(w*tp).^2);
 
% Window for EI
A_plusei=0.8;
A_minusei=-1;
 H0ei=(A_plusei + A_minusei)*tp; 
H1ei=(A_plusei - A_minusei)*tp;

 Hwei=(H0ei+1i*w.*tp*H1ei)./(1+(w*tp).^2);
 
% Window for EE
A_plusee=1;
A_minusee=-0.8;
 H0ee=(A_plusee + A_minusee)*tp; 
H1ee=(A_plusee - A_minusee)*tp;

 Hwee=(H0ee+1i*w.*tp*H1ee)./(1+(w*tp).^2); 
 
 % Window for ES
A_pluses=1;
A_minuses=-0.8;
 H0es=-(A_pluses + A_minuses)*tp; 
H1es=(A_pluses - A_minuses)*tp;

 Hwes=(H0es+1i*w.*tp*H1es)./(1+(w*tp).^2); 

  % Window for RE
A_plusre=1;
A_minusre=-0.8;
 H0re=(A_plusre+ A_minusre)*tp; 
H1re=(A_plusre - A_minusre)*tp;

 Hwre=(H0re+1i*w.*tp*H1re)./(1+(w*tp).^2); 
 
 dw=w(2)-w(1);
 
 %------------------------------------------------------------------------
 
 %parameters:
 
 
%initial nus at wake = nuw1
%-----------------------------
% get G_ab for initial point:

p=model.params();
%p.nus=nuw1;
p.nus=nu(726018,:);
gab1=p.gab(1,:); 

gab=zeros(10,8);
gab(1,:)=gab1;
nus(1,:)=p.nus(1,:);%nuw1;

for c=1:10
    


G_ee=gab(c,1);
G_ei=gab(c,2);
G_es=gab(c,3);
G_re=gab(c,7);
G_rs=gab(c,8);
G_sr=gab(c,5);
G_sn=gab(c,6);
G_se=gab(c,4);

%nus:
nu_ee=nus(c,1);
nu_ei=nus(c,2);
nu_es=nus(c,3);
nu_re=nus(c,7);
nu_rs=nus(c,8);
nu_sr=nus(c,5);
nu_sn=nus(c,6);
nu_se=nus(c,4);



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

%sleep
%  alpha=45;
%  beta=185;


gamma_ee=116;
L=((1-1i*w./alpha).*(1-1i*w./beta)).^-1;

t_0=0.085;
tau_re=t_0/2;
tau_se=t_0/2;
tau_es=t_0/2;

r_ee=0.086;
Gamma_e=(1-1i*w/gamma_ee).^-2;


phi_n=0.001;

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
 
for i=2:length(w)
   
    %_see
    dsee(1)=(1/2*pi)*conj(Hwee(1))*Qe(1)*conj(Gamma_e(1))*conj(Qe(1));
    dsee(i)=(1/2*pi)*conj(Hwee(i))*Qe(i)*conj(Gamma_e(i))*conj(Qe(i));
    Iee(1)=dsee(1)*dw;
    Iee(i)=Iee(i-1)+dsee(i)*dw;
    %s_ei
    dsei(1)=(1/2*pi)*conj(Hwei(1))*Qe(1)*conj(Qi(1));
    dsei(i)=(1/2*pi)*conj(Hwei(i))*Qe(i)*conj(Qi(i));
    Iei(1)=dsei(1)*dw;
    Iei(i)=Iei(i-1)+dsei(i)*dw;
    %s_es
    dses(1)=(1/2*pi)*conj(Hwes(1))*Qe(1)*conj(Qs(1));
    dses(i)=(1/2*pi)*conj(Hwes(i))*Qe(i)*conj(Qs(i));
    Ies(1)=dses(1)*dw;
    Ies(i)=Ies(i-1)+dses(i)*dw;
    %s_ie
    dsie(1)=0; %(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qe(1));
    dsie(i)=0; %(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qe(i));
    Iie(1)=dsie(1)*dw;
    Iie(i)=Iie(i-1)+dsie(i)*dw;
    %s_ii
    dsii(1)=0; %(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qi(1));
    dsii(i)=0; %(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qi(i));
    Iii(1)=dsii(1)*dw;
    Iii(i)=Iii(i-1)+dsii(i)*dw;
    %s_is
    dsis(1)=(1/2*pi)*conj(Hw(1))*Qi(1)*conj(Qs(1));
    dsis(i)=(1/2*pi)*conj(Hw(i))*Qi(i)*conj(Qs(i));
    Iis(1)=dsis(1)*dw;
    Iis(i)=Iis(i-1)+dsis(i)*dw;
    %s_re
    dsre(1)=(1/2*pi)*conj(Hwre(1))*Qr(1)*conj(Qe(1));
    dsre(i)=(1/2*pi)*conj(Hwre(i))*Qr(i)*conj(Qe(i));
    Ire(1)=dsre(1)*dw;
    Ire(i)=Ire(i-1)+dsre(i)*dw;
    %s_rs
    dsrs(1)=(1/2*pi)*conj(Hwrs(1))*Qr(1)*conj(Qs(1));
    dsrs(i)=(1/2*pi)*conj(Hwrs(i))*Qr(i)*conj(Qs(i));
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
last=length(w);
Ilastp(c,1:11)=[Iee(last),Iei(last),Ies(last),Iie(last),Iii(last),Iis(last),Ire(last),Irs(last),Ise(last),Isr(last),Isn(last)];


Ilastpr(c,:)=[Ilastp(c,1) Ilastp(c,2) Ilastp(c,3) Ilastp(c,9) Ilastp(c,10) Ilastp(c,11) Ilastp(c,7) Ilastp(c,8)];


nus(c+1,:)=nus(c,:)+real(Ilastpr(c,:)).*Nr(1,:)/100000;

    p=model.params();
    p.nus=nus(c+1,:);
    gab(c+1,:)=p.gab(1,:);

end

Xgab=gab(:,1)./(1-gab(:,2));
Ygab=(gab(:,3).*gab(:,4)+gab(:,3).*gab(:,5).*gab(:,7))./((1-gab(:,5).*gab(:,8)).*(1-gab(:,2)));
scatter(Xgab(1),Ygab(1),'filled');
hold on
plot(Xgab,Ygab);
hold on
scatter(Xgab,Ygab);
hold on
tent.draw_blobs({'ec','n3','n1'},0.1)


