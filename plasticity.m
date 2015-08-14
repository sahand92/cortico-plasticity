
clear
load('/suphys/sahanda/neurofield/corticothalamic-model/example_parameters.mat');
wgab = [2.074,-4.110,0.772,7.768,-3.301,8.097,0.656,0.196];

G=zeros(5,8);
%wake
for j=1:8
    G(1,j)=wgab(j); %wake 

    G(2,j)=n1.gab(j); %N1 sleep 

    G(3,j)=n2.gab(j); %N2

    G(4,j)=n3.gab(j); %N3
    G(5,j)=n2s.gab(j); %N2s
end
for c=1:5
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


t_0=0.085;
tau_re=t_0/2;
tau_se=t_0/2;
tau_es=t_0/2;

r_ee=0.086;
k0=10;
Lx=0.5;
Ly=0.5;
phi_n=0.001;
%frequency range and resolution
nw=1000;
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

loglog(f,P/P(1))
set(gca,'XLim',[1 100]);
%plot(k_x,k_y,log(P_j))
hold on
end
        

