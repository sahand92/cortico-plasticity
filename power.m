clear


for c=1:2
%min radis
 G(1,1:8)=[1.460054689899175,-6.433444266744750,3.416884723262748 ,2.163056508698735, -0.145951540377034, 0.006282944866352,4.162976108515869,3.462423553355372];
%max radius    
 G(2,1:8)=[9.174490544063291,-14.512659326324547,4.754422470509468,5.332531962375183,-1.161263288195806,7.313578776683109,3.252438507388800,0.181085848746451];

 
 
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

t_0=0.085;
tau_re=t_0/2;
tau_se=t_0/2;
tau_es=t_0/2;

r_ee=0.086;
phi_n=0.001;

k0=10;
Lx=0.5;
Ly=0.5;
nw=5000;
alpha=83.3;
beta=770;
gamma_ee=116;

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
loglog(w,P)
hold on
end