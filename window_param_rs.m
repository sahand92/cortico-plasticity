clear

tp=0.01; %plasticity timescale

A_plus=1;
A_minus=-0.8;


 fmax=45;
 nw=1000;
 f = linspace(0,fmax,nw);
 w=2*pi*f;
 
H0=(A_plus + A_minus)*tp; 
H1=(A_plus - A_minus)*tp;
Hw=(H0+1i*w.*tp*H1)./(1+(w*tp).^2);
 
 
 
 %Window function for RS-SR population------------------------
 
A_plusrs=1;
%A_minusrs=-0.2;



 

for A_minusrs=-0.6
    
    
H0rs=-(A_plusrs + A_minusrs)*tp; 
H1rs=(A_plusrs - A_minusrs)*tp;    
Hwrs=(H0rs+1i*w.*tp*H1rs)./(1+(w*tp).^2); 

    


% H(w) ---------------

t1=0:0.001:0.05;
Ht1=A_plus*exp(-t1/tp);


t2=-0.05:0.001:0;
Ht2=A_minus*exp(t2/tp);


%H(w)-RS---------------

Ht1rs=A_plusrs*exp(-t1/tp);
Ht2rs=A_minusrs*exp(t2/tp);





subplot(2,2,1)
plot(t1,Ht1,t2,Ht2,'blue')
title('H(\tau)')
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',0.1) %x-axis
line([0 0], yL,'color','k','linewidth',0.1) %y-axis
xlabel('\tau')


subplot(2,2,2)
plot(t1,Ht1rs,t2,Ht2rs,'blue')
title('H_{rs,sr}(\tau)')
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',0.1) %x-axis
line([0 0], yL,'color','k','linewidth',0.1) %y-axis
xlabel('\tau')

subplot(2,2,3)
plot(w,Hw,w,imag(Hw))
title('H(\omega)')
xlabel('\omega')

subplot(2,2,4)
plot(w,Hwrs,w,imag(Hwrs))
title('H_{rs,sr}(\omega)')
xlabel('\omega')


end