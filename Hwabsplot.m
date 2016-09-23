clear
w=-1000:0.01:1000;
f=w/(2*pi);
a=1;
b=-0.3;
tstpd=0.01;
tcdp=0.01;
Hw=a.*(2.*1i.*w.*tstpd.^2)./(1+(w.*tstpd).^2)+ (b.*(2*tcdp)./(1+(w.*tcdp).^2));
%f1=linspace(-200,200,length(w));
figure
plot(w,Hw,w,imag(Hw),'red','LineWidth',2)
hold on
plot(w,(real(Hw)+imag(Hw)),'--','color','black')
xlabel('f(Hz)','FontSize',14)
ylabel('H','FontSize',14)
legend('H_{CDP}(\omega)','H_{STDP}(\omega)','\Sigma H(\omega)')
 xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',1) %x-axis
line([0 0], yL,'color','k','linewidth',1) %y-axis

% hold on
% plot(8.9*ones(10),linspace(-10*10^-3,8*10^-3,10),'--','color','green')

box off
%tri-phasic

% 
figure
 t=0:0.0001:0.07;
 Hts=0.1*exp(-t/0.010);
 plot(t,Hts)
 hold on
Apos=0.25;
Aminus=0.1;
alpha=0.005;
Ht=Apos*exp(-((t-alpha).^2)./0.00002)-Aminus*exp(-((t-alpha).^2)./0.0002);

plot(t,Ht)
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',1) %x-axis
line([0 0], yL,'color','k','linewidth',1) %y-axis

figure
a=0.00002;
b=0.0002;
Aminus=-0.1;
Hwtri=sqrt(pi).*exp(-1i*w.*alpha).*(Apos.*sqrt(a).*exp(-0.25.*w.^2 .*a) + Aminus.*sqrt(b).*exp(-0.25.*w.^2 .*b));
plot(f,Hwtri)
hold on
plot(f,imag(Hwtri),'red')