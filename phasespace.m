

%wake--------------------------------------------



for k=1:2
G_minw=G(6009,k);
G_maxw=G_fw(k);
theta=[0:0.01:2*pi] ;
G1w=linspace(G_minw,G_maxw,1000);

Tw=16;
nw=2*pi./Tw;
bw=(G_minw+G_maxw)./2;
%dGeew=sqrt(nw^2 .*(1-(G1w-bw).^2 ./bw.^2));
dGeew=-(G_maxw-G_minw)/Tw;
%sleep--------------------------------------------

Gmini=[16.611 -17.4237];
Gmaxi=[16.599 -17.4237];
G_mins=Gmini(k);
G_maxs=Gmaxi(k);

theta=[0:0.01:2*pi] ;
G1s=linspace(G_mins,G_maxs,1000);



Ts=8;
ns=2*pi./Ts;
bs=(G_mins+G_maxs)./2;
%dGees=-sqrt(ns^2 .*(1-(G1s-bs).^2 ./bs.^2));
dGees=(G_maxs-G_mins)/Ts;



subplot(2,2,[1,3])
plot(G1w,dGeew.*(G1w./G1w),'blue',G1s,dGees.*(G1w./G1w),'red','LineWidth',3)

hold on
line([G_minw,G_mins],[min(dGeew),(dGees(end))],'LineStyle','--','color','black')
hold on
line([G_maxw,G_maxs],[(dGeew(end)),dGees(end)],'LineStyle','--','color','black')
hold on

end
hold on
h_legend=legend('Wake','Sleep');
set(h_legend,'FontSize',20);
xlabel('$G$','FontSize',30,'interpreter','latex')
ylabel('$\dot{G}$','FontSize',30,'interpreter','latex')
set(gca,'FontSize',16)

%Window plots
t2=-0.05:0.001:0;

t1=0:0.001:0.05;

%EE
A_plusrs=(H0rs+H1rs)./(2*tp);
A_minusrs=(H0rs-H1rs)./(2*tp);
Ht1rs=A_plusrs*exp(-t1/tp);
Ht2rs=A_minusrs*exp(t2/tp);
subplot(2,2,2)
plot(t1,Ht1rs,t2,Ht2rs,'blue')
ylabel('H_{rs}(\tau)','FontSize',20)
xlabel('\tau','FontSize',20)
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',0.1) %x-axis
line([0 0], yL,'color','k','linewidth',0.1) %y-axis
%EI
A_plus=(H0+H1)./(2*tp);
A_minus=(H0-H1)./(2*tp);
Ht1=A_plus*exp(-t1/tp);
Ht2=A_minus*exp(t2/tp);
set(gca,'FontSize',16)
subplot(2,2,4)
plot(t1,Ht1,t2,Ht2,'blue')
ylabel('H_{sr}(\tau)','FontSize',20)
xlabel('\tau','FontSize',20)
xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',0.1) %x-axis
line([0 0], yL,'color','k','linewidth',0.1) %y-axis
set(gca,'FontSize',16)
