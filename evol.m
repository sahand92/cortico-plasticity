

%load('700kresults.mat');
%load('/suphys/sahanda/cortico plasticity/ind700k.mat');
%load('/import/ghrian1/sahanda/romesh-large-files/pdb_wake.mat');
% X=xyz_final(index,1);
% Y=xyz_final(index,2);
% Z=xyz_final(index,3);
% G=gab_final(index,:);
% nu=nus_final(index,:);

figure
ep=0.02;
dt=0.02;
%Starting poit:
clear lw
%lw(:,1)=find(X>0.1 & X<0.2 & Y>0.6 & Y<0.75);
lw(:,1)=find(X>0.15 & X<0.2 & Y>0.55 & Y<0.65);
lw1=lw;

     %quiver(X(lw(:,1)),Y(lw(:,1)),transpose(dUdt(lw(:,1))),transpose(dVdt(lw(:,1))),0.5,'red')
     hold on
     quiver(mean(X(lw(:,1))),mean(Y(lw(:,1))),mean(dUdt(lw(:,1))),mean(dVdt(lw(:,1))),0.1,'LineWidth',1,'MaxHeadSize',1,'color','red')
     hold on
     
     X1=mean(X(lw(:,1))+dt*mean(dUdt(lw(:,1))));
     Y1=mean(Y(lw(:,1))+dt*mean(dVdt(lw(:,1))));
     %scatter(X1,Y1)
     mean(G(lw,:));
     nuw1=mean(nu(lw,:));
     
%Iteration for next points:
for i=1:13

clear lw     
lw(:,i+1)=find(X>X1-ep & X<X1+ep & Y>Y1-ep & Y<ep+Y1);
mean(G(lw(:,i+1),:));

G_trackw(i,:)=mean(G(lw(:,i+1),:));

%quiver(X(lw(:,i+1)),Y(lw(:,i+1)),transpose(dUdt(lw(:,i+1))),transpose(dVdt(lw(:,i+1))),0.5,'red')
     hold on
     quiver(mean(X(lw(:,i+1))),mean(Y(lw(:,i+1))),mean(dUdt(lw(:,i+1))),mean(dVdt(lw(:,i+1))),0.1,'LineWidth',1,'MaxHeadSize',1,'color','red')
     hold on
     
     X1=mean(X(lw(:,i+1))+dt*mean(dUdt(lw(:,i+1))));
     Y1=mean(Y(lw(:,i+1))+dt*mean(dVdt(lw(:,i+1))));
     %scatter(X1,Y1)
     hold on

end
     

     
  
clear ls
%ls(:,1)=find(X>0.8 & X<1 & Y>-0.15 & Y<-0.0);
ls(:,1)=find(X>0.75 & X<0.85 & Y>0 & Y<0.1);
ls1=ls;
%quiver(X(ls),Y(ls),transpose(dUdt(ls)),transpose(dVdt(ls)),0.5,'blue')
hold on
quiver(mean(X(ls)),mean(Y(ls)),mean(dUdt(ls)),mean(dVdt(ls)),0.1,'LineWidth',2,'MaxHeadSize',2,'color','blue')
hold on

X1=mean(X(ls(:,1))+dt*mean(dUdt(ls(:,1))));
     Y1=mean(Y(ls(:,1))+dt*mean(dVdt(ls(:,1))));
     %scatter(X1,Y1)
%Iteration for next points:
dt1=0.01;
for i=1:15

clear ls     
ls(:,i+1)=find(X>X1-ep & X<X1+ep & Y>Y1-ep & Y<ep+Y1);

G_tracks(i,:)=mean(G(ls(:,i+1),:));

%quiver(X(ls(:,i+1)),Y(ls(:,i+1)),transpose(dUdt(ls(:,i+1))),transpose(dVdt(ls(:,i+1))),0.1,'blue')
     hold on
     quiver(mean(X(ls(:,i+1))),mean(Y(ls(:,i+1))),mean(dUdt(ls(:,i+1))),mean(dVdt(ls(:,i+1))),0.05,'LineWidth',1,'MaxHeadSize',1,'color','blue')
     hold on
     
     X1=mean(X(ls(:,i+1))+dt1*mean(dUdt(ls(:,i+1))));
     Y1=mean(Y(ls(:,i+1))+dt1*mean(dVdt(ls(:,i+1))));
     %scatter(X1,Y1)
     hold on

end
hold on


%Sleep wake switch-------constant velocity

alpha=83.3; 
beta=770;

Gw2=G(704714,:);%704715,G_trackw(end,:);
Gs1=G(355953,:);%356829,mean(G(ls1(:,1),:));%G_tracks(1,:);
ws=ndlinspace(Gw2,Gs1,100);

Xws=ws(:,1)./(1-ws(:,2));
Yws=(ws(:,3).*ws(:,4)+ws(:,3).*ws(:,5).*ws(:,7))./((1-ws(:,5).*ws(:,8)).*(1-ws(:,2)));
Zws=(ws(:,5).*ws(:,8)).*(alpha*beta)/(alpha+beta)^2;
%Y=(G_ese+G_erse)/((1-G_srs)*(1-G_ei));
%Z=-G_srs*(alpha*beta)/(alpha+beta)^2;

Gw1=G(666269,:);%mean(G(lw1(:,1),:));%G_trackw(1,:);
Gs2=G(206279,:);%G_tracks(end,:);
sw=ndlinspace(Gs2,Gw1,100);

Xsw=sw(:,1)./(1-sw(:,2));
Ysw=(sw(:,3).*sw(:,4)+sw(:,3).*sw(:,5).*sw(:,7))./((1-sw(:,5).*sw(:,8)).*(1-sw(:,2)));
Zsw=(sw(:,5).*sw(:,8)).*(alpha*beta)/(alpha+beta)^2;

plot(Xws,Yws,'black');
hold on
plot(Xsw,Ysw,'black');

%----------------------------------------


hold on
%quiver(X,Y,transpose(dUdt),transpose(dVdt),0.3,'Color',[0.7,0.7,0.7])
%quiver(X,Y,transpose(dXdt_sm),transpose(dYdt_sm),0.3,'color',[0.7,0.7,0.7])
% tent.draw_blobs({'ec','n1','n3'},0.1)
% tent.surface
% xlabel('X')
% ylabel('Y')
% view(2)

figure
subplot(1,2,1)
for i=1:3 
fitx=linspace(1,13,100);    
fity=interp1([1:13],G_trackw(:,i),fitx,'spline');
fitx2=linspace(0,16,100);
h=plot(fitx2,(fity),'LineWidth',1.5,'LineSmoothing','on');
c=h.Color;

hold on
line([16,17],[fity(end),G_tracks(1,i)],'LineStyle','--','color',c)
%scatter([1:15],G_trackw(:,i));
hold on
line([24,25],[G_tracks(13,i),G_trackw(1,i)],'LineStyle','--','color',c)
hold on

fitx3=linspace(25,43,100);
h3=plot(fitx3,fity,'Color',c,'LineWidth',1.5);
%end
%legend('G_{wee}','G_{ei}','G_{es}','G_{se}','G_{sr}','G_{sn}','G_{re}','G_{rs}')

hold on
% figure
%for i=1:1
fitx=linspace(1,15,100);    
fity=interp1([1:15],G_tracks(:,i),fitx,'spline');
fitx1=linspace(17,24,100);
plot(fitx1,fity,'Color',c,'LineWidth',1.5)
hold on
%scatter([1:15],G_tracks(:,i));
hold on
axis([0,42,-20,15]);
end


for i=8 
fitx=linspace(1,13,100);    
fity=interp1([1:13],G_trackw(:,i),fitx,'spline');
fitx2=linspace(0,16,100);


h1=plot(fitx2,(fity),'LineWidth',1.5,'LineSmoothing','on');
c=h1.Color;



hold on
fitx3=linspace(25,42,100);
h3=plot(fitx3,fity,'Color',c,'LineWidth',1.5);

%scatter([1:15],G_trackw(:,i));
hold on

%end
%legend('G_{wee}','G_{ei}','G_{es}','G_{se}','G_{sr}','G_{sn}','G_{re}','G_{rs}')

hold on
% figure
%for i=1:1
fitx=linspace(1,15,100);    
fity=interp1([1:15],G_tracks(:,i),fitx,'spline');
fitx1=linspace(17,24,100);
h2=plot(fitx1,fity,'Color',c,'LineWidth',1.5);
hold on
%scatter([1:15],G_tracks(:,i));

line([24,25],[fity(end),0.78],'LineStyle','-','color',c)
hold on
line([16,17],[G_trackw(end),G_tracks(1,i)],'LineStyle','-','color',c)
hold on
%h4=plot(fitx3,fity,'Color',c,'LineWidth',1.5);

end

hold on

patch([0 16 16 0],[-20 -20 15 15],'white','FaceAlpha',0.05,'EdgeColor','none')
hold on
patch([17 24 24 17],[-20 -20 15 15],'blue','FaceAlpha',0.05,'EdgeColor','none')
hold on
patch([25 43 43 25],[-20 -20 15 15],'white','FaceAlpha',0.05,'EdgeColor','none')
xlabel('Time (h)','FontSize',20)
ylabel('G_{ab}','FontSize',22)
l1=legend('G_{ee}','G_{ei}','G_{es}','G_{se}','G_{sr}','G_{sn}','G_{re}','G_{rs}','','','','','','');
l1.FontSize=12;
set(gca,'FontSize',16)

subplot(1,2,2)
for i=4:7 
fitx=linspace(1,13,100);    
fity=interp1([1:13],G_trackw(:,i),fitx,'spline');
fitx2=linspace(0,16,100);
fitx4=linspace(25,43,100);
h=plot(fitx2,(fity),'LineWidth',1.5,'LineSmoothing','on','LineStyle','-.');

c=h.Color;
hold on
h4=plot(fitx4,fity,'Color',c,'LineWidth',1.5,'LineStyle','-.');
hold on
l1=line([16,17],[G_trackw(13,i),G_tracks(1,i)],'LineStyle','--','color',c);
hold on
l2=line([24,25],[G_tracks(13,i),G_trackw(1,i)],'LineStyle','--','color',c);
%scatter([1:15],G_trackw(:,i));
hold on

%end
%legend('G_{wee}','G_{ei}','G_{es}','G_{se}','G_{sr}','G_{sn}','G_{re}','G_{rs}')

hold on
% figure
%for i=1:1
fitx=linspace(1,15,100);    
fity=interp1([1:15],G_tracks(:,i),fitx,'spline');
fitx1=linspace(17,24,100);
plot(fitx1,fity,'Color',c,'LineWidth',1.5,'LineStyle','-.')

%scatter([1:15],G_tracks(:,i));
hold on

end
hold on

patch([0 16 16 0],[-10 -10 20 20],'white','FaceAlpha',0.05,'EdgeColor','none')
hold on
patch([17 24 24 17],[-10 -10 20 20],'blue','FaceAlpha',0.05,'EdgeColor','none')
hold on
patch([25 43 43 25],[-10 -10 20 20],'white','FaceAlpha',0.05,'EdgeColor','none')
axis([0,42,-10,20]);
xlabel('Time (h)','FontSize',20)
ylabel('G_{ab}','FontSize',22)
set(gca,'FontSize',16)
l2=legend('G_{ee}','G_{ei}','G_{es}','G_{se}','G_{sr}','G_{sn}','G_{re}','G_{rs}','','','','','','');
l2.FontSize=12;



