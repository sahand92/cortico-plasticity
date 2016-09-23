
load('/suphys/sahanda/cortico plasticity/data/pdb_wake.mat');
Xw=xyz_final(1:100000,1);
Yw=xyz_final(1:100000,2);
Zw=xyz_final(1:100000,3);
Gw=gab_final(1:100000,:);
nuw=nus_final(1:100000,:);






lw1=find(Xw>0.10 & Xw<0.2 & Yw>0.6 & Yw < 0.75);
lw2=find(Xw>0.15 & Xw<0.45 & Yw>0.7 & Yw<0.85);
Gw1m=[mean(Gw(lw1,1)),mean(Gw(lw1,2)),mean(Gw(lw1,3)),mean(Gw(lw1,4)),mean(Gw(lw1,5)),mean(Gw(lw1,6)),mean(Gw(lw1,7)),mean(Gw(lw1,8))];
Gw2m=[mean(Gw(lw2,1)),mean(Gw(lw2,2)),mean(Gw(lw2,3)),mean(Gw(lw2,4)),mean(Gw(lw2,5)),mean(Gw(lw2,6)),mean(Gw(lw2,7)),mean(Gw(lw2,8))];

X1w=Gw1m(:,1)./(1-Gw1m(:,2));
Y1w=(Gw1m(:,3).*Gw1m(:,4) + Gw1m(:,3).*Gw1m(:,5).*Gw1m(:,7))./((1-Gw1m(:,5).*Gw1m(:,8)).*(1-Gw1m(:,2)));

X2w=Gw2m(:,1)./(1-Gw2m(:,2));
Y2w=(Gw2m(:,3).*Gw2m(:,4) + Gw2m(:,3).*Gw2m(:,5).*Gw2m(:,7))./((1-Gw2m(:,5).*Gw2m(:,8)).*(1-Gw2m(:,2)));



load('/suphys/sahanda/cortico plasticity/data/pdb_sleep.mat');
Xs=xyz_final(1:100000,1);
Ys=xyz_final(1:100000,2);
Zs=xyz_final(1:100000,3);
Gs=gab_final(1:100000,:);
nus=nus_final(1:100000,:);



ls1=find(Xs>0.8 & Xs<1 & Ys>-0.15 & Ys<-0.0);
ls2=find(Xs>0.75 & Xs<0.85 & Ys>-0.2 & Ys<-0.15);
Gs1m=[mean(Gs(ls1,1)),mean(Gs(ls1,2)),mean(Gs(ls1,3)),mean(Gs(ls1,4)),mean(Gs(ls1,5)),mean(Gs(lw1,6)),mean(Gs(ls1,7)),mean(Gs(ls1,8))];
Gs2m=[mean(Gs(ls2,1)),mean(Gs(ls2,2)),mean(Gs(ls2,3)),mean(Gs(ls2,4)),mean(Gs(ls2,5)),mean(Gs(lw2,6)),mean(Gs(ls2,7)),mean(Gs(ls2,8))];

X1s=Gs1m(:,1)./(1-Gs1m(:,2));
Y1s=(Gs1m(:,3).*Gs1m(:,4) + Gs1m(:,3).*Gs1m(:,5).*Gs1m(:,7))./((1-Gs1m(:,5).*Gs1m(:,8)).*(1-Gs1m(:,2)));

X2s=Gs2m(:,1)./(1-Gs2m(:,2));
Y2s=(Gs2m(:,3).*Gs2m(:,4) + Gs2m(:,3).*Gs2m(:,5).*Gs2m(:,7))./((1-Gs2m(:,5).*Gs2m(:,8)).*(1-Gs2m(:,2)));

scatter(X1w,Y1w,'red')
hold on
scatter(X2w,Y2w,'blue')


hold on


scatter(X1s,Y1s,'red')
hold on
scatter(X2s,Y2s,'blue')
hold on


scatter(Xw(lw1),Yw(lw1),'red')
hold on
scatter(Xw(lw2),Yw(lw2),'blue')
hold on
scatter(Xs(ls1),Ys(ls1),'red')
hold on
scatter(Xs(ls2),Ys(ls2),'blue')

hold on
tent.draw_blobs({'ec','n1'},0.1)