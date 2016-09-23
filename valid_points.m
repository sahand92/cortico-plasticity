% clear ls
% clear ls2
% ls=find(X>0.8 & X<1 & Y>-0.2 & Y<0.1 & dGredt>0 & dGeedt<0 & dGeidt<0 & dGeedt<0 & dGesdt>0 & dGsedt<0 & dGsrdt<0 & dGrsdt>0);
% for k=1:length(ls)
% 
% ls2=find(G(ls,7)>G(ls(k),7) & G(ls,1)<G(ls(k),1) & G(ls,2)<G(ls(k),2) & G(ls,3)>G(ls(k),3) & G(ls,4)<G(ls(k),4) & G(ls,5)<G(ls(k),5) & G(ls,8)>G(ls(k),8));
% 
% scatter(X(ls2),Y(ls2),'blue')
% hold on
% scatter(X(ls(k)),Y(ls(k)),'red')
% hold on
% end

%k=6009;
k=6340;
G_fw=G(k,:)+0.0005*real(dGdtr(k,:));

nu_fw=Nr(k,:).*(real(Ilastr(k,:)))*0.0005+nu(k,:);



X_fw=G_fw(:,1)./(1-G_fw(:,2));
Y_fw=(G_fw(:,3).*G_fw(:,4) + G_fw(:,3).*G_fw(:,5).*G_fw(:,7))./((1-G_fw(:,5).*G_fw(:,8)).*(1-G_fw(:,2)));

scatter(X_fw,Y_fw,'blue')
hold on
scatter(X(k),Y(k),'red')
hold on
quiver(X(k),Y(k),transpose(dUdt(k)),transpose(dVdt(k)),0.3,'Color',[0.7,0.7,0.7])