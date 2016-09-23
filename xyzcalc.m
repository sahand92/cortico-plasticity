
clear
 wgab = [2.074,-4.110,0.772,7.768,-3.301,8.097,0.656,0.196];
 Gint(5,:)=wgab;
 
 
 
 Gint(5,:)=Gint(5,:)+rand(size(Gint(5,:)))/5;
 Gint(2,:)=Gint(5,:)-rand(size(Gint(5,:)))/5;
 Gint(3,:)=Gint(5,:)+rand(size(Gint(5,:)))/5;
 Gint(4,:)=Gint(5,:)-rand(size(Gint(5,:)))/5;
 Gint(1,:)=Gint(5,:);
Gint(6,:)=Gint(5,:)+rand(size(Gint(5,:)))/5;
Gint(7,:)=Gint(5,:)-rand(size(Gint(5,:)))/5;
Gint(8,:)=Gint(5,:)+rand(size(Gint(5,:)))/5;
Gint(9,:)=Gint(5,:)-rand(size(Gint(5,:)))/5;
Gint(10,:)=Gint(5,:)+rand(size(Gint(5,:)))/5;


for p=1:1
%initial G
G=Gint(p,:);
G_ee=G(1);
G_ei=G(2);
G_es=G(3);
G_re=G(7);
G_rs=G(8);
G_sr=G(5);
G_sn=G(6);
G_se=G(4);
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
X=G_ee/(1-G_ei);
Y=(G_ese+G_erse)/((1-G_srs)*(1-G_ei));
Z=-G_srs*(alpha*beta)/(alpha+beta)^2;

scatter(X,Y)
hold on
end
hold on
patch([1.2 1.2 0],[-0.2 1 1],[0.9 0.9 0.9],'EdgeAlpha',0.2)