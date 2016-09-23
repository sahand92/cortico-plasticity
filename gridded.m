% X1=linspace(0,1,1259);
% Y1=linspace(-1,1,1259);
% F=griddedInterpolant(X1,Y1);
% x=F(X);
% y=F(Y);
% u=F(real(dUdt));
% v=F(real(dVdt));
% 
% quiver(x,y,u,v,0.5,'black')

X1=linspace(0,1,1000);
Y1=linspace(0,1,1000);
F=griddedInterpolant(X1,Y1,linear);
x=F(X);
y=F(Y);
u=F(real(dUdt));
v=F(real(dVdt));
%scatter(x,y,'.')
%r=rand(length(radius),3);

[xx, yy]=meshgrid(x, y);
[uu, vv]=meshgrid(u,v);
% for i=1:length(X)
% quiver(x(i),y(i),u(i),v(i),0.1,'color',[1,radius(i)/82,0])
% %quiver(x(i),y(i),u(i),v(i),0.1,'color',[r(i,1),r(i,2),r(i,3)])
% hold on
% end
% % hold on

verts=stream2(xx,yy,uu,vv,0:0.01:1,-1:0.02:1);
iverts=interpstreamspeed(xx,yy,uu,vv,verts);
streamline(iverts)

