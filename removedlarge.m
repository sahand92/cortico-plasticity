A=(find(radius<0.1));
XX=X(A);
YY=Y(A);
UU=dUdt(A);
VV=dVdt(A);
dXXdt=dXdt(A);
dYYdt=dYdt(A);
rr=sqrt(real(dYYdt).^2+real(dXXdt).^2);
c=rr*10;
for i=1:length(XX)
quiver(XX(i),YY(i),UU(i),VV(i),0.02,'color',[1-c(i),c(i),0.5],'AutoScale','off','LineWidth',2,'MaxHeadSize',1)
%quiver(x(i),y(i),u(i),v(i),0.1,'color',[r(i,1),r(i,2),r(i,3)])
hold on
end
% %quiver(XX,YY,UU,VV,0.5)