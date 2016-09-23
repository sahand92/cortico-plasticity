% %do a smoothing first
% 

X1=ceil(10*(Xsort(1:1225)));
Y1=ceil(10*Ysort(1:1225));
dUdt1=dXdts(1:1225);
dVdt1=dYdts(1:1225);
% 
for n=1:35
X2(n,:)=[X1((1+35*(n-1):35*n))];
end
for n=1:35
Y2(n,:)=[Y1((1+35*(n-1):35*n))];
end
for n=1:35
dUdt2(n,:)=[dUdt1((1+35*(n-1):35*n))];
end
for n=1:35
dVdt2(n,:)=[dVdt1((1+35*(n-1):35*n))];
end
quiver(X2,Y2,dUdt2,dVdt2,0.5)
hold on
[sx,sy]=meshgrid(0:0.1:1,0:0.1:1);
fuck=stream2(X2,Y2,dUdt2,dVdt2,sx,sy);

hold on
streamline(fuck)

% 
% X1=Xsort(1:9);
% Y1=Ysort(1:9);
% dUdt1=dXdts(1:9);
% dVdt1=dYdts(1:9);
% 
% for n=1:3
% X2(n,:)=[X1((1+3*(n-1):3*n))];
% end
% for n=1:3
% Y2(n,:)=[Y1((1+3*(n-1):3*n))];
% end
% for n=1:3
% dUdt2(n,:)=[dUdt1((1+3*(n-1):3*n))];
% end
% for n=1:3
% dVdt2(n,:)=[dVdt1((1+3*(n-1):3*n))];
% end
% quiver(X2,Y2,dUdt2,dVdt2,0.2);
% [sx,sy]=meshgrid(0.1,-1:0.1:1);
% fuck=stream2(X2,Y2,dUdt2,dVdt2,sx,sy);
% streamline(fuck)
% hold on
% scatter(sx,sy)
% grid square
