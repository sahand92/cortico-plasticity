

%after loading feather object
Y=combined(:,2);
X=combined(:,1);
values=hist3([Y X],[50 50]);
values1=values;
values1(size(values,1)+1,size(values,2)+1)=0;

xb=linspace(min(X),max(X),size(values,1)+1);
yb=linspace(min(Y),max(Y),size(values,1)+1);

surf(xb,yb,values1);
shading interp
view(2)
%colormap('bone')
colormap('jet')
%colormap(flipud(colormap))
grid off
hold on

% hold on
% tent.draw_blobs({'ec','n1'},0.05)
% hold on
% tent.surface
% view(2)

% values=hist3([Y(ind1) X(ind1)],[50 50]);
% values1=values;
% values1(size(values,1)+1,size(values,2)+1)=0;
% 
% xb=linspace(min(X(ind1)),max(X(ind1)),size(values,1)+1);
% yb=linspace(min(Y(ind1)),max(Y(ind1)),size(values,1)+1);
% 
% surf(xb,yb,values1);
% shading interp
% grid off