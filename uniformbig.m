clear
load('/suphys/sahanda/cortico plasticity/data/pdb_wake.mat');
X=xyz_final(:,1);
Y=xyz_final(:,2);
Z=xyz_final(:,3);
phie=phi_final(:,1);
% U=transpose(S.u);
% V=transpose(S.v);
% W=transpose(S.w);
%scatter(X,Y,'.')
% newX=X(X<0.2&X>0);
% newY=Y(Y>0.4&Y<0.8);
% M=[X Y];
% 
% newX=X(M(:,1)<0.2 & M(:,2)>0.4)
% newY=Y(M(:,2)>0.4 & M(:,1)<0.2) 
% 
% 
% 
% %blocks (squares) are M
% M=zeros(size(M))
% M=[newX newY];
n=0.001;  %square size is n x n
for i=0:ceil(0.9/n)
    for j=0:ceil(1.8/n)
   
l=find(X>0+n*i&X<0+n*(i+1) & Y>-1+n*j&Y<-1+n*(j+1),1);
if isempty(l)==1;
     l=0;  
 end 
 %disp(l)
 
 index(i+1,j+1)=l;

% M2=[X(l) Y(l)];
% scatter(X(l),Y(l),'.')
% hold on
  end  

end
a=numel(index);
index=reshape(index,1,a);
index=index(index~=0);
clear i
clear j