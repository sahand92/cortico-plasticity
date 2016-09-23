clear

%  S=load('/import/ghrian1/sahanda/uniform_wake');
%  X=S.X;
%  Y=S.Y;
%  Z=S.Z;
%  G=S.G;
%  nu=S.nu;

load('/suphys/sahanda/cortico plasticity/data/pdb_wake.mat');
X=xyz_final(1:30000,1);
Y=xyz_final(1:30000,2);
Z=xyz_final(1:30000,3);
G=gab_final(1:30000,:);
nu=nus_final(1:30000,:);

n=0.01;  


for i=0:ceil(0.9/n)
    for j=0:ceil(1.8/n)
   
l=find(X>0+n*i&X<0+n*(i+1) & Y>-1+n*j&Y<-1+n*(j+1),1);


    

if isempty(l)==1;
     l=0;  
 end 
 %disp(l)
 
%  if l~=[];
%      x=l;
%  else
%      x=0;
%  end
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
scatter(X(index),Y(index),'.','black')

