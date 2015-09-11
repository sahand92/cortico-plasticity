clear
load('/suphys/sahanda/cortico plasticity/points.mat');  %data of 3000 points containing positions, velocities, Gabs and Nus
X=S.xyz(:,1);
Y=S.xyz(:,2);
Z=S.xyz(:,3);
U=transpose(S.u);
V=transpose(S.v);
W=transpose(S.w);
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
n=0.05;  %square size is n x n
%index=zeros(1,27);



for i=0:ceil(0.9/n)
    for j=0:ceil(1.8/n)
   
l=find(X>0+n*i&X<0+n*(i+1) & Y>-1+n*j&Y<-1+n*(j+1),1);


    

%  radius=sqrt(real(U).^2+real(V).^2);
%     dUdt=(real(U))./real(radius);
%      dVdt=(real(V))./real(radius);    
%      
% 
%   

% 
% M2=[X(l) Y(l)];
scatter(X(l),Y(l),'.')
% quiver(X(l),Y(l),dUdt(l),dVdt(l),0.05)
 hold on
 if isempty(l)==1;
     l=nan;  
 end 
 disp(l)
 
 
 
    end
    
%hold on

end


