traj=open('/suphys/sahanda/fitted_trajectories.mat');

dGee=diff(traj.T.Gee);
dGee(dGee==0) = [];
Gee=traj.T.Gee;
Gee(diff(traj.T.Gee)==0)=[];
Gee(end)=[];

%hist(dGee,1000);

dGei=diff(traj.T.Gei);
dGei(dGei==0) = [];
Gei=traj.T.Gei;
Gei(diff(traj.T.Gei)==0)=[];
Gei(end)=[];



%hist(dGei,1000);


 %dX=((1-Gei).*dGee+Gee.*dGei)./(1-Gei).^2;
 
 dX=diff(traj.T.X);
 dX(dX==0)=[];
 X=traj.T.X;
 X(diff(traj.T.X)==0)=[];
 X(end)=[];
 
  dY=diff(traj.T.Y);
 dY(dY==0)=[];
 Y=traj.T.Y;
 Y(diff(traj.T.Y)==0)=[];
Y(end)=[];
 
 
 
 dZ=diff(traj.T.Z);
 dZ(dZ==0)=[];
 Z=traj.T.Z;
 Z(diff(traj.T.Z)==0)=[];
 Z(end)=[]; 
 for i=1:1000:length(X)
     
 comet(X(i),Y(i));
 hold on
 pause(0.05)
 end
 %comet(X,)