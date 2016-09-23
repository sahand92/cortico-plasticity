

% scatter(nus_data(:,1),nus_data(:,2),'.')
% hold on
% 
% for m=1:100
% 
% u(m)=find(dist==inde(m));
% 
% plot(nus_data(u,1),nus_data(u,2),'red')
% 
% hold on
% 
% end



% 
% 
scatter(xyz_data(:,1),xyz_data(:,2),'.')
hold on

for m=1:100

u(m)=find(dist==inde(m));

plot(xyz_data(u,1),xyz_data(u,2),'red')

hold on

end