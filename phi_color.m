
X=xyz_final(1:10000,1);
Y=xyz_final(1:10000,2);
Z=xyz_final(1:10000,3);
phie=phi_final(1:10000,1);
% X=X(index);
% Y=Y(index);
% phie=phie(index);
for i=1:length(phie)
    scatter3(X(i),Y(i),phie(i),50,'MarkerFaceColor',[1-phie(i)/max(phie),1-phie(i)/max(phie),phie(i)/max(phie)],'MarkerEdgeColor',[1-phie(i)/max(phie),1-phie(i)/max(phie),phie(i)/max(phie)]);
hold on
end
