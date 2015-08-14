%plots

set(gcf, 'renderer', 'painters')
plot(w/(2*pi),I,'LineWidth',1,'linesmoothing','on')
xlabel('f(Hz)')
ylabel('Integral')
%legend('','','','','','','','','','')
