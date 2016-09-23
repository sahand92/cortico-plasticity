
function plot_integrand(states)

  
   
	states_all = {'eo','ec','rem','n1','n2','n3','n2s'};
  
    
    if nargin < 1 || isempty(states)
		states = states_all;
	end
    
    
    ind=find(ismember(states_all,states));
    
 
  %Load Data---------------
  data=load('data.mat');
  I=data.I;
  f=data.f;
  P=data.Pow;
  X=data.X;
  Y=data.Y;
  Z=data.Z;
  In=data.In;
  %------------------------
  
  alpha_freq=data.alpha_freq;
  
  for j=1:length(ind)
      for i=1:8
      Id=reshape((I(ind(j),i,:)),[1 length(f)]);
      ID(i,:)=Id;
      end
      
 figure('Position',[100 100,1500,600])
 subplot(1,2,1)
    % ti=linspace(sinh(-10)/10^9,sinh(10)/10^9,1000);
  ti=[10^6 10^5 10^4 10^3 10^2 10^1 10^0 0 -10^0 -10^1 -10^2 -10^3 -10^4 -10^5 -10^6] ;
     for k=1:13
     y1=asinh(ti);
     line([0,1],[y1(k),y1(k)],'color','black')
     tis=ti(k).*(1:10);
     y2=asinh(tis);
     for m=1:10
     line([0,0.2],[y2(m),y2(m)],'color','black')
     end
     %line([0,80],[0,0],'color','black')
     %y2(k);
    
     end
     size(y1);
     
     hold on
  %h=plot(f,10^6.*ID);
%  hold on
 h=plot(f,asinh(10^10*ID));
 set(h,'LineWidth',1.3);
 set(h(8),'LineStyle','--');
 set(h(7),'LineStyle','--');
 set(h(6),'LineStyle','-.');
 %l=legend('I_{ee}','I_{ei}','I_{es}','I_{se}','I_{sr}','I_{sn}','I_{re}','I_{rs}','Location','bestoutside');
 %l.FontSize=18;
 
 dim=[0.2 0.6 0.3 0.3];
 annotation('textbox',dim,'String',states(j),'FitBoxToText','on','FontSize',20);
 %Tickmarks=scale
 
  ylim([-6 6]);   
   
 xlabel('f (Hz)','FontSize',20)
 ylabel('Integrand','FontSize',22)
 set(gca,'ytick',[]);
 set(gca,'yticklabel',[]);
 set(gca,'FontSize',16)
 set(gca,'XTick',0:5:50);
 set(gca,'XTickLabel',{'0',[] ,'10',[], '20',[], '30',[], '40',[],'50'});
 xlim([0 50]);
 box on
 
 hold all
 
 handaxes2 = axes('Position', [0.5 0.3 0.2 0.1]);

  loglog(f,P(j,:))
 xlabel('f (Hz)')
 ylabel('P (s^{-1})')
 hold on
 plot(alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
 hold on
 plot(2*alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
 hold on
 plot(3*alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
 text(alpha_freq+1,10^-8*1.1,'f_\alpha')
 text(2*alpha_freq+1,10^-8*1.1,'2f_\alpha')
  text(3*alpha_freq+1,10^-8*1.1,'3f_\alpha')
 set(gca,'XLim',[1 50],'YLim',[min(P((ind(j)),:)) max(P(ind(j),:))]);
 
  subplot(1,2,2)
   for i=1:8
      Id=reshape((In(ind(j),i,:)),[1 length(f)]);
      ID(i,:)=Id;
      end
 g=plot(f,asinh(ID*10^10));
 set(g,'LineWidth',1.3);
 set(g(8),'LineStyle','--');
 set(g(7),'LineStyle','--');
 set(g(6),'LineStyle','-.');
 l=legend('I_{ee}','I_{ei}','I_{es}','I_{se}','I_{sr}','I_{sn}','I_{re}','I_{rs}','Location','bestoutside');
 l.FontSize=18;
 dim=[0.4 0.6 0.3 0.3];
 %annotation('textbox',dim,'String',states(j),'FitBoxToText','on','FontSize',20);
         
  ti3=[10^6 10^5 10^4 10^3 10^2 10^1 10^0 0 -10^0 -10^1 -10^2 -10^3 -10^4 -10^5 -10^6] ;
     for k=1:13
     y3=asinh(ti3);
     line([0,1],[y3(k),y3(k)],'color','black')
     tis=ti3(k).*(1:10);
     y4=asinh(tis);
     for m=1:10
     line([0,0.5],[y4(m),y4(m)],'color','black')
     end
 xlabel('f_{max} (Hz)','FontSize',20)
 ylabel('Integral','FontSize',22)
 set(gca,'XTick',0:5:50);
 set(gca,'XTickLabel',{'0',[] ,'10',[], '20',[], '30',[], '40',[],'50'});
  xlim([0 45]);
  ylim([-11 11]);
   set(gca,'ytick',[]);
 set(gca,'yticklabel',[]);


 set(gca,'FontSize',16)
     end
  

 
 figure
 tent.draw_blobs({'n3','ec'},'0.1')
 hold on
 scatter(X,Y,'filled','red')
 t=text(X,Y+0.02,states);
 set(t,'FontSize',12);
 tent.surface
 xlabel('X')
 ylabel('Y')
 view(2)
 
 
 
% hold on
% 
% xlabel('f (Hz)','FontSize',20)
% ylabel('Integrand','FontSize',22)
% dim=[0.4 0.6 0.3 0.3];
% 
%  str='EC';
% 
% annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',20);
% set(gca,'FontSize',16)
% legend('I_{ee}','I_{ei}','I_{es}','I_{re}','I_{rs}','I_{se}','I_{sr}','I_{sn}','Location','bestoutside')
% 
% %legend('ds_{ee}/dt','ds_{ei}/dt','ds_{es}/dt','ds_{re}/dt','ds_{rs}/dt','ds_{se}/dt','ds_{sr}/dt','ds_{sn}/dt')
% box off
% handaxes2 = axes('Position', [0.5 0.3 0.2 0.1]);
% 
%   loglog(f,P)
%  xlabel('f (Hz)')
%  ylabel('P (s^{-1})')
%  hold on
%  plot(alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
%  hold on
%  plot(2*alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
%  hold on
%  plot(3*alpha_freq*ones(10),linspace(10^-8,10^-2,10),'--','color','green')
%  text(alpha_freq+1,10^-8*1.1,'f_\alpha')
%  text(2*alpha_freq+1,10^-8*1.1,'2f_\alpha')
%   text(3*alpha_freq+1,10^-8*1.1,'3f_\alpha')
%  set(gca,'XLim',[1 50],'YLim',[min(P) max(P)]);
 
end