function [] = plotter(ii,t,Ej,Yj,Uj,Rj)
  figure(1)
%  Wykres b��du Ej dla danej iteracji
  subplot(1,3,1);
  plot(t,Ej,'LineWidth',2);
  title('B��d, Ej','FontSize',12);
  ylabel('Przemieszczenie [m]','FontSize',12);
%  ylim([-5 5])
%  Wykres sygna�u steruj�cego Uj dla danej iteracji
  subplot(1,3,2);
  plot(t,Uj,'-k','LineWidth',2);
  title({['Iteracja: ', num2str(ii)],'Sygna� steruj�cy, Uj'},'FontSize',12);
  xlabel('Czas (s)','FontSize',12);
 ylabel('Nat�enie [A]','FontSize',12);
%  ylim([0 7])
%  Wykres sygna�u wyj�ciowego Yj dla danej iteracji
  subplot(1,3,3);
  plot(t,Yj,t,Rj,'-k','LineWidth',2);
  title('Sygna� wyj�ciowy, Yj','FontSize',12);
 ylabel('Przemieszczenie [m]','FontSize',12);
 % ylim([0 7])
   pause(0.02);
%   figure(76)
%   plot(t,Yj,t,Rj,'-k','LineWidth',2);
%   title('Sygna� wyj�ciowy','FontSize',12);
%   ylabel('Przemieszczenie wirnika [m]','FontSize',12);
%   set(gcf,'Color',([1 1 1]));
%   set(findall(gcf,'-property','FontSize'),'FontSize',18);
% 
%  xlim([0 0.01])
end