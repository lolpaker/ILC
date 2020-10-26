function [Ej_max] = plotter_Ej_max(ii,Ej)
  figure(235)
%  Wykres b³êdu Ej max dla danej iteracji
    Ej_max=max(abs(Ej));
%     Ej_max=max(Ej);
    plot(ii,Ej_max,'*','LineWidth',1.5);
    title('B³¹d, Ej max','FontSize',12);
    xlabel('Iteracja','FontSize',12);
    ylabel('B³¹d maksymalny[\mum]','FontSize',12);
    set(gca, 'YScale', 'log'); 
    set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'-property','FontSize'),'FontSize',18);

    hold on
  
end