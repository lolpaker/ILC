function [Ej_rmse] = plotter_Ej(ii,t,Ej)
  figure(234)
%  Wykres b³êdu Ej dla danej iteracji
  
    Ej_sum=0;
    for i2=1:1:length(Ej)
        Ej_pow= power(Ej(i2,1),2);
        Ej_sum = Ej_sum+Ej_pow;
    end
    
     Ej_mse=Ej_sum/length(Ej);
    Ej_rmse= sqrt(Ej_mse);
    
  plot(ii,Ej_rmse,'*','LineWidth',8);
  title('B³¹d RMSE','FontSize',12);
  xlabel('Iteracje','FontSize',12);
  ylabel('Ej_RMS','FontSize',12);
  set(gca, 'YScale', 'log'); 
  set(gcf,'Color',([1 1 1]));
  set(findall(gcf,'-property','FontSize'),'FontSize',18)

   hold on
%    grid on
  
end