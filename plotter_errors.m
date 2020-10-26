function [Ej_max,Ej_mse,Ej_rmse] = plotter_errors(ii,Ej)
  figure(236)
    Ej_sum=0;
    for i2=1:1:length(Ej)
        Ej_pow= power(Ej(i2,1),2);
        Ej_sum = Ej_sum+Ej_pow;
    end
    
    Ej_mse=Ej_sum/length(Ej);
    Ej_rmse= sqrt(Ej_mse);
    Ej_max=max(Ej);

    plot(ii,Ej_max,'b*','LineWidth',2.5);
    hold on
%     plot(ii,Ej_mse,'go','LineWidth',2.5);
    plot(ii,Ej_rmse,'r+','LineWidth',2.5);
    legend('B≥πd MAX','B≥πd RMSE')%'B≥πd MSE',
    title('Odpowiedü uk≥adu na wymuszenie');
    xlabel('Iteracje');
    ylabel('Przemieszczenie wirnika [m]');
    set(gca, 'YScale', 'log');
    set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
%     grid on
end

 