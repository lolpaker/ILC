clc, clear all, close all
%transmitancja obiektu AMB
% num = [85.0824]; 
% den = [3 0 -40.2436];
% num = [435440]; 
% den = [3 0 -40.2436];
num = [40.2436]; 
den = [3 0 -435440];
Ts = 0.00001;
sysc = tf(num,den); %postac ciagla
sysd = c2d(sysc,Ts,'ZOH'); %postac dyskretna

%parametry regulatora PID
% kp=5; 
% ki=10; 
% kd=0.7; 

% sysd_PID=pid(kp,ki,kd);
sysd_PID=pidtune(sysd,'pidf');
sysd_PID.Ts=0.00001;

sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zast?pczej uk?adu z regulatorem PID
sysd_f=feedback(sysd_se,1); 
% figure(2)
% step(sysd_f);
%%
    [Ad Bd Cd Dd] = ssdata(sysd_f); %zamiana transmitancji na model przestrzeni stanu
    [Ad_PID Bd_PID Cd_PID Dd_PID] = dssdata(sysd_PID); %zamiana transmitancji na model przestrzeni stanu regulatora PID
     
% Warunki poczatkowe x0, zakres czasu t, opoznienie n0, r, macierz N
    
    x0 = 0;
    t = 0:Ts:0.05;
    n0 = 0;
    r = 1;
    N = length(t);
    % Definiowanie wejsciowego wektora U i sygna?u wartosci zadanej Rj
    rj_step=-0.0002;
%     rj_step=0;
    xx=t;
%    Rj = (-9.477*xx.^6+28.43*power(xx,5)-29.34*power(xx,4)+11.29*power(xx,3)-0.926*power(xx,2)+0.016*power(xx,1))';

  % Rj=[0  1.99999999994649e-06*power(xx,6) 1.09999999999832e-05*power(xx,5)  0.000105999999999995*power(xx,4) 0.000241999999999964*power(xx,3)  0.00563699999999995*power(xx,2) 0.035937*power(xx,1) 0]';

     Rj = [rj_step*ones(1,200) 0*ones(1,4801)]';

    % Sformu?owanie macierzy obiektu
    Gvec = zeros(N,1);
    rVec = ((r-1):(N-n0-1))';
    
    for ii = 1:length(rVec)
      ApowVec = Ad^rVec(ii);
      Gvec(ii) = Cd*ApowVec*Bd;
    end
    G = tril(toeplitz(Gvec)); %utworzenie macierzy toeplitza i dolnotrójk?tnej obiektu G
    
    for ii_PID = 1:length(rVec)
      ApowVec_PID = Ad_PID^rVec(ii_PID);
      Gvec_PID(ii_PID) = Cd_PID*ApowVec_PID*Bd_PID;
    end
    G_PID = tril(toeplitz(Gvec_PID)); %utworzenie macierzy toeplitza i dolnotrójk?tnej 
                                      %transmitancji regulatora GPID

                                      
    % czestotliwosci filtrów
    f1=400;
    f2=500;
    f3=700;
                                      
    jmax = 50;
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,f1*(2/Fs),'low'); 
        %Yj - wyj?cie z uk?adu regulacji
        %Uj - sygna? wyj?ciowy z ILC
        %Ej - sygna? b??du
        %Rj - warto?? zadana na wymuszeniu uk?adu
        %G - transmitancja obiektu
        %Q - macierz, filtr
        %L - macierz ucz?ca PD
    % Run ILC and plot the response for each iteration
    for ii = 1:jmax
      Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn? próbk? czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
    Ej_F1(ii,1)=plotter_Ej(ii,t,Ej);
%     Ej_F1_max(ii,1)=plotter_Ej_max(ii,Ej);
%     plotter(ii,t,Ej,Yj,Uj,Rj); 
%         plotter_errors(ii,Ej);
    end
    Yj_w1=Yj;
%%     
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    [bB,aB] = butter(4,f2*(2/Fs),'low'); 
        %Yj - wyj?cie z uk?adu regulacji
        %Uj - sygna? wyj?ciowy z ILC
        %Ej - sygna? b??du
        %Rj - warto?? zadana na wymuszeniu uk?adu
        %G - transmitancja obiektu
        %Q - macierz, filtr
        %L - macierz ucz?ca PD
    % Run ILC and plot the response for each iteration
    for ii = 1:jmax
      Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn? próbk? czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
       Ej_F2(ii,1)=plotter_Ej(ii,t,Ej);  
%        Ej_F2_max(ii,1)=plotter_Ej_max(ii,Ej);
%     plotter(ii,t,Ej,Yj,Uj,Rj); 

%        Ej_F2(ii,1)=plotter(ii,t,Ej,Yj,Uj,Rj);
%         plotter_errors(ii,Ej);
    end
    Yj_w2=Yj;
%%     
        Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    [bB,aB] = butter(4,f3*(2/Fs),'low'); 
        %Yj - wyj?cie z uk?adu regulacji
        %Uj - sygna? wyj?ciowy z ILC
        %Ej - sygna? b??du
        %Rj - warto?? zadana na wymuszeniu uk?adu
        %G - transmitancja obiektu
        %Q - macierz, filtr
        %L - macierz ucz?ca PD
    % Run ILC and plot the response for each iteration
    for ii = 1:jmax
      Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn? próbk? czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
    Ej_F3(ii,1)=plotter_Ej(ii,t,Ej);
%     Ej_F3_max(ii,1)=plotter_Ej_max(ii,Ej);
%       Ej_F3(ii,1)=plotter(ii,t,Ej,Yj,Uj,Rj);
      
%          plotter(ii,t,Ej,Yj,Uj,Rj); 

%         plotter_errors(ii,Ej);
    end
    Yj_w3=Yj;

%%    
figure()
plot(t,Yj_w1,'b',t,Yj_w2,'k',t,Yj_w3,'r');
grid on
% legend(strcat('f1= Hz',num2str(f1)),strcat('f2=',num2str(f2)),strcat('f3=',num2str(f3)))
legend('f1=400Hz', 'f2=500Hz', 'f3=700Hz');
title('OdpowiedŸ uk³adu na wymuszenie pr¹dowe');
xlabel('Czas [s]');
ylabel('Przemieszczenie wirnika [m]');
xlim([0 0.03]);
set(gcf,'Color',([1 1 1]));
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(findall(gcf,'type','line'),'linewidth',2);

%%
    figure(922)
    iteration = 1:1:length(Ej_F1);
     plot(iteration,Ej_F1,'b');
    hold on
%     plot(iteration,Ej_F1_max);
    plot(iteration,Ej_F2,'k');  
%     plot(iteration,Ej_F2_max);
    plot(iteration,Ej_F3,'r');
%     plot(iteration,Ej_F3_max);
%     plot(iteration,Ej_FPoiprzed,'c+')
%     legend('Filtr 1', 'Filtr 1max','2','2max','3','3max');
%      title('B³¹d RMSE')
    legend('f1=400Hz', 'f2=500Hz', 'f3=700Hz');
    xlabel('Iteracje');
%     ylabel('Przemieszczenie wirnika [m]');
    ylabel('B³¹d RMSE [m]');
    set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
    set(gca, 'YScale', 'log');
    set(findall(gcf,'type','line'),'linewidth',2)
    grid on