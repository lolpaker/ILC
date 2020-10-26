clc, clear all, close all
%transmitancja obiektu AMB
% num = [85.0824]; 
% den = [3 0 -40.2436];
% num = [435440]; 
% den = [3 0 -40.2436];
num = [40.2436]; 
den = [3 0 -435440];
Ts = 0.00001;
sysc = tf(num,den); %postaæ ci¹g³a
sysd = c2d(sysc,Ts,'ZOH'); %postaæ dyskretna

%parametry regulatora PID
kp=5; 
ki=10; 
kd=0.7; 
sysd_PID=pidtune(sysd,'pidf');
% sysd_PID=pid(kp,ki,kd);
sysd_PID.Ts=0.00001;

sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zastêpczej uk³adu z regulatorem PID
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
    % Definiowanie wejsciowego wektora U i sygna³u wartosci zadanej Rj
    rj_step=-0.0002;
    Rj = [rj_step*ones(1,200) 0*ones(1,4801)]';% tu mo¿na zmieniæ wartoœæ zadan¹ 0*ones(1,20) 

    % Sformu³owanie macierzy obiektu
    Gvec = zeros(N,1);
    rVec = ((r-1):(N-n0-1))';
    
    for ii = 1:length(rVec)
      ApowVec = Ad^rVec(ii);
      Gvec(ii) = Cd*ApowVec*Bd;
    end
    G = tril(toeplitz(Gvec)); %utworzenie macierzy toeplitza i dolnotrójk¹tnej obiektu G
    
    for ii_PID = 1:length(rVec)
      ApowVec_PID = Ad_PID^rVec(ii_PID);
      Gvec_PID(ii_PID) = Cd_PID*ApowVec_PID*Bd_PID;
    end
    G_PID = tril(toeplitz(Gvec_PID)); %utworzenie macierzy toeplitza i dolnotrójk¹tnej 
                                      %transmitancji regulatora GPID

    jmax = 200;
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,500*(2/Fs),'low'); 
        %Yj - wyjœcie z uk³adu regulacji
        %Uj - sygna³ wyjœciowy z ILC
        %Ej - sygna³ b³êdu
        %Rj - wartoœæ zadana na wymuszeniu uk³adu
        %G - transmitancja obiektu
        %Q - macierz, filtr
        %L - macierz ucz¹ca PD
    % Run ILC and plot the response for each iteration
    for ii = 1:jmax
      Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn¹ próbkê czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
        Ej_FPrzed(ii,1)=plotter_Ej(ii,t,Ej);
    end
    
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,400*(2/Fs),'low'); 

    for ii = 1:jmax
      Uj_1=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      Uj=filtfilt(bB,aB,Uj_1);
      Yj = G*Uj+rj_step;
      Ej = Rj - Yj; Ej(1) = 0;
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn¹ próbkê czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
    Ej_FPo(ii,1)=plotter_Ej(ii,t,Ej);
      
    end
%%     
    figure(921)
    iteration = 1:1:length(Ej_FPo);
    plot(iteration,Ej_FPrzed,'r-');
    hold on
    plot(iteration,Ej_FPo,'b-');
    legend('Filtr przed ILC', 'Filtr po ILC');
    title('OdpowiedŸ uk³adu na wymuszenie');%B³¹d RMSE')
    xlabel('Iteracja');
    ylabel('Przemieszczenie wirnika [m]');%B³¹d RMSE
    set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
    set(gca, 'YScale', 'log'); 
    set(findall(gcf,'type','line'),'linewidth',2)
    ylim([10e-9 3e-4]);
    %% dodatkowo bez filtra
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,400*(2/Fs),'low'); 

    for ii = 1:jmax
      Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn¹ próbkê czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
        Ej_FPrzed(ii,1)=plotter_Ej(ii,t,Ej);
    end
    
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,500*(2/Fs),'low'); 

    for ii = 1:jmax
      Uj_1=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      Uj=filtfilt(bB,aB,Uj_1);
      Yj = G*Uj+rj_step;
      Ej = Rj - Yj; Ej(1) = 0;
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn¹ próbkê czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
    Ej_FPo(ii,1)=plotter_Ej(ii,t,Ej);
      
    end
    

%%     

    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,500*(2/Fs),'low'); 
    
    for ii = 1:jmax
      Uj=(Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold));
     
      Yj = G*Uj+rj_step;
      Ej = Rj - Yj; Ej(1) = 0;
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
      Ej_BF(ii,1)=plotter_Ej(ii,t,Ej);
      Ej_F3_max(ii,1)=plotter_Ej_max(ii,Ej); 
    end
%     
%        Uj = zeros(N,1); Ujold = Uj;
%     Ej = zeros(N,1); Ejold = Ej;
%     
%     Ejold_add_1=zeros(N,1);
%     kp1=5;
%     kd1=5;
%     Fs=10000;
%     [bB,aB] = butter(4,1500*(2/Fs),'low'); 
% 
%     for ii = 1:jmax
%       Uj_1=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
%       Uj=filtfilt(bB,aB,Uj_1);
%       Yj = G*Uj+rj_step;
%       Ej_1 = Rj - Yj; Ej(1) = 0;
%       Ej=filtfilt(bB,aB,Ej_1);
%       Ejold = Ej;
%       Ujold = Uj;
%         for i=1:1:N-1 %przesuwa wektor Ej o jedn¹ próbkê czasu Ej(q+1)
%             Ejold_add_1(i+1,1)=Ejold(i,1);
%         end
%     Ej_FPoiprzed(ii,1)=plotter_Ej(ii,t,Ej);
%       
%     end
%%     
    figure(922)
    iteration = 1:1:length(Ej_BF);
%     plot(iteration,Ej_FPrzed,'r-');
    hold on
%     plot(iteration,Ej_FPo,'b');  
%      plot(iteration,Ej_F3_max,'b');
    plot(iteration,Ej_BF,'r');
    
%     plot(iteration,Ej_FPoiprzed,'c+')
    legend('ILC bez filtra');%'Filtr przed ILC', 'Filtr po ILC',
    xlabel('Iteracje');
    ylabel('B³¹d RMSE [m]');
    set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
    set(gca, 'YScale', 'log');
    set(findall(gcf,'type','line'),'linewidth',2)
    grid on