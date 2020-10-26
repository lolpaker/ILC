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
% kp=10; 
% ki=5; 
% kd=0.7; 

%  sysd_PID=pid(kp,ki,kd);
sysd_PID=pidtune(sysd,'pidf');%Tworzenie transmitancji zast?pczej uk?adu z regulatorem PID
sysd_PID.Ts=0.00001;

sysd_se=series(sysd,sysd_PID);
% sysd_f=feedback(sysd_se,1); 
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
    Rj = [rj_step*ones(1,200) 0*ones(1,4801)]';
        xx=t;
%     Rj = (784313.725*power(xx,6)- 117647.0588*power(xx,5)+ 5867.269981*power(xx,4)- 96.53091996*power(xx,3)- 0.062021803*power(xx,2)+0.001215721*power(xx,1)+ 6.06746E-07)';
    % sin
%    Rj=(0.0001*sin(20*2*pi*t))';    
    
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
%% 
                                      
    % cz?stotliwo?ci filtrów
    f1=400;
    f2=f1;
    f3=f1;
    f4=f1;                                  
    jmax = 50;
    
    kp1=1; %0.5
    kd1=10;
    
    kp2=0.3; %5
    kd2=0;
    
    kp3=0.6; %10
    kd3=0;
   
    kp4=0.9;
    kd4=0;
    
    %1
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);

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
        
%      Ej_F1(ii,1)=plotter_Ej(ii,t,Ej);
%     Ej_F1_max(ii,1)=plotter_Ej_max(ii,Ej);       
%         plotter_Ej(ii,t,Ej);
       plotter(ii,t,Ej,Yj,Uj,Rj);
%         plotter_errors(ii,Ej);
    end
    Yj_w1=Yj;
 %%
 %2   
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    Fs=100000;
    [bB,aB] = butter(4,f2*(2/Fs),'low'); 

    for ii = 1:jmax
     Uj=Ujold+kp2*Ejold+kd2*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);

      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn? próbk? czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
%     Ej_F2(ii,1)=plotter_Ej(ii,t,Ej);
%     Ej_F2_max(ii,1)=plotter_Ej_max(ii,Ej);        
%         plotter_Ej(ii,t,Ej);
%         plotter_errors(ii,Ej);
    end
    Yj_w2=Yj;
%3    
        Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    Fs=100000;
    [bB,aB] = butter(4,f3*(2/Fs),'low'); 

    for ii = 1:jmax
     Uj=Ujold+kp3*Ejold+kd3*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn? próbk? czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
%     Ej_F3(ii,1)=plotter_Ej(ii,t,Ej);
%     Ej_F3_max(ii,1)=plotter_Ej_max(ii,Ej);       
%         plotter_Ej(ii,t,Ej);
%         plotter_errors(ii,Ej);
    end
    Yj_w3=Yj;
  
%4
        Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    Fs=100000;
    [bB,aB] = butter(4,f4*(2/Fs),'low'); 

    for ii = 1:jmax
     Uj=Ujold+kp4*Ejold+kd4*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn? próbk? czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
%     Ej_F4(ii,1)=plotter_Ej(ii,t,Ej);
%     Ej_F4_max(ii,1)=plotter_Ej_max(ii,Ej);       
%         plotter_Ej(ii,t,Ej);
%         plotter_errors(ii,Ej);
    end
    Yj_w4=Yj;

%% wykres    
figure()
plot(t,Yj_w1,'b',t,Yj_w2,'k',t,Rj,'k--')%t,Yj_w3,'r',,t,Yj_w4,'g'
grid on
%     legend('kp=0.1, kd=1','kp=1, kd=10','kp=1, kd=30','kp=1.6, kd=10');
legend(strcat('kp=',num2str(kp1)),strcat('kp=',num2str(kp2)),strcat('kp=',num2str(kp3)));%,strcat('kp=',num2str(kp4))
% legend(strcat('kd=',num2str(kd1)),strcat('kd=',num2str(kd2)));%,strcat('kd=',num2str(kd3))
% title('OdpowiedŸ uk³adu na wymuszenie');
xlabel('Czas [s]');
ylabel('Przemieszczenie wirnika [m]');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(gcf,'Color',([1 1 1]));
set(findall(gcf,'type','line'),'linewidth',2);
hold on
grid on
%  xlim([0 0.03]);
 s1=stepinfo(Yj_w1,t)
 s2=stepinfo(Yj_w2,t)
%  s3=stepinfo(Yj_w3,t)
%  s4=stepinfo(Yj_w4,t)
 %% wykresy b³êdu
     figure(922)
    iteration = 1:1:length(Ej_F1_max);
%     plot(iteration,Ej_F1,'blue');
    hold on
    plot(iteration,Ej_F1_max,'blue--');
%     plot(iteration,Ej_F2,'black');  
    plot(iteration,Ej_F2_max,'black--');
%    plot(iteration,Ej_F3,'r');
%     plot(iteration,Ej_F3_max),'red';
%     plot(iteration,Ej_F4,'g');
%     plot(iteration,Ej_F4_max,'g');
%     plot(iteration,Ej_FPoiprzed,'c+')
%     title('B³¹d RMSE');
%     legend('kp=0.1, kd=1','kp=1, kd=10','kp=1, kd=30','kp=1.6, kd=10');
%        legend('kp=0.1','kp=0.3','kp=0.6','kp=0.9');
    legend('PD MAX','P MAX')%,'kp=0.9')'PD MAX',,'P MAX'
% legend('B³¹d RMSE', 'B³¹d MAX');
    xlabel('Iteracje');
    ylabel('B³¹d [m]');
    set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
    set(gca, 'YScale', 'log');
    set(findall(gcf,'type','line'),'linewidth',2);
    grid on
     ylim([0.0000001 0.0003])