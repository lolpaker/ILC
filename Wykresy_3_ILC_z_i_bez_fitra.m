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
% kp=1; 
% ki=1; 
% kd=1; 

% sysd_PID=pid(kp,ki,kd);
sysd_PID=pidtune(sysd,'pidf');
sysd_PID.Ts=0.00001;
sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zastêpczej uk³adu z regulatorem PID
sysd_f=feedback(sysd_se,1); 
figure(2)
opt=stepDataOptions('InputOffset',-0.0002,'StepAmplitude',0.0002);
step(sysd_f,opt);
% [A,B,C,D]=ssdata(ss(sysd(1,1)));
% p=[-134 -200];
% K=place(A,B,p);
% sys_K=ss(A-B(:,1)*K,B,C,D);
% figure()
% pzmap(sys_K)
% step(sys_K);
% KDC=dcgain(sys_K);
% Kr=1/KDC;
% sys_Kskala=ss(A-B(:,1)*K,B*Kr,C,D);
% sysd_f=feedback(sys_Kskala,1);

% figure();initial(sysd_f,[-0.0002,0]);
% figure();initial(sys_K,[-0.0002,0]);

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



   %% ILC z filtrem 
       jmax = 50;
       
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    
    for ii = 1:jmax
    czesto=500;    
    [bB,aB] = butter(4,czesto*(2/Fs),'low'); 
    transfiltra=tf(bB,aB,Ts)
    
      Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
     plotter(ii,t,Ej,Yj,Uj,Rj)
%       Ej_ZF(ii,1)=plotter_Ej(ii,t,Ej);
    end
    
    %% ILC bez filtra

    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    
    for ii = 1:jmax
      Uj=(Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold));
     
      Yj = G*Uj+rj_step;
      Ej = Rj - Yj; Ej(1) = 0;
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
%         figure()
%         plot(t,Rj);
%       plotter(ii,t,Ej,Yj,Uj,Rj)
      Ej_BF(ii,1)=plotter_Ej(ii,t,Ej);
    end
    %% wykres Ej z filtrem i bez filtra 
    
    iteration = 1:1:length(Ej_ZF);
    
    figure(190)
    plot(iteration,Ej_BF,'r-');
%     hold on
%     plot(iteration,Ej_ZF,'b-');
     legend('ILC bez filtra')%, 'ILC z filtrem');
    title('OdpowiedŸ uk³adu na wymuszenie pr¹dowe')
    set(gca, 'YScale', 'log');
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
    xlabel('Iteracja');
    ylabel('Przemieszczenie wirnika [m]');
    set(gcf,'Color',([1 1 1]));

       ylim([10e-10 10e30]);
%      xlim([1 jmax]);
    set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'type','line'),'linewidth',2)
 
 
    
    
    
    