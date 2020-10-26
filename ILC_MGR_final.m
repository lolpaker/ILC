clc, clear all, close all
%transmitancja obiektu AMB
% num = [85.0824]; 
% den = [3 0 -40.2436];
% num = [435440]; 
% den = [3 0 -40.2436];
num = [40.2436]; 
den = [3 0 -435440];
Ts = 0.00001;
sysc = tf(num,den); %posta� ci�g�a
sysd = c2d(sysc,Ts,'ZOH'); %posta� dyskretna

%parametry regulatora PID
% kp=50; 
% ki=100; 
% kd=40; 

% sysd_PID=pid(kp,ki,kd);
sysd_PID=pidtune(sysd,'pidf');
sysd_PID.Ts=0.00001;

sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zast�pczej uk�adu z regulatorem PID
sysd_f=feedback(sysd_se,1); 
% figure(2)
% step(sysd_f);
%%
    [Ad Bd Cd Dd] = ssdata(sysd_f); %zamiana transmitancji na model przestrzeni stanu
    [Ad_PID Bd_PID Cd_PID Dd_PID] = dssdata(sysd_PID); %zamiana transmitancji na model przestrzeni stanu regulatora PID
     
% Warunki poczatkowe x0, zakres czasu t, opoznienie n0, r, macierz N
    
    
    %t = 0:Ts:0.01;
    t = 0:Ts:0.05;
    n0 = 0;
    r = 1;
    N = length(t);
    % Definiowanie wejsciowego wektora U i sygna�u wartosci zadanej Rj
    rj_step=-0.0002;
    rj_step=0;
    %Rj = [rj_step*ones(1,5) 0*ones(1,96)]';
%     Rj = [rj_step*ones(1,200) 0*ones(1,4801)]';
    xx=t;
    Rj = (784313.725*power(xx,6)- 117647.0588*power(xx,5)+ 5867.269981*power(xx,4)- 96.53091996*power(xx,3)- 0.062021803*power(xx,2)+0.001215721*power(xx,1)+ 6.06746E-07)';
% sin
%    Rj=(0.0001*sin(20*2*pi*t))';    
%     grid on
    % Sformu�owanie macierzy obiektu
    Gvec = zeros(N,1);
    rVec = ((r-1):(N-n0-1))';
    
    for ii = 1:length(rVec)
      ApowVec = Ad^rVec(ii);
      Gvec(ii) = Cd*ApowVec*Bd;
    end
    G = tril(toeplitz(Gvec)); %utworzenie macierzy toeplitza i dolnotr�jk�tnej obiektu G
    
    for ii_PID = 1:length(rVec)
      ApowVec_PID = Ad_PID^rVec(ii_PID);
      Gvec_PID(ii_PID) = Cd_PID*ApowVec_PID*Bd_PID;
    end
    G_PID = tril(toeplitz(Gvec_PID)); %utworzenie macierzy toeplitza i dolnotr�jk�tnej 
                                      %transmitancji regulatora GPID

    jmax = 50;
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,400*(2/Fs),'low'); 
        %Yj - wyj�cie z uk�adu regulacji
        %Uj - sygna� wyj�ciowy z ILC
        %Ej - sygna� b��du
        %Rj - warto�� zadana na wymuszeniu uk�adu
        %G - transmitancja obiektu
        %Q - macierz, filtr
        %L - macierz ucz�ca PD
    % Run ILC and plot the response for each iteration
    for ii = 1:jmax
      Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      Yj = G*Uj+rj_step;
%       Yj = G*Uj;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn� pr�bk� czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
              plotter(ii,t,Ej,Yj,Uj,Rj)
%        plotter_Ej(ii,t,Ej)
%      plotter_errors(ii,Ej);
    end
    
%% wykres wartosci zadanej

    figure()
    plot(t,Rj)
%     legend('')%, 'ILC z filtrem');
%     title('Wymuszenie skokowe');
%     title('Wymuszenie nieliniowe');    ,
%     title('Wymuszenie sinusoidalne'); 
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
    xlabel('Czas[s]');
    ylabel('Przemieszczenie wirnika [m]');
    set(gcf,'Color',([1 1 1]));

%        ylim([-0.00015 0.00015]);
%      xlim([1 jmax]);
    set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'type','line'),'linewidth',2)
    grid on
    
