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
kp=50; 
ki=100; 
kd=40; 

sysd_PID=pid(kp,ki,kd);
%  sysd_PID=pidtune(sysd,'pidf');
sysd_PID.Ts=0.00001;

 
%%
    [Ad Bd Cd Dd] = ssdata(sysd); %zamiana transmitancji na model przestrzeni stanu
    [Ad_PID Bd_PID Cd_PID Dd_PID] = dssdata(sysd_PID); %zamiana transmitancji na model przestrzeni stanu regulatora PID
     
% Warunki poczatkowe x0, zakres czasu t, opoznienie n0, r, macierz N
    
    
    %t = 0:Ts:0.01;
    t = 0:Ts:0.05;
    n0 = 0;
    r = 1;
    N = length(t);
    % Definiowanie wejsciowego wektora U i sygna³u wartosci zadanej Rj
    rj_step=-0.0002;

    %Rj = [rj_step*ones(1,5) 0*ones(1,96)]';
    Rj = [rj_step*ones(1,200) 0*ones(1,4801)]';

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

    jmax = 20;
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,500*(2/Fs),'low'); 

    for ii = 1:jmax
       % ILC 
      Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      % wyjœcie z C
      Uc=G_PID*Ej;
      % wyjœcie z sumatora 2
      Ujx=Uc+Uj;
      % wyjœcie z uk³adu  
      Yj = G*Ujx+rj_step;
      % wyjœcie z sumatora 1 (yd-yj)
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      
      Ej=filtfilt(bB,aB,Ej_1);
      
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn¹ próbkê czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
              plotter(ii,t,Ej,Yj,Uj,Rj)
%        plotter_Ej(ii,t,Ej)
%      plotter_errors(ii,Ej);
    end
    