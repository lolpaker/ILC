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

    jmax = 50;
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,400*(2/Fs),'low');  
%     sysfiltr=tf(bB,aB);
%     sysfiltrd=c2d(sysfiltr,Ts)
    for ii = 1:jmax
     Uj=Ujold+kp1*Ejold+kd1*(Ejold_add_1-Ejold);
      
      Yj = G*Uj+rj_step;
      Ej = Rj - Yj; Ej(1) = 0;
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
    end
     

    % charakterystyki dla sygna³u Uj
    figure(8) %charakterystyka amplitudowa zak³óconego
Y2 = fft(Ejold); %szybka transformata fouriera
q22=length(Ejold);
P22 = abs(Y2/q22);  %elementy zaczerpniête z pomocy programu Matlab
P12 = P22(1:q22/2+1);
P12(2:end-1) = 2*P12(2:end-1);
FX2 = Fs*(0:(q22/2))/q22;
plot(FX2,P12)
title('Charakterystyka amplitudowa sygna³u b³êdu')
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Amplituda');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(gcf,'Color',([1 1 1]));
set(findall(gcf,'type','line'),'linewidth',2);
xlim([0 5000]);
ylim([0 0.0004]);
grid on
hold on
%% 
figure(71)%charakterystyka fazowa zak³óconego
theta22 = unwrap(angle(Y2/q22));  
theta11 = theta22(1:q22/2+1);
theta11(2:end-1) = theta11(2:end-1);
plot(FX2,theta11)
title('Charakterystyka fazowa sygna³u Uj')
xlabel('Czêstotliwoœæ [Hz]');  
ylabel('Faza [rad]');
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(gcf,'Color',([1 1 1]));

  
%% filtr dolnoprzepustowy butterworda
% [bB,aB] = butter(4,2000*(2/Fs),'low'); % 1 rz¹d,górna czêstotliwoœæ filtrowania 0.1, dolna czêstotliwoœæ filtrowania 60
[HB,wB] = freqz(bB,aB,N); %zwraca odpowiedŸ czêstotliwoœciow¹ dla systemu filtracji z zastosowaniem 1001 próbek

figure(113)%charakterystyka amplitudowa(liniowa,decybelowa) i fazowa filtru
subplot(3,1,1)
plot(wB*Fs/(2*pi),abs(HB),'r','LineWidth',2)
title('Charakterystyka amplitudowa filtru Butterwortha dolnoprzepustowego (skala liniowa)')
xlabel('Czêstotliwoœæ [Hz]');
ylabel('Amplituda');
xlim([0 5000]);
grid on
subplot(3,1,2)
plot(wB*Fs/(2*pi),mag2db(abs(HB)),'r','LineWidth',2)
title('Charakterystyka amplitudowa filtru Butterwortha  dolnoprzepustowego (skala decybelowa)')
xlabel('Czêstotliwoœæ [Hz]');
ylabel('Amplituda [dB]');
xlim([0 5000]);
ylim([-100 20]);
grid on
subplot(3,1,3)
plot(wB*Fs/(2*pi),unwrap(angle(HB)),'b','LineWidth',2)
title('Charakterystyka fazowa filtru Butterwortha dolnoprzepustowego');
xlabel('Czêstotliwoœæ [Hz]');
ylabel('Faza [rad]');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(gcf,'Color',([1 1 1]));
xlim([0 5000]);
grid on
%%   
figure(13)
BUT_X2 = filtfilt(bB,aB,Ejold);
Y_BUT_X2 = fft(BUT_X2); %szybka transformata fouriera
q2_BUT_X2=length(BUT_X2);
P2_BUT_X2 = abs(Y_BUT_X2/q2_BUT_X2);  %elementy zaczerpniête z pomocy programu Matlab
P1_BUT_X2 = P2_BUT_X2(1:q2_BUT_X2/2+1);
P1_BUT_X2(2:end-1) = 2*P1_BUT_X2(2:end-1);
FX_BUT_X2 = Fs*(0:(q2_BUT_X2/2))/q2_BUT_X2;
subplot(2,1,1)% charakterystyka amplitudowa sygna³u liniowa
plot(FX2,P12,'b','LineWidth',1)
hold on
plot(FX_BUT_X2,P1_BUT_X2,'r','LineWidth',1)
title('Charakterystyka amplitudowa b³êdu z filtrem Butterwortha (skala liniowa)')
xlabel('Czêstotliwoœæ [Hz]');  
ylabel('Amplituda');
legend({'Sygna³ b³êdu przed filtracj¹','Sygna³ b³êdu po filtracji'});
set(findall(gcf,'-property','FontSize'),'FontSize',18);
xlim([0 5000]);
grid on

subplot(2,1,2)% charakterystyka amplitudowa sygna³u decybelowa
plot(FX2,mag2db(abs(P12)),'b','LineWidth',1)
hold on
plot(FX_BUT_X2,mag2db(abs(P1_BUT_X2)),'r','LineWidth',1)
title('Charakterystyka amplitudowa b³êdu z filtrem Butterwortha (skala decybelowa)')
xlabel('Czêstotliwoœæ [Hz]');  
ylabel('Amplituda [dB]');
legend({'Sygna³ b³êdu przed filtracj¹','Sygna³ b³êdu po filtracji'});
set(findall(gcf,'-property','FontSize'),'FontSize',18);
grid on
set(gcf,'Color',([1 1 1]));
xlim([0 5000]);
figure(73)%charakterystyka fazowa zak³óconego
thetaB2 = unwrap(angle(Y_BUT_X2/q2_BUT_X2));  
thetaB1 = thetaB2(1:q2_BUT_X2/2+1);
thetaB1(2:end-1) = thetaB1(2:end-1);
plot(FX2,theta11,'b','LineWidth',1)
hold on
plot(FX_BUT_X2,thetaB1,'r','LineWidth',1)
title('Charakterystyka fazowa sygna³u zak³óconego')
xlabel('Czêstotliwoœæ [Hz]');  
ylabel('Faza [rad]');
legend({'Sygna³ Uj przed filtracj¹','Sygna³ Uj po filtracji'});
set(findall(gcf,'-property','FontSize'),'FontSize',18);
grid on
set(gcf,'Color',([1 1 1]));
xlim([0 5000]);
