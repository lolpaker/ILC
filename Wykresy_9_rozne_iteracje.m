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
% kp=50; 
% ki=100; 
% kd=40; 

%sysd_PID=pid(kp,ki,kd);
 sysd_PID=pidtune(sysd,'pidf');
sysd_PID.Ts=0.00001;

sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zastêpczej uk³adu z regulatorem PID
sysd_f=feedback(sysd_se,1); 

% figure(1)
% pzmap(sysd_f)
% set(gcf,'Color',([1 1 1]));
% ylabel('Oœ liczb urojonych');
% xlabel('Oœ liczb rzeczywistych');
% set(findall(gcf,'-property','FontSize'),'FontSize',14)
% figure(2)
% pzmap(sysd)
% set(gcf,'Color',([1 1 1]));
% ylabel('Oœ liczb urojonych');
% xlabel('Oœ liczb rzeczywistych');
% set(findall(gcf,'-property','FontSize'),'FontSize',14)

%%
    [Ad Bd Cd Dd] = ssdata(sysd_f); %zamiana transmitancji na model przestrzeni stanu
   
     
% Warunki poczatkowe x0, zakres czasu t, opoznienie n0, r, macierz N
    
    
    %t = 0:Ts:0.01;
    t = 0:Ts:0.05;
    n0 = 0;
    r = 1;
    N = length(t);
    % Definiowanie wejsciowego wektora U i sygna³u wartosci zadanej Rj
%     rj_step=-0.0002;
    rj_step=0;
   
%    Rj = [rj_step*ones(1,200) 0*ones(1,4801)]';

    
   xx=t;
   % pseudo gauss
%    Rj = (784313.725*power(xx,6)- 117647.0588*power(xx,5)+ 5867.269981*power(xx,4)- 96.53091996*power(xx,3)- 0.062021803*power(xx,2)+0.001215721*power(xx,1)+ 6.06746E-07)';
   
   % sin
   Rj=(0.0001*sin(20*2*pi*t))';
   
%    figure()
%     plot(t,Rj)
%     grid on
    
    %%
    % Sformu³owanie macierzy obiektu
    Gvec = zeros(N,1);
    rVec = ((r-1):(N-n0-1))';
    
    for ii = 1:length(rVec)
      ApowVec = Ad^rVec(ii);
      Gvec(ii) = Cd*ApowVec*Bd;
    end
    G = tril(toeplitz(Gvec)); %utworzenie macierzy toeplitza i dolnotrójk¹tnej obiektu G
    
    jmax = 52;
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=1;
    kd1=10;
    Fs=100000;
    [bB,aB] = butter(4,400*(2/Fs),'low'); 
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
%       Yj = G*Uj;
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn¹ próbkê czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
            %  plotter(ii,t,Ej,Yj,Uj,Rj)
                 ri=1;
      if ii> ri 
          for xi=1:1:N
              Ej_xi(xi,ii-ri)=Ej(xi,1);
              Uj_xi(xi,ii-ri)=Uj(xi,1);
              Yj_xi(xi,ii-ri)=Yj(xi,1);
          end
      end     
              
%        plotter_Ej(ii,t,Ej)
%      plotter_errors(ii,Ej);
    end
    
    % ustawienie iteracji
    it1=1;
    it2=5;
    it3=10;
    it4=50;
    
    %% b³¹d Ej
    figure(764)
    plot(t,Ej_xi(:,it1),'b', t,Ej_xi(:,it2),'k', t,Ej_xi(:,it3),'r', t,Ej_xi(:,it4),'g','LineWidth',2) 
    grid on
   legend(strcat('iteracja =',num2str(it1)),strcat('iteracja = ',num2str(it2)),strcat('iteracja = ',num2str(it3)),strcat('iteracja = ',num2str(it4)));
     set(gcf,'Color',([1 1 1]));
     xlabel('Czas [s]');
     ylabel('Przemieszczenie wirnika [m]');
     set(findall(gcf,'-property','FontSize'),'FontSize',18);
     xlim([0 0.015])
     %sterowanie Uj
       figure(765)
    plot(t,Uj_xi(:,it1),'b', t,Uj_xi(:,it2),'k', t,Uj_xi(:,it3),'r', t,Uj_xi(:,it4),'g','LineWidth',2) 
    grid on
   legend(strcat('iteracja =',num2str(it1)),strcat('iteracja = ',num2str(it2)),strcat('iteracja = ',num2str(it3)),strcat('iteracja = ',num2str(it4)));
     set(gcf,'Color',([1 1 1]));
     xlabel('Czas [s]');
     ylabel('Natê¿enie [A]');
     set(findall(gcf,'-property','FontSize'),'FontSize',18); 
     xlim([0 0.015])
     % Wyjœcie Yj
        figure(766)
    plot(t,Yj_xi(:,it1),'b', t,Yj_xi(:,it2),'k', t,Yj_xi(:,it3),'r', t,Yj_xi(:,it4),'g',t,Rj,'k--','LineWidth',2) 
    grid on
   legend(strcat('iteracja =',num2str(it1)),strcat('iteracja = ',num2str(it2)),strcat('iteracja = ',num2str(it3)),strcat('iteracja = ',num2str(it4)),['wymuszenie']);
     set(gcf,'Color',([1 1 1]));
     xlabel('Czas [s]');
     ylabel('Przemieszczenie wirnika [m]');
     set(findall(gcf,'-property','FontSize'),'FontSize',18);
%     xlim([0 0.015]);
    ylim([-0.000125 0.000125])