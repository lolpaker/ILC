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
% kp=5;%5.626; 
% ki=10;%10.11; 
% kd=0.7; 

% sysd_PID=pid(kp,ki,kd);
sysd_PID=pidtune(sysd,'pidf');
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
%     rj_step=-0.0002;
    rj_step=0;
%     Rj = [rj_step*ones(1,200) 0*ones(1,4801)]'; 
    xx=t;
%         Rj = (784313.725*power(xx,6)- 117647.0588*power(xx,5)+ 5867.269981*power(xx,4)- 96.53091996*power(xx,3)- 0.062021803*power(xx,2)+0.001215721*power(xx,1)+ 6.06746E-07)';
    %sin
       Rj=(0.0002*sin(20*2*pi*t))';    
    
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

    jmax = 30;
    
    Uj = zeros(N,1); Ujold = Uj;
    Ej = zeros(N,1); Ejold = Ej;
    
    Ejold_add_1=zeros(N,1);
    kp1=0.3;
    kd1=0;
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
      Ej_1 = Rj - Yj; Ej_1(1) = 0;
      Ej=filtfilt(bB,aB,Ej_1);
      Ejold = Ej;
      Ujold = Uj;
        for i=1:1:N-1 %przesuwa wektor Ej o jedn¹ próbkê czasu Ej(q+1)
            Ejold_add_1(i+1,1)=Ejold(i,1);
        end
        
        ri=1;
      if ii> ri 
          for xi=1:1:N
              Ej_xi(xi,ii-ri)=Ej(xi,1);
              Uj_xi(xi,ii-ri)=Uj(xi,1);
              Yj_xi(xi,ii-ri)=Yj(xi,1);
          end
      end
      
%        plotter(ii,t,Ej,Yj,Uj,Rj)
%      plotter_Ej(ii,t,Ej);
    end
    

     %%
figure(457)
x=linspace(ri+1,jmax-ri,jmax-ri); % wyœwietlanie wykresu od 2 iteracji, mo¿na napisaæ, 
% ¿e poczatkowe iteracje s¹ zerowe i Ÿle wp³ywaj¹ na czytelnoœæ wykresu 
y=t;
Z=Ej_xi;

[X,Y]=meshgrid(x,y);

%  surf(X,Y,Z,'FaceAlpha',1,'EdgeAlpha',0.1);
 surf(X,Y,Z,'FaceAlpha',1,'EdgeAlpha',0.1,'MeshStyle','column');

 title('Sygna³ b³êdu')
 set(gca, 'YDir','reverse')
 ylabel('Czas [s]');
 zlabel('Przemieszczenie [m]');
 xlabel('Iteracje');
  xlim([0 jmax]);
    ylim([0 0.05]);
%     zlim([-0.000001 0.0000045]);

 colormap jet;%
set(gcf,'Color',([1 1 1]));
  set(findall(gcf,'-property','FontSize'),'FontSize',18)
     %%
figure(458)

x=linspace(ri+1,jmax-ri,jmax-ri);
y=t;
Z=Uj_xi;

[X,Y]=meshgrid(x,y);

%  surf(X,Y,Z,'FaceAlpha',1,'EdgeAlpha',0.1);
 surf(X,Y,Z,'FaceAlpha',1,'EdgeAlpha',0.1,'MeshStyle','column');

 title('Sygna³ steruj¹cy')
 set(gca, 'YDir','reverse')
 ylabel('Czas [s]');
 zlabel('Natê¿enie [A]');
 xlabel('Iteracje');
 xlim([0 jmax]);
%  zlim([-0.0001 0.00025]);
  ylim([0 0.05]);
colormap(jet); % 
  set(gcf,'Color',([1 1 1]));
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
     %%
figure(459)
x=linspace(ri+1,jmax-ri,jmax-ri);
y=t;
Z=Yj_xi;

[X,Y]=meshgrid(x,y);

%  surf(X,Y,Z,'FaceAlpha',1,'EdgeAlpha',0.1);
 surf(X,Y,Z,'FaceAlpha',1,'EdgeAlpha',0.1,'MeshStyle','column');

 title('Sygna³ wyjœciowy')
 set(gca, 'YDir','reverse')
 ylabel('Czas [s]');
 zlabel('Przemieszczenie [m]');
 xlabel('Iteracje');
  xlim([0 jmax]);
%   zlim([-0.00021 0.0001]);
   ylim([0 0.05])
  colormap(jet); %
  set(gcf,'Color',([1 1 1]));
  set(findall(gcf,'-property','FontSize'),'FontSize',18)
%%  
% figure(100)
% subplot(1,3,1)
% [X,Y]=meshgrid(x,y);
% 
%  surf(X,Y,Z,'FaceAlpha',1,'EdgeAlpha',0.2);
% subplot(1,3,2)
% 
% subplot(1,3,3)
      