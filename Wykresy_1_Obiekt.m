clc, clear all, close all
%DEFINIOWANIE OBIEKTU
%transmitancja obiektu AMB
% num = [40.2]; 
% den = [3 0 -40.2436];
% num = [435440]; 
% den = [3 0 -40.2436];
num = [40.2436]; 
den = [3 0 -435440];
Ts = 0.00001;
sysc = tf(num,den); %postaæ ci¹g³a
sysd = c2d(sysc,Ts,'ZOH'); %postaæ dyskretna

%% mapa biegunów i zer
figure()
% pole(sysd)
pzmap(sysd);
% grid on
%% charakterystyki czasowe
figure(3)
step(sysd);
title('OdpowiedŸ skokowa uk³adu na wymuszenie');
xlabel('Czas [s]');
ylabel('Przemieszczenie wirnika [mm]');
labels = findall(figure(3),'Type','Text');
label = labels(strncmp(get(labels,'String'),'Czas [s]',8));
set(label,'String',regexprep(get(label,'String'),'(\(\w*)\)',''));
set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'-property','FontSize'),'FontSize',18);
set(gcf,'Color',([1 1 1]));

%% charakterystyki czêstotliwoœciowe
figure(4)
str=nyquistoptions;
str.ShowFullContour='off';
str.FreqUnits='Hz';
nyquist(sysd,str);

figure(5)
str1=bodeoptions;
str1.FreqUnits='Hz';
bode(sysd,str1);
set(findall(gcf,'-property','FontSize'),'FontSize',14);
title('Charakterystyka Bodego','FontWeight', 'bold','FontSize',18);
xlabel('Czêstotliwoœæ','FontSize', 18);
axes=findobj('type','axes');
h_magnitude=get(axes(2),'YLabel');
h_phase=get(axes(1),'YLabel');
set(h_magnitude,'String','Modu³ (dB)');
set(h_magnitude,'FontSize', 18);
set(h_phase,'String','Faza (deg)');
set(h_phase,'FontSize', 18);
set(gcf,'Color',([1 1 1]));
set(findall(gcf,'type','line'),'linewidth',2);
% xlabel('FontSize', 18);
 xlim([10e-3 50000]);

grid on
% (Do Charakterystyki Bodego) Uk³ad bêdzie stabilny, gdy:
% - logarytmiczna charakterystyka fazowa osi¹gnie –?, to logarytmiczna charakterystyka amplitudowa bêdzie ujemna,
% - logarytmiczna charakterystyka amplitudowa jest równa 0, to logarytmiczna charakterystyka fazowa le¿y powy¿ej –?.

%%
%parametry regulatora PID
kp=5; 
ki=10; 
kd=0.7; 
% Tf=0.006437;
% Ts=0.001;
%sysd_PID=pid(kp,ki,kd);
sysd_PID=pidtune(sysd,'pidf');
sysd_PID.Ts=0.00001;

sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zastêpczej uk³adu z regulatorem PID
sysd_f=feedback(sysd_se,1);

figure(6)
% step(sysd_f);
opt=stepDataOptions('InputOffset',-0.0002,'StepAmplitude',0.0002);
step(sysd_f,opt);
title('OdpowiedŸ skokowa uk³adu na wymuszenie');
xlabel('Czas [s]');
ylabel('Przemieszczenie wirnika [mm]');
set(gcf,'Color',([1 1 1]));
set(findall(gcf,'-property','FontSize'),'FontSize',18);
labels = findall(figure(6),'Type','Text');
label = labels(strncmp(get(labels,'String'),'Czas [s]',8));
set(label,'String',regexprep(get(label,'String'),'(\(\w*)\)',''));
set(findall(gcf,'type','line'),'linewidth',2)
ylim([-0.0002 0.00005]);
%% charakterystyki czêstotliwoœciowe
figure(7)
str=nyquistoptions;
str.ShowFullContour='off';
str.FreqUnits='Hz';
nyquist(sysd_f,str);
title('Charakterystyka Nyquista');
xlabel('Oœ rzeczywista');
ylabel('Oœ urojona');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
%% 
figure(8)
str1=bodeoptions;
str1.FreqUnits='Hz';
bode(sysd_f,str1);
set(findall(gcf,'-property','FontSize'),'FontSize',14);
title('Charakterystyka Bodego','FontWeight', 'bold','FontSize',18');
xlabel('Czêstotliwoœæ','FontSize', 18);
axes=findobj('type','axes');
h_magnitude=get(axes(2),'YLabel');
h_phase=get(axes(1),'YLabel');
set(h_magnitude,'String','Zapas modu³u (dB)');
set(h_magnitude,'FontSize', 18);
set(h_phase,'String','Zapas fazy (deg)');
set(h_phase,'FontSize', 18);
set(gcf,'Color',([1 1 1]));

set(findall(gcf,'type','line'),'linewidth',2);
xlim([10e-3 50000]);

grid on
%% POLE PLACEMENT
[A,B,C,D]=ssdata(ss(sysd(1,1)));
p=[-1+1i -1-1i];
% p=[-134 -2.4];
K=place(A,B,p);
sys_K=ss(A-B(:,1)*K,B,C,D);
figure()
pzmap(sys_K);
% step(sys_K);
KDC=dcgain(sys_K)
Kr=1/KDC
sys_Kskala=ss(A-B(:,1)*K,B*Kr,C,D);
opt=stepDataOptions('InputOffset',-0.0002,'StepAmplitude',0.0002);
% step(sys_Kskala,opt);
%%
figure()
margin(sysd_f);
grid on
%%
figure()
% pole(sysd)
pzmap(sysd_f);
% (Do Charakterystyki Bodego) Uk³ad bêdzie stabilny, gdy:
% - logarytmiczna charakterystyka fazowa osi¹gnie –?, to logarytmiczna charakterystyka amplitudowa bêdzie ujemna,
% - logarytmiczna charakterystyka amplitudowa jest równa 0, to logarytmiczna charakterystyka fazowa le¿y powy¿ej –?.