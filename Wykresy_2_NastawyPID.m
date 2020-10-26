clc, clear all, close all
%DEFINIOWANIE OBIEKTU, REGULATORA I FILTRA
%transmitancja obiektu AMB
% num = [435440]; 
% den = [3 0 -40.2436];
% num = [85.0824]; 
% den = [3 0 -40.2436];
num = [40.2436]; 
den = [3 0 -435440];
Ts = 0.00001;
sysc = tf(num,den); %postaæ ci¹g³a
sysd = c2d(sysc,Ts,'ZOH'); %postaæ dyskretna


%parametry regulatora PID
% kp=0.002; 
% ki=0.004; 
% kd=0.09; 
% kp=5;
% ki=10;
% kd=50;
% sysd_PID=pid(kp,ki,kd);
sysd_PID=pidtune(sysd,'pidf')
sysd_PID.Ts=0.00001;

sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zastêpczej uk³adu z regulatorem PID
sysd_f=feedback(sysd_se,1); 
fig1=figure(1)
% opt=stepDataOptions('InputOffset',-0.0002,'StepAmplitude',0.0002);
step(sysd_f);
title('OdpowiedŸ skokowa uk³adu na wymuszenie');
xlabel('Czas [s]');
ylabel('Amplituda');
% ylabel('Przemieszczenie wirnika [m]');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
labels = findall(fig1,'Type','Text');
label = labels(strncmp(get(labels,'String'),'Czas [s]',8));
set(label,'String',regexprep(get(label,'String'),'(\(\w*)\)',''));
set(findall(gcf,'type','line'),'linewidth',2);
%  ylim([-0.0002 0.00005]);
set(gcf,'Color',([1 1 1]));
% grid on;
%% zmiana kp
fig2=figure(2)
for p=200000:5000:205000
% Kp = 0.00221, Ki = 0.00425, Kd = 0.000282, Tf = 0.00445, Ts = 0.0001
%     ki=10.11; 
%     kd=0.7828; 
%     sysd_PID=pid(kp,ki,kd);
    sysd_PID=pidtune(sysd,'pidf');
    sysd_PID.kp=p; 
    sysd_PID.Ts=0.00001;

    sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zastêpczej uk³adu z regulatorem PID
    sysd_f=feedback(sysd_se,1);
    step(sysd_f,opt);
    hold on;
    grid on;
end
title('OdpowiedŸ skokowa uk³adu na wymuszenie');
legend('kp=0.002','kp=0.042','kp=0.082','kp=0.122');
set(gcf,'Color',([1 1 1]));
xlabel('Czas [s]');
ylabel('Przemieszczenie wirnika [m]');
set(findall(gcf,'-property','FontSize'),'FontSize',18);

labels = findall(fig2,'Type','Text');
label = labels(strncmp(get(labels,'String'),'Czas [s]',8));
set(label,'String',regexprep(get(label,'String'),'(\(\w*)\)',''));
set(findall(gcf,'type','line'),'linewidth',2);
 ylim([-0.0002 0.0002]);
%% zmiana ki
fig3=figure(3)
for i=0.0002:0.01:0.04
%     kp=5.626; 
%     ki=i; 
%     kd=0.7828; 
%     sysd_PID=pid(kp,ki,kd);
%     sysd_PID.Ts=0.0001;
    sysd_PID=pidtune(sysd,'pidf');
    sysd_PID.ki=i; 
    sysd_PID.Ts=0.00001;
    
    sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zastêpczej uk³adu z regulatorem PID
    sysd_f=feedback(sysd_se,1);
    step(sysd_f,opt);
    hold on;
    grid on;
end
title('OdpowiedŸ skokowa uk³adu na wymuszenie')
legend('ki=0.002','ki=0.022','ki=0.042','ki=0.062');
set(gcf,'Color',([1 1 1]));
xlabel('Czas [s]');
ylabel('Przemieszczenie wirnika [m]');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
labels = findall(fig3,'Type','Text');
label = labels(strncmp(get(labels,'String'),'Czas [s]',8));
set(label,'String',regexprep(get(label,'String'),'(\(\w*)\)',''));
set(findall(gcf,'type','line'),'linewidth',2);
% ylim([-0.0002 0.0001]);
%% zmiana kd
fig4=figure(4);
for d=0.0002:0.0004:0.0014
%     kp=5.626; 
%     ki=10.11; 
%     kd=d; 
%     sysd_PID=pid(kp,ki,kd);
%     sysd_PID.Ts=0.0001;
    sysd_PID=pidtune(sysd,'pidf');
    sysd_PID.kd=d; 
    sysd_PID.Ts=0.00001;

    sysd_se=series(sysd,sysd_PID); %Tworzenie transmitancji zastêpczej uk³adu z regulatorem PID
    sysd_f=feedback(sysd_se,1);
    step(sysd_f,opt);
    hold on;
    grid on;
end
title('OdpowiedŸ skokowa uk³adu na wymuszenie')
xlabel('Czas [s]');
ylabel('Przemieszczenie wirnika [m]');
legend('kd=0.0002','kd=0.0006','kd=0.0010','kd=0.0014');
set(gcf,'Color',([1 1 1]));
set(findall(gcf,'-property','FontSize'),'FontSize',18);
labels = findall(fig4,'Type','Text');
label = labels(strncmp(get(labels,'String'),'Czas [s]',8));
set(label,'String',regexprep(get(label,'String'),'(\(\w*)\)',''));
set(findall(gcf,'type','line'),'linewidth',2);
% ylim([-0.0002 0.0001]);
% xlim([0 0.6]);