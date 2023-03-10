%% General
%clear all; clc; close all
% Load path
path = 'C:\Users\JeppeLaptop\OneDrive - Danmarks Tekniske Universitet\6. semester\Diplomingeniørprojekt\Målinger\Bode 100\';
%% PI Open, LQR closed 4
% Load
log_name = 'xtf_0dbm_pi_open_lqr_closed_4ohm';
pi_open_lqr_closed_4ohm_table = readtable([path log_name '.log']);

pi_open_lqr_closed_4ohm_data = table2array(pi_open_lqr_closed_4ohm_table);
figure(1)
subplot(2,1,1)
semilogx(pi_open_lqr_closed_4ohm_data(:,1), pi_open_lqr_closed_4ohm_data(:,2), 'k');
hold on
grid on
xlim([min(pi_open_lqr_closed_4ohm_data(:,1)) max(pi_open_lqr_closed_4ohm_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold off
subplot(2,1,2)
semilogx(pi_open_lqr_closed_4ohm_data(:,1), pi_open_lqr_closed_4ohm_data(:,3),'--','color','black');
hold on
grid on
xlim([min(pi_open_lqr_closed_4ohm_data(:,1)) max(pi_open_lqr_closed_4ohm_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
hold off
%sgt = sgtitle('PI controller frequency response','Color','black');
%sgt.FontSize = 20;

%% PI Open, LQR closed
log_name = 'tf_0dbm_pi_open_lqr_closed';
pi_open_lqr_closed_table = readtable([path log_name '.log']);

pi_open_lqr_closed_data = table2array(pi_open_lqr_closed_table);
figure(2)
subplot(2,1,1)
semilogx(pi_open_lqr_closed_data(:,1), pi_open_lqr_closed_data(:,2), 'k');
hold on
grid on
xlim([min(pi_open_lqr_closed_data(:,1)) max(pi_open_lqr_closed_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold off
subplot(2,1,2)
semilogx(pi_open_lqr_closed_data(:,1), pi_open_lqr_closed_data(:,3),'','color','black');
hold on
grid on
xlim([min(pi_open_lqr_closed_data(:,1)) max(pi_open_lqr_closed_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
hold off

%% PI open, LQR closed 300mohm ESR
log_name = 'tf_0dbm_pi_open_lqr_closed_lout_300mohmESR_both';
pi_open_lqr_closed_L_ESR_table = readtable([path log_name '.log']);
pi_open_lqr_closed_L_ESR_data = table2array(pi_open_lqr_closed_L_ESR_table);

figure(2)

subplot(2,1,1)
semilogx(pi_open_lqr_closed_data(:,1), pi_open_lqr_closed_data(:,2),'-','color','black');
hold on
semilogx(pi_open_lqr_closed_L_ESR_data(:,1), pi_open_lqr_closed_L_ESR_data(:,2), '--','color','black');
grid on
%xlim([min(pi_open_lqr_closed_L_ESR_data(:,1)) max(pi_open_lqr_closed_L_ESR_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend('PI control loop with closed LQR loop','Added Inductor ESR 300 m\Omega','location','best');
hold off
subplot(2,1,2)
semilogx(pi_open_lqr_closed_data(:,1), pi_open_lqr_closed_data(:,3),'-','color','black');
hold on
semilogx(pi_open_lqr_closed_L_ESR_data(:,1), pi_open_lqr_closed_L_ESR_data(:,3),'--','color','black');
grid on
%xlim([min(pi_open_lqr_closed_L_ESR_data(:,1)) max(pi_open_lqr_closed_L_ESR_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
legend('PI control loop with closed LQR loop','Added Inductor ESR 300 m\Omega','location','best');
hold off

%% PI open, LQR closed, ESR 300 mOhm, 240pF ceramic
log_name = 'tf_0dbm_pi_open_lqr_closed_lout_300mohmESR_240pF_both';
pi_open_lqr_closed_L_ESR_Cparallel120pF_table = readtable([path log_name '.log']);
pi_open_lqr_closed_L_ESR_Cparallel120pF_data = table2array(pi_open_lqr_closed_L_ESR_Cparallel120pF_table);

figure(3)
subplot(2,1,1)
semilogx(pi_open_lqr_closed_data(:,1), pi_open_lqr_closed_data(:,2),'-','color','black');
hold on
semilogx(pi_open_lqr_closed_L_ESR_data(:,1), pi_open_lqr_closed_L_ESR_data(:,2), '--','color','black');
semilogx(pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,1), pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,2), ':','color','black');
grid on
%xlim([min(pi_open_lqr_closed_L_ESR_data(:,1)) max(pi_open_lqr_closed_L_ESR_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
%legend('PI control loop with closed LQR loop','Inductor ESR 300 m\Omega','Inductor parallel 120pF','location','best');
hold off
subplot(2,1,2)
semilogx(pi_open_lqr_closed_data(:,1), pi_open_lqr_closed_data(:,3),'-','color','black');
hold on
semilogx(pi_open_lqr_closed_L_ESR_data(:,1), pi_open_lqr_closed_L_ESR_data(:,3),'--','color','black');
semilogx(pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,1), pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,3),':','color','black');
grid on
%xlim([min(pi_open_lqr_closed_L_ESR_data(:,1)) max(pi_open_lqr_closed_L_ESR_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
legend('PI control loop with closed LQR loop','Inductor ESR 300 m\Omega','Inductor parallel 120pF','location','best');
hold off

%% General
%clear all; clc; close all
% Load path
%path = 'C:\Users\JeppeLaptop\OneDrive - Danmarks Tekniske Universitet\6. semester\Diplomingeniørprojekt\Målinger\Bode 100\';

%% PI, LQR open
log_name = 'tf_0dbm_pi_open_lqr_open';
pi_open_lqr_open_table = readtable([path log_name '.log']);

pi_open_lqr_open_data = table2array(pi_open_lqr_open_table);
figure(1)
subplot(2,1,1)
semilogx(pi_open_lqr_open_data(:,1), pi_open_lqr_open_data(:,2), 'k');
hold on
grid on
xlim([min(pi_open_lqr_open_data(:,1)) max(pi_open_lqr_open_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold off
subplot(2,1,2)
semilogx(pi_open_lqr_open_data(:,1), pi_open_lqr_open_data(:,3),'','color','black');
hold on
grid on
xlim([min(pi_open_lqr_open_data(:,1)) max(pi_open_lqr_open_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
hold off
%% PI, LQR open 300mohm ESR
log_name = 'tf_0dbm_pi_open_lqr_open_lout_300mohmESR_both';
pi_open_lqr_open_L_ESR_table = readtable([path log_name '.log']);
pi_open_lqr_open_L_ESR_data = table2array(pi_open_lqr_open_L_ESR_table);

figure(2)

subplot(2,1,1)
semilogx(pi_open_lqr_open_data(:,1), pi_open_lqr_open_data(:,2),'-','color','black');
hold on
semilogx(pi_open_lqr_open_L_ESR_data(:,1), pi_open_lqr_open_L_ESR_data(:,2), '--','color','black');
grid on
%xlim([min(pi_open_lqr_closed_L_ESR_data(:,1)) max(pi_open_lqr_closed_L_ESR_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend('PI control loop with closed LQR loop','Added Inductor ESR 300 m\Omega','location','best');
hold off
subplot(2,1,2)
semilogx(pi_open_lqr_open_data(:,1), pi_open_lqr_open_data(:,3),'-','color','black');
hold on
semilogx(pi_open_lqr_open_L_ESR_data(:,1), pi_open_lqr_open_L_ESR_data(:,3),'--','color','black');
grid on
%xlim([min(pi_open_lqr_closed_L_ESR_data(:,1)) max(pi_open_lqr_closed_L_ESR_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
legend('PI control loop with open LQR loop','Inductor ESR 300 m\Omega','location','best');
hold off

%% PI open, LQR open ESR 300 mOhm, 240pF ceramic
log_name = 'tf_0dbm_pi_open_lqr_open_lout_300mohmESR_240pF_both';
pi_open_lqr_open_L_ESR_Cparallel120pF_table = readtable([path log_name '.log']);
pi_open_lqr_open_L_ESR_Cparallel120pF_data = table2array(pi_open_lqr_open_L_ESR_Cparallel120pF_table);

figure(3)
subplot(2,1,1)
semilogx(pi_open_lqr_open_data(:,1), pi_open_lqr_open_data(:,2),'-','color','black');
hold on
semilogx(pi_open_lqr_open_L_ESR_data(:,1), pi_open_lqr_open_L_ESR_data(:,2), '--','color','black');
semilogx(pi_open_lqr_open_L_ESR_Cparallel120pF_data(:,1), pi_open_lqr_open_L_ESR_Cparallel120pF_data(:,2), ':','color','black');
grid on
%xlim([min(pi_open_lqr_closed_L_ESR_data(:,1)) max(pi_open_lqr_closed_L_ESR_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
%legend('PI control loop with closed LQR loop','Inductor ESR 300 m\Omega','Inductor parallel 120pF','location','best');
hold off
subplot(2,1,2)
semilogx(pi_open_lqr_open_data(:,1), pi_open_lqr_open_data(:,3),'-','color','black');
hold on
semilogx(pi_open_lqr_open_L_ESR_data(:,1), pi_open_lqr_open_L_ESR_data(:,3),'--','color','black');
semilogx(pi_open_lqr_open_L_ESR_Cparallel120pF_data(:,1), pi_open_lqr_open_L_ESR_Cparallel120pF_data(:,3),':','color','black');
grid on
%xlim([min(pi_open_lqr_closed_L_ESR_data(:,1)) max(pi_open_lqr_closed_L_ESR_data(:,1))])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
legend('PI control loop with open LQR loop','Inductor ESR 300 m\Omega','Inductor parallel 120pF','location','best');
hold off

%% Comparison of parasitic elements open and closed LQR

figure(4)
subplot(2,1,1)
semilogx(pi_open_lqr_open_L_ESR_Cparallel120pF_data(:,1), pi_open_lqr_open_L_ESR_Cparallel120pF_data(:,2), '--','color','black');
hold on
semilogx(pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,1), pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,2),':','color','black');
grid on
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold off
subplot(2,1,2)
semilogx(pi_open_lqr_open_data(:,1), pi_open_lqr_open_data(:,3),'--','color','black');
hold on
semilogx(pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,1), pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,3), ':','color','black');
grid on
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
legend('PI control loop with closed LQR loop and parasitic elements','PI control loop with open LQR loop and parasitic elements','location','best');
hold off

%% Transfer function estimation
%% First
%Freq = pi_open_lqr_closed_4ohm_data(:,1);
%Mag_4ohm = 10.^(pi_open_lqr_closed_4ohm_data(:,2)./20);
%Pha_4ohm = pi_open_lqr_closed_4ohm_data(:,3);
%zfr_4ohm = Mag_4ohm.*exp(1i*Pha_4ohm*pi/180);
%W = Freq*2*pi;
%Ts = 1/(2*max(Freq)); % your sampling time
%gfr_4ohm = frd(zfr_4ohm,W,Ts);
%sys=tfest(gfr_4ohm,4);
%[mag,ph]=bode(sys,W);
%[poles,zeros] = pzmap(sys);
% figure(5)
% compare(gfr_4ohm,sys)
% grid on
% figure(6)
% pzmap(sys)
% grid on
%% Load data

% Load first system
Freq_closed = pi_open_lqr_closed_data(:,1);
Mag_closed = 10.^(pi_open_lqr_closed_data(:,2)./20);
Pha_closed = pi_open_lqr_closed_data(:,3);
zfr_closed = Mag_closed.*exp(1i*Pha_closed*pi/180);
W = Freq_closed*2*pi;
Ts_closed = 1/(2*max(Freq_closed)); % your sampling time
gfr_closed = frd(zfr_closed,W,Ts_closed);
sys_closed=tfest(gfr_closed,8);
% Parasites 1
Mag_closed_para = 10.^(pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,2)./20);
Pha_closed_para = pi_open_lqr_closed_L_ESR_Cparallel120pF_data(:,3);
zfr_closed_para = Mag_closed_para.*exp(1i*Pha_closed_para*pi/180);
gfr_closed_para = frd(zfr_closed_para,W,Ts_closed);
sys_closed_para=tfest(gfr_closed_para,5);

% Load second system
Freq_open = pi_open_lqr_open_data(:,1);
Mag_open = 10.^(pi_open_lqr_open_data(:,2)./20);
Pha_open = pi_open_lqr_open_data(:,3);
zfr_open = Mag_open.*exp(1i*Pha_open*pi/180);
W2 = Freq_open*2*pi;
Ts_open = 1/(2*max(Freq_open)); % your sampling time
gfr_open = frd(zfr_open,W2,Ts_open);
sys_open = tfest(gfr_open,8);
% Parasites 2
Mag_open_para = 10.^(pi_open_lqr_open_L_ESR_Cparallel120pF_data(:,2)./20);
Pha_open_para = pi_open_lqr_open_L_ESR_Cparallel120pF_data(:,3);
zfr_open_para = Mag_open_para.*exp(1i*Pha_open_para*pi/180);
gfr_open_para = frd(zfr_open_para,W2,Ts_open);
sys_open_para=tfest(gfr_open_para,5);

%%
% Display
% Wrap phase to [-180,180]
opts = bodeoptions('cstprefs');
opts.PhaseWrapping='on';
%opts.PhaseWrappingBranch=-180;

[mag2,ph2]=bode(sys_closed,W,opts);
[mag2_para,ph2_para]=bode(sys_closed_para,W);
[poles,zeros] = pzmap(sys_closed);
%figure(5)
%pzmap(sys_closed)
%grid on
figure(6)
%subplot(2,1,1)
%compare(gfr,sys2)
%compare(gfr_closed_para,sys2_para)
%grid on
%figure(7)
%subplot(2,1,2)
hold on
bode(sys_closed,W,opts,'k');
bode(sys_closed_para,W,opts,'--k');
%[mag2_para,ph2_para]=bode(sys2_para,W);
grid on
legend('Closed loop without inductor parasitic elements', 'Closed loop with inductor parasitic elements', 'location','best');
hold off
figure(8)
hold on
tsum_closed=sys_closed_para/sys_closed;
tsum_open=sys_open_para/sys_open;
bode(tsum_closed,W,opts,'-k');
bode(tsum_open,W2,opts,'--k');
legend('LQR closed, est. TF of R_{L} and C_{L}', 'LQR open, est. TF of R_{L} and C_{L}', 'location','best');
grid on
hold off
%w=pi_open_lqr_closed_4ohm_data(:,1)*2*pi; %convert Hz to rad/sec 
%gfr = idfrd(response,w,Ts);
%sys=tfest(gfr,4);
%tfest = tfest(pi_open_lqr_closed_4ohm_data,4)


%% TF estimation
opts = bodeoptions('cstprefs');
opts.PhaseWrapping='on';
opts.PhaseWrappingBranch=-180;
%lintf = tf([4.285e12,8.578e16],[1,1.263e5,1.428e11,0]);
lintf_real = tf([1.428e11,4.761e15],[1,1.263e5,1.428e11,0]);
lintf_rl = tf([2.142e12,4.289e16],[1,2.959e5,1.643e11,0]);
lintf_cs = tf([4.159e12,8.326e16],[1,1.225e5,1.386e11,0]);
%[lintf_mag,lintf_pha,lintf_w] = bode(lintf);
%lintf_mag=squeeze(lintf_mag);lintf_pha=squeeze(lintf_pha);
lintf_rl2 = lintf_rl/lintf_real;
lintf_cs2 = lintf_cs/lintf_real;
lintf_real2 = lintf_real*lintf_rl2*lintf_cs2;
figure(9)
hold on
bode(sys_open,Freq_open,opts,'-k')
bode(lintf_real,Freq_open,opts,'--k')
grid on
legend('Open loop measurement','Theoretical function','location','best');
hold off

figure(10)
hold on
bode(sys_open_para,Freq_open,opts,'-k')
bode(linsys1_tf,Freq_open,opts,'--k')
grid on
legend('Open loop measurement','Theoretical function','location','best');
hold off