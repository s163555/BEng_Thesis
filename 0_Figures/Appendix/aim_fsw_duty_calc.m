%% Clear workspace
clc
clear all
close all
exportDir = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\TEX\0_Figures\Synthesis\';
set(gcf,'defaultFigurePaperPositionMode','auto')
%% --------------- DESIGN SPECIFICATIONS --------------- %%
% Initial values
V_ref   = 2.5;          % Reference voltage
V_span  = 2;            % Exprected input voltage span
V_hw    = 0.5;          % Hysteresis window
V_out   = 5;            % Output voltage span (V_H - V_L)
f_idle  = 600e3;        % Idle switching frequency
C       = 1.5e-9;       % Chosen capacitor C

%% ----------------- DESIGN CALCULATION ---------------- %%
% HIGH and LOW output voltage (PWM)
V_H     = V_ref + V_out/2;
V_L     = V_ref - V_out/2;

% Hysteresis window upper and lower threshold voltage
V_thh   = V_hw * (V_H-V_ref)/(V_H-V_L) + V_ref;
V_thl   = V_hw * (V_L-V_ref)/(V_H-V_L) + V_ref;

% Calculate Thevenin resistance (evaluate at V_in=V_ref)
V_CH0   = V_ref + (V_span*(V_H-V_ref) + V_hw*(V_H-V_ref)) / (V_H-V_L+V_span);
V_CL0   = V_ref + (V_span*(V_L-V_ref) + V_hw*(V_L-V_ref)) / (V_H-V_L+V_span);
R_t     = 1 / (f_idle*C*log( ((V_CH0-V_thl)*(V_CL0-V_thh)) / ((V_CH0-V_thh)*(V_CL0-V_thl)) ));

% Input voltage range
V_in    = (V_ref-V_span/2):0.001:(V_ref+V_span/2);

% Carrier waveform upper and lower limit
V_CH    = V_in + (V_span*(V_H-V_in) + V_hw*(V_H-V_in)) / (V_H-V_L+V_span);
V_CL    = V_in + (V_span*(V_L-V_in) + V_hw*(V_L-V_in)) / (V_H-V_L+V_span);

% Modulator timings
tau     = R_t*C;
t_high  = tau * log((V_CH-V_thl) ./ (V_CH-V_thh));
t_low   = -tau * log((V_CL-V_thl) ./ (V_CL-V_thh));

% Modulator frequency
f_sw    = 1 ./ (t_high + t_low);
D       = t_high ./ (t_high + t_low);

% R_in & R_fb voltage divider
k2      = (V_span + V_hw) / (V_span + V_out);
% Calculate R_in & R_fb
R_fb    = (V_span + V_out) / (V_span + V_hw) * R_t;
R_in    = R_fb*R_t / (R_fb - R_t);

%% ----------------- PLOT RESULTS ---------------- %%
figure(1)
colororder({'k','k'})
%subplot(1,2,1)
% Switching frequency
%yyaxis left
%subplot(2,1,1)
plot(V_in,f_sw*1e-3,'k')
hold on
plot(sim_fosc_x,sim_fosc_y,'--k');
xlabel('V_{in} [V]')
ylabel('f_{sw} [kHz]');
xlim([(V_ref-V_span/2) (V_ref+V_span/2)]);
ylim([min(f_sw) max(f_sw*1e-3)])
grid on
hold off
legend('Calculated','Simulated')
%saveas(gcf,[exportDir 'modulator_fsw.eps'],'eps')

figure(2)
%subplot(1,2,2)
% Duty cycle
%yyaxis right
%subplot(2,1,2)
plot(V_in,D,'k','linestyle','-')
hold on
plot(sim_duty_x,sim_duty_y,'--k')
xlabel('V_{in} [V]')
ylabel('D');
grid on
legend('Calculated','Simulated')
xlim([(V_ref-V_span/2) (V_ref+V_span/2)]);
ylim([0 1]);
hold off
%saveas(gcf,[exportDir 'modulator_duty.eps'],'eps')

% Export
%fig_width   = 20;
%fig_height  = 8;
%set(gcf,'units','centimeters','Position',[0 2 fig_width fig_height]);
%saveas(gcf,[exportDir 'modulator_duty.eps'],'eps')