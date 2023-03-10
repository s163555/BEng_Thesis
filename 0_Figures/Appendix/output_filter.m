%% Analytical calculation of transfer function in output filter with parasitic elements
close all; clc; clear all;
opts = bodeoptions('cstprefs');
opts.PhaseWrapping='on'; % Bode plot phase wrap enable
%opts.PhaseWrappingBranch=-180;
opts.Xlim=[1e4, 1e6]; % Axis limitation
Rbtl = 4; Cbtl = 1.98e-6; % Differential filter values
R = Rbtl/2; C = Cbtl*2; L = 1.768e-6; R_L = 300e-3; C_L = 120e-9; % Single-ended filter values
sys = tf([1],[L*C, L/R, 1]);
sys_rl = tf([R],[R*L*C, L+R*R_L*C, R_L+R]);
sys_cl = tf([R*L*C_L, 0, R],[R*L*C_L+R*L*C, L, R]);
[mag,phase,wout] = bode(sys); % Get Plot Data
mag = squeeze(mag); % Reduce (1x1xN) Matrix To (1xN)
phase=squeeze(phase); % Reduce (1x1xN) Matrix To (1xN)
magr2 = (mag/max(mag)).^2; % Calculate Power Of Ratio Of mag/max(mag)
dB3 = interp1(magr2, [wout phase mag], 0.5, 'spline'); % Find corner frequency

% Display
figure(9)
hold on
bode(sys,opts,'-k');
bode(sys_rl,opts,'--k');
bode(sys_cl,opts,':k');
grid on
legend('Ideal', 'ESR 300 m\Omega', 'Parallel capacitance 120 pF','location','best')
hold off
