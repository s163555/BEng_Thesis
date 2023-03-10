%% Clear all
clear all
close all
clc
%% Set up
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\';
circuit_name = '2020-11-24 - self oscillating class-d cleanup test filter_v2';

raw_data = LTspice2Matlab([circuit_path circuit_name '.raw']);

%% Plot results
% Display plot variables
fprintf('PLOT VARIABLES\nvar \t name\n')
for i=1:raw_data.num_variables
    fprintf('%d \t\t %s \n', raw_data.selected_vars(i), raw_data.variable_name_list{i})
end

% Define plot variables
var_V_signal_in    = 1;    % Input voltage
var_v_SPKP    = 2;    % Switch voltage
var_v_SPKM   = 3;    % Resonant tank voltage
var_V_rec   = 4;    % Diode voltage
var_V_o     = 5;    % Output voltage
var_I_D     = 6;    % Diode current
var_I_L_R   = 7;    % Resonant tank current
var_I_Lf    = 8;    % Input current
var_I_o     = 9;    % Output current
var_i_SW    = 10;   % Switch current

%% waveforms+calculations
% Calculate output power and resonant component voltages
%V_CR        = raw_data.variable_mat(var_v_SPKP,:)-raw_data.variable_mat(var_v_SPKM,:);
%V_LR        = raw_data.variable_mat(var_v_SPKM,:)-raw_data.variable_mat(var_V_rec,:);
V_R_BTL      = raw_data.variable_mat(var_v_SPKP,:)-raw_data.variable_mat(var_v_SPKM,:);

% Plots
figure(1)
% V_IN and v_SW
subplot(2,1,1)
plot(raw_data.time_vect, raw_data.variable_mat(var_V_signal_in,:), 'k');
hold on
%plot(raw_data.time_vect, raw_data.variable_mat(var_v_SPKP,:), 'k', 'linestyle', '--');
hold off
%legend('Input voltage', 'Switch voltage')
grid on
xlim([min(raw_data.time_vect) max(raw_data.time_vect)])
title('Input voltage V_{IN}');
ylabel('Voltage [V]');
xlabel('Time [s]');
% Output voltage and diode voltage
subplot(2,1,2)
%plot(raw_data.time_vect, raw_data.variable_mat(var_V_o,:), 'k', 'linestyle', '--');
plot(raw_data.time_vect, V_R_BTL, 'k', 'linestyle', '-');
hold on
%plot(raw_data.time_vect, raw_data.variable_mat(var_V_rec,:), 'k');
hold off
%legend('Output voltage', 'Diode voltage')
grid on
xlim([min(raw_data.time_vect) max(raw_data.time_vect)])
title('Output voltage across bridge-tied-load R_{L}');
ylabel('Voltage [V]');
xlabel('Time [s]');

% Export
fig_width   = 25;
fig_height  = 25;
%set(gcf,'paperunits','centimeters','Paperposition',[0 0 fig_width fig_height]);
%set(gcf,'units','centimeters','Position',[0 1 fig_width fig_height]);
%saveas(gcf, 'figures_con/waveforms.eps', 'epsc')
%saveas(gcf, 'figures_con/waveforms_tuned.eps', 'epsc')

%% Plot preamp input filter ac
% load
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\Simuleringsarkiv\';
circuit_name = 'input_filter';
raw_data = LTspice2Matlab([circuit_path circuit_name '.raw']);
% Display variables
fprintf('var \t name\n')
for i=1:raw_data.num_variables
    fprintf('%d \t\t %s \n', raw_data.selected_vars(i), raw_data.variable_name_list{i})
end
var_out = 3;
Log_Magnitude_dB = 20*log10(abs(raw_data.variable_mat(var_out,:)));
Norm_Phase_Degrees = angle(raw_data.variable_mat(var_out,:))*180/pi;
Freq_Vect = raw_data.freq_vect;
figure(2)
subplot(2,1,1)
semilogx(Freq_Vect, Log_Magnitude_dB, 'k');
hold on
grid on
xlim([min(raw_data.freq_vect) max(raw_data.freq_vect)])
%ylim([-20 8])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold off
subplot(2,1,2)
semilogx(Freq_Vect, Norm_Phase_Degrees, 'k');
hold on
grid on
xlim([min(raw_data.freq_vect) max(raw_data.freq_vect)])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
hold off
%sgt = sgtitle('PI controller frequency response','Color','black');
%sgt.FontSize = 20;
%% Plot PI ac analysis
% load
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\Simuleringsarkiv\';
circuit_name = 'PI_regulator_ac';
raw_data = LTspice2Matlab([circuit_path circuit_name '.raw']);
% Display variables
fprintf('var \t name\n')
for i=1:raw_data.num_variables
    fprintf('%d \t\t %s \n', raw_data.selected_vars(i), raw_data.variable_name_list{i})
end
var_out = 1;
Log_Magnitude_dB = 20*log10(abs(raw_data.variable_mat(var_out,:)));
Norm_Phase_Degrees = angle(raw_data.variable_mat(var_out,:))*180/pi;
Freq_Vect = raw_data.freq_vect;
figure(2)
subplot(2,1,1)
semilogx(Freq_Vect, Log_Magnitude_dB, 'k');
hold on
grid on
xlim([min(raw_data.freq_vect) max(raw_data.freq_vect)])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
hold off
subplot(2,1,2)
semilogx(Freq_Vect, Norm_Phase_Degrees, 'k');
hold on
grid on
xlim([min(raw_data.freq_vect) max(raw_data.freq_vect)])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
hold off
%sgt = sgtitle('PI controller frequency response','Color','black');
%sgt.FontSize = 20;
%% Plot output filter frequency response
% load
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\Simuleringsarkiv\';
circuit_name = 'output_filter_smallsignal';
raw_data = LTspice2Matlab([circuit_path circuit_name '.raw']);
% Display variables
fprintf('var \t name\n')
for i=1:raw_data.num_variables
    fprintf('%d \t\t %s \n', raw_data.selected_vars(i), raw_data.variable_name_list{i})
end
var_out = 2;
var_out_l_esr = 4;
var_out_l_c = 6;
var_out_c_esr = 8;
var_out_c_esl = 10;

Log_Magnitude_dB = 20*log10(abs(raw_data.variable_mat(var_out,:)));
Log_Magnitude_dB_l_esr = 20*log10(abs(raw_data.variable_mat(var_out_l_esr,:)));
Log_Magnitude_dB_l_c = 20*log10(abs(raw_data.variable_mat(var_out_l_c,:)));
Log_Magnitude_dB_c_esr = 20*log10(abs(raw_data.variable_mat(var_out_c_esr,:)));
Log_Magnitude_dB_c_esl = 20*log10(abs(raw_data.variable_mat(var_out_c_esl,:)));
Norm_Phase_Degrees = angle(raw_data.variable_mat(var_out,:))*180/pi;
Norm_Phase_Degrees_l_esr = angle(raw_data.variable_mat(var_out_l_esr,:))*180/pi;
Norm_Phase_Degrees_l_c = angle(raw_data.variable_mat(var_out_l_c,:))*180/pi;
Norm_Phase_Degrees_c_esr = angle(raw_data.variable_mat(var_out_c_esr,:))*180/pi;
Norm_Phase_Degrees_c_esl = angle(raw_data.variable_mat(var_out_c_esl,:))*180/pi;
Freq_Vect = raw_data.freq_vect(1,:);
figure(2)
subplot(2,1,1)
hold on
semilogx(Freq_Vect, Log_Magnitude_dB, 'color', 'black');
semilogx(Freq_Vect, Log_Magnitude_dB_l_esr, '--', 'color', 'black');
semilogx(Freq_Vect, Log_Magnitude_dB_l_c, ':', 'color', 'black');
grid on
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xlim([1e4 1e6])
set(gca, 'XScale', 'log');
legend('Ideal', 'Inductor ESR 300 m\Omega', 'Inductor parallel capacitance 120 pF', 'location', 'best');
hold off
subplot(2,1,2)
hold on
semilogx(Freq_Vect, Norm_Phase_Degrees, 'color', 'black');
semilogx(Freq_Vect, Norm_Phase_Degrees_l_esr, '--', 'color', 'black');
semilogx(Freq_Vect, Norm_Phase_Degrees_l_c, ':', 'color', 'black');
grid on
xlim([1e4 1e6])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
legend('Ideal', 'Inductor ESR 300 m\Omega', 'Inductor parallel capacitance 120 pF', 'location', 'best');
hold off
%sgt = sgtitle('PI controller frequency response','Color','black');
%sgt.FontSize = 20;

figure(3)
subplot(2,1,1)
semilogx(Freq_Vect, Log_Magnitude_dB, 'color', 'black');
hold on
semilogx(Freq_Vect, Log_Magnitude_dB_c_esr, '--', 'color', 'black');
semilogx(Freq_Vect, Log_Magnitude_dB_c_esl, ':', 'color', 'black');
grid on
xlim([1e4 1e6])
ylim([-60 20])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
%legend('Ideal', 'Inductor ESR 300 m\Omega, C_{s} 120 pF', 'Capacitor ESR 300 m\Omega, ESL 100 nH', 'location', 'southwest');
hold off
subplot(2,1,2)
semilogx(Freq_Vect, Norm_Phase_Degrees, 'color', 'black');
hold on
semilogx(Freq_Vect, Norm_Phase_Degrees_l_esr, '--', 'color', 'black');
semilogx(Freq_Vect, Norm_Phase_Degrees_l_c, ':', 'color', 'black');
grid on
xlim([1e4 1e6])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
legend('Ideal', 'Capacitor ESR 300 m\Omega', 'Capacitor ESL 100nH', 'location', 'best');
hold off


figure(4)
subplot(2,1,1)
semilogx(Freq_Vect, Log_Magnitude_dB, 'color', 'black');
hold on
semilogx(Freq_Vect, Log_Magnitude_dB_l_c, '--', 'color', 'black');
semilogx(Freq_Vect, Log_Magnitude_dB_c_esl, ':', 'color', 'black');
grid on
xlim([1e4 1e6])
ylim([-60 20])
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend('Ideal', 'Inductor ESR 300m\Omega, Parallel capacitance 120 pF', 'Capacitor ESR 300 m\Omega, ESL 100nH', 'location', 'best');
hold off
subplot(2,1,2)
semilogx(Freq_Vect, Norm_Phase_Degrees, 'color', 'black');
hold on
semilogx(Freq_Vect, Norm_Phase_Degrees_l_c, '--', 'color', 'black');
semilogx(Freq_Vect, Norm_Phase_Degrees_c_esl, ':', 'color', 'black');
grid on
xlim([1e4 1e6])
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
hold off
%% Plot AIM
% load
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\Simuleringsarkiv\';
circuit_name = 'aim_modulator_duty_full';
raw_data = LTspice2Matlab([circuit_path circuit_name '.log.raw'], []);
% Display variables
disp( sprintf('\n\nThis file contains %.0f variables:\n', raw_data.num_variables) );
fprintf('var \t name\n')
for i=1:raw_data.num_variables
    fprintf('%d \t\t %s \t\t %s \n', i, char(raw_data.variable_type_list{i}), raw_data.variable_name_list{i})
end
var_v_spk_p = 1;
var_v_spk_m = 2;
var_v_pos_out = 5;
var_v_neg_out = 7;
var_v_out = 8;
var_v_in = 11;
var_i_e1 = 41;
var_i_e2 = 40;
var_i_v1 = 42;
%raw_data = LTspice2Matlab([circuit_path circuit_name '.raw'], [var_v_spk_p,var_v_spk_m,var_v_pos_out,var_v_neg_out,var_v_out,var_v_in,var_i_e1,var_i_e2,var_i_v1]);
raw_data = LTspice2Matlab([circuit_path circuit_name '.raw'],[1,2,5,7,8,11,40,41,42],10);

%% Plot eff
% load
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\Simuleringsarkiv\';
circuit_name = 'aim_modulator_duty_full';
%para_data = readtable([circuit_path circuit_name '.eff_para.txt']);
%para_data = readtable([circuit_path circuit_name '.eff_esr.txt']);
eff_data = readtable([circuit_path circuit_name '.eff.bak.txt']);
eff_data = readtable([circuit_path circuit_name '.eff.txt']);
%para_array = table2array(para_data);
eff_array = table2array(eff_data);
%eff_array(:,1) = [];
%para_array(:,1) = [];
%s = 20; 
%for c = 1:s
%   eff_array(c,:) = []; 
%end
eff_array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],:) = [];
figure(4)
%semilogx(para_array(:,2),para_array(:,1),'k');
hold on
semilogx(eff_array(:,3),eff_array(:,2),'k');
grid on
xlim([min(eff_array(:,3)) max(eff_array(:,3))]);
ylim([0 100]);
New_XTickLabel = get(gca,'xtick');
set(gca,'XTickLabel',New_XTickLabel);
xlabel('Output Power [W]');
ylabel('Efficiency [%]');
%legend('With Inductor ESR','Without Inductor ESR','location','best');
hold off

%% Plot aim fsw duty
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\Simuleringsarkiv\';
circuit_name = 'aim_modulator_duty_full';
open_data = readtable([circuit_path circuit_name '.fosc.txt']);
open_array = table2array(open_data);

figure(5)
%subplot(1,2,1)
plot(open_array(:,1),open_array(:,3)*1e-3,'k')
hold on
%yticks([0 100.0e3 200e3 300e4 400e3 500e3 600e3])
grid on
xlabel('V_{in} [V]');
ylabel('f_{sw} [kHz]');
hold off
%subplot(1,2,2)
figure(6)
plot(open_array(:,1),open_array(:,2),'k')
hold on
grid on
xlabel('V_{in} [V]');
ylabel('D [%]');
ylim([0 100]);
hold off
%% Plot aim duty proper
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\Simuleringsarkiv\';
circuit_name = 'aim_modulator_duty_full';
open_data = readtable([circuit_path circuit_name '.duty.txt']);
open_array = table2array(open_data);
figure(6)
plot(open_array(:,1),open_array(:,2),'k')
hold on
grid on
xlabel('V_{in} [V]');
ylabel('D [%]');
ylim([0 1]);
xlim([1.3 3.7]);
hold off
%% Plot THD
circuit_path = 'C:\Users\JeppeLaptop\Danmarks Tekniske Universitet\BEng Jeppe audio amplifier control loop - General\Simuleringsfiler\Simuleringsarkiv\';
circuit_name = 'aim_modulator_duty_full';
data = readtable([circuit_path circuit_name '.thd.log']);
data_array = table2array(data);
figure(7)
plot(data_array(:,1)*1e-3,data_array(:,2),'k')
hold on
grid on
xlabel('f_{in} [kHz]');
ylabel('THD [%]');
xlim([min(data_array(:,1))*1e-3 max(data_array(:,1))*1e-3])
ylim([0 max(data_array(:,2))])
hold off
