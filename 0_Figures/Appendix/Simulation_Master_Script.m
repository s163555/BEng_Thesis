%% Clear all
clear all
close all
clc
%% Set up parameters in circuit
% Supply parameters
Vs = 5;
Vref = 2.5;
% Input signal parameters
InDC = 2.5;
InAmpl = 0.9;
InFreq = 10e3;
InDelay = 200e-6;
% Modulator parameters
Rin = 1.5e3;
Rfb = 2.7e3;
C_aim = 1.5e-9;
Cin = 1e-15;
R1 = 2.21e3;
% Notat: Open loop fungerer bedst med 2.2k og closed loop fungerer bedst
% med 1k (?)
R2 = 20e3;
% Output filter parameters
Lind = 1.768e-6;
Cbtl = 1.98e-6;
Rbtl = 4;
shunt = 15e-3;
% Control loop parameters
R_cur_sense_fb = 16e3;
R_spkv_fb = 8.3e3;
R_voltacq = 60e3;
C_cur_ref = 0.1e-6;
% PI controller
R_pic = 40e3;
C_pic = 1e-6;
% Simulation control
SimStop = 600e-6;
SimDelay = 200e-6;
MaxDeltaT = 5e-9;
% Schmitt trigger hysteresis
Vt = 2.5;
Vh = 0.05;

%% CIRCUIT: CLASSD OPEN LOOP
%% Export data to SPICE directive files
% Add circuit models
%FID_models = fopen('../LTspice/models.txt', 'w');
%fprintf(FID_models, '.model SW SW(Vt=%d)\n', V_t);
%fprintf(FID_models, '.model D_SW D(Ron=%d Roff=%d Vfwd=%d Vrev=%d)\n', R_on_SW, R_off_SW, V_fwd_SW, V_rev_SW);
%fprintf(FID_models, '.model D D(Ron=%d Roff=%d Vfwd=%d Vrev=%d)\n', R_on_D, R_off_D, V_fwd_D, V_rev_D);
%fclose(FID_models);

% Circuit: Class D open loop
path_LTSpice = '"C:/Program Files/LTC/LTspiceXVII/XVIIx64.exe"';
circuit_path = 'C:/Users/JeppeLaptop/Danmarks Tekniske Universitet/BEng Jeppe audio amplifier control loop - General/Simuleringsfiler/classd_openloop/';
circuit_name = 'classd_open';
command = ' -Run -ascii -b ';


% Add circuit parameters
FID_param = fopen(append(circuit_path,'param.txt'), 'w');
% Define components and parameters
param_nam = ["Vs", "Vref", "InDC", "InAmpl", "InFreq", "InDelay", "Rin", "Rfb", "C_aim", "Cin", "R1", "R2", "Cbtl", "Rbtl", "shunt", "R_cur_sense_fb", "R_spkv_fb", "R_voltacq", "C_cur_ref", "R_pic", "C_pic", "SimStop", "SimDelay", "MaxDeltaT"];
param_val = [Vs, Vref, InDC, InAmpl, InFreq, InDelay, Rin, Rfb, C_aim, Cin, R1, R2, Cbtl, Rbtl, shunt, R_cur_sense_fb, R_spkv_fb, R_voltacq, C_cur_ref, R_pic, C_pic, SimStop, SimDelay, MaxDeltaT];
% Define symmetrical components (will be printed for pos and neg circuit)
param_nam_sym = ["Lind"];
param_val_sym = [Lind];
% Export parameters to external file
for i=1:length(param_nam)
    fprintf(FID_param, '.param %s = %d\n', param_nam(i), param_val(i));
end
% Export symmetrical parameters to external file
for i=1:length(param_nam_sym)
    fprintf(FID_param, '.param %s_P = %d\n', param_nam_sym(i), param_val_sym(i));
    fprintf(FID_param, '.param %s_N = %d\n', param_nam_sym(i), param_val_sym(i));
end
fclose(FID_param);
%% Simulate and import data
% Run LTSpice simulation
command = [path_LTSpice command '"' circuit_path circuit_name '.asc"'];
dos(command);

% Import data
raw_data = LTspice2Matlab([circuit_path circuit_name '.raw']);

% Display variables
fprintf('var \t name\n')
for i=1:raw_data.num_variables
    fprintf('%d \t\t %s \n', raw_data.selected_vars(i), raw_data.variable_name_list{i})
end
%% Plot data
% Trace numbers
var_input = 1;
var_spk_p = 2;
var_spk_m = 3;

Open_TimeVect = raw_data.time_vect;
Open_varspk_p = raw_data.variable_mat(var_spk_p,:);
Open_varspk_m = raw_data.variable_mat(var_spk_m,:);
Open_input = raw_data.variable_mat(var_input,:);

figure(1)
subplot(2,1,1)
hold on
plot(raw_data.time_vect, raw_data.variable_mat(var_spk_p,:), '--k');
plot(raw_data.time_vect, raw_data.variable_mat(var_spk_m,:), ':k');
title(sprintf('Waveform %s and %s', raw_data.variable_name_list{var_spk_p},raw_data.variable_name_list{var_spk_m}));
grid on
ylabel(raw_data.variable_type_list{var_input});
xlabel('Time [sec]' );
xlim([SimDelay,SimStop]);
%ylim([min(raw_data.variable_mat(var_input,:)),max(raw_data.variable_mat(var_input,:))])
legend(raw_data.variable_name_list{var_spk_p}, raw_data.variable_name_list{var_spk_m},'location','best');
hold off
subplot(2,1,2)
hold on
plot(raw_data.time_vect, (raw_data.variable_mat(var_spk_p,:)-raw_data.variable_mat(var_spk_m,:)), 'k');
title(sprintf('Waveform %s-%s', raw_data.variable_name_list{var_spk_p},raw_data.variable_name_list{var_spk_m}));
grid on
ylabel(raw_data.variable_type_list{var_spk_p});
xlabel('Time [sec]' );
xlim([SimDelay,SimStop]);
ylim([min(raw_data.variable_mat(var_spk_p,:)-raw_data.variable_mat(var_spk_m,:)),max(raw_data.variable_mat(var_spk_p,:)-raw_data.variable_mat(var_spk_m,:))])
legend(sprintf('%s-%s', raw_data.variable_name_list{var_spk_p},raw_data.variable_name_list{var_spk_m}));
hold off
%sgt = sgtitle('Open loop','Color','black');
%sgt.FontSize = 20;

%% CIRCUIT: CLASSD REGULATOR
%% Export data to SPICE directive files
% Add circuit models
%FID_models = fopen('../LTspice/models.txt', 'w');
%fprintf(FID_models, '.model SW SW(Vt=%d)\n', V_t);
%fprintf(FID_models, '.model D_SW D(Ron=%d Roff=%d Vfwd=%d Vrev=%d)\n', R_on_SW, R_off_SW, V_fwd_SW, V_rev_SW);
%fprintf(FID_models, '.model D D(Ron=%d Roff=%d Vfwd=%d Vrev=%d)\n', R_on_D, R_off_D, V_fwd_D, V_rev_D);
%fclose(FID_models);

% Circuit: Class D open loop
path_LTSpice = '"C:/Program Files/LTC/LTspiceXVII/XVIIx64.exe"';
circuit_path = 'C:/Users/JeppeLaptop/Danmarks Tekniske Universitet/BEng Jeppe audio amplifier control loop - General/Simuleringsfiler/classd_PI_LQR/';
circuit_name = 'classd_PI_LQR';
command = ' -Run -ascii -b ';

% Add circuit parameters
FID_param = fopen(append(circuit_path,'param.txt'), 'w');
% Define components and parameters
param_nam = ["Vs", "Vref", "InDC", "InAmpl", "InFreq", "InDelay", "Rin", "Rfb", "C_aim", "Cin", "R1", "R2", "Cbtl", "Rbtl", "shunt", "R_cur_sense_fb", "R_spkv_fb", "C_cur_ref", "R_pic", "C_pic", "SimStop", "SimDelay", "MaxDeltaT"];
param_val = [Vs, Vref, InDC, InAmpl, InFreq, InDelay, Rin, Rfb, C_aim, Cin, R1, R2, Cbtl, Rbtl, shunt, R_cur_sense_fb, R_spkv_fb, C_cur_ref, R_pic, C_pic, SimStop, SimDelay, MaxDeltaT];
% Define symmetrical components (will be printed for pos and neg circuit)
param_nam_sym = ["Lind", "R_voltacq"];
param_val_sym = [Lind, R_voltacq];
% Export parameters to external file
for i=1:length(param_nam)
    fprintf(FID_param, '.param %s = %d\n', param_nam(i), param_val(i));
end
% Export symmetrical parameters to external file
for i=1:length(param_nam_sym)
    fprintf(FID_param, '.param %s_P = %d\n', param_nam_sym(i), param_val_sym(i));
    fprintf(FID_param, '.param %s_N = %d\n', param_nam_sym(i), param_val_sym(i));
end
fclose(FID_param);
%% Simulate and import data
% Run LTSpice simulation
command = [path_LTSpice command '"' circuit_path circuit_name '.asc"'];
dos(command);

% Import data
raw_data = LTspice2Matlab([circuit_path circuit_name '.raw']);

%Print variables
fprintf('var \t name\n')
for i=1:raw_data.num_variables
    fprintf('%d \t\t %s \n', raw_data.selected_vars(i), raw_data.variable_name_list{i})
end
%% Plot data
% Trace numbers
var_spk_p = 1;
var_spk_m = 2;
var_input = 20;
var_aim_n = 8;
var_aim_p = 10;

Regulator_TimeVect = raw_data.time_vect;
Regulator_varspk_p = raw_data.variable_mat(var_spk_p,:);
Regulator_varspk_m = raw_data.variable_mat(var_spk_m,:);

figure(2)
subplot(2,1,1)
hold on
plot(raw_data.time_vect, raw_data.variable_mat(var_spk_p,:), '--k');
plot(raw_data.time_vect, raw_data.variable_mat(var_spk_m,:), ':k');
title(sprintf('Waveform %s and %s', raw_data.variable_name_list{var_spk_p},raw_data.variable_name_list{var_spk_m}));
grid on
ylabel(raw_data.variable_type_list{var_spk_p});
xlabel('Time [sec]' );
xlim([SimDelay,SimStop]);
%ylim([min(raw_data.variable_mat(var_input,:)),max(raw_data.variable_mat(var_input,:))])
legend(raw_data.variable_name_list{var_spk_p}, raw_data.variable_name_list{var_spk_m},'location','best');
hold off
subplot(2,1,2)
hold on
plot(raw_data.time_vect, (raw_data.variable_mat(var_spk_p,:)-raw_data.variable_mat(var_spk_m,:)), 'k');
title(sprintf('Waveform %s-%s', raw_data.variable_name_list{var_spk_p},raw_data.variable_name_list{var_spk_m}));
grid on
ylabel(raw_data.variable_type_list{var_spk_p});
xlabel('Time [sec]' );
xlim([SimDelay,SimStop]);
ylim([min(raw_data.variable_mat(var_spk_p,:)-raw_data.variable_mat(var_spk_m,:)),max(raw_data.variable_mat(var_spk_p,:)-raw_data.variable_mat(var_spk_m,:))])
legend(sprintf('%s-%s', raw_data.variable_name_list{var_spk_p},raw_data.variable_name_list{var_spk_m}));
hold off
%sgt = sgtitle('PI+LQR closed loop','Color','black');
%sgt.FontSize = 20;

figure(3)
plot(raw_data.time_vect, raw_data.variable_mat(var_aim_n,:), '--k');
hold on
plot(raw_data.time_vect, raw_data.variable_mat(var_aim_p,:), '-k');
title(sprintf('Waveform %s and %s', raw_data.variable_name_list{var_aim_p},raw_data.variable_name_list{var_aim_n}));
grid on
xlabel('Time [sec]' );
ylabel(raw_data.variable_type_list{var_aim_p});
xlim([SimDelay,SimDelay+5e-6]);
ylim([2,3]);
legend(raw_data.variable_name_list{var_aim_p}, raw_data.variable_name_list{var_aim_n},'location','best');
hold off
%% Plot both open and regulator
figure(3)
hold on
plot(Regulator_TimeVect, (Regulator_varspk_p-Regulator_varspk_m),'-','color','black');
plot(Open_TimeVect, (Open_varspk_p-Open_varspk_m),'-','color','blue');
ylabel(raw_data.variable_type_list{var_spk_p});
xlabel('Time [sec]' );
xlim([SimDelay,SimStop]);
legend({'PI+LQR','Open loop'},'FontSize',12,'Location','best')
hold off