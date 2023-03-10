% Clear workspace
clc
clear all
close all

%% --------------- DESIGN SPECIFICATIONS --------------- %%
f_hp        = 10;       % High pass filter cutoff frequency
f_lp        = 80e3;     % Low pass filter cutoff frequency
G           = 2;        % Amplification
C_1         = 1e-6;     % High pass filter capacitor
C_2         = 1e-9;     % Low pass filter capacitor

%% --------------- CALCULATIONS --------------- %%
% Calculate ideal values
R_1         = 1/(2*pi*f_hp*C_1);    % High pass filter resistor
R_2         = 1/(2*pi*f_lp*C_2);    % Low pass filter resistor

% Chose component values
R_1         = 16.2e3;
R_2         = 2e3;

% Calculate ideal R_3
R_3         = R_2/(G-1);

f_hp = 1/(2*pi*R_1*C_1);
f_lp = 1/(2*pi*R_2*C_2);

%% --------------- PRINT RESULTS --------------- %%
fprintf('DESIGN SPECIFICATIONS (Input Filter)\n')
fprintf(['f_hp   = ' num2str(f_hp) ' Hz \n'])
fprintf(['f_lp   = ' num2str(f_lp*1e-3) ' kHz \n'])
fprintf(['G      = ' num2str(G) ' V/V \n'])
fprintf(['C_1    = ' num2str(C_1*1e6) ' uF\n'])
fprintf(['C_2    = ' num2str(C_2*1e9) ' nF\n'])

fprintf('\nCOMPONENT VALUES\n')
fprintf(['R_1  = ' num2str(R_1*1e-3,4) ' kOhm \n'])
fprintf(['R_2  = ' num2str(R_2*1e-3,4) ' kOhm \n'])
fprintf(['R_3  = ' num2str(R_3*1e-3,4) ' kOhm \n'])

