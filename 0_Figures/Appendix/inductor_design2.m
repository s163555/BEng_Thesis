%% Clear workspace
clc
clear all
close all

%% Design specifications
% Amplifier
L           = 1.768e-6;            % Inductance                        [H]
V_cc        = 30;                   % Rail-to-rail voltage              [V]
R_spk       = 4;                    % Speaker load                      [ohm]
I_max       = V_cc/R_spk;           % Max output current                [A]

% Core(s)
core        = ['T80-2'; 'T94-2'];   % Core type
B_sat       = [400e-3, 400e-3];     % Max allowable flux density        [T]
A_ef        = [.231e-4, .362e-4];   % Effective cross-sectional area    [m^2]
A_L         = [5.5e-9, 8.4e-9];     % Inductance rating                 [H/N^2]
% Core physical dimensions
OD          = [20.2e-3, 23.9e-3];   % Outer diameter                    [m]
ID          = [12.6e-3, 14.2e-3];   % Inner diameter                    [m]
Ht          = [6.35e-3, 7.92e-3];   % Height/thickness                  [m]
l           = [5.14e-2, 5.97e-3];   % 

% Wire specifications
D_wire      = 0.95e-3;              % Wire diameter

%% Calculations
vectorSize  = size(core,1);
N           = zeros(1,vectorSize);
B           = zeros(1,vectorSize);

% Calculate required number of turns N and peak core flux B
for i=1:vectorSize
    N(i)    = (L/A_L(i))^(1/2);
    B(i)    = (L*I_max) / (A_ef(i) * N(i));
end

% Calculate minimum required wire length
l_wire      = zeros(1,vectorSize);
for i=1:vectorSize
    r_mean  = OD(i)/2 - (OD(i)-ID(i))/4;
    l_wire(i) = N(i) * (2*Ht(i)+ (OD(i)-ID(i)) + 4*D_wire) + 2*pi*r_mean;
end

%% Print results
fprintf('\nDESIGN SPECIFICATIONS\n')
fprintf(['L        = ' num2str(L*1e6) ' uH\n'])
fprintf(['V_cc     = ' num2str(V_cc) ' V\n'])
fprintf(['R_spk    = ' num2str(R_spk) ' ohm\n'])
fprintf(['I_max    = ' num2str(I_max) ' A\n'])

for i=1:vectorSize
    fprintf('--------------------------------------------------')
    fprintf(['\nCORE ' core(i,:) '\n'])
    fprintf('parameters:\n')
    fprintf(['B_sat    = ' num2str(B_sat(i)*1e3) ' mT\n'])
    fprintf(['A_ef     = ' num2str(A_ef(i)*1e4) ' cm^2\n'])
    fprintf(['A_L      = ' num2str(A_L(i)*1e9) ' nH/N^2\n'])
    fprintf('results:\n')
    fprintf(['N        = ' num2str(N(i)) ' turns\n'])
    fprintf(['B        = ' num2str(B(i)*1e3) ' mT\n'])
    fprintf(['l_wire   = >' num2str(l_wire(i)*1e2) ' cm\n'])
    fprintf('Evaluation: ')
    if (B(i) < B_sat(i))
        fprintf('Valid design (B < B_sat)\n')
    else
        fprintf('Invalid design (B > B_sat)\n')
    end
end
