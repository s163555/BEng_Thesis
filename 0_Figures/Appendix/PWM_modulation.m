clc; clear all; close all;
Simulation_Time = 1;
input = 1.5; % Signal frekvens

sim('PWM_modulation1',Simulation_Time)


% Plot switching frequency
figure(1)
hold on
plot(get(ans,'pwm'),'-','color','black')
plot(get(ans,'input'),'--','color','black')
plot(get(ans,'carrier'),':','color','black')
xlabel('Normalised time [s]')
ylabel('Normalised amplitude [V]');
legend({'PWM','Input','Carrier'},'FontSize',12,'Location','northeast')
grid on
hold off