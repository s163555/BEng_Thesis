%clc; clear all; close all;

Cbtl = 1.98e-6; Rbtl = 4; L = 1.768e-6;
Cf = Cbtl*2; Rf = Rbtl/2;

C_pi = 1.5e-9; R_pi = 33e3;
tau_i = C_pi*R_pi;
gain = 30;

% Parasitic elements
Rl = 300e-3; Cs = 120e-9;

s = tf('s');

sim('control_loop')

%[b,a] = ss2tf(linsys1.A,linsys1.B,linsys1.C,linsys1.D)
%[b,a] = ss2tf(linsys1.A,linsys1.B,linsys1.C,linsys1.D); linsys1_tf = tf(b,a)

lintf = tf([4.285e12,8.578e16],[1,1.263e5,1.428e11,0]);
lintf_rl = tf([2.142e12,4.289e16],[1,2.959e5,1.643e11,0]);
lintf_cs = tf([4.159e12,8.326e16],[1,1.225e5,1.386e11,0]);
