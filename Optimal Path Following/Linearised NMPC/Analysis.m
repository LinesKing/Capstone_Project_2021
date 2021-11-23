clc;
clear all;
close all;

load('record_q20_r5')
u20_5 = record.uHis;
load('record_q20_r01')
u20_01 = record.uHis;
load('record_q20_r20')
u20_20 = record.uHis;
times = record.times;

subplot(2,1,1)
stairs(times(1:length(times)), u20_5(1,:),'r','LineWidth',1)
hold on 
stairs(times(1:length(times)),  u20_01(1,:),'k','LineWidth',1)
hold on 
stairs(times(1:length(times)),  u20_20(1,:),'b','LineWidth',1)
title('Inputs')
ylabel("\omega [rad/s]")
legend('r=5','r=0.1','r=20')
grid on
% Plot a history
subplot(2,1,2)
stairs(times(1:length(times)), u20_5(2,:),'r','LineWidth',1)
hold on 
stairs(times(1:length(times)),  u20_01(2,:),'k','LineWidth',1)
hold on 
stairs(times(1:length(times)),  u20_20(2,:),'b','LineWidth',1)
grid on
ylabel("a [m/s^2]")
legend('r=5','r=0.1','r=20')
xlabel("Time (s)")
