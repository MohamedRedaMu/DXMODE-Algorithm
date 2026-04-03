clc;
clear;
close all;

[goalReached, GlobalBest, countFE] = DXMODE_algorithm();

disp('==============================');
disp(['Goal reached: ', num2str(goalReached)]);
disp(['Best cost: ', num2str(GlobalBest.Cost)]);
disp(['Function evaluations: ', num2str(countFE)]);
disp('Best position:');
disp(GlobalBest.Position);
disp('==============================');
