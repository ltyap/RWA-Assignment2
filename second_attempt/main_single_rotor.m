%% single rotor calculations
clc
clear all
close all

[CT,CP,CQ] = lifting_line_loop(20,1,0,0);

plot_BEM_vs_LLT;