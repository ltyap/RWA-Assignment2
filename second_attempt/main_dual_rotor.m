%% all dual rotor calculations
clc
clear all
close all

Radius = 50;
D = 2*Radius;
L = [1,2,5,1000]*D; % distance between rotors

for i=1:length(L)
    [CT(i,:),CP(i,:),CQ(i,:)] = lifting_line_loop(20,1,L(i),1)
end