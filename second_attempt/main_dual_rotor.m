%% all dual rotor calculations
clc
clear all
close all

Radius = 50;
D = 2*Radius;
N=20;
spacing=1; % 1- uniform, 0-cosine 
sec_rot=1;  % is there a second rotor
a_wake_bem = 0.2602;
%L = [1,2,5,1000]*D; % distance between rotors
L=2*D;
Nrotations = 10; %number of rotations done by the rotor (length of wake)
for i=1:length(L)
    [CT(i,:),CP(i,:),CQ(i,:)] = lifting_line_loop(N,spacing,L(i),sec_rot, a_wake_bem, Nrotations);
end