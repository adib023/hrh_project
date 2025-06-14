clc
clear all
close all

addpath("../pre")

fileName =  "../file/hrh.txt";
lattice =  read_lattice(fileName);
%strip_disp(lattice);
strip_modeshape(lattice,0);