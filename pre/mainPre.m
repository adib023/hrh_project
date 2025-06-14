currentFolder = pwd;
if currentFolder ~= "../pre"
    cd ../pre
end

close all
clear all
clc


N1 = 4;

theta =  -15;
pattern = 'hrh';
N2_list =  [8 1 8];


k_list = [1 1 1];
kt_list = [0.01 0.01 0.01];
zeta = 0.1;
TB_disp = [0 2; 0 -2];

build_lattice(theta, pattern,N1, N2_list, k_list, kt_list, zeta, TB_disp)




