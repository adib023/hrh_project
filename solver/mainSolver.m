%function mainSolver(start_w,final_w)
currentFolder = pwd;
if currentFolder ~= "../solver"
    cd ../solver
end

clc
clear all
close all
addpath("../pre")


w_ext = 0;
filename = "../file/hrh.txt";
nDofNode = 2;
lattice = read_lattice(filename);
A = 0.0001;
startTime = 0;
finalTime = 300;
dt =  0.001;


eigenSolver(lattice) 
%defrmdConfig(lattice)
% 
% def_xCoord = read('../result/NRresult/def_xCoord.txt');
% lattice.xCoord = def_xCoord;

%[strpID,strpElem,strpfix] = read('../file/strpIDforDispy.txt');
% % 
%disp_anal_strip(lattice, strpID,strpElem,strpfix)



% 
%rk_static(lattice,bndryLoadDur,dt,resultFileName_stat)

%count = 0;

% v = [0.025 0.05 0.1 0.15 0.2];
% %
% parfor loop = 1:length(v)
%      rk_static(lattice,200,dt,0.1);
% end

% v_pres = 0.02;
%rk_static(lattice,startTime,finalTime,dt,v_pres);
%parallel_rk(lattice,startTime,dt,finalTime,v_pres,A)


%freq_res_anal(lattice)

%end


%rk4(lattice,startTime,dt,finalTime,w_ext,0.01,1.5)
