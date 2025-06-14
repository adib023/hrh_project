%function mainSolver(start_w,final_w)

currentFolder = pwd;
if currentFolder ~= "../solver"
    cd ../solver
end

clc
clear all
close all
addpath("../pre")

filename = "../file/re_hexa_pbc.txt";
nDofNode = 2;
lattice = read_lattice(filename);
force = force_info;
force.id = lattice.frcID;
%force.w =  0.37227;
force.N =  10;
force.A = 0.001;
startTime = 0;
finalTime = 1000;
dt =  0.01;
bndryLoadStart = 500;
bndryLoadDur = 50;

RKresultFile = 'rk_all';
% resultFileName_stat = 'rk_def_stat.txt';


%  
% defrmdConfig(lattice)
% 
% def_xCoord = read('../result/NRresult/def_xCoord.txt');
% lattice.xCoord = def_xCoord;

%[strpID,strpElem,strpfix] = read('../file/strpIDforDispy.txt');
% % 
%disp_anal_strip(lattice, strpID,strpElem,strpfix)

% eigenSolver(lattice)


% 
%rk_static(lattice,bndryLoadDur,dt,resultFileName_stat)

%count = 0;

v = 0.01:0.01:0.05;
%
parfor loop = 1:length(v)
    parallel_rk(lattice,force,startTime,dt,finalTime,v(loop))
end


% freq_res_anal(lattice)

%end


function parallel_rk(lattice,force,startTime,dt,finalTime,v_pres)
range =  linspace(0.3,0.6,100);
parfor loop =  1:length(range)
    w = range(loop)
    rk4(lattice,force,startTime,dt,finalTime,w,loop,v_pres)
end
end