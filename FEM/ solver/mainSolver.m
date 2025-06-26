clear all
clc

format long


global DIM left_disp right_disp nMatType


%problemType = 1; %%% newton raphson
%problemType = 2; %%% modeshape of deformed lattice
%problemType = 3; %%% modeshape of undeformed lattice
%problemType = 4; %%% freq_response analysis of undeformed lattice
problemType = 5; %%% freq_response analysis of deformed lattice
%problemType = 6; %%% post of freq_response analysis for undeformed lattice
%problemType = 7; %%% post of freq_response analysis for deformed lattice

pres_disp = 0.020;
initial_configuration(pres_disp)

isSparse = 1 ;

frcCoord_and_fctr = [0.8660602627 0 0.005 0 -1 0;
    0.8660602627 0 0 0 -1 0;
    0.8368745646 0 0.005 0 1 0;
    0.8368745646 0 0 0 1 0];


res_anal_node_coord =  [0.8660602627 0 0.005;
                        0.8368745646 0 0.005;
                        0.789044056, 0, 0.005
                        0.724713929, 0, 0.005;
                        0.5317235479, 0, 0.005;
                       ]; 


f_range =  5000:0.05:5200;

% quadratic hexahedral elements:
%meshFile = '../hrh_tet2_fem_long.txt';
meshFile = '../hrh_new_long_fem.txt';


% TODO: make sure that fixed node is not a trailing node.
%fixedNode(1) = 18 ;


[elemType,nNode, nElem, xNode, indx ,matProp,elemMatIndx,pbcNode, left_node, right_node]...
    = preProcess_new(meshFile);


[nDof,eqNum, matStiff,Wmat, forceNode, resNode] =...
    initialize(nNode,xNode,pbcNode,matProp,left_node, right_node, ...
    frcCoord_and_fctr, res_anal_node_coord);





if problemType ==  1

    outputFile = '../Results/output.txt';
    fid = fopen(outputFile,'w');


    [fid,strainList] =  NRsolver(isSparse, elemType,nNode, nElem ,...
        xNode, indx, elemMatIndx, matStiff,...
        eqNum , Wmat, DIM,...
        pbcNode,left_node, right_node, fid);

    fclose(fid);


elseif problemType == 2

    alpha = 0.01;
    UUR =  load('../Results/final_UUR.txt');

    M_global = mass_matrix(elemType,nElem, ...
        xNode, indx, elemMatIndx, eqNum, UUR, Wmat);


    [~, ~, K_global] = ...
        total_assembly(isSparse , elemType,nElem, xNode, indx,...
        elemMatIndx, matStiff, eqNum, Wmat , UUR, ...
        false,true);

    [freq, modeshape] = determineMode(M_global, K_global, 131, 20, Wmat, ...
        'defFreq.txt','defMode.txt');

    writeMode(alpha, UUR,freq,modeshape,'vtkDefMode.txt')


elseif problemType == 3


    UUR =  zeros(nNode,3);
    alpha =  0.001;

    M_global = mass_matrix(elemType,nElem, ...
        xNode, indx, elemMatIndx, eqNum, UUR, Wmat);


    [~, ~, K_global] = ...
        total_assembly(isSparse , elemType,nElem, xNode, indx,...
        elemMatIndx, matStiff, eqNum, Wmat , UUR, ...
        false,true);

    [freq, modeshape] = determineMode(M_global, K_global,1000 , 1000, Wmat, ...
        'freq.txt','mode.txt');


    writeMode(alpha,UUR,freq,modeshape,'vtkMode.txt')


elseif problemType == 4
    %%%% case1 - decay
    %%%% case2 - leaky

    UUR =  zeros(nNode,3);

    M_global = mass_matrix(elemType,nElem, ...
        xNode, indx, elemMatIndx, eqNum, UUR, Wmat);

    [~, ~, K_global] = ...
        total_assembly(isSparse , elemType,nElem, xNode, indx,...
        elemMatIndx, matStiff, eqNum, Wmat , UUR, ...
        false,true);

    freq_response(M_global, K_global,forceNode, f_range, resNode, ...
    Wmat, 'undefFreqList.txt', 'undefRes.txt');

elseif problemType == 5
    UUR =  load('../Results/final_UUR.txt');

    M_global = mass_matrix(elemType,nElem, ...
        xNode, indx, elemMatIndx, eqNum, UUR, Wmat);


    [~, ~, K_global] = ...
        total_assembly(isSparse , elemType,nElem, xNode, indx,...
        elemMatIndx, matStiff, eqNum, Wmat , UUR, ...
        false,true);

    freq_response(M_global, K_global,forceNode, f_range, resNode, ...
    Wmat, 'defFreqList.txt', 'defRes.txt');

elseif problemType == 6

    freq= load("../Results/undefFreqList.txt");
    res = load("../Results/undefRes.txt");
size(freq)
size(res)
    markers = ["g-","r-", "k-","m-", "b-"];
    figure
    for p = 1:size(res,2)
        plot(freq, res(:,p),markers(p), 'linewidth', 2)
        hold on
    end

elseif problemType == 7
    freq= load("../Results/defFreqList.txt");
    res = load("../Results/defRes.txt");

    markers = ["g--","r--", "k--","m--", "b--"];
    figure
    for p = 1:size(res,2)
        plot(freq, res(:,p),markers(p),'linewidth', 2)
        hold on
    end
    
else
    error('the problem type is incorrectly defined')

end



exit
%

