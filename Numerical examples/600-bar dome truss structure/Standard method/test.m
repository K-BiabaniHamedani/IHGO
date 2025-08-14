%-----------------------------------------------------------------
% The MATLAB source code of the 600-bar dome truss structure
% Written by Kiarash Biabani Hamedani - Postdoctoral researcher at Iran University of Science and Technology
% Supervisor: Distinguished professor Ali Kaveh
%-----------------------------------------------------------------

clc;
clear;
close all;
format short e
rng('shuffle')

%% Data of the problem
NUM_S = 24;
NUM_NODE = 9;
THETA = 2*pi/NUM_S;
RADIUS = [1 1 3 5 7 9 11 13 14];
Z_COORDS = [7 7.5 7.25 6.75 6 5 3.5 1.5 0];

nodes = [];
elements = [];
Constr = zeros(3,216);

for ns = 1:NUM_S
    nodes = [nodes;[RADIUS*cos((ns-1)*THETA);RADIUS*sin((ns-1)*THETA);Z_COORDS]'];
end

for ns = 1:NUM_S
    firstNodeCur = (ns-1)*NUM_NODE+1;
    if ns ~= NUM_S
        firstNodeNext = firstNodeCur+NUM_NODE;
    else
        firstNodeNext = 1;
    end
    elements = [elements;[firstNodeCur:firstNodeCur+NUM_NODE-2;firstNodeCur+1:firstNodeCur+NUM_NODE-1]'];
    elements = [elements;[firstNodeCur,firstNodeCur+2]];
    elements = [elements;[firstNodeCur:firstNodeCur+NUM_NODE-2;firstNodeNext:firstNodeNext+NUM_NODE-2]'];
    elements = [elements;[firstNodeCur+2:firstNodeCur+NUM_NODE-1;firstNodeNext+1:firstNodeNext+NUM_NODE-2]'];
    elements = [elements;[firstNodeCur,firstNodeNext+1]];
end

for j = 1:24
    Constr(:,9*j) = 1;
end

Ur = 1-Constr;
f = find(Ur);

xmin = 0.0001;            % Lower bound of design variables
xmax = 0.01;              % Upper bound of design variables
dimension = 25;           % Problem dimension

%% Parameter setting
popsize = 30;             % Population size (20 for GO and 30 for IHGO)
MaxFEs = 30000;           % Maximum number of objective function evaluations
Runs = 20;                % Number of runs

%% Record the best results
B_Designs = zeros(Runs,dimension);
B_Weights = zeros(Runs,1);
B_P_Weights = zeros(Runs,1);
B_Frequencies = zeros(Runs,5);
Weights_History = zeros(MaxFEs,Runs);
P_Weights_History = zeros(MaxFEs,Runs);
NFEs_Runs = zeros(Runs,1);
run_id = 0;
runs = 0;

while run_id < Runs
    
    run_id = run_id + 1;
    runs = runs + 1;
    
    %[gbestX,gbestweight,gbestfitness,gbestfrequencies,gbesthistory_weight,gbesthistory_fitness,FEs,e1,e2,e3,P1,P3] = GO(nodes,elements,f,popsize,dimension,xmin,xmax,MaxFEs,run_id);
    [gbestX,gbestweight,gbestfitness,gbestfrequencies,gbesthistory_weight,gbesthistory_fitness,FEs,e1,e2,e3,P1,P3] = IHGO(nodes,elements,f,popsize,dimension,xmin,xmax,MaxFEs,run_id);
    
    B_Designs(run_id,:) = gbestX;
    B_Weights(run_id,1) = gbestweight;
    B_P_Weights(run_id,1) = gbestfitness;
    B_Frequencies(run_id,:) = gbestfrequencies;
    Weights_History(:,run_id) = gbesthistory_weight;
    P_Weights_History(:,run_id) = gbesthistory_fitness;
    NFEs_Runs(run_id,1) = FEs;
    
    if B_Weights(run_id,1) ~= B_P_Weights(run_id,1)
        run_id = run_id-1;
    end
    
end