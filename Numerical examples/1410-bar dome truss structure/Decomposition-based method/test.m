%-----------------------------------------------------------------
% The MATLAB source code of the 1410-bar dome truss structure
% Written by Kiarash Biabani Hamedani - Postdoctoral researcher at Iran University of Science and Technology
% Supervisor: Distinguished professor Ali Kaveh
%-----------------------------------------------------------------

clc;
clear;
close all;
format short e
rng('shuffle')

%% Data of the problem
NUM_S = 30;
NUM_NODE = 13;
THETA = 2*pi/NUM_S;
X_COORDS = [1 3 5 7 9 11 13 1.989 3.978 5.967 7.956 9.945 11.934];
Y_COORDS = [0 0 0 0 0 0 0 0.209 0.418 0.627 0.836 1.0453 1.2543];
Z_COORDS = [4 3.75 3.25 2.75 2 1.25 0 3 2.75 2.25 1.75 1 -0.5];
coorda = [X_COORDS;Y_COORDS];

nodes = [];
elements = [];
Constr = zeros(3,39);

for ns = 1:3
    trans = [cos((ns-1)*THETA) -sin((ns-1)*THETA);sin((ns-1)*THETA) cos((ns-1)*THETA)];
    nodes = [nodes;[trans*coorda;Z_COORDS]'];
end

for ns = 1:2
    firstNodeCur = (ns-1)*NUM_NODE+1;
    firstNodeNext = firstNodeCur+NUM_NODE;
    elements = [elements;[firstNodeCur:firstNodeCur+5;firstNodeCur+1:firstNodeCur+6]'];
    elements = [elements;[firstNodeCur+7:firstNodeCur+11;firstNodeCur+8:firstNodeCur+12]'];
    elements = [elements;[firstNodeCur+1:firstNodeCur+6;firstNodeCur+7:firstNodeCur+12]'];
    elements = [elements;[firstNodeCur+1:firstNodeCur+5;firstNodeCur+8:firstNodeCur+12]'];
    elements = [elements;[firstNodeCur,firstNodeCur+7]];
    elements = [elements;[firstNodeCur:firstNodeCur+5;firstNodeNext:firstNodeNext+5]'];
    elements = [elements;[firstNodeCur+7:firstNodeCur+12;firstNodeNext+7:firstNodeNext+12]'];
    elements = [elements;[firstNodeCur+7:firstNodeCur+12;firstNodeNext:firstNodeNext+5]'];
    elements = [elements;[firstNodeCur+7:firstNodeCur+12;firstNodeNext+1:firstNodeNext+6]'];
end

for j = 1:3
    Constr(:,13*(j-1)+7) = 1;
end

Ur = 1-Constr;
f = find(Ur);

xmin = 0.0001;            % Lower bound of design variables
xmax = 0.01;              % Upper bound of design variables
dimension = 47;           % Problem dimension

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