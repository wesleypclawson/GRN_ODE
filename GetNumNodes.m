% Get the number of nodes in each network
clear;
clc;
numNodes = 35;
Bin_numNodes = zeros(1,numNodes);

for i = 1:numNodes 
    ModelId = i;
    addpath(sprintf('Biomodels/Model%d',ModelId)); % Connect the model to current path

    % Simulate the model without intervention with number matching onset and withdrawal of drug
    x0 = Initialize_ODE();
    Bin_numNodes(i) = numel(x0);
    addpath(sprintf('Biomodels/Model%d',ModelId)); % Connect the model to current path
end