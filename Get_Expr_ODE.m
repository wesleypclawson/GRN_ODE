clear;
clc;
NumNet = 3;
Bin_Expr = cell(1,NumNet);
for j = 1:NumNet
    % Simulate the network multiple times to settle in a steady state
    addpath(sprintf('Biomodels1/Model%d',j));
    Expressions = Get_Expressions();
    Bin_Expr(j) = {Expressions};
    rmpath(sprintf('Biomodels/Model%d',j));
end