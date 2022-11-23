% This program finds the memroy of R (each gene in GRN) in a GRN, in terms of
% association between stimuli (US and CS) for all possible combination of US and R  
%clearvars -except Bin_Mem_Prop Bin_Memory_Type Bin_US_R_CS Bin_US_R_CS_Genes;
clear;
clc;

load Expr/RandAsODE500 Bin_Expr;
Z = length(Bin_Expr);

% Variables
    Bin_US_R_CS_Genes1 = cell(1,Z);
    Bin_US_R_CS1 = cell(1,Z);
    Bin_Memory_Type1 = zeros(Z,8);
    Bin_Mem_Prop1 = zeros(Z,8);
    Bin_Time1 = zeros(1,Z);
% Parameters
    simCnt = 500;     % Length of time administering stimuli
    numStimuli = 2;     % Number of stimuli instance for long relax used as reference
    USscaleUp = 100;    % Upregulation of US above maximum point
    RscaleUp = 2;       % Upregulation of R above maximum point
    totLenmultiRelax = (2*numStimuli+1)*simCnt;

    % Perfor memory evaluation of each network
for i = 1:Z
    tic
    % preceeded by relax
        Expressions = Bin_Expr{i};
        N = size(Expressions,1);
        Genes  = zeros(1,N);
        [ Ref ] = multiRelax(Expressions, Genes, totLenmultiRelax);
        X1 = Ref(1:simCnt,:);
        Genes_SS = (X1(end,:)).';

    % Get maximum amplitude of different species/node
        MaxNd = zeros(1,size(X1,2));
        MinNd = zeros(1,size(X1,2));
        for JJ = 1:size(X1,2)
            MaxNd(JJ) = max(X1(:,JJ));
            MinNd(JJ) = min(X1(:,JJ));
        end

    % UP/Down regulation applied over maximum point
        FdUp = MaxNd.*USscaleUp; % Up-regulating UCS
        FdDn = MinNd./USscaleUp; % Down-regulating UCS
        Fd = [FdUp;FdDn]; % Populating fold-changes for up/down regulation experiment

    % Get all UCS and CS for each R in both up-regulated and down-regulated form
        [ Bin_US_List, Bin_CS_List ] = Get_R_US_NS_Exhaustive( Expressions, X1, Ref, Fd, RscaleUp, simCnt ); % Function call

    % Create network specific variables to store memory evaluation results
            Mem_Prop = zeros(1,8);
            Mem_Type = zeros(1,8);
            Mem_US_R_CS = [];
            US_R_CS_Genes = [];

    % Perform memory evaluation
        for j = 1:length(Genes_SS)            % Estimate memory in each R
            R = j;
            US_List = Bin_US_List{j};   % Get all US, their activation status (upDn) and activation status of current R
            CS_List = Bin_CS_List{j};   % Get all CS, and their activation status (upDn)

            if ~isempty(US_List)
                [ Temp_US_R_CS_Genes ] = Mem_Eval_US_R( Expressions, X1, Ref, US_List, CS_List, Fd, RscaleUp, simCnt );
                US_R_CS_Genes = [ US_R_CS_Genes; Temp_US_R_CS_Genes];
            end
        end
        % Clean results
        US_R_CS_Genes( ~any(US_R_CS_Genes,2), : ) = [];  % Delete all 0 rows
        US_R_CS_Genes = unique(US_R_CS_Genes,'rows');
        if ~isempty(US_R_CS_Genes)
            [Mem_Type, Mem_Prop, Mem_US_R_CS] = Clean_Results(US_R_CS_Genes);
        end
        % Assignments to repositories of all networks
            Bin_Memory_Type1(i,:) = Mem_Type;
            Bin_Mem_Prop1(i,:) = Mem_Prop;
            Bin_Time1(i) = toc;
            Bin_US_R_CS1(i) = {Mem_US_R_CS};
            Bin_US_R_CS_Genes1(i) = {US_R_CS_Genes};
           disp(i);
end