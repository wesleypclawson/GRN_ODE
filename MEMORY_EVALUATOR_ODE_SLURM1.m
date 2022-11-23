%function [ ] = MEMORY_EVALUATOR_ODE_SLURM1(I,J)
   clear;
   clc;
   % Parameters
        ModelId = 4;        % Enter the network model id 
        simCnt = 5000;     % Length of time administering stimuli
        numStimuli = 2;     % Number of stimuli instance for long relax used as reference
        USscaleUp = 100;    % Upregulation of US above maximum point
        RscaleUp = 2;       % Upregulation of R above maximum point
        totLenmultiRelax = (2*numStimuli+1)*simCnt;

    % Simulate the network multiple times to settle in a steady state
        addpath(sprintf('Biomodels/Model%d',ModelId));
        x0 = Initialize_ODE();
    
    % Simulate the network to get a stable value upon 2 stimuli Onset 
    %(highest number of stimuli onset upon any memory evaluation)
    % preceeded by relax
        [ T, Ref ] = multiRelax(x0,totLenmultiRelax);
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
        [ Bin_US_List, Bin_CS_List ] = Get_R_US_NS_Exhaustive( X1, Ref((simCnt+1):(2*simCnt),:), Fd, RscaleUp, simCnt ); % Function call
    
    % Create network specific variables to store memory evaluation results
        Mem_Prop = zeros(1,8);
        Mem_Type = zeros(1,8);
        Mem_US_R_CS = [];
        US_R_CS_Genes = [];
        
    % Perform memory evaluation
        for i = 1:length(x0)            % Estimate memory in each R
            R = i;
            US_List = Bin_US_List{i};   % Get all US, their activation status (upDn) and activation status of current R
            CS_List = Bin_CS_List{i};   % Get all CS, and their activation status (upDn)
            
            if ~isempty(US_List)
                [ Temp_US_R_CS_Genes ] = Mem_Eval_US_R( X1, Ref, US_List, CS_List, Fd, RscaleUp, simCnt );
                US_R_CS_Genes = [ US_R_CS_Genes; Temp_US_R_CS_Genes];
            end
        end
    % Clean results
        US_R_CS_Genes( ~any(US_R_CS_Genes,2), : ) = [];  % Delete all 0 rows 
        US_R_CS_Genes = unique(US_R_CS_Genes,'rows');
        if ~isempty(US_R_CS_Genes)
            [Mem_Type, Mem_Prop, Mem_US_R_CS] = Clean_Results(US_R_CS_Genes);
        end
    % Remove current folder from path
        rmpath(sprintf('Biomodels/Model%d',ModelId)); % Remove the model from current path
    % Save results
%         fileName = sprintf('Out_%s.mat',int2str(ModelId));
%         cd OUT2;
%         save(fileName,'Mem_Type','US_R_CS_Genes', 'Mem_Prop', 'Mem_US_R_CS');
%end 