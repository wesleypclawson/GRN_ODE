function [ Ref, RscaleUp, Fd, simCnt ] = Initialize_Recheck(x0)
 % Parameters
        simCnt = 5000;      % Length of time administering stimuli
        numStimuli = 2;     % Number of stimuli instance for long relax used as reference
        USscaleUp = 100;    % Upregulation of US above maximum point
        RscaleUp = 2;       % Upregulation of R above maximum point
        totLenmultiRelax = (2*numStimuli+1)*simCnt;
        
    
    % Simulate the network to get a stable value upon 2 stimuli Onset 
    %(highest number of stimuli onset upon any memory evaluation)
    % preceeded by relax
        [ ~, Ref ] = multiRelax( x0, totLenmultiRelax );
        X1 = Ref(1:simCnt,:);
    
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
end

