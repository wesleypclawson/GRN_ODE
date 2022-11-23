function [ MemUS, Genes_Leave ] = Test_US_Memory( X1, US_R_Info, Fd, simCnt )
%   Given the Boolean expressions, US and steady state in Genes   
%   This function Tests if there is US-memory
    
    % Stimulate the state of US and Generate states flipping US
    % Clamp stimuli in the same state and simulate
    
        US = US_R_Info(1);
        UpDnUS = US_R_Info(2);
        R = US_R_Info(3);
        UpDnR = US_R_Info(4);
        Genes_SS = X1(end,:);
    
    % Stimulate US according to original activation status and simulate the network
        [ E1 ] = stimulateRelaxSingle( Genes_SS.', US, UpDnUS,  Fd, simCnt );
        Genes_US = (E1(end,:));
    
    % Relax the network  
        US = [];
        UpDnUS = [];
        [ E2 ] = stimulateRelaxSingle( Genes_US.', US, UpDnUS,  Fd, simCnt );
        Genes_Leave = E2(end,:);
    
    % Detect US Memory
        MemUS  = DetectMem( X1, E1, E2, R, UpDnR );
 end
 