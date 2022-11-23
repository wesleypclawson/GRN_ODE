%function [ MemUS, Genes_Leave ] = Test_US_Memory( X1, US_R_Info, Fd, simCnt )
%   Given the Boolean expressions, US and steady state in Genes   
%   This function Tests if there is US-memory
    
    close all;
    clear;
    clc;
    load Result/ODE_New Bin_US_R_CS; 
    
    % Parameters
        ModelId = 19;       % Enter the network model id 
        Num = 4;           % Memory idx within the network i.e. row number within network's memory instance
        
    % Add model to current path and initialize parameters
        addpath(sprintf('Biomodels/Model%d',ModelId));
        x0 = Initialize_ODE();
        [ Ref, RscaleUp, Fd, simCnt ] = Initialize_Recheck(x0);
        X1 = Ref(1:simCnt,:);
        Genes_SS = (X1(end,:)).';
        
    % Obtain US, R, CS and their Up-Dn status for current memory check
        A = Bin_US_R_CS{ModelId};
        curRow = A(Num,:);
        US = curRow(2);
        UpDnUS = curRow(3);
        R = curRow(4);
        UpDnR = curRow(5);
        CS = curRow(6);
        UpDnCS = curRow(7);
    
    % Stimulate US according to original activation status and simulate the network
        [ E1 ] = stimulateRelaxSingle( Genes_SS.', US, UpDnUS,  Fd, simCnt );
        Genes_US = (E1(end,:));
    
    % Relax the network  
        tempUS = [];
        tempUpDnUS = [];
        [ E2 ] = stimulateRelaxSingle( Genes_US.', tempUS, tempUpDnUS,  Fd, simCnt );
        Genes_Leave = E2(end,:);
        
        Z1 = [X1;E1;E2];
        seq = [1, 3];
        Plot_Timeseries(Z1, US, R, seq);
        
    % Detect US Memory
        MemUS  = DetectMem( X1, E1, E2, R, UpDnR );
        
    rmpath(sprintf('Biomodels/Model%d',ModelId));
 %end
 