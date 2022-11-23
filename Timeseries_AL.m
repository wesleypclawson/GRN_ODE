%function [ L1_US_R_CS_Genes ] = AL( X1, Ref, US_R_Info, CS_List, Fd, RscaleUp, simCnt )
close all;
clear;
clc;
    load Result/ODE Bin_US_R_CS; 
    
    % Parameters
        ModelId = 15;       % Enter the network model id 
        Num = 145;           % Memory idx within the network i.e. row number within network's memory instance
        
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
    
    
        L1_US_R_CS_Genes = [];
        Temp = [];
        Bool1 = 999;
        
   
 
    %% Check for Transfer memory
    ST = [];
    ST = US;
    UpDnST = UpDnUS;
    [ Bool1, ~, Genes_Final, ~, ~, ~ ] = TimeseriesUCS_NS_Asso( X1, Ref, ST, UpDnST, [CS, UpDnCS], R, Fd, RscaleUp, simCnt );
    
    if Bool1 == 1 % If UCS alone can turn candidate CS to actual CS
        Temp = [ 3, curRow(2:5), curRow(6:7), Genes_Final]; % Assigning memory Type 3 as Transfer Memory
        L1_US_R_CS_Genes = [L1_US_R_CS_Genes; Temp];
        Temp = [];
        
        % Check for Associative Memory
    else % Check US paired with candidate CS can turn candidate CS to actual CS
        ST = [];
        ST = [ US, CS ];
        UpDnST = [ UpDnUS, UpDnCS ];
        Bool2 = 999; % Resetting Boolean variable
        [ Bool2, Bool_Pair, Genes_Final, E2, E3, UpDnR ] = TimeseriesUCS_NS_Asso( X1, Ref, ST, UpDnST, [CS, UpDnCS], R, Fd, RscaleUp, simCnt );
        if Bool2 == 1 % If actual association established between {US,CS}
            Temp = [ 4, curRow(2:5), curRow(6:7), Genes_Final]; % Assigning memory Type 4 as actual Associative Memory
            L1_US_R_CS_Genes = [L1_US_R_CS_Genes; Temp];
            Temp = [];
            
            % Leave the network for N times and let it go into an attract
            ST = [];
            UpDnST = [];
            E4 = stimulateRelaxSingle( Genes_Final.', ST, UpDnST,  Fd, simCnt );
            Genes_Leave = E4(end,:);
            
            % Compare and detect if R has flipped (Up/Down regulated)
            [ MemLRAM ] = DetectMem( E2, E3, E4, R, UpDnR );
            
            %  Check long Recall
            if (MemLRAM ~= 0) % No Memory
                Temp = [5, curRow(2:5), curRow(6:7), Genes_Leave]; % Assigning memory Type 5 as Long Recall
                L1_US_R_CS_Genes = [L1_US_R_CS_Genes; Temp];
                Temp = [];
                
            else % Short recall
                Temp = [6, curRow(2:5), curRow(6:7), Genes_Final]; % Assigning memory Type 6 as Short Recall
                L1_US_R_CS_Genes = [L1_US_R_CS_Genes; Temp];
                Temp = [];
            end
        elseif Bool_Pair == 1 % If Pairing passed but R did not flip with CS
            % Leave the network for N times and let it go into an attractor
            ST = [];
            UpDnST = [];
            E5 = stimulateRelaxSingle( Genes_Final.', ST, UpDnST,  Fd, simCnt );
            Genes_Leave = E5(end,:);
            
            % Compare and detect if R has flipped (Up/Down regulated)
            UpDnR = DetectRegR( E2, E5, Ref((4*simCnt+1):(5*simCnt),:), R, RscaleUp );
            
            % Check for Consolidation Memory
            if UpDnR ~= 0 % R flipped (Up/Down regulated)
                Temp = [7, curRow(2:5), curRow(6:7), Genes_Leave]; % Assigning memory Type 7 as consolidation memory
                L1_US_R_CS_Genes = [L1_US_R_CS_Genes; Temp];
                Temp = [];
                
            else % If no memory
                Temp = [8, curRow(2:5), curRow(6:7), Genes_Final]; % Assigning memory Type 8 as no memory
                L1_US_R_CS_Genes = [L1_US_R_CS_Genes; Temp];
                Temp = [];
            end
        end
    end
rmpath(sprintf('Biomodels/Model%d',ModelId));