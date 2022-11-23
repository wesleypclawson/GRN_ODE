%function [ US_R_CS_Genes ] = Test_Pairing_Memory( X1, Ref, US_R_Info, CS_List, Fd, RscaleUp, simCnt )
%   Given Expressions, US, CS and R
%   This function test whether there is a Pairing Memory
        close all;
        clear;
        clc;
        load Result/ODE Bin_US_R_CS; 
    
    % Parameters
        ModelId = 7;       % Enter the network model id 
        Num = 64;           % Memory idx within the network i.e. row number within network's memory instance
        
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
        
        ST = [US, CS];
        UpDnST = [ UpDnUS, UpDnCS ];
        
        % Perform pairing
        [ E1 ] = stimulateRelaxSingle( Genes_SS.', ST, UpDnST,  Fd, simCnt );
        Genes_Pairing = (E1(end,:));
        
        % Check regulation status of R through Pairing
        [ UpDnR ] = DetectRegR( X1, E1, Ref, R, RscaleUp );
        
        % Test if memory
        if UpDnR ~= 0 % If pairing passes
            
            % Leave/Relax the network for N times
            TempST = [];
            TempUpDnST = [];
            E2 = stimulateRelaxSingle( Genes_Pairing.', TempST, TempUpDnST,  Fd, simCnt ); % Simulate network without intervention
            Genes_Leave = E2(end,:);
            
            % Detect Pairing Memory
            MemPair  = DetectMem( X1, E1, E2, R, UpDnR );
            
            if (MemPair ~= 0 ) % If there is pairing memory
                US_R_CS_Genes(:) = [ 2, curRow(2:5), curRow(6:7), Genes_Leave];
                Z1 = [X1;E1;E2];
                seq = [1, 2, 3];
                Plot_Timeseries(Z1, ST, R, seq);
            else
                US_R_CS_Genes(KK,:) = [ 8, curRow(2:5), curRow(6:7), Genes_Leave];
            end
        else
            US_R_CS_Genes(KK,:) = [ 8, US_R_Info, CS_List(KK,:), Genes_Pairing];
        end
    %end

