function [ US_R_CS_Genes ] = Test_Pairing_Memory( X1, Ref, US_R_Info, CS_List, Fd, RscaleUp, simCnt )
%   Given Expressions, US, CS and R
%   This function test whether there is a Pairing Memory
        US = US_R_Info(1);
        UpDnUS = US_R_Info(2);
        R = US_R_Info(3);
        Genes_SS = (X1(end,:));
        N = length(Genes_SS);
    % Test Pairing memory by pairing US and each CS in CS List
        C = size(CS_List,1);
        US_R_CS_Genes = zeros(C, N+7);
        for KK = 1:C  % For each CS do
            Temp_CS = CS_List(KK,1);
            UpDnCS =  CS_List(KK,2);
            ST = [ US, Temp_CS ];
            UpDnST = [ UpDnUS, UpDnCS ];

            % Perform pairing
                [ E1 ] = stimulateRelaxSingle( Genes_SS.', ST, UpDnST,  Fd, simCnt );
                Genes_Pairing = (E1(end,:));
                
            % Check regulation status of R through Pairing
               [ UpDnR ] = DetectRegR( X1, E1, Ref, R, RscaleUp );

            % Test if memory
            if UpDnR ~= 0 % If pairing passes

                % Leave/Relax the network for N times 
                    ST = [];
                    UpDnST = [];
                    E2 = stimulateRelaxSingle( Genes_Pairing.', ST, UpDnST,  Fd, simCnt ); % Simulate network without intervention
                    Genes_Leave = E2(end,:);
                    
                % Detect Pairing Memory
                    MemPair  = DetectMem( X1, E1, E2, R, UpDnR );

                if (MemPair ~= 0 ) % If there is pairing memory
                    if US_R_Info(1) ~= CS_List(KK,1)
                        US_R_CS_Genes(KK,:) = [ 2, US_R_Info, CS_List(KK,:), Genes_Leave];
                    end
                else
                    US_R_CS_Genes(KK,:) = [ 8, US_R_Info, CS_List(KK,:), Genes_Leave];
                end
            else
                US_R_CS_Genes(KK,:) = [ 8, US_R_Info, CS_List(KK,:), Genes_Pairing];
            end
        end
    end

