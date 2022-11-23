function [ US_R_CS_Genes ] = Mem_Eval_US_R( X1, Ref, US_List, CS_List , Fd, RscaleUp, simCnt )
    % This function performs memory evaluation of each row in US_List 
    % Create categories
        US_R_CS_Genes = [];
        U = size(US_List,1);            % Number of US
        N = size(X1,2);
        
        % Evaluate memory in UCS [UCS, Activation status] and R [R, Regulation status]
            for JJ = 1:U                    % for each row of US List do
                US_R_Info = US_List(JJ,:);  % Current US, Activation status, R, Regulation staus
                
                % Test UCS-based memory
                [ MemUS, Genes_USMem ] = Test_US_Memory( X1, US_R_Info, Fd, simCnt );
                
                if MemUS == 1 % If there is US Memory
                    if ~isempty(CS_List) % If there are CS of current R
                        C = size(CS_List,1);
                        Temp = zeros(C, N+7);
                        for KK = 1:C
                            if US_R_Info(1) ~= CS_List(KK,1)
                                Temp(KK,:) = [ 1, US_R_Info, CS_List(KK,:), Genes_USMem];
                            end
                        end
                        US_R_CS_Genes = Temp;
                    else
                        US_R_CS_Genes = [ 1, US_R_Info, [0 0], Genes_USMem];
                    end
                else % If no UCS Memory found, check for Pairing Memory
                    if ~isempty(CS_List) % If there are CS of current R check for Pairing Memory
                        [ L1_US_R_CS_Genes ] = Test_Pairing_Memory( X1, Ref, US_R_Info, CS_List, Fd, RscaleUp, simCnt );
                        US_R_CS_Genes = L1_US_R_CS_Genes;
                    else %if there is neither UCS nor Pairing Memory
                        US_R_CS_Genes = [ 8, US_R_Info, [0 0], Genes_USMem];
                    end
                end
            end
            RscaleUp = 1.5;
           [ L2_US_R_CS_Genes ] = AL( X1, Ref, US_R_Info, CS_List, Fd, RscaleUp, simCnt );
           US_R_CS_Genes = [US_R_CS_Genes; L2_US_R_CS_Genes];
end


