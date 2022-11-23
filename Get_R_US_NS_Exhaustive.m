function [ Bin_US_List, Bin_CS_List ] = Get_R_US_NS_Exhaustive( X1, Ref, Fd, RscaleUp, simCnt )
    % This function helps receive the set of UCS and NS for each R
    
    % The function takes stable states of all genes (Genes_SS) and 
    % up/down regulations to be applied (Fd) to each gene as stimuli
    
    % Declare output variables
    Genes_SS = (X1(end,:)).';
    n = length(Genes_SS);
    Bin_US_List = cell(1,n);
    Bin_CS_List = cell(1,n);
    
    % Get US and CS lists for each R
    for R = 1:n % treat each gene as R
        US_List = [];
        CS_List = [];
        for ST = 1:n % Treat each gene as stimulus
            if ST ~= R
                for UpDn = 1:2 % Up (1) or Down (2) stimulated stimulus
                    [ X2 ] = stimulateRelaxSingle( Genes_SS, ST, UpDn,  Fd, simCnt );
                    if (mean(X2(:,R)) >= RscaleUp*mean(X1(:,R)))&& (mean(X2(:,R)) > RscaleUp*mean(Ref(:,R)))% If R Up-Regulated
                        US_List(end+1,:) = [ST,UpDn,R,1];
                    elseif (mean(X2(:,R)) <= (1/RscaleUp)*mean(X1(:,R)))&& (mean(X2(:,R)) < (1/RscaleUp)*mean(Ref(:,R)))% If R Down-Regulated
                        US_List(end+1,:) = [ST,UpDn,R,2];
                    else
                        CS_List(end+1,:) = [ST,UpDn];
                    end
                end
            end
        end
        Bin_US_List(R) = {US_List};
        Bin_CS_List(R) = {CS_List};
    end
end

