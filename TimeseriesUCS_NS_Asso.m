function [ Bool, Bool_Pair, Genes_Final, E2, E3, UpDnR ] = TimeseriesUCS_NS_Asso( X1, Ref, ST, UpDnST, CSInfo, R, Fd, RscaleUp, simCnt )
%   Giving current UCS and CS given
%   This function checks whether UCS alone can turn candidate CS 
%   to actual CS
        Bool = 999;
        Bool_Pair = 999;
        E2 = [];
        E3 = [];
        curSTLen = 2;
    %% Stimulate US to see whether R is regulated 
        Genes_SS = X1(end,:);
        E1 = stimulateRelaxSingle( Genes_SS.', ST, UpDnST,  Fd, simCnt );
        Genes_Alone_Pair = E1(end,:);
        
    % Plot 
        if length(ST) == curSTLen
        % Plot UCS controls R
            Z1 = stimulateRelaxSingle( Genes_SS.', ST(1), UpDnST(1),  Fd, simCnt );
            Z1 = [X1;Z1];
            seq = [1, 3];
            Plot_Timeseries(Z1, ST(1), R, seq);
        % Plot NS does not control R
            Z2 = stimulateRelaxSingle( Genes_SS.', CSInfo(1), CSInfo(2),  Fd, simCnt );
            Z2 = [X1;Z2];
            seq = [2, 3];
            Plot_Timeseries(Z2, CSInfo(1), R, seq);
        % Plot Training
            Z1 = [X1;E1];
            seq = [1, 2, 3]; % In case of Consoloidation memory
            Plot_Timeseries(Z1, ST, R, seq);
        end
        
    
        
    % Compare and detect if R has Up/Down regulated through stimulated simulation
        RscaleUp = 1.5;
        UpDnR = DetectRegR( X1, E1, Ref((simCnt+1):(2*simCnt),:), R, RscaleUp );

    % Check if {US} alone/ {US, NS} pair can convert NS into CS
    %% Pairing successful?
    if UpDnR ~= 0 % If R being regulated by {US}/{US,NS}
        
        % Pairing / clamping experiment (UCS clampped) passed
            Bool_Pair = 1;
        %% Relax network for sufficient time
            tempST = [];
            tempUpDnST = [];
            E2 = stimulateRelaxSingle( Genes_Alone_Pair.', tempST, tempUpDnST,  Fd, simCnt );
            Genes_Leave = E2(end,:);
         
        %% Confirm if candidate CS  has turned to actual CS
        % First - simulate the network stimulating CS
            E3 = stimulateRelaxSingle( Genes_Leave.', CSInfo(1), CSInfo(2),  Fd, simCnt );
            Genes_CS = E3(end,:);
        % Second - Compare and detect if R has Up/Down regulated through stimulated simulation
            UpDnR = DetectRegR( E2, E3, Ref((3*simCnt+1):(4*simCnt),:), R, RscaleUp );
            
        % If NS has turned into CS
            if UpDnR ~= 0 
                disp('I am here');
                Bool = 1; % Candidate CS has become actual CS (Associative memory)
                Genes_Final = Genes_CS; 
                % Plot 
                if length(ST) == curSTLen
                % Plot UCS controls R
                    
                    Z1 = [E2;E3];
                    seq = [2, 3];
                    Plot_Timeseries(Z1, CSInfo(1), R, seq);
                end
            else
                Bool = 999; % Test fails
                Genes_Final = Genes_CS;
                % For plotting known consolidation memory only
                    tempST = [];
                    tempUpDnST = [];
                    E4 = stimulateRelaxSingle( Genes_Final.', tempST, tempUpDnST,  Fd, simCnt );
                    if length(ST) == curSTLen
                    % Plot UCS controls R
                        Z1 = [E2;E3;E4];
                        seq = [2, 3];
                        Plot_Timeseries(Z1, CSInfo(1), R, seq);
                    end
            end
     else
         Bool = 999;
         Genes_Final = Genes_Alone_Pair;
     end
end