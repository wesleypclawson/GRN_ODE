function [ Bool, Bool_Pair, Genes_Final, E2, E3, UpDnR ] = UCS_NS_Asso( X1, Ref, ST, UpDnST, CSInfo, R, Fd, RscaleUp, simCnt )
%   Giving current UCS and CS given
%   This function checks whether UCS alone can turn candidate CS 
%   to actual CS
        Bool = 999;
        Bool_Pair = 999;
        E2 = [];
        E3 = [];
        UpDnR = 0;
    %% Stimulate US to see whether R is regulated 
        Genes_SS = X1(end,:);
        E1 = stimulateRelaxSingle( Genes_SS.', ST, UpDnST,  Fd, simCnt );
        Genes_Alone_Pair = E1(end,:);
    
    % Compare and detect if R has Up/Down regulated through stimulated simulation
        UpDnR = DetectRegR( X1, E1, Ref((simCnt+1):(2*simCnt),:), R, RscaleUp );

    % Check if {US} alone/ {US, NS} pair can convert NS into CS
    %% Pairing successful?
    if UpDnR ~= 0 % If R being regulated by {US}/{US,NS}
        
        % Pairing / clamping experiment (UCS clampped) passed
            Bool_Pair = 1;
        %% Relax network for sufficient time
            ST = [];
            UpDnST = [];
            E2 = stimulateRelaxSingle( Genes_Alone_Pair.', ST, UpDnST,  Fd, simCnt );
            Genes_Leave = E2(end,:);
         
        %% Confirm if candidate CS  has turned to actual CS
        % First - simulate the network stimulating CS
            E3 = stimulateRelaxSingle( Genes_Leave.', CSInfo(1), CSInfo(2),  Fd, simCnt );
            Genes_CS = E3(end,:);
        % Second - Compare and detect if R has Up/Down regulated through stimulated simulation
            UpDnR = DetectRegR( E2, E3, Ref((3*simCnt+1):(4*simCnt),:), R, RscaleUp );
            
        % If NS has turned into CS
            if UpDnR ~= 0 
                Bool = 1; % Candidate CS has become actual CS (Associative memory)
                Genes_Final = Genes_CS;  
            else
                Bool = 999; % Test fails
                Genes_Final = Genes_CS;
            end
     else
         Bool = 999;
         Genes_Final = Genes_Alone_Pair;
     end
end