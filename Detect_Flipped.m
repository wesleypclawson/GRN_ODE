function [Genes_B] = Detect_Flipped(Start_Stop, Fd)
% This function detects wheather R is flipped
    Sp = size(Start_Stop,2);
    Genes_B = zeros(1,Sp);
    for i = 1:Sp
        if Start_Stop(2,i)/Start_Stop(1,i) >= Fd(1)
            Genes_B(i) = 1; % Upregulated
        elseif Start_Stop(2,i)/Start_Stop(1,i) <= Fd(2)
            Genes_B(i) = -1; % Downregulated
        else
            Genes_B(i) = 0; % No flip
        end
    end
end

