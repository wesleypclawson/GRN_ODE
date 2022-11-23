function [Genes] = Flip_Stimuli( Genes, ST, Fd )
    %UNTITLED Summary of this function goes here
    if length(ST) == 1 % Flipping single stimulus
        %x = randi([1,2], 1,1); % Choose up/down regulation
        Genes(ST) = Genes(ST)*Fd(1);
    elseif length(ST) == 2 % In case of pairing (flipping two stimuli together)
        %x = randi([1,2], 1,2); % Choose up/down regulation
        Genes(ST(1)) = Genes(ST(1))*Fd(1);
        Genes(ST(2)) = Genes(ST(2))*Fd(1);
    end
end

