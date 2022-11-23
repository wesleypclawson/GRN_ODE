function [ Bool ] = Recheck_Mem_PharmacoResist( Expressions, US, R  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Fold change
    FdUp = 2;
    FdDn = 0.5;
    Fd = [FdUp,FdDn];
    
    % Get initial states of GRN
    Genes = Initialize_ODE(); % Put initial values in all genes
    ST = [];
    States = Flip_Clamp_Simulate( Expressions, Genes, ST, Fd ); 
    Genes_SS = States(end,:); 
    
   
    [ Genes_B, Genes_SS ] = Test_US_Memory( Expressions, Genes_SS, US, Fd );
    if Genes_B(R) == 0
        N_M = 1;
    else
        N_M = 8;
    end
        
    if N_M ~= M % if memory vanishes
        Bool = 1;
    else
        Bool = 0;
    end
end

