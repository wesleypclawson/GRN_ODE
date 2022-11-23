function [ Genes_Pairing ] = Pairing_ODE( Expressions, US, CS, N, Genes_SS, Sim_Cnt )
%   Given Expressions, US, CS and R
%   This function pair the stimuli (US & CS) and let the GRN memorize this

    % Freeze the network 
    newG = Genes_SS; % Steady state 
 
    % Start pairing
    k1 = newG(US);
    newG(US) = ~k1;    % Flip US
    k2 = newG(CS);
    newG(CS) = ~k2;    % Flip CS
    
    % Simulate expressions of down-stream nodes 
    % Leaving the UCS-CS pair 
    % Peform pairing experiment
    for Count = 1:Sim_Cnt
        Genes = newG;
        for j = 1:N
            if (j ~= US) && (j ~= CS)
                newG(j) = eval(Expressions(j));
            end
        end
    end
    
    Genes_Pairing = newG;
end

