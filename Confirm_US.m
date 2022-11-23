function [ Genes_B, Genes_CS ] = Confirm_US( Expressions, Genes, US, Fd )  

    % Generate states flipping UCS  
    [ States ] = Flip_Clamp_Simulate( Expressions, Genes, US, Fd ); 
    Genes_CS = States(end,:);
    
    % Compare and detect if R has flipped (Up/Down regulated)
    Start_Stop = [Genes;Genes_CS];
    Genes_B = Detect_Flipped(Start_Stop,Fd);
 end
 