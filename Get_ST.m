clear;
clc;

Num_BioNet = 35;
Bin_Expr = cell(1,Num_BioNet);
Bin_Reac = cell(1,Num_BioNet);
Bin_ReacName = cell(1,Num_BioNet);
for k = 1:Num_BioNet
    NId = k; % Current model
    addpath(sprintf('Biomodels/Model%d',NId));
    x0 = Initialize_ODE();
    % Get reactions and expressions
    ReacName = "";
    Reactions = "";
    Expressions = "";
    fid = fopen('f.m');
    o = 0;
    p = 0;
    
    % Treat every line of simulation file to get reactions and expressions
    while ~feof(fid)
        tline = fgetl(fid);
        str = tline;
        str = strtrim(str);
        Expr = '\=';
        
        % Check and store the names of reactions
        pattern1 = "reaction_";
        TF1 = startsWith(str,pattern1);
        
        if TF1 == 1 % If line starts with "reaction_", it is a reaction
            % Split the reaction in occurance of '='
            splitStr1 = regexp(str,Expr,'split');
            o = o+1;
            ReacName(o) = splitStr1{1};
            Reactions(o) = splitStr1{2}; % Right part goes to Expression
        end
        
        % Check and store the names of expression
        pattern2 = "xdot(";
        TF2 = startsWith(str,pattern2);
        if TF2 == 1 % If line starts with "xdot(", it is an expression
            % Split the expression in occurance of '='
            splitStr2 = regexp(str,Expr,'split');
            p = p+1;
            Expressions(p) = splitStr2{2}; % Right part goes to Expression
        end
    end
    fclose(fid);
    rmpath(sprintf('Biomodels/Model%d',NId)); 
    
    % Remove anything other than reaction from the reaction list
    Reactions(contains(Reactions,pattern2)) = [];
    Bin_Expr(k) = {Expressions};
    Bin_Reac(k) = {Reactions};
    Bin_ReacName(k) = {ReacName};
end