clear;
clc;

% Global variables

Num_BioNet = 35;
Num_ConfigNet = 500;
Bin_Expr = cell(1,Num_BioNet*Num_ConfigNet);

% For each Biomodel do
for k = 1:Num_BioNet
    NId = k; % Current model
    addpath(sprintf('Biomodels/Model%d',NId));
    % Get reactions and expressions
    Reactions = "";
    Expressions = "";
    fid = fopen('f.m','r');
    o = 0;
    p = 0;
    
    % Treat every line of simulation file to get reactions and expressions
    while ~feof(fid)
        tline = fgetl(fid);
        if tline ~= -1
            str = tline;
            str = strtrim(str);
            Expr = '\=';

            % Check and store the names of reactions
            pattern1 = "reaction_";
            TF1 = startsWith(str,pattern1);

            if TF1 == 1
                % Split the reaction in occurance of '='
                splitStr1 = regexp(str,Expr,'split');
                o = o+1;
                Reactions(o) = splitStr1{1}; % Right part goes to Expression
            end

            pattern2 = "xdot(";
            TF2 = startsWith(str,pattern2);
            if TF2 == 1
                % Split the expression in occurance of '='
                splitStr2 = regexp(str,Expr,'split');
                p = p+1;
                Expressions(p) = splitStr2{2}; % Right part goes to Expression
            end
        end
    end
    fclose(fid);
    
    % Remove anything other than reaction from the reaction list
    Reactions(contains(Reactions,pattern2)) = [];

    % For each configuration model or randomized model do
    for m = 1:Num_ConfigNet
        % Randomize each expressions
        Expr = "";
        for i =1:length(Expressions)
            my_str = Expressions(i);
            my_str = insertBefore(my_str, ')' , ' ');
            my_str = insertBefore(my_str, '(' , ' ');
            % Replace each reaction with a random reaction
            newStr = strrep(my_str, extractBetween(my_str, 'reaction_', ' ', 'Boundaries','inclusive'), Reactions(randi(length(Reactions),1,1)));
            if isempty(newStr)
                Expr(i) = my_str;
            else
                Expr(i) = newStr(1); 
            end
        end
        Bin_Expr((k-1)*Num_ConfigNet+m) = {Expr};
    end
    rmpath(sprintf('Biomodels/Model%d',NId));
end
