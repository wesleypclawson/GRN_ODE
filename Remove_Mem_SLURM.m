function [] = Remove_Mem_SLURM(I,J)
    load Result/ODE Bin_US_R_CS;
    
    % Prepare to be executed in SLURM    
    if J == 1  % If the network is indexed less than 2000
        NId = I;
        A = Bin_US_R_CS{I}; % Assign result of current GRN
    else
        K = 2000+I;
        A = Bin_US_R_CS{K}; % Assign result of current GRN
        NId = ceil(K/100);
    end
    
    addpath(sprintf('Biomodels/Model%d',NId));

    % Get reactions and expressions
    Reactions = "";
    Expressions = "";
    fid = fopen('Simulate_ODE.m');
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
            Reactions(o) = splitStr1{1}; % Right part goes to Expression
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

    % Remove anything other than reaction from the reaction list
    Reactions(contains(Reactions,pattern2)) = [];


    % Find Robastness
    
    A((A(:,1)==5|A(:,1)==6|A(:,1)==8),:) = []; % Delete Subcategories and no memory from result
    Num = size(A,1); % Number of relevent memory cases in current network
    Robustness = zeros(1,5); % For each of five types of memory
    M_Rob = zeros(Num,2);
    if round(((length(Reactions)*length(Expressions))*10/100))> 15
        cMax = round(((length(Reactions)*length(Expressions))*10/100));
    else
        cMax = 15;
    end
    % For each S-R combination do
    for i = 1:Num
        M = A(i,1); % Memory type
        US = A(i,2);
        R = A(i,3);
        CS = A(i,4);
        Expr = Expressions;

        % Attempts to forget memory
        s = 0;
        for j = 1:cMax
            idx1 = randi(length(Expr),1,1);
            my_str = Expr(idx1);
            my_str = insertBefore(my_str, ')' , ' ');
            my_str = insertBefore(my_str, '(' , ' ');

            % Replace just one reaction with a random reaction
            idx2 = regexp(my_str,'reaction_'); % use regexp to find the number of matches
            % syntax: regexprep(str, expression, replace, n) where n = #match
            newStr = regexprep(my_str, extractBetween(my_str, 'reaction_', ' ', 'Boundaries','inclusive'), Reactions(randi(length(Reactions),1,1)), randi(length(idx2),1,1));

            Expr(idx1) = newStr;
            Bool = Recheck_Mem( Expr, M, US, R, CS );
            if Bool ~= 1
                if M == 7
                    Robustness(end) = Robustness(end)+1;
                else
                    Robustness(M) = Robustness(M)+1;
                end
                s = s+1;
            end
        end
        M_Rob(i,:) = [M,s];
    end
    
    rmpath(sprintf('Biomodels/Model%d',NId));
    
    if J == 1
        fileName = sprintf('Out_%s',int2str(I));
    else
        K = 2000+I;
        fileName = sprintf('Out_%s',int2str(K));
    end
    
    cd OUT;
    save(fileName,'M_Rob','Robustness');
end
