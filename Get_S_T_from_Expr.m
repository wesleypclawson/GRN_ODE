%function [S,T] = Get_S_T_from_Expr(Expressions)
% Given the Boolean equations 
%   this function produce the source target node pairs of the underlying
%   graph
clear
clc;
    load forST Bin_Expr Bin_Reac_name Bin_Reac;
    Z = length(Bin_Expr);
    %Z= 6;
    Bin_ST = cell(1,Z);
    Bin_Degree = cell(1,Z);
    Connectivity = zeros(Z,5);
    for j = 1:Z  % for each network do
        S = [];
        T = [];
        Expressions = Bin_Expr{j};
        Reac_name = Bin_Reac_name{j};
        Reactions = Bin_Reac{j};
        N = length(Expressions);
        startpat = 'x(';
        endpat = ')';
        for i = 1:N % for each expression do
            A = Expressions(i);  % Get an expression
            idx1 = unique(str2double(extractBetween(A,startpat,endpat))); % Assign all s if any x(s) present in current expression
            idx2 = []; 
            for k = 1:length(Reac_name) % Search for each reaction name if exists in current expression
                pat = Reac_name(k); % Current reaction name
                % If current expression have the reaction
                TF = contains(A,pat); 
                if TF == 1 % Current reaction if exist in expression
                    idx2(end+1) = k; % Include index of current rection from reaction name list
                end
            end
            idx2 = unique(idx2);
            S2 = Reactions(idx2); % Brings all reactions present in current expression
            if ~isempty(S2)
                for l = 1:length(S2)
                    idx3 = (unique(str2double(extractBetween(S2(l),startpat,endpat)))); % Obtain predecessor for node i from the expression
                    if iscolumn(idx3) == 1
                        idx3 = idx3.';
                    end
                    idx1 = [idx1,idx3];
                end
            end
            idx1 = unique(idx1);
            if iscolumn(S) == 1
                S = S.';
            end
            if iscolumn(idx1) == 1
                idx1 = idx1.';
            end
            S = [S,idx1]; % All regulators of the current node
            temp = [];
            temp = ones(1,length(idx1))*i; % Current expression id is the source
            T = [T,temp];
            idx1 = [];
        end
        S = S.';
        T = T.';
        G = graph(S,T);
        D = degree(G);
        Bin_Degree(j) = {D};
        S_T = [S,T];
        S_T = unique(S_T, 'rows');
        Connectivity(j,:) = [j,N,N^2,size(S_T,1),(size(S_T,1)/N^2)*100];
        Bin_ST(j) = {S_T};
    end
    AvgConnec = mean(Connectivity(:,5));
%end

