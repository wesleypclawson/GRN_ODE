function [ CS ] = Get_probable_CS( Results, N, R )
%   Given the Results obtained from treating every node as a stimuli
%   This function finds the list of probable CS
j = 1;
CS = [];
for i = 1:N
    A = Results(i,:);
    if (i~=R)  % If current node has less than 1/5th of the total genes as down-stram genes activated
        if ~ismember(A, R)  % Current node is not R and R is not activated as a downstream of current node
            CS(j) = i;
            j = j+1;
        end
    end 
end


