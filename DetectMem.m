function [ Mem ] = DetectMem( X1, E1, E2, R, UpDnR )
    %  This function detects memory
    %  It takes initial, regulating and relaxed simulation results - X, E1
    %  and E2
    if (UpDnR == 1) % If R Up-Regulated
        if mean(E2(:,R)) >= (mean(X1(:,R))+(mean(E1(:,R))-mean(X1(:,R)))/2) % If R retains 50% of its regulated value after relaxation
            Mem = 1;
        else
            Mem = 9999;
        end
    elseif (UpDnR == 2) % If R Down-Regulated
        if mean(E2(:,R)) <= (mean(X1(:,R))-(mean(X1(:,R))-mean(E1(:,R)))/2) % If R retains 50% of its regulated value after relaxation
            Mem = 1;
        else
            Mem = 9999;
        end
    end
end

