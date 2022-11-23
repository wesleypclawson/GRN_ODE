function [totRobust] = Box_Robust_Prep1( Bin_M_Rob )
   
    x1 = [];
    x2 = [];
    x3 = [];
    x4 = [];
    x5 = [];
    totRobust = zeros(5,6);
    for i = 1:length(Bin_M_Rob) % for each GRN do
        B = Bin_M_Rob{i}; % Assign robustness information of current GRN
        if ~isempty(B)
            Mem = [1 2 3 4 7];
            for j= 1:length(Mem) % for each memory type do
                X = B(B(:,1)==Mem(j),:); % memory cases with current memory type
                if ~isempty(X)
                    Y = X(:,2);
                    switch j
                        case 1
                            %x1(end+1,1) = round(mean(Y));
                            x1 = [x1;Y];
                        case 2
                            %x2(end+1,1) = round(mean(Y));
                            x2 = [x2;Y];
                        case 3
                            %x3(end+1,1) = round(mean(Y));
                            x3 = [x3;Y];
                        case 4
                            %x4(end+1,1) = round(mean(Y));
                            x4 = [x4;Y];
                        otherwise
                            %x5(end+1,1) = round(mean(Y));
                            x5 = [x5;Y];
                    end
                    Y = [];
                end
            end
        end
    end
    for k = 1:length(Mem)
        switch k
            case 1
                totRobust(k,1) = length(x1);
                temp = x1;
            case 2
                totRobust(k,1) = length(x2);
                temp = x2;
            case 3
                totRobust(k,1) = length(x3);
                temp = x3;
            case 4
                totRobust(k,1) = length(x4);
                temp = x4;
            case 5
                totRobust(k,1) = length(x5);
                temp = x5;
        end
        
        for j = 1:length(Mem)  
            t1 = numel(find(temp==j));
            if ~isempty(t1)
                totRobust(k,j+1) = totRobust(k,j)-t1;
            else
                totRobust(k,j+1) = totRobust(k,j);
            end
        end
        temp = [];
        Ref = totRobust(k,1);
        totRobust(k,1) = 100;
        for j = 2:length(Mem)+1 
            totRobust(k,j) = round((totRobust(k,j)/Ref)*100);
        end
    end
end

