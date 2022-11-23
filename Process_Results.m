function [Mem_Type, Mem_Prop, Mem_US_R_CS] = Process_Results(US_R_CS_Genes)
    % PROCESS_RESULTS takes US_R_CS_Genes
    % And returns Mem_Type, Mem_Prop, Mem_US_R_CS
    A = US_R_CS_Genes;
    Mem_Type = zeros(1,8);
    Mem_Prop = zeros(1,8);
    Z = [];
    if ~isempty(A)
        A = A(:,1:4);
        A = unique(A, 'rows');
        L = size(A,1);
        for i = 1:L % 
            a = A(i,2:4);
            for j = 1:L
                b = A(j,2:4);
                if i ~= j
                    if (a == b) % Compare for duplicate entries 
                        if A(j,1) == 8 % Discard no memory if comes with any other type of memory
                            Z = [Z; A(j,:)];
                        end
                    end
                end
            end
        end
        if ~isempty(Z) 
            A = setdiff(A, Z, 'rows');
        end
        if ~isempty(A)
            for o = 1:8
                Mem_Type(o) = sum(A(:,1) == o);
            end
        end
        Mem_US_R_CS = A;
        sum_except4 = sum(Mem_Type)-Mem_Type(4);
        Mem_Prop =(Mem_Type./sum_except4).*100;
    end
end

