function [ Expr ] = Expr_2_comp_ODE(X,Y)
%  Given two variables 
%  This function generates Boolean Expression
    R = randsrc(1,1,[1 2; 0.87 0.13]); %AND/OR in 80% probability, XOR in 20% probability
    r1 = randi([0  1], 1,1); % State of first input (complemented or not)
    r2 = randi([0  1], 1,1); % State of second input
    r3 = randi([0  1], 1,1); % Whether ANDing or ORing
    r4 = randi([0  1], 1,1); % State of final outcome
    if R == 1
        if r1 == 0 && r2 == 0
            if r3 == 1
                if r4 == 1
                    Expr = convertCharsToStrings(sprintf('(~(%s) & ~(%s))', X,Y));
                else
                    Expr = convertCharsToStrings(sprintf('(~(~(%s) & ~(%s)))', X,Y));
                end
            else
                if r4 == 1
                    Expr = convertCharsToStrings(sprintf('(~(%s) | (%s))', X,Y));
                else
                    Expr = convertCharsToStrings(sprintf('(~(~(%s) | (%s)))', X,Y));
                end
            end

        elseif r1 == 0 && r2 == 1
            if r3 == 1
                if r4 == 1
                    Expr = convertCharsToStrings(sprintf('(~(%s) & (%s))', X,Y));
                else
                    Expr = convertCharsToStrings(sprintf('(~(~(%s) & (%s)))', X,Y));
                end
            else
                if r4 == 1
                    Expr = convertCharsToStrings(sprintf('(~(%s) | (%s))', X,Y));
                else
                    Expr = convertCharsToStrings(sprintf('(~(~(%s) | (%s)))', X,Y));
                end
            end

        elseif r1 == 1 && r2 == 0
            if r3 == 1
                if r4 == 1
                    Expr = convertCharsToStrings(sprintf('((%s) & ~(%s))', X,Y));
                else
                    Expr = convertCharsToStrings(sprintf('(~((%s) & ~(%s)))', X,Y));
                end
            else
                if r4 == 1
                    Expr = convertCharsToStrings(sprintf('((%s) | ~(%s))', X,Y));
                else
                    Expr = convertCharsToStrings(sprintf('(~((%s) | ~(%s)))', X,Y));
                end
            end

        else
            if r3 == 1
                if r4 == 1
                    Expr = convertCharsToStrings(sprintf('((%s) & (%s))', X,Y));
                else
                    Expr = convertCharsToStrings(sprintf('(~((%s) & (%s)))', X,Y));
                end
            else
                if r4 == 1
                    Expr = convertCharsToStrings(sprintf('((%s) | (%s))', X,Y));
                else
                    Expr = convertCharsToStrings(sprintf('(~((%s) | (%s)))', X,Y));
                end
            end
        end
    else
        if r1 == 0 && r2 == 0
            if r4 == 1
                Expr = convertCharsToStrings(sprintf('(xor(~(%s),~(%s)))', X,Y));
            else
                Expr = convertCharsToStrings(sprintf('(~(xor(~(%s),~(%s))))', X,Y));
            end
        elseif r1 == 0 && r2 == 1
            if r4 == 1
                Expr = convertCharsToStrings(sprintf('(xor(~(%s),(%s)))', X,Y));
            else
                Expr = convertCharsToStrings(sprintf('(~(xor(~(%s),~(%s))))', X,Y));
            end
        elseif r1 == 1 && r2 == 0
            if r4 == 1
                Expr = convertCharsToStrings(sprintf('(xor((%s),~(%s)))', X,Y));
            else
                Expr = convertCharsToStrings(sprintf('(~(xor((%s),~(%s))))', X,Y));
            end
        else
            if r4 == 1
                Expr = convertCharsToStrings(sprintf('(xor((%s),(%s)))', X,Y));
            else
                Expr = convertCharsToStrings(sprintf('(~(xor((%s),(%s))))', X,Y));
            end
        end
    end
end

