function [Expressions] = Get_Reactions()   
    fid = fopen('Simulate_ODE.m');
    Reactions = "";
    i = 0;
    while ~feof(fid)
        tline = fgetl(fid);
        str = tline;
        pattern = "xdot(";
        TF = contains(str,pattern);
        if TF == 1
            Expr = '\=';
            % Split the expression in occurance of '='
            splitStr = regexp(str,Expr,'split');
            i = i+1;
            Expressions(i) = splitStr{2}; % Right part goes to Expression
        end
    end
    fclose(fid);
end