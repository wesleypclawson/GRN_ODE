function [BinT,BinUpX,BinDnX] = multiDrug(x0, lenStimuli, Fd)

    % This function performs the experiments administering multiple stimuli
    % in sequence
   
    % Set total length of conducting experiment
    totLenStimuli = lenStimuli(end);
    BinUpX = {};
    BinDnX = {};
    BinT = {};
    for II = 1:2 % Apply Up/Down regulation - 1 for up and 2 for down regulation
        for JJ = 1:length(x0) % Treat each node as UCS
            ST = JJ;
            toUp = Fd(II,JJ);
            X = zeros((totLenStimuli*10)+1,length(x0));
            T = zeros((totLenStimuli*10)+1,1);
            t0 = 0;
            h = 0.01;
            x1 = x0;
            k = 1;
            % Do for each ODE solver
            for i = 1:totLenStimuli
                tf = t0+h*10;
                tspan = t0:h:tf;
                [t,x]=ode23tb(@f,tspan,x1);
                x1 = x(end,:).';
                if (i == lenStimuli(1)) || (i == lenStimuli(3)) || (i == lenStimuli(5)) || (i == lenStimuli(7))|| (i == lenStimuli(9)) || (i == lenStimuli(11))
                    x1(ST) = toUp; % UpReg Fixed
                    toClamp = x1(ST);
                end
                if (i>lenStimuli(1) && i<=lenStimuli(2)) || (i>lenStimuli(3) && i<=lenStimuli(4)) || (i>lenStimuli(5) && i<=lenStimuli(6)) || (i>lenStimuli(7) && i<=lenStimuli(8)) || (i>lenStimuli(9) && i<=lenStimuli(10)) || (i>lenStimuli(11) && i<=lenStimuli(12))
                    x1(ST) = toClamp;
                end
                X(k:k+10,:) = x;
                T(k:k+10,:) = t;
                
                k = k+10;
                t0 = tf;
            end
            
            if II == 1
                BinUpX(end+1) = {X};
                BinT(end+1) = {T};
            else
                BinDnX(end+1) = {X};
            end
        end
    end
end