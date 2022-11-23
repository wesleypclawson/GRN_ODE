function [BinT,BinX,BinST_Stat] = multiDrugNew(x0, lenStimuli, Fd, Bin_US_List)

    % This function performs the experiments administering multiple stimuli
    % in sequence
   
    % Set total length of conducting experiment
    totLenStimuli = lenStimuli(end);
    BinX = {};
    BinT = {};
    BinST_Stat = [];
    for KK = 1:length(x0) % Treat each node as R
        temp = Bin_US_List{KK};
        if isempty(temp) == 0 % If temp not empty i.e. there are UCS of current R
            US_List = temp(:,1:2);
            numUS = size(US_List,1);
            for JJ = 1:numUS % Do for eacb UCS
                ST = US_List(JJ,1); % Current UCS
                II = US_List(JJ,2); % Up-Dn status of US
                toUp = Fd(II,ST); % Stimulation level
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

                BinX(end+1) = {X};
                BinT(end+1) = {T};
                BinST_Stat(end+1,:) = [II,JJ];
            end
        end
    end
end