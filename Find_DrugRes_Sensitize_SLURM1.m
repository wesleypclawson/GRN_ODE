%function [] = Find_DrugRes_Sensitize_SLURM1(I, J)
    % This program finds PharmacoResistance and Sensitization in Biomodels
    % UCS Up/Down scaling - 100 multiplied with the maximum/minimum point of the node during normal simulation
    % Parameters
    
    ModelId = 4; % Enter the network model id
    startLenStimuli = 5000; % Length of time administering stimuli
    numStimuli = 6; % Number of stimuli instance for sensitization/habituation determination
    lenInc = 1000; % Incremental period length upon each stimulionset
    USscaleUp = 100; % Upregulation of US above maximum point
    RscaleUp = 1.5; % Upregulation of R above maximum point

    % Determining clamping/un-clamping points of stimuli
    lenStimuli = zeros(1,((numStimuli*2)+1)); % create list of stimuli clamping/unclamping points
    c = startLenStimuli;
    lenStimuli(1) = c;
    for i = 2:((numStimuli*2)+1)
        if rem(i,2) ~= 0
            c = c+lenInc; % Incremental spacing for stimuli application
        end
        lenStimuli(i) = lenStimuli(i-1)+c;
    end

    totLenStimuli = lenStimuli(end)/10; % Total time steps conducting the experiment

    addpath(sprintf('Biomodels/Model%d',ModelId)); % Connect the model to current path

    % Simulate the model without intervention with number matching onset and withdrawal of drug
    % onset
    x0 = Initialize_ODE();

    [T1,X1] = multiRelax(x0,totLenStimuli);

    % Get maximum amplitude of different species/node
    MaxNd = zeros(1,size(X1,2));
    MinNd = zeros(1,size(X1,2));
    for JJ = 1:size(X1,2) % Find maximum/minimum of
        nZCur = nonzeros(X1(:,JJ));
        if ~isempty(nZCur)
            MaxNd(JJ) = max(nZCur);
            MinNd(JJ) = min(nZCur);
        else
            MaxNd(JJ) = 0;
            MinNd(JJ) = 0;
        end
    end

    % UP/Down regulation applied over maximum point
    FdUp = MaxNd.*USscaleUp; % Up-regulating UCS
    FdDn = MinNd./USscaleUp; % Down-regulating UCS
    Fd = [FdUp;FdDn]; % Populating fold-changes for up/down regulation experiments

    % Conduct Multi-Drug experiment stimulating (Upregulating and Down regulating) each node with multiple Stimuli
    [BinT,BinUpX,BinDnX] = multiDrug(x0,(lenStimuli./10), Fd);


    %% Check Sensitization and/or Drug Resistance

    US_R_Reg_Status = []; %Up-Down status of US, US index, R, up-Down status of R upon 6 US onset
    Break_DrugRes = []; % Up-Down status of US, US index, R index, First and last UpDns of R, stimuli breaks resistance, stimuli UpDn
    Break_Sensit = [];
    cntSensit = 0;
    cntDrugRes = 0;
    FailedBreakSensit = 0;
    FailedBreakDrugRes = 0;
    lock1 = 0;
    lock2 = 0;
    for II = 1:2 % For up/Down regulation of UCS
        for JJ = 1:size(X1,2) % for each node as UCS do
            if Fd(II,JJ) ~= 0
                if II == 1
                    X2 = BinUpX{JJ}; % Assign experimental results upon current UCS up-regulation
                else
                    X2 = BinDnX{JJ}; % Assign experimental results upon current UCS down-regulation
                end
                % Create label for stimuli vector
                Lb = NaN(lenStimuli(end)+1,1);
                for t = 1:numStimuli
                    Lb(lenStimuli(2*t-1):(lenStimuli(2*t))) = Fd(II,JJ);
                end

                for KK = 1:size(X2,2) % for each node treated as R do
                    R_On_US_temp = zeros(1,numStimuli);

                    if (KK ~= JJ && range(X2(:,JJ)) ~= 0) && range(X2(:,KK)) ~= 0% If UCS and R are different then proceed
                        if (mean(X2(lenStimuli(1):lenStimuli(2)-1,KK))> RscaleUp*mean(X2(1:lenStimuli(1)-1,KK)))  % if real R up-regulated
                            R_On_US_temp(1) = 1; % R is up-regulated upon first UCS
                            for t = 1:(numStimuli-1) % for each UCS up-regulation counting from 2
                                if (mean(X2(lenStimuli((2*(t+1)-1)):(lenStimuli(2*(t+1))),KK)) > RscaleUp*mean(X2(lenStimuli(2*t-1):(lenStimuli(2*t)),KK)))
                                    R_On_US_temp(t+1) = 1;
                                elseif mean(X2(lenStimuli((2*(t+1)-1)):(lenStimuli(2*(t+1))),KK)) <= (1/RscaleUp)*mean(X2(lenStimuli(2*t-1):(lenStimuli(2*t)),KK))
                                    R_On_US_temp(t+1) = 2;
                                end
                            end
                        elseif (mean(X2(lenStimuli(1):lenStimuli(2)-1,KK))<= (1/RscaleUp)*mean(X2(1:lenStimuli(1)-1,KK))) % if real R down-regulated
                            R_On_US_temp(1) = 2; % R up-regulated upon first UCS
                            for t = 1:(numStimuli-1) % for each UCS up-regulation counting from 2
                                if mean(X2(lenStimuli((2*(t+1)-1)):(lenStimuli(2*(t+1))),KK)) <= (1/RscaleUp)*mean(X2(lenStimuli(2*t-1):(lenStimuli(2*t)),KK))
                                    R_On_US_temp(t+1) = 1;
                                elseif mean(X2(lenStimuli((2*(t+1)-1)):(lenStimuli(2*(t+1))),KK)) > RscaleUp*mean(X2(lenStimuli(2*t-1):(lenStimuli(2*t)),KK))
                                    R_On_US_temp(t+1) = 2;
                                end
                            end
                        end
                        US_R_Reg_Status(end+1,:) = [II, JJ, KK, R_On_US_temp];

                        %% Create and save figures for Drug resistance and sensitization

                        %                 if numel(find(R_On_US_temp)) >= 3
                        %                     FigH = Plot_US_R(T1,X1,X2,Lb1, JJ, KK,lenStimuli);
                        %                     fpath = strcat('Result\MultiDrugResistSensit\',sprintf('Model%d',ModelId),'\');
                        %                     filename = sprintf('Net = %d_Up_DOWN = %d_ST-%d_R-%d',ModelId, II, JJ,KK);
                        %                     saveas(FigH, fullfile(fpath, filename), 'png');
                        %                     pause(2)
                        %                 end
                        %% Count Drug resistace and sensitization
                        temp1 = nonzeros(R_On_US_temp);

                        if length(temp1)>1
                            if (temp1(1) ~= temp1(end)) % If a case of drug resistance
                                cntDrugRes = cntDrugRes+1;
                            end
                            if range(temp1) == 0 % If all are same positive integer
                                cntSensit = cntSensit+1;
                            end
                        end

                        %% Important variables for drug resistance and sensitization
                        simCnt2 = 10000;
                        XX1 = X2(lenStimuli(end-5):lenStimuli(end),:); % Last 2 stimulation of muliti-drug experiment
                        TT1 = (0:(size(XX1,1)-1))*.01;
                        Lb1 = Lb(lenStimuli(end-5):lenStimuli(end)); % Labels of last 2 stimulation of multiDrug experiment
                        lenStimuli1 = (lenStimuli((end-5):(end)));
                        %% Break Drug resistance if any

                        % If the first and last non-zero UpDn status for R in
                        % multi-drug experiment are different
                        if (length(temp1) > 1) && (temp1(1) ~= temp1(end))
                            GenesDrugResis = X2(end,:);
                            for LL = 1:size(X2,2) % For each other genes treated as stimuli do
                                if (LL ~= JJ) && (LL ~= KK) % Check if the gene not equals current UCS or R
                                    for MM = 1:2 % Up and down regulation of new stimulus
                                        if Fd(MM,LL) ~= 0
                                            [ X3 ] = stimulateRelaxSingle( GenesDrugResis, LL, MM,  Fd, simCnt2 ); % Stimulate with current new stimuli 1st time
                                            GenesDrugResis1 = X3(end,:);
                                            [ X4 ] = stimulateRelaxSingle( GenesDrugResis1, [], [],  Fd, simCnt2 ); % Relax
                                            GenesDrugResis2 = X4(end,:);
                                            [ X5 ] = stimulateRelaxSingle( GenesDrugResis2, JJ, II,  Fd, simCnt2 ); % Stimulate with UCS
                                            GenesDrugResis3 = X5(end,:);
                                            [ X6 ] = stimulateRelaxSingle( GenesDrugResis3, [], [],  Fd, simCnt2 ); % Relax
                                            % If R (KK) scaled Up on UCS followed by new Stimuli & last regulation type of R was = 2
                                            if (mean(X5(:,KK)) >= RscaleUp*mean(X4(:,KK))) && (temp1(end) ~= 1)
                                                Break_DrugRes(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,0,0,0,0,1]; % Case of Breaking drug resistance
                                                Resist_Sensit = 1;
                                                % Start plot
                                                XX2 = vertcat(X2((end-simCnt2):end,:),X3,X4); % Stimulations breaking resistance
                                                TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                XX3 = vertcat(X4,X5,X6); % Resistance went
                                                TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                for t = 1:length(TT2)/(simCnt2+1)
                                                    if rem(t,2) == 0
                                                        if t == 4 || t == 8
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                        elseif t == 2 || t == 6
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                        elseif t == 10
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                        end
                                                    end
                                                end
                                                DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3, II, JJ, KK, MM,LL,0,0,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                % End plot

                                            elseif (mean(X5(:,KK)) <= (1/RscaleUp)*mean(X4(:,KK))) && (temp1(end) ~= 2)
                                                Break_DrugRes(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,0,0,0,0,2]; % Case of Breaking drug resistance
                                                Resist_Sensit = 1;
                                                % Start plot
                                                XX2 = vertcat(X2((end-simCnt2):end,:),X3,X4); % Stimulations breaking resistance
                                                TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                XX3 = vertcat(X4,X5,X6); % Resistance went
                                                TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                for t = 1:length(TT2)/(simCnt2+1)
                                                    if rem(t,2) == 0
                                                        if t == 4 || t == 8
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                        elseif t == 2 || t == 6
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                        elseif t == 10
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                        end
                                                    end
                                                end
                                                DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3, II, JJ, KK, MM,LL,0,0,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                % End plot
                                            else
                                                GenesDrugResis4 = X6(end,:);
                                                [ X7 ] = stimulateRelaxSingle( GenesDrugResis4, LL, MM,  Fd, simCnt2 ); % Stimulate with current new stimuli 2nd time
                                                GenesDrugResis5 = X7(end,:);
                                                [ X8 ] = stimulateRelaxSingle( GenesDrugResis5, [], [],  Fd, simCnt2 ); % Relax
                                                GenesDrugResis6 = X8(end,:);
                                                [ X9 ] = stimulateRelaxSingle( GenesDrugResis6, JJ, II,  Fd, simCnt2 ); % Stimulate with UCS
                                                GenesDrugResis7 = X8(end,:);
                                                [ X10 ] = stimulateRelaxSingle( GenesDrugResis7, [], [],  Fd, simCnt2 ); % Relax

                                                if (mean(X9(:,KK)) >= RscaleUp*mean(X8(:,KK))) && (temp1(end) ~= 1)
                                                    Break_DrugRes(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,MM,LL,0,0,1]; % Case of Breaking drug resistance
                                                    Resist_Sensit = 1;
                                                    % Start plot
                                                    XX2 = vertcat(X2((end-simCnt2):end,:),X3,X4,X5,X6,X7,X8); % Stimulations breaking resistance
                                                    TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                    Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                    XX3 = vertcat(X8,X9,X10); % Resistance went
                                                    TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                    Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                    Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                    for t = 1:length(TT2)/(simCnt2+1)
                                                        if rem(t,2) == 0
                                                            if t == 4 || t == 8
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                            elseif t == 2 || t == 6
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                            elseif t == 10
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                            end
                                                        end
                                                    end
                                                    DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,0,0,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                    % End plot

                                                elseif (mean(X9(:,KK)) <= (1/RscaleUp)*mean(X8(:,KK))) && (temp1(end) ~= 2)
                                                    Break_DrugRes(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,MM,LL,0,0,2]; % Case of Breaking drug resistance
                                                    Resist_Sensit = 1;
                                                    % Start plot
                                                    XX2 = vertcat(X2((end-simCnt2):end,:),X3,X4,X5,X6,X7,X8); % Stimulations breaking resistance
                                                    TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                    Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                    XX3 = vertcat(X8,X9,X10); % Resistance went
                                                    TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                    Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                    Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                    for t = 1:length(TT2)/(simCnt2+1)
                                                        if rem(t,2) == 0
                                                            if t == 4 || t == 8
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                            elseif t == 2 || t == 6
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                            elseif t == 10
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                            end
                                                        end
                                                    end
                                                    DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,0,0,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                    % End plot

                                                else % Choose a further new stimuli
                                                    All = setdiff(1:size(X2,2),[JJ,KK,LL]); % All genes exept current R and UCS
                                                    for newST2 = 1:length(All)
                                                        NN = All(newST2);
                                                        for OO = 1:2
                                                            if Fd(OO,NN) ~= 0
                                                                GenesDrugResis8 = X10(end,:);
                                                                [ X11 ] = stimulateRelaxSingle( GenesDrugResis8, NN, OO,  Fd, simCnt2 ); % Stimulate with current further new stimulus
                                                                GenesDrugResis9 = X11(end,:);
                                                                [ X12 ] = stimulateRelaxSingle( GenesDrugResis9, [], [],  Fd, simCnt2 ); % Relax
                                                                GenesDrugResis10 = X12(end,:);
                                                                [ X13 ] = stimulateRelaxSingle( GenesDrugResis10, JJ, II,  Fd, simCnt2 ); % Stimulate with UCS
                                                                GenesDrugResis11 = X13(end,:);
                                                                [ X14 ] = stimulateRelaxSingle( GenesDrugResis11, [], [],  Fd, simCnt2 ); % Relax

                                                                if (mean(X13(:,KK)) >= RscaleUp*mean(X12(:,KK))) && (temp1(end) ~= 1)
                                                                    Break_DrugRes(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,MM,LL,OO,NN,1]; % Case of Breaking drug resistance
                                                                    Resist_Sensit = 1;
                                                                    % Start plot
                                                                    XX2 = vertcat(X2((end-simCnt2):end,:),X3,X4,X5,X6,X7,X8,X9,X10,X11,X12); % Stimulations breaking resistance
                                                                    TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                                    Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                                    XX3 = vertcat(X12,X13,X14); % Resistance went
                                                                    TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                                    Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                                    Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                                    for t = 1:length(TT2)/(simCnt2+1)
                                                                        if rem(t,2) == 0
                                                                            if t == 4 || t == 8
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                                            elseif t == 2 || t == 6
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                                            elseif t == 10
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                                            end
                                                                        end
                                                                    end
                                                                    DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,OO,NN,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                                    % End plot

                                                                elseif (mean(X13(:,KK)) <= (1/RscaleUp)*mean(X12(:,KK))) && (temp1(end) ~= 2)
                                                                    Break_DrugRes(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,MM,LL,OO,NN,2]; % Case of Breaking drug resistance
                                                                    Resist_Sensit = 1;
                                                                    % Start plot
                                                                    XX2 = vertcat(X2((end-simCnt2):end,:),X3,X4,X5,X6,X7,X8,X9,X10,X11,X12); % Stimulations breaking resistance
                                                                    TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                                    Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                                    XX3 = vertcat(X12,X13,X14); % Resistance went
                                                                    TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                                    Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                                    Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                                    for t = 1:length(TT2)/(simCnt2+1)
                                                                        if rem(t,2) == 0
                                                                            if t == 4 || t == 8
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                                            elseif t == 2 || t == 6
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                                            elseif t == 10
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                                            end
                                                                        end
                                                                    end
                                                                    DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,OO,NN,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                                    % End plot

                                                                else
                                                                    FailedBreakDrugRes = FailedBreakDrugRes+1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        %% Break Sensitization instances if any
                        % If the first and last non-zero UpDn status for R in
                        % multi-drug experiment are different
                        if (length(temp1) > 1) && (range(temp1) == 0)
                            GenesSensit = X2(end,:); % Get last states of current Sensitization instance
                            simCnt2 = 10000;
                            for LL = 1:size(X2,2) % For each other genes treated as stimuli do
                                if (LL ~= JJ) && (LL ~= KK) % Check if the gene not equals current UCS or R
                                    for MM = 1:2 % Up and down regulation of new stimulus
                                        if Fd(MM,LL) ~= 0
                                            [ Y3 ] = stimulateRelaxSingle( GenesSensit, LL, MM,  Fd, simCnt2 ); % Stimulate with current new stimuli 1st time
                                            GenesSensit1 = Y3(end,:);
                                            [ Y4 ] = stimulateRelaxSingle( GenesSensit1, [], [],  Fd, simCnt2 ); % Relax
                                            GenesSensit2 = Y4(end,:);
                                            [ Y5 ] = stimulateRelaxSingle( GenesSensit2, JJ, II,  Fd, simCnt2 ); % Stimulate with UCS
                                            GenesSensit3 = Y5(end,:);
                                            [ Y6 ] = stimulateRelaxSingle( GenesSensit3, [], [],  Fd, simCnt2 ); % Relax

                                            % If R (KK) scaled Up on UCS followed by new Stimuli & last regulation type of R was = 2
                                            if (mean(Y5(:,KK)) >= RscaleUp*mean(Y4(:,KK))) && (temp1(end) ~= 1)
                                                Break_Sensit(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,0,0,0,0,1]; % Case of Breaking drug resistance
                                                Resist_Sensit = 2;
                                                % Start plot
                                                XX2 = vertcat(X2((end-simCnt2):end,:),Y3,Y4); % Stimulations breaking resistance
                                                TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                XX3 = vertcat(Y4,Y5,Y6); % Resistance went
                                                TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                for t = 1:length(TT2)/(simCnt2+1)
                                                    if rem(t,2) == 0
                                                        if t == 4 || t == 8
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                        elseif t == 2 || t == 6
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                        elseif t == 10
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                        end
                                                    end
                                                end
                                                DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,0,0,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                % End plot

                                            elseif (mean(Y5(:,KK)) <= (1/RscaleUp)*mean(Y4(:,KK))) && (temp1(end) ~= 2)
                                                Break_Sensit(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,0,0,0,0,2]; % Case of Breaking drug resistance
                                                Resist_Sensit = 2;
                                                % Start plot
                                                XX2 = vertcat(X2((end-simCnt2):end,:),Y3,Y4); % Stimulations breaking resistance
                                                TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                XX3 = vertcat(Y4,Y5,Y6); % Resistance went
                                                TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                for t = 1:length(TT2)/(simCnt2+1)
                                                    if rem(t,2) == 0
                                                        if t == 4 || t == 8
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                        elseif t == 2 || t == 6
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                        elseif t == 10
                                                            Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                        end
                                                    end
                                                end
                                                DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,0,0,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                % End plot

                                            else
                                                GenesSensit4 = Y6(end,:);
                                                [ Y7 ] = stimulateRelaxSingle( GenesSensit4, LL, MM,  Fd, simCnt2 ); % Stimulate with current new stimuli 2nd time
                                                GenesSensit5 = Y7(end,:);
                                                [ Y8 ] = stimulateRelaxSingle( GenesSensit5, [], [],  Fd, simCnt2 ); % Relax
                                                GenesSensit6 = Y8(end,:);
                                                [ Y9 ] = stimulateRelaxSingle( GenesSensit6, JJ, II,  Fd, simCnt2 ); % Stimulate with UCS
                                                GenesSensit7 = Y9(end,:);
                                                [ Y10 ] = stimulateRelaxSingle( GenesSensit7, [], [],  Fd, simCnt2 ); % Relax

                                                if (mean(Y9(:,KK)) >= RscaleUp*mean(Y8(:,KK))) && (temp1(end) ~= 1)
                                                    Break_Sensit(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,MM,LL,0,0,1]; % Case of Breaking drug resistance
                                                    Resist_Sensit = 2;
                                                    % Start plot
                                                    XX2 = vertcat(X2((end-simCnt2):end,:),Y3,Y4,Y5,Y6,Y7,Y8); % Stimulations breaking resistance
                                                    TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                    Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                    XX3 = vertcat(Y8,Y9,Y10); % Resistance went
                                                    TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                    Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                    Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                    for t = 1:length(TT2)/(simCnt2+1)
                                                        if rem(t,2) == 0
                                                            if t == 4 || t == 8
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                            elseif t == 2 || t == 6
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                            elseif t == 10
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                            end
                                                        end
                                                    end
                                                    DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,0,0,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                    % End plot

                                                elseif (mean(Y9(:,KK)) <= (1/RscaleUp)*mean(Y8(:,KK))) && (temp1(end) ~= 2)
                                                    Break_Sensit(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,MM,LL,0,0,2]; % Case of Breaking drug resistance
                                                    Resist_Sensit = 2;
                                                    % Start plot
                                                    XX2 = vertcat(X2((end-simCnt2):end,:),Y3,Y4,Y5,Y6,Y7,Y8); % Stimulations breaking resistance
                                                    TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                    Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                    XX3 = vertcat(Y8,Y9,Y10); % Resistance went
                                                    TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                    Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                    Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                    for t = 1:length(TT2)/(simCnt2+1)
                                                        if rem(t,2) == 0
                                                            if t == 4 || t == 8
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                            elseif t == 2 || t == 6
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                            elseif t == 10
                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                            end
                                                        end
                                                    end
                                                    DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,0,0,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                    % End plot

                                                else % Choose a further new stimuli
                                                    All = setdiff(1:size(X2,2),[JJ,KK,LL]); % All genes exept current R and UCS
                                                    for newST2 = 1:length(All)
                                                        NN = All(newST2);
                                                        for OO = 1:2
                                                            if Fd(OO,NN) ~= 0
                                                                GenesSensit8 = Y10(end,:);
                                                                [ Y11 ] = stimulateRelaxSingle( GenesSensit8, NN, OO,  Fd, simCnt2 ); % Stimulate with current further new stimulus
                                                                GenesSensit9 = Y11(end,:);
                                                                [ Y12 ] = stimulateRelaxSingle( GenesSensit9, [], [],  Fd, simCnt2 ); % Relax
                                                                GenesSensit10 = Y12(end,:);
                                                                [ Y13 ] = stimulateRelaxSingle( GenesSensit10, JJ, II,  Fd, simCnt2 ); % Stimulate with UCS
                                                                GenesSensit11 = Y13(end,:);
                                                                [ Y14 ] = stimulateRelaxSingle( GenesSensit11, [], [],  Fd, simCnt2 ); % Relax

                                                                if (mean(Y13(:,KK)) >= RscaleUp*mean(Y12(:,KK))) && (temp1(end) ~= 1)
                                                                    Break_Sensit(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,MM,LL,OO,NN,1]; % Case of Breaking drug resistance
                                                                    Resist_Sensit = 2;
                                                                    % Start plot
                                                                    XX2 = vertcat(X2((end-simCnt2):end,:),Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12); % Stimulations breaking resistance
                                                                    TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                                    Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                                    XX3 = vertcat(Y12,Y13,Y14); % Resistance went
                                                                    TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                                    Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                                    Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                                    for t = 1:length(TT2)/(simCnt2+1)
                                                                        if rem(t,2) == 0
                                                                            if t == 4 || t == 8
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                                            elseif t == 2 || t == 6
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                                            elseif t == 10
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                                            end
                                                                        end
                                                                    end
                                                                    DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,OO,NN,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                                    % End plot

                                                                elseif (mean(Y13(:,KK)) <= (1/RscaleUp)*mean(Y12(:,KK))) && (temp1(end) ~= 2)
                                                                    Break_Sensit(end+1,:) = [II,JJ,KK,temp1(1),temp1(end),MM,LL,MM,LL,OO,NN,2]; % Case of Breaking drug resistance
                                                                    Resist_Sensit = 2;
                                                                    % Start plot
                                                                    XX2 = vertcat(X2((end-simCnt2):end,:),Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12); % Stimulations breaking resistance
                                                                    TT2 = (0:(size(XX2,1)-1))*.01; % Total time breaking resistance
                                                                    Lb2 = NaN(length(TT2),1); % Stimulation labels breaking resistance
                                                                    XX3 = vertcat(Y12,Y13,Y14); % Resistance went
                                                                    TT3 = (0:(size(XX3,1)-1))*.01; % Total time breaking resistance
                                                                    Lb3 = NaN(length(TT3),1); % Stimulation labels breaking resistance
                                                                    Lb3((simCnt2+1)+1:2*(simCnt2+1)) = Fd(II,JJ);
                                                                    for t = 1:length(TT2)/(simCnt2+1)
                                                                        if rem(t,2) == 0
                                                                            if t == 4 || t == 8
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(II,JJ);
                                                                            elseif t == 2 || t == 6
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(MM,LL);
                                                                            elseif t == 10
                                                                                Lb2(((t-1)*(simCnt2+1))+1:(t*(simCnt2+1))) = Fd(OO,NN);
                                                                            end
                                                                        end
                                                                    end
                                                                    DrawBreakResisSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, KK, MM,LL,OO,NN,lenStimuli1,ModelId,Fd,Resist_Sensit)
                                                                    % End plot

                                                                else
                                                                    FailedBreakSensit = FailedBreakSensit+1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    rmpath(sprintf('Biomodels/Model%d',ModelId)); % Remove the model from current path

    % Important results
        numSp = length(x0);
        PrcntCntSensit = (cntSensit/size(US_R_Reg_Status,1))*100;
        PrcntCntDrugRes = (cntDrugRes/size(US_R_Reg_Status,1))*100;
        PrcntBreakDrugRes = (size(Break_DrugRes,1)/(size(Break_DrugRes,1)+FailedBreakDrugRes))*100;
        PrcntBreakSensit = (size(Break_Sensit,1)/(size(Break_Sensit,1)+FailedBreakSensit))*100;

    %   Save results in files
     %   fileName = sprintf('Out_%s',int2str(I));
      %  cd OUT;
       % save(fileName,'US_R_Reg_Status','numSp','Break_DrugRes','Break_Sensit','PrcntCntSensit','PrcntCntDrugRes','PrcntBreakDrugRes','PrcntBreakSensit');
%end
