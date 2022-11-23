function [FigH] = Plot_US_R_BreakResSensit1(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, R, MM,LL,OO,NN,lenStimuli1,Fd,Resist_Sensit,ModelId)
    
%% Initial setting    
    
    FigH = figure('Position', get(0, 'Screensize'));
    if Resist_Sensit == 1
        titleText = 'Breaking Drug Resistance';
    else
        titleText = 'Breaking Sensitization';
    end
    subTitleText = sprintf('  UCS-Status = %d-%d, R = %d, NewST1-Status = %d-%d, NewST2-Status = %d-%d',JJ,II,R,LL,MM,NN,OO);

    % Resetting indices - last 2 stimulation instance of UCS developing drug resis 
        lenStimuli = zeros(1,length(lenStimuli1)-1);
        s = 0;
        for i = 2:length(lenStimuli1)
            lenStimuli(i-1) = s+(lenStimuli1(i)-lenStimuli1(i-1));
            s = lenStimuli(i-1);
        end
        
    % Stimulation duration simCnt+1 during trial breaking resistance    
        simCnt = 10000;
        
    % Determining number of trials
        sz = length(TT2)/(simCnt+1); % How many stimulations for breaking drug resistance/sensitization
        if sz == 3
            prt = 1;
        elseif sz == 7
            prt = 3;
        elseif sz == 11
            prt = 5;
        else
            disp('Invalid selection');
        end
%% Developing Resistance/growth of response
            % 4th subplot demonstrating stimulation of UCS in making resis/sensit
                h1 = subplot(2,2+prt,1);
                plot(TT1,XX1(:,JJ),'LineWidth',2,'Color','r');
                set(gca,'FontSize',20)
                xlim([0 500])
                xtk1 = get(h1, 'xtick');
                xtkEnd1 = xtk1(end);
                xticklabels(xtk1)
                grid on

            % 7th subplot demonstrating regulation of R in making resist/sensit
                subplot(2,2+prt,(2+prt)+1);
                plot(TT1,XX1(:,R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xlim([0 500])
                xticklabels(xtk1)
                grid on        

        
%% Trials breaking resistance
        if   prt == 1  % If 1 instance of trial with new stimulus1 breaks resistance/sensitization
            
            
            % 5th subplot demonstrating trials breaking resist/sensit with new stimulus1
                prt1 = subplot(2,2+prt,2);
                plot(TT2,XX2(:,LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                xtk2 = get(prt1, 'xtick');
                xtk2 = xtk2+xtkEnd1;
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on
                
            % 8th subplot demonstrating regulation of R in trials breaking resist/sensit
                subplot(2,2+prt,(2+prt)+2);               
                plot(TT2,XX2(:,R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on 
                
                 
                
        elseif prt == 3 % If 2 instances of new stimulus1 (1st instance follows UCS stimulation breaks resistance/sensit
           
            % 7th Subplot demonstrating 1st trial breaking resist/sensit with newST1
                prt1 = subplot(2,2+prt,2);
                plot(TT2(1:3*(simCnt+1)),XX2((1:3*(simCnt+1)),LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                xtk2 = get(prt1, 'xtick');
                xtk2 = xtk2+xtkEnd1;
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on                
                
           % 12th subplot breaking resist/sensit showing R regulation
                subplot(2,2+prt,(2+prt)+2);
                plot(TT2(1:3*(simCnt+1)),XX2((1:3*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on
              
                              
            % 8th subplot demonstrating stimulation with UCS after 1st
            % trial with newST1 to break resist/sensit
                prt2 = subplot(2,2+prt,3);
                plot(TT2((2*(simCnt+1)+1):5*(simCnt+1)),XX2(((2*(simCnt+1)+1):5*(simCnt+1)),JJ),'LineWidth',2,'Color','r');
                set(gca,'FontSize',20)
                xlim([0 300])
                xtk2 = get(prt2, 'xtick');
                xtk2 = xtk2+(xtkEnd2-xtk2(1));
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on
                
            % 13th subplot breaking resist/sensit showing R regulation
                subplot(2,2+prt,(2+prt)+3);
                plot(TT2((2*(simCnt+1)+1):5*(simCnt+1)),XX2(((2*(simCnt+1)+1):5*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on 
                
                
            % 9th subplot demonstrating 2nd trial breaking resist/sensit with newST1
                prt3 = subplot(2,2+prt,4);
                plot(TT2((4*(simCnt+1)+1):7*(simCnt+1)),XX2(((4*(simCnt+1)+1):7*(simCnt+1)),LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                xlim([0 300])
                xtk2 = get(prt3, 'xtick');
                xtk2 = xtk2+(xtkEnd2-xtk2(1));
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on
                
            % 14th subplot breaking resist/sensit showing R regulation
                subplot(2,2+prt,(2+prt)+4);
                plot(TT2((4*(simCnt+1)+1):7*(simCnt+1)),XX2(((4*(simCnt+1)+1):7*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on   
                                
        elseif prt == 5 % if 2 instances new stimulus1 (follows UCS) plus 1 instance of new stimulus2 breaks resis/sensit            
            
            % 9th Subplot demonstrating 1st trial breaking resist/sensit with newST1
                prt1 = subplot(2,2+prt,2);
                plot(TT2(1:3*(simCnt+1)),XX2((1:3*(simCnt+1)),LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                xtk2 = get(prt1, 'xtick');
                xtk2 = xtk2+xtkEnd1;
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on                
                
           % 16th subplot breaking resist/sensit showing R regulation
                subplot(2,2+prt,(2+prt)+2);
                plot(TT2(1:3*(simCnt+1)),XX2((1:3*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on
                              
            % 10th subplot demonstrating stimulation with UCS after 1st
            % trial with newST1 to break resist/sensit
                prt2 = subplot(2,2+prt,3);
                plot(TT2((2*(simCnt+1)+1):5*(simCnt+1)),XX2(((2*(simCnt+1)+1):5*(simCnt+1)),JJ),'LineWidth',2,'Color','r');
                set(gca,'FontSize',20)
                %xlim([0 300])
                xtk2 = get(prt2, 'xtick');
                xtk2 = xtk2+(xtkEnd2-xtk2(1));
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on
                
            % 17th subplot breaking resist/sensit showing R regulation
                subplot(2,2+prt,(2+prt)+3);
                plot(TT2((2*(simCnt+1)+1):5*(simCnt+1)),XX2(((2*(simCnt+1)+1):5*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on 
                
                
            % 11th subplot demonstrating 2nd trial breaking resist/sensit with newST1
                prt3 = subplot(2,2+prt,4);
                plot(TT2((4*(simCnt+1)+1):7*(simCnt+1)),XX2(((4*(simCnt+1)+1):7*(simCnt+1)),LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                %xlim([0 300])
                xtk2 = get(prt3, 'xtick');
                %xtk2 = xtk2+xtkEnd2;
                xtk2 = xtk2+(xtkEnd2-xtk2(1));
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on
                
            % 18th subplot breaking resist/sensit showing R regulation
                subplot(2,2+prt,(2+prt)+4);
                plot(TT2((4*(simCnt+1)+1):7*(simCnt+1)),XX2(((4*(simCnt+1)+1):7*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on   
                                
                
           % 12th subplot demonstrating 2nd trial breaking resist/sensit with newST1
                prt4 = subplot(2,2+prt,5);
                plot(TT2((6*(simCnt+1)+1):9*(simCnt+1)),XX2(((6*(simCnt+1)+1):9*(simCnt+1)),JJ),'LineWidth',2,'Color','r');
                set(gca,'FontSize',20)
                %xlim([0 300])
                xtk2 = get(prt4, 'xtick');
                %xtk2 = xtk2+xtkEnd2;
                xtk2 = xtk2+(xtkEnd2-xtk2(1));
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on
                
            % 18th subplot breaking resist/sensit showing R regulation
                subplot(2,2+prt,(2+prt)+5);
                plot(TT2((6*(simCnt+1)+1):9*(simCnt+1)),XX2(((6*(simCnt+1)+1):9*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on   
                                
            
           % 13th subplot demonstrating 2nd trial breaking resist/sensit with newST1
                prt5 = subplot(2,2+prt,6);
                plot(TT2((8*(simCnt+1)+1):11*(simCnt+1)),XX2(((8*(simCnt+1)+1):11*(simCnt+1)),NN),'LineWidth',2,'Color','m');
                set(gca,'FontSize',20)
                %xlim([0 300])
                xtk2 = get(prt5, 'xtick');
                %xtk2 = xtk2+xtkEnd2;
                xtk2 = xtk2+(xtkEnd2-xtk2(1));
                xtk2(1) = xtk2(1)+1;
                xticklabels(xtk2)
                grid on
                
            % 19th subplot breaking resist/sensit showing R regulation
                subplot(2,2+prt,(2+prt)+6);
                plot(TT2((8*(simCnt+1)+1):11*(simCnt+1)),XX2(((8*(simCnt+1)+1):11*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                xticklabels(xtk2)
                xtkEnd2 = xtk2(end);
                grid on   
                                
        else
            disp('invalid selection');
        end
 %% Resistance/growth of response overrided        
            % 6th subplot demonstrating stimulation of UCS in breaking resist/sensit
                h6 = subplot(2,2+prt,(2+prt));
                plot(TT3,XX3(:,JJ),'LineWidth',2,'Color','r');
                set(gca,'FontSize',20)
                %xlim([0 300])
                xtk3 = get(h6, 'xtick');
                xtk3 = xtk3+xtkEnd2;
                xtk3(1) = xtk3(1)+1;
                xticklabels(xtk3)
                grid on

            % 9th subplot demonstrating regulation of R in breaking resist/sensit
                subplot(2,2+prt,2*(2+prt));
                plot(TT3,XX3(:,R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                %xlim([0 300])
                xticklabels(xtk3)
                grid on  
        
    sgtitle(strcat(titleText,subTitleText))
        
    fpath = strcat('BreakResSensit\',sprintf('Model%d',ModelId),'\');
    filename = sprintf('Net = %d_Up_DOWN = %d_ST-%d_R-%dDrugRes',ModelId, II, JJ,R);
    %saveas(FigH, fullfile(fpath, filename), 'png');

end

