function [FigH] = Plot_US_R_BreakResSensit(XX1,TT1,XX2,TT2,XX3,TT3,Lb1,Lb2,Lb3,II, JJ, R, MM,LL,OO,NN,lenStimuli1)
    
    
    FigH = figure('Position', get(0, 'Screensize'));
    
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
        
        % 4th subplot demonstrating stimulation of UCS in making resis/sensit
        h4 = subplot(3,2+prt,(2+prt)+1);
        LbPos1 = get(h4,'position');
        plot(TT1,XX1(:,JJ),'LineWidth',2,'Color','r');
        set(gca,'FontSize',20)
        grid on
        
    % 6th subplot demonstrating stimulation of UCS in breaking resist/sensit
        h6 = subplot(3,2+prt,2*(2+prt));
        LbPos3 = get(h6,'position');
        plot(TT3,XX3(:,JJ),'LineWidth',2,'Color','r');
        set(gca,'FontSize',20)
        grid on
    
    % 7th subplot demonstrating regulation of R in making resist/sensit
        subplot(3,2+prt,(2*(2+prt))+1);
        plot(TT1,XX1(:,R),'LineWidth',2,'Color','g');
        set(gca,'FontSize',20)
        grid on
        
    % 9th subplot demonstrating regulation of R in breaking resist/sensit
        subplot(3,2+prt,3*(2+prt));
        plot(TT3,XX3(:,R),'LineWidth',2,'Color','g');
        set(gca,'FontSize',20)
        grid on
            
        if   prt == 1  % If 1 instance of trial with new stimulus1 breaks resistance/sensitization
                 
             % 2nd subplot demonstrating UCS-labels of trials breaking drug-res/sensit   
                h2 = subplot(3,2+prt,2);
                h2.Position = [0.4108  0.7093    0.2134    0.08];
                hold on
%                 plot(TT2(1 :(simCnt+1)),Lb2(1 :(simCnt+1)),'Color','w');
%                 hold on
                plot(TT2((simCnt+1)+1 :2*(simCnt+1)),Lb2((simCnt+1)+1 :2*(simCnt+1)),'LineWidth',5,'Color',[91, 207, 244] / 255);
                hold on
%                 plot(TT2((2*(simCnt+1)+1):end),Lb2((2*(simCnt+1)+1):end),'Color','w');
%                 hold on
                set(gca,'FontSize',20)
                set(gca,'xtick',[]); 
                yMax = max(Lb2);    
                yticks([yMax]);
                yticklabels({yMax});
            
            % 5th subplot demonstrating trials breaking resist/sensit with new stimulus1
                subplot(3,2+prt,5);
                plot(TT2,XX2(:,LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                grid on
                
            % 8th subplot demonstrating regulation of R in trials breaking resist/sensit
                h2 = subplot(3,2+prt,8);
                LbPos2 = get(h2,'position');
                plot(TT2,XX2(:,R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                grid on
                
            % title({'Stimulated newStimulus1'});
        elseif prt == 3 % If 2 instances of new stimulus1 (1st instance follows UCS stimulation breaks resistance/sensit

            % 2nd subplot demonstrating UCS-labels of trials breaking drug-res/sensit   
                h2 = subplot(3,2+prt,[2,3,4]);
                h2.Position = [0.4108    0.7093    0.2134    0.08];
                hold on
                plot(TT2(1 :(simCnt+1)),Lb2(1 :(simCnt+1)),'Color','w'); % No stimulation
                hold on
                plot(TT2((simCnt+1)+1 :2*(simCnt+1)),Lb2((simCnt+1)+1 :2*(simCnt+1)),'LineWidth',5,'Color',[91, 207, 244] / 255); % 1st stimulation of new stimulus1
                hold on
                plot(TT2((2*(simCnt+1)+1):3*(simCnt+1)),Lb2((2*(simCnt+1)+1):3*(simCnt+1)),'Color','w'); % No stimulation
                hold on
                plot(TT2(3*(simCnt+1)+1 :4*(simCnt+1)),Lb2(3*(simCnt+1)+1:4*(simCnt+1)),'LineWidth',5,'Color','y'); % Stimulation with UCS
                hold on
                plot(TT2((4*(simCnt+1)+1):5*(simCnt+1)),Lb2((4*(simCnt+1)+1):5*(simCnt+1)),'Color','w'); % No stimulation
                hold on
                plot(TT2(5*(simCnt+1)+1 :6*(simCnt+1)),Lb2(5*(simCnt+1)+1:6*(simCnt+1)),'LineWidth',5,'Color',[91, 207, 244] / 255); % 2nd stimulation of new stimulus1
                hold on
                plot(TT2((6*(simCnt+1)+1):7*(simCnt+1)),Lb2((6*(simCnt+1)+1):7*(simCnt+1)),'Color','w'); % No stimulation
                hold on
                set(gca,'FontSize',20)
                set(gca,'xtick',[]); 
                yMax1 = max(Lb2(1:3*(simCnt+1))); 
                yMax2 = 3*(simCnt+1)+1 :5*(simCnt+1);
                yticks([yMax1,yMax2]);
                yticklabels({yMax1,yMax2});
            
            % 7th subplot demonstrating 1st trial breaking resist/sensit with new stimulus1
                subplot(3,2+prt,7);
                plot(TT2(1:3*(simCnt+1)),XX2((1:3*(simCnt+1)),LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                grid on
                
            % 8th subplot demonstrating stimulation with UCS after 1st
            % trial with new stimulus1 to break resist/sensit
                subplot(3,2+prt,8);
                plot(TT2((2*(simCnt+1)+1):5*(simCnt+1)),XX2(((2*(simCnt+1)+1):5*(simCnt+1)),JJ),'LineWidth',2,'Color','r');
                set(gca,'FontSize',20)
                grid on
                
            % 9th subplot demonstrating 2nd trial breaking resist/sensit with new stimulus1
                subplot(3,2+prt,9);
                plot(TT2((4*(simCnt+1)+1):7*(simCnt+1)),XX2(((4*(simCnt+1)+1):7*(simCnt+1)),LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                grid on
                
            % 12th subplot breaking resist/sensit showing R regulation
                h2 = subplot(3,2+prt,[12,13,14]);
                LbPos2 = get(h2,'position');
                plot(TT2(1 :7*(simCnt+1)),XX2((1 :7*(simCnt+1)),LL),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                grid on
            % title({'Stimulated newStimulus1'});
        elseif prt == 5 % if 2 instances new stimulus1 (follows UCS) plus 1 instance of new stimulus2 breaks resis/sensit 
            
           
            
            % 9th subplot demonstrating 1st trial breaking resist/sensit with new stimulus1
                subplot(3,2+prt,9);
                plot(TT2(1:3*(simCnt+1)),XX2((1:3*(simCnt+1)),LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                grid on
                
            % 10th subplot demonstrating stimulation with UCS after 1st
            % trial with new stimulus1 to break resist/sensit
                subplot(3,2+prt,10);
                plot(TT2((2*(simCnt+1)+1):5*(simCnt+1)),XX2(((2*(simCnt+1)+1):5*(simCnt+1)),JJ),'LineWidth',2,'Color','r');
                set(gca,'FontSize',20)
                grid on
                
            % 11th subplot demonstrating 2nd trial breaking resist/sensit with new stimulus1
                subplot(3,2+prt,11);
                plot(TT2((4*(simCnt+1)+1):7*(simCnt+1)),XX2(((4*(simCnt+1)+1):7*(simCnt+1)),LL),'LineWidth',2,'Color','b');
                set(gca,'FontSize',20)
                grid on
                
            % 12th subplot demonstrating stimulation with UCS after 2nd
            % trial with new stimulus1 to break resist/sensit
                subplot(3,2+prt,12);
                plot(TT2((6*(simCnt+1)+1):9*(simCnt+1)),XX2(((6*(simCnt+1)+1):9*(simCnt+1)),JJ),'LineWidth',2,'Color','r');
                set(gca,'FontSize',20)
                grid on
                
            % 13th subplot demonstrating trial breaking resist/sensit with new stimulus2
                subplot(3,2+prt,13);
                plot(TT2((8*(simCnt+1)+1):11*(simCnt+1)),XX2(((8*(simCnt+1)+1):11*(simCnt+1)),NN),'LineWidth',2,'Color','m');
                set(gca,'FontSize',20)
                grid on    
                
            % 16th subplot breaking resist/sensit showing R regulation
                h2 = subplot(3,2+prt,[16,17,18,19,20]);
                LbPos2 = get(h2,'position');
                plot(TT2(1 :7*(simCnt+1)),XX2((1 :7*(simCnt+1)),R),'LineWidth',2,'Color','g');
                set(gca,'FontSize',20)
                grid on
                
                 % 2nd subplot demonstrating UCS-labels of trials breaking drug-res/sensit   
                h2 = subplot(3,2+prt,[2,3,4,5,6]);
                h2.Position = [LbPos2(1)    0.7093    LbPos2(3)    0.08];
                hold on
%                 plot(TT2(1 :(simCnt+1)),Lb2(1 :(simCnt+1)),'Color','w'); % No stimulation
%                 hold on
                plot(TT2((simCnt+1)+1 :2*(simCnt+1)),Lb2((simCnt+1)+1 :2*(simCnt+1)),'LineWidth',5,'Color',[91, 207, 244] / 255); % 1st stimulation of new stimulus1
                hold on
%                 plot(TT2((2*(simCnt+1)+1):3*(simCnt+1)),Lb2((2*(simCnt+1)+1):3*(simCnt+1)),'Color','w'); % No stimulation
%                 hold on
                plot(TT2(3*(simCnt+1)+1 :4*(simCnt+1)),Lb2(3*(simCnt+1)+1:4*(simCnt+1)),'LineWidth',5,'Color','y'); % Stimulation with UCS
                hold on
%                 plot(TT2((4*(simCnt+1)+1):5*(simCnt+1)),Lb2((4*(simCnt+1)+1):5*(simCnt+1)),'Color','w'); % No stimulation
%                 hold on
                plot(TT2(5*(simCnt+1)+1 :6*(simCnt+1)),Lb2(5*(simCnt+1)+1:6*(simCnt+1)),'LineWidth',5,'Color',[91, 207, 244] / 255); % 2nd stimulation of new stimulus1
                hold on
%                 plot(TT2((6*(simCnt+1)+1):7*(simCnt+1)),Lb2((6*(simCnt+1)+1):7*(simCnt+1)),'Color','w'); % No stimulation
%                 hold on
                plot(TT2(7*(simCnt+1)+1 :8*(simCnt+1)),Lb2(7*(simCnt+1)+1:8*(simCnt+1)),'LineWidth',5,'Color','y'); % Stimulation with UCS
                hold on
%                 plot(TT2((8*(simCnt+1)+1):9*(simCnt+1)),Lb2((8*(simCnt+1)+1):9*(simCnt+1)),'Color','w'); % No stimulation
%                 hold on
                plot(TT2(9*(simCnt+1)+1 :10*(simCnt+1)),Lb2(9*(simCnt+1)+1:10*(simCnt+1)),'LineWidth',5,'Color','m'); % Stimulation of new stimulus2
                hold on
%                 plot(TT2((10*(simCnt+1)+1):11*(simCnt+1)),Lb2((10*(simCnt+1)+1):11*(simCnt+1)),'Color','w'); % No stimulation
%                 hold on
                set(gca,'FontSize',20)
                set(gca,'xtick',[]); 
                yMax1 = max(Lb2(1:3*(simCnt+1))); 
                yMax2 = max(Lb2(3*(simCnt+1)+1 :5*(simCnt+1)));
                yMax3 = max(Lb2(8*(simCnt+1)+1 :11*(simCnt+1)));
                A = [yMax1,yMax2,yMax3];
                A = sort(A);
                yticks(A);
                yticklabels(A);
            
        else
            disp('invalid selection');
        end
    
    % 1st subplot demonstrating stimulus-labels of making drug-res/sensit
        h1=subplot(3,2+prt,1);
        h1.Position = [0.1300    0.7093    0.2134    0.08];
        hold on
%         plot(TT1(1 :lenStimuli(1)),Lb1(1 :lenStimuli(1)),'Color','w');        
%         hold on
        plot(TT1(lenStimuli(1)+1 :lenStimuli(2)),Lb1(lenStimuli(1)+1 :lenStimuli(2)),'LineWidth',5,'Color','y');
        hold on
%         plot(TT1(lenStimuli(3)+1 :lenStimuli(4)),Lb1(lenStimuli(3)+1 :lenStimuli(4)),'LineWidth',5,'Color','y')
%         hold on
%         plot(TT1(lenStimuli(4)+1 :lenStimuli(5)),Lb1(lenStimuli(4)+1 :lenStimuli(5)));
        set(gca,'FontSize',20)
        set(gca,'xtick',[]); 
        yMax = max(Lb1);    
        yticks([yMax]);
        yticklabels({yMax});
        
     % 3rd subplot demonstrating UCS-labels of breaking drug-res/sensit   
        h3 = subplot(3,2+prt,2+prt);
        h3.Position = [0.6916    0.7093    0.2134    0.08];
        hold on
        plot(TT3(1 :(simCnt+1)),Lb3(1 :(simCnt+1)),'Color','w');
        hold on
        plot(TT3((simCnt+1)+1 :2*(simCnt+1)),Lb3((simCnt+1)+1 :2*(simCnt+1)),'LineWidth',5,'Color','y');
        hold on
        plot(TT3((2*(simCnt+1)+1):end),Lb3((2*(simCnt+1)+1):end),'Color','w');
        hold on
        set(gca,'FontSize',20)
        set(gca,'xtick',[]); 
        yMax = max(Lb3);    
        yticks([yMax]);
        yticklabels({yMax});
        
    
        
%     % Construct a Legend with the data from the sub-plots
%     hL = legend([l1,l2,l3],{'UCS','R','UCS Stimulation'});
%     %Programatically move the Legend
%     newPosition = [0.25 0.73 0.1 0.05];
%     newUnits = 'normalized';
%     set(hL,'Position', newPosition,'Units', newUnits);

end

