function [] = Plot_Timeseries(States, ST, R, seq)
% This function makes the timeseries plots of US,R and or CS
    X = 1:size(States,1);
    if length(seq) == 3
        A = States(:,ST(1));
        A = A.';
        B = States(:,ST(2));
        B = B.';
        C = States(:,R);
        C = C.';
        
        % Plot the graph
        figure
%         h = semilogy(X,A, X,B, X,C); % Plot UCS, CS, R
%         set(h(1),'Color','r','LineWidth',4); % Set UCS color
%         set(h(2),'Color',[.5 0 .5],'LineWidth',4); % Set CS color
%         set(h(3),'Color','g','LineWidth',4); % Set R color

        subplot(3,1,1);
        plot(X,A,'Color','r','LineWidth',4);
        subplot(3,1,2);
        plot(X,B,'Color',[.5 0 .5],'LineWidth',4);
        subplot(3,1,3);
        plot(X,C,'Color','g','LineWidth',4);
    else
        A = States(:,ST);
        A = A.';
        C = States(:,R);
        C = C.';
        D = States(:,1);
        D = D.';
        % Plot the graph
        figure
%         h = semilogy(X,A, X,C); % Plot UCS/CS and R
%         if seq(1) == 2 % If current stimuli is CS
%             set(h(1),'Color',[.5 0 .5],'LineWidth',4);
%         else
%             set(h(1),'Color','r','LineWidth',4);
%         end
%         set(h(2),'Color','g','LineWidth',4);    
        subplot(3,1,1);
        if seq(1) == 1
            plot(X,A,'Color','r','LineWidth',4);
        else
            plot(X,A,'Color',[.5 0 .5],'LineWidth',4);
        end
        subplot(3,1,2);
        plot(X,C,'Color','g','LineWidth',4);
        subplot(3,1,3);
        plot(X,D,'Color','y','LineWidth',4);
     end   

%     xticks([0 2000 4000 6000 8000 10000]);
set(findall(gcf,'-property','FontSize'),'FontSize',20)
%     set(h,'linewidth',8);
%     set(gca, 'linewidth',4,'FontSize',60)
%     xlabel('Time')
%     ylabel('Expressions')
end


