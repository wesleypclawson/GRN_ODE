function [FigH] = Plot_US_R(T,X1,X2,X3, ST, R,lenStimuli)
    %   UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    FigH = figure('Position', get(0, 'Screensize'));

    h1=subplot(3,2,2);
    %h1.Position = [0.5703    0.7093    0.3347    0.08];
    h1.Position = [0.5703    0.72    0.3347    0.08];
    hold on
    l3=plot(T(lenStimuli(1)+1 :lenStimuli(2)-1),X3(lenStimuli(1)+1 :lenStimuli(2)-1),'LineWidth',5,'Color','y');
    hold on
    plot(T(lenStimuli(3)+1 :lenStimuli(4)-1),X3(lenStimuli(3)+1 :lenStimuli(4)-1),'LineWidth',5,'Color','y')
    hold on
    plot(T(lenStimuli(5)+1 :lenStimuli(6)-1),X3(lenStimuli(5)+1 :lenStimuli(6)-1),'LineWidth',5,'Color','y')
    hold on
    plot(T(lenStimuli(7)+1 :lenStimuli(8)-1),X3(lenStimuli(7)+1 :lenStimuli(8)-1),'LineWidth',5,'Color','y')
    hold on
    plot(T(lenStimuli(9)+1 :lenStimuli(10)-1),X3(lenStimuli(9)+1 :lenStimuli(10)-1),'LineWidth',5,'Color','y')
    hold on
    plot(T(lenStimuli(11)+1 :lenStimuli(12)-1),X3(lenStimuli(11)+1 :lenStimuli(12)-1),'LineWidth',5,'Color','y')
    set(gca,'FontSize',20)
    set(gca,'xtick',[]);
    %yMin = min(X2(:,ST));
    
    yMax = max(X3);
    %ylim([yMin, yMax]);
    yticks([yMax])
    yticklabels({yMax})

    subplot(3,2,3);
    l1 = plot(T,X1(:,ST),'LineWidth',2,'Color','r');
    set(gca,'FontSize',20)
    grid on
    title({'Normal UCS'});
    xlim([0 1000])

    subplot(3,2,4);
    plot(T,X2(:,ST),'LineWidth',2,'Color','r');
    set(gca,'FontSize',20)
    grid on
    title({'Stimulated UCS'});
    %ylim([0 10])
    xlim([0 1000])
    
    subplot(3,2,5);
    l2=plot(T,X1(:,R),'LineWidth',2,'Color','g');
    set(gca,'FontSize',20)
    grid on
    title({'Normal R'});
    xlim([0 1000])
    
    subplot(3,2,6)
    plot(T,X2(:,R),'LineWidth',2,'Color','g');
    set(gca,'FontSize',20)
    grid on
    title({'R upon UCS Onset'}); 
    xlim([0 1000])

    % Construct a Legend with the data from the sub-plots
    hL = legend([l1,l2,l3],{'UCS','R','UCS Stimulation'});
    %Programatically move the Legend
    newPosition = [0.25 0.73 0.1 0.05];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits);


end

