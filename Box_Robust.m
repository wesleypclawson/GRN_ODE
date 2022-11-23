clear;
clc;
close all;

load Result/ODERobustnessNewAlt.mat;
%[x1, x2, x3, x4, x5] = Box_Robust_Prep1( Bin_M_Rob );
[totRobust1] = Box_Robust_Prep1( Bin_M_Rob );


%load Result/RobustDiscDiff14_Alt.mat;
load ResultsRandDiscDiff500_Size10/RobustnessDiscDiffRandAsODE500.mat;
%[y1, y2, y3, y4, y5] = Box_Robust_Prep1( Bin_M_Rob );
[totRobust2] = Box_Robust_Prep1( Bin_M_Rob );

% ttest 2 sample for each memory type
numMem = size(totRobust1,1);
h1 = zeros(numMem,1);
p = zeros(numMem,1);
for i = 1:numMem
    [h1(i),p(i)] = ttest2(totRobust1(i,:),totRobust2(i,:));
end

h = subplot(5,2,1);
bar(totRobust1(1,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,2);
bar(totRobust2(1,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,3);
bar(totRobust1(2,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,4);
bar(totRobust2(2,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,5);
bar(totRobust1(3,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,6);
bar(totRobust2(3,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,7);
bar(totRobust1(4,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,8);
bar(totRobust2(4,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,9);
bar(totRobust1(5,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);

h = subplot(5,2,10);
bar(totRobust2(5,:),'BarWidth', 1);
set(h, 'xTickLabel', 0:5, 'FontSize', 14);


