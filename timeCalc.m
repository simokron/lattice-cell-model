clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

directory = ['frames' '/'];
frame = importdata([directory 'frame-' num2str(1) '.dat']);

numIters = 2^23;

a = dir([directory '*.dat']);
b = numel(a)

for n = 1:b
    frame = importdata([directory 'frame-' num2str(n) '.dat']);
    c0(n) = 1 - nnz(frame)/numel(frame);
    MCS(n) = numIters*(n-1)/(size(frame,1)*size(frame,2));
end

f = fit(c0',MCS','linear')

100*MCS(end)/f(0.1)

% Plotting
%figure
h1 = axes;
%set(gca,'FontSize',12)
hold on
plot(c0,MCS,'.k')
plot(f)
%plot(MCS,y,'.k')
%plot(xp,f(xp),'-m')
% plot([0 1.4],[42 42], '-.', 'Color', [0 0 0] + 0.5) % y = 25
% plot([0 1.4],[17 17], '-.', 'Color', [0 0 0] + 0.5) % y = 140
% plot([0.27 0.27],[0 160], '-.', 'Color', [0 0 0] + 0.5) % x = 0.27
% plot([0.73 0.73],[0 160], '-.', 'Color', [0 0 0] + 0.5) % x = 0.75
hold off
set(h1, 'Xdir', 'reverse')

% Cosmetic plot stuff.
xlabel('Concentration')
ylabel('MCS')
%title('Line profiles')
legend('Data points','Exponential fit','Location','northwest')
box on

xlim([0.1, max(c0)]);
%ylim([0, 30]);
% 
%xticks([29.66,50,65,80])
%xticklabels({'30','50','65','80'})
% yticks([0, 23, 143])
%yticklabels({'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'})

%Rotate ylabel, taking into account its size/centre relation.
% ylh = get(gca,'ylabel');
% gyl = get(ylh);
% ylp = get(ylh, 'Position');
% set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
%tightfig;