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

cT = [0 0 25];
currentTime = cT(1)*60*60 + cT(2)*60 + cT(3);
numIters = 2^18;
cutoffConc = 0.0;

a = dir([directory '*.dat']);
b = numel(a);

fprintf(['Number of frames:              ' num2str(b)])

for n = 1:b
    frame = importdata([directory 'frame-' num2str(n) '.dat']);
    c0(n) = 1 - nnz(frame)/numel(frame);
    MCS(n) = numIters*(n-1)/(size(frame,1)*size(frame,2));
    if mod(n, 10) == 0 || n == 2
        if n + 1 <= b
            m = n;
        else
            m = n - 1;
        end
    end
end

%% Test
clc;
clf;
fprintf(['Number of frames:              ' num2str(b)])

f = fit(c0',MCS','exp1');

percCompActual = 100*(cutoffConc/c0(end));
percComp = 100*MCS(end)/f(cutoffConc);
percCompCalc = 100*MCS(m+1)/f(cutoffConc);
timeLeft = datestr((currentTime/percCompCalc)*(100-percCompCalc)/(24*60*60), 'HH:MM:SS');

fprintf(['\nActual percentage complete:    ' num2str(round(percCompActual,0)) '%%\n'])

if percCompCalc <= 100
    fprintf(['\nEstimated percentage complete:  ' num2str(round(percComp,0)) '%%'])
    fprintf(['\nEstimated time left:           ' timeLeft '\n'])
end

% Plotting
%figure
h1 = axes;
%set(gca,'FontSize',12)
hold on
plot(c0,MCS,'.k')
plot([cutoffConc:0.01:max(c0)],f([cutoffConc:0.01:max(c0)]),'-m')
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
legend('Data points','Fit','Location','northwest')
box on

xlim([cutoffConc, max(c0)]);
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