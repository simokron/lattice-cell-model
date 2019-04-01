clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% x = 1024./[2 4 8 16 32 64];
% y = [29 26 24 22 20 18]+5;

%x = 512./[2 4 8 16 32 64];
%y = [29 26 24 22 20 18];

% x = 256./[2 4 8 16 32 64];
% y = [29 26 24 22 20 18]-4;

x = 128./[2 4 8 16 32 64];
y = [29 26 24 22 20 18]-10;

f=fit(x',y','power1');
coeffvalues(f)

%f = fit(log(x)',log(y)','poly1');
c = coeffvalues(f);
%fprintf(['k = ' num2str(exp(c(2))) '*exp(' num2str(c(1)) '*lambda)\n']);

% Plotting
%figure
%h1 = axes;
%set(gca,'FontSize',12)
hold on
plot(x,y,'.k')
%plot(exp([-1:0.01:max(x)]),exp(f([-1:0.01:max(x)])),'-m')
plot([-1:0.01:max(x)],f([-1:0.01:max(x)]),'-m')
%plot(MCS,y,'.k')
%plot(xp,f(xp),'-m')
% plot([0 1.4],[42 42], '-.', 'Color', [0 0 0] + 0.5) % y = 25
% plot([0 1.4],[17 17], '-.', 'Color', [0 0 0] + 0.5) % y = 140
% plot([0.27 0.27],[0 160], '-.', 'Color', [0 0 0] + 0.5) % x = 0.27
% plot([0.73 0.73],[0 160], '-.', 'Color', [0 0 0] + 0.5) % x = 0.75
hold off
%set(h1, 'Xdir', 'reverse')
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')

% Cosmetic plot stuff.
xlabel('$L/\lambda$')
ylabel('Scaling factor')
%title('Line profiles')
legend('Data points','Fit','Location','northwest')
box on

%xlim([min(x) - 0.2*min(x), max(x) + 0.3*max(x)]);
%ylim([0 2]);
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
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