clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

x = [4, 8, 16, 32, 64];
yF = [1.05727875, 1.13801157, 1.28465378, 1.49189711, 1.77337289]; %FBC; L = 1024
yP = [1.06241822, 1.14819252, 1.30451047, 1.53615487, 1.85011041]; %PBC; L = 1024

yF2 = [1.21179819, 1.27475214, 1.43546546, 1.71248770, 1.98802304]; %FBC; L = 512

%f=fit(x',y','linear')

%f = fit(log(x)',log(y)','poly1');
%c = coeffvalues(f);
%fprintf(['k = ' num2str(exp(c(2))) '*exp(' num2str(c(1)) '*lambda)\n']);

% Plotting
%figure
%h1 = axes;
set(gca,'FontSize',14)
grid on
hold on
plot(1024./x,yP,'.k', 'MarkerSize',20)
plot(1024./x,yF,'.m', 'MarkerSize',20)
%plot(512./x,yF2,'.m', 'MarkerSize',20)
%plot(exp([-1:0.01:max(x)]),exp(f([-1:0.01:max(x)])),'-m')
%plot([-1:0.01:max(x)],f([-1:0.01:max(x)]),'-m')
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
xlabel('$L/\lambda$ [$\log_2$]')
ylabel('$\langle \mathcal{H}_{\textnormal{int}} \rangle$')
%title('Line profiles')
legend('PBC','FBC','Location','northeast')
%legend('$L = 1024$','$L = 512$','Location','northeast')
box on

%xlim([min(x) - 1*min(x), max(x) + 0.1*max(x)]);
%ylim([1 2]);
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
%
%xticks(unique([sort(512./x) sort(1024./x)]))
%xticklabels(split(num2str(log2(unique([sort(512./x) sort(1024./x)])))))
xticks(unique([sort(1024./x)]))
xticklabels(split(num2str(log2(unique([sort(1024./x)])))))
yticks([0:0.2:100])
%yticklabels({'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'})

%Rotate ylabel, taking into account its size/centre relation.
% ylh = get(gca,'ylabel');
% gyl = get(ylh);
% ylp = get(ylh, 'Position');
% set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
tightfig;