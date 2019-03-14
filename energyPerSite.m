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

export = true;

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
plot(x,yF,'.k', 'MarkerSize',20)
plot(x,yP,'.m', 'MarkerSize',20)

%plot(x,yF,'.k', 'MarkerSize',20)
%plot(x,yF2,'.m', 'MarkerSize',20)

hold off

% Cosmetic plot stuff.
xlabel('$\lambda$ [$\log_2$]')
ylabel('$\langle \mathcal{H}_{\textnormal{int}} \rangle$')
%title('Line profiles')
legend('FBC','PBC','Location','southeast')
%legend('$L = 1024$','$L = 512$','Location','southeast')
box on

xlim([min(x) - 1*min(x), max(x) + 0.1*max(x)]);
ylim([min([yF yF2]) - 0.1*min([yF yF2]), max([yF yF2]) + 0.1*max([yF yF2])]);
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
%
%xticks(unique([sort(512./x) sort(1024./x)]))
%xticklabels(split(num2str(log2(unique([sort(512./x) sort(1024./x)])))))
xticks(unique([sort(x)]))
xticklabels(split(num2str(log2(unique([sort(x)])))))
yticks([0:0.2:100])
%yticklabels({'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'})

set(gcf,'Units','pixels');
set(gcf,'Position', [0 0 550 400])
%set(gca,'Position', [0 0 1 1])
%pbaspect([1.5 1 1])
tightfig;

if export == true
    fig = gcf;
    filename = 'farts';
    
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig,filename,'-dpdf','-r0')
end

%Rotate ylabel, taking into account its size/centre relation.
% ylh = get(gca,'ylabel');
% gyl = get(ylh);
% ylp = get(ylh, 'Position');
% set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
%tightfig;