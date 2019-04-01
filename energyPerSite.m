clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

x = [2, 4, 8, 16, 32, 64]; leng = 1:size(x,2);

yF = [1.08747876, 1.01214349, 1.14117622, 1.23938262, 1.97696221]; %FBC; L = 128

yF = [yF, 1.02918613, 0.935339808, 1.01214576, 1.20282924, 1.25657368, 1.97355664]; %FBC; L = 256
%yP = [0.928037643, 1.14819252, 1.30451047, 1.53615487, 1.85011041]; %PBC; L = 256

yF = [yF, 0.989822507, 0.859253883, 0.926216125, 1.00834835, 1.20077586, 1.25075674]; %FBC; L = 512

export = false;
bol = exist('yP');

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
if exist('yP') == 1
    plot(x, yF(leng),'.k', 'MarkerSize',20)
    plot(x, yP(leng),'.m', 'MarkerSize',20)
else
    leng2 = 1:size(x,2)-1;
    plot(x(leng2), yF(leng2),'.k', 'MarkerSize',20)
    plot(x(leng), yF(leng + max(leng2)),'.m', 'MarkerSize',20)
    plot(x(leng), yF(leng + max(leng2) + size(x,2)),'.b', 'MarkerSize',20)
end

hold off

% Cosmetic plot stuff.
xlabel('$\lambda$ [$\log_2$]')
ylabel('$\langle \mathcal{H}_{\textnormal{int}} \rangle$')
%title('Line profiles')
if exist('yP') == 1
    legend('FBC','PBC','Location','southeast')
else
    legend('$L = 128$','$L = 256$','$L = 512$','Location','northwest')
end
box on

xlim([min(x) - 1*min(x), max(x) + 0.1*max(x)]);
ylim([min(yF) - 0.1*min(yF), max(yF) + 0.1*max(yF)]);
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
%
%xticks(unique([sort(512./x) sort(1024./x)]))
%xticklabels(split(num2str(log2(unique([sort(512./x) sort(1024./x)])))))
xticks(unique([sort(x)]))
xticklabels(split(num2str(log2(unique([sort(x)])))))
yticks([0:0.2:100])
%yticklabels({'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'})

if export ~= true
    set(gcf,'Units','pixels');
    set(gcf,'Position', [0 0 550 400]*1.5)
    set(gcf,'color','w');
    tightfig;
else
    set(gcf,'Units','pixels');
    set(gcf,'Position', [0 0 550 400])
    set(gcf,'color','w');
    tightfig;
    fig = gcf;
    filename = 'energyPerSite';
    
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