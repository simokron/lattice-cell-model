clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%prefix = 'squareScaling/';
%prefix = 'PBCvsFBC/';
prefix = 'automatedRun/';
%prefix = '';
if isempty(prefix) == true
    %directory = [prefix 'frames'];
    directory = [prefix 'lambda_4-L_512-J_0.0000_0.7500_1.2500-numIters_2-24-initialDist_60_30_10-FBC'];
else
    directory = [prefix 'lambda_2-L_512-J_0.0000_0.7500_1.2500-numIters_2-29-initialDist_60_20_20-FBC'];
end

skipFrames = 1;
evapFront = 0.25;
polDeg = 15;
lateralView = false;

F = 'pdf'; %pdf or png!
export = true; %Turns on the frame export! For GIF exporting, use exporGIF below. DO NOT USE BOTH!
exportGIF = false;
pauseTime = 0.2; %The time between each frame in the GIF.
timeDep = false; %Shows the time dependence after completion (and exports if export = true).

if exportGIF == true
    fprintf('Exporting frames as gif...\n')
elseif export == true
    fprintf(['Exporting frames as ' F '...\n'])
end

current = 1.0; %relates to exporting - largest concentration to save (usually best keps at 1.0)

ipos = strfind(directory,'lambda_') + strlength("lambda_");
iposLim = strfind(directory,'-L_') - 1;
lambda = str2num(directory(ipos:iposLim));

ipos = strfind(directory,'numIters_2-') + strlength("numIters_2-");
iposLim = strfind(directory,'-initialDist_') - 1;
exponent = str2num(directory(ipos:iposLim)); numIters = 2^exponent;

a = dir([directory '/*.dat']);
b = numel(a);

h = 1;
for n = 1:b-1:b
    frame = importdata([directory '/frame-' num2str(n) '.dat']);
    k = 1;
    for i = 1:lambda:size(frame,1)
        cTemp = 0;
        for x1 = i:1:i+lambda-1
            cTemp = cTemp + size(frame,1) - nnz(frame(x1,:));
        end
        cTemp = cTemp/(lambda^2 * size(frame,1)/lambda);
        c0(k) = cTemp;
        k = k+1;
    end
    theLims(h) = sum(c0)/size(c0,2);
    h = h+1;
end
kLim = (0.8/theLims(1))*(theLims(1)-theLims(end))/b;

for n = 1:skipFrames:b
    
    frame = importdata([directory '/frame-' num2str(n) '.dat']);
    MCSTemp = numIters*(n-1)/(size(frame,1)*size(frame,2));
    
    k = 1;
    if lateralView == false
        for i = 1:lambda:size(frame,1)
            cTemp = 0;cTemp1 = 0;cTemp2 = 0;
            for x1 = i:1:i+lambda-1
                cTemp = cTemp + size(frame,1) - nnz(frame(x1,:));
                cTemp1 = cTemp1 + sum(frame(x1,:) == 1);
                cTemp2 = cTemp2 + sum(frame(x1,:) == -1);
            end
            cTemp = cTemp/(lambda^2 * size(frame,1)/lambda);
            cTemp1 = cTemp1/(lambda^2 * size(frame,1)/lambda);
            cTemp2 = cTemp2/(lambda^2 * size(frame,1)/lambda);
            c0(k) = cTemp;
            c1(k) = cTemp1;
            c2(k) = cTemp2;
            k = k+1;
        end
    else
        for j = 1:lambda:size(frame,2)
            cTemp = 0;cTemp1 = 0;cTemp2 = 0;
            for x2 = j:1:j+lambda-1
                cTemp = cTemp + size(frame,1) - nnz(frame(:,x2));
                cTemp1 = cTemp1 + sum(frame(:,x2) == 1);
                cTemp2 = cTemp2 + sum(frame(:,x2) == -1);
            end
            cTemp = cTemp/(lambda^2 * size(frame,1)/lambda);
            cTemp1 = cTemp1/(lambda^2 * size(frame,1)/lambda);
            cTemp2 = cTemp2/(lambda^2 * size(frame,1)/lambda);
            c0(k) = cTemp;
            c1(k) = cTemp1;
            c2(k) = cTemp2;
            k = k+1;
        end
    end
    
    X = [1:1:size(frame,1)/lambda]';
    
    if lateralView == false
        ws = warning('off','all');  % Turn off warnings
        f = polyfit(X,c0',polDeg);
        warning(ws)  % Turn it back on.
        
        fun = @(x)polyval(f,x)-evapFront; rootExists = false;
        if fun(X(1)) < 0 && fun(X(end-1)) > 0
            x0(n) = fzero(fun,X([1,end-1])); rootExists = true;
            MCS(n) = numIters*(n-1)/(size(frame,1)*size(frame,2));
        elseif fun(X(end-1)) < 0 && fun(X(1)) > 0
            x0(n) = fzero(fun,X([1,end-1])); rootExists = true;
            MCS(n) = numIters*(n-1)/(size(frame,1)*size(frame,2));
        end
    end
    
    clf;
    % Plotting
    %figure
    h1 = axes;
    set(gca,'FontSize',12)
    hold on
    if timeDep == false
        %plot(X,c2,'ok', 'MarkerSize',5)
        plot(X,c0,'.r', 'MarkerSize',20)
        %plot(X,c1,'.k', 'MarkerSize',20)
    else
        plot(X,c0,'.r', 'MarkerSize',20)
    end
    if timeDep == true
        if lateralView == false
            plot(X,polyval(f,X),'-k')
            if rootExists == true
                plot(x0(n),evapFront,'d', 'MarkerSize',10, 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
                plot([x0(n) x0(n)],[0 evapFront], '--', 'Color', [0 0 0] + 0.5) % x = x0
                plot([0 x0(n)],[evapFront evapFront], '--', 'Color', [0 0 0] + 0.5) % y = evapFront
                plot(x0(n),evapFront,'d', 'MarkerSize',10, 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
            end
        end
    end
    hold off
    
    % Cosmetic plot stuff.
    if lateralView == false
        xlabel('$i$')
    else
        xlabel('$j$')
    end
    if export ~= true
        title(['Concentration of zeros = ' num2str(round(sum(c0)/size(c0,2),2)) '; MCS = ' num2str(round(MCSTemp,0))])
    end
    if timeDep == false
        %ylabel('Concentrations')
        ylabel('Solvent concentration')
        %legend('$c_{-1}$','$c_{0}$','$c_{+1}$','Location','northwest')
    else
        ylabel('Solvent concentration')
        legend('$c_{0}$','Polynomial fit','Location','northwest')
    end
    box on
    grid on
    
    xlim([1, size(frame,1)/lambda]);
    %ylim([0 1.0-kLim*n]);
    ylim([0 1.0]);
    xticks([0:size(c0,2)/8:size(c0,2)])
    
    % yticks([0, 23, 143])
    %yticklabels({'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'})
    
    %     %Rotate ylabel, taking into account its size/centre relation.
    %     ylh = get(gca,'ylabel');
    %     gyl = get(ylh);
    %     ylp = get(ylh, 'Position');
    %     set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
    %tightfig;
    
    set(gca,'FontSize',14)
    if export ~= true
        set(gcf,'Units','pixels');
        set(gcf,'Position', [0 0 550 400]*1.5)
        set(gcf,'color','w');
        tightfig;
    elseif export == true
        set(gcf,'Units','pixels');
        set(gcf,'Position', [0 0 550 400])
        set(gcf,'color','w');
        tightfig;
    end
    if export == true
        for k = 1:9
            if c0 <= 0.1
                k = 0.1;
            end
            if round(sum(c0)/size(c0,2),2) == k/10 && k/10 < current
                current = k/10;
                fig = gcf;
                filename = sprintf([directory '_MCS_' num2str(round(MCSTemp,0)) '_c0_0%d-solventDistribution.' F],str2num(strrep(num2str(round(sum(c0)/size(c0,2),2)),'.','')));
                if sum(F == 'png') == 3
                    frame = getframe(fig);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    imwrite(imind,cm,filename,F);
                else
                    set(fig,'Units','Inches');
                    pos = get(fig,'Position');
                    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
                    print(fig,filename,'-dpdf','-r0')
                end
            end
        end
    elseif exportGIF == true
        h = gcf;
        filename = [directory '-solventDistribution.gif'];
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',2);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',pauseTime);
        end
    end
    %pause(1)
    %pause(0.5)
    pause(0.0333); % 0.0167 for 60 FPS, 0.0333 for 30 FPS
    %pause(0.0167);
end
if exportGIF == true
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',5);
end

%%
if timeDep == true && lateralView == false
    frame = importdata([directory '/frame-' num2str(n-skipFrames) '.dat']);
    MCS = MCS(MCS~=0); x0 = x0(x0~=0); %Remove zeros.
    X = 0:0.1:MCS(end);
    fun = fit(MCS',x0','power1'); coeffs = coeffvalues(fun);
    
    clf;
    % Plotting
    %figure
    h1 = axes;
    set(gca,'FontSize',12)
    hold on
    plot(MCS,x0,'.k', 'MarkerSize',20)
    ws = warning('off','all');  % Turn off warnings.
    plot(X,fun(X),'-m')
    warning(ws)  % Turn them back on.
    hold off
    
    s2 = ['Power law fit: $y = ' num2str(round(coeffs(1),2)) '\cdot x^{' num2str(round(coeffs(2),2)) '}$'];
    % Cosmetic plot stuff.
    xlabel('MCS')
    ylabel('$i_{\textnormal{evap}}$')
    %title(['Concentration of zeros = ' num2str(round(sum(c0)/size(c0,2),2))])
    legend('Data points',s2,'Location','northwest')
    %legend('Data points','Exponential fit','Linear fit','Location','northwest')
    box on
    grid on
    
    xlim([0, max(MCS)]);
    ylim([1, size(frame,1)/lambda]);
    yticks([0:size(c0,2)/8:size(c0,2)])
    % yticks([0, 23, 143])
    %yticklabels({'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'})
    
    %     %Rotate ylabel, taking into account its size/centre relation.
    %     ylh = get(gca,'ylabel');
    %     gyl = get(ylh);
    %     ylp = get(ylh, 'Position');
    %     set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
    %tightfig;
    
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
        filename = sprintf([directory 'timeDependence.' F]);
        if sum(F == 'png') == 3
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,filename,F);
        else
            set(fig,'Units','Inches');
            pos = get(fig,'Position');
            set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(fig,filename,'-dpdf','-r0')
        end
    end
end
