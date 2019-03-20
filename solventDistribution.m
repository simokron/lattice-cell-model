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
    directory = [prefix 'lambda_8-L_512-J_0.0000_0.7500_1.2500-numIters_2-24-initialDist_60_20_20-FBC'];
end

skipFrames = 1;
polDeg = 15;
lateralView = false;

export = false;
exportGIF = false;
pauseTime = 0.2;
timeDep = true;

if exportGIF == true
    fprintf('Exporting frames as GIF...\n')
end

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
        
        F = @(x)polyval(f,x)-0.35; rootExists = false;
        if F(X(1)) < 0 && F(X(end)) > 0
            x0(n) = fzero(F,X([1,end])); rootExists = true;
            MCS(n) = numIters*(n-1)/(size(frame,1)*size(frame,2));
        elseif F(X(end)) < 0 && F(X(1)) > 0
            x0(n) = fzero(F,X([1,end])); rootExists = true;
            MCS(n) = numIters*(n-1)/(size(frame,1)*size(frame,2));
        end
    end
    
    clf;
    % Plotting
    %figure
    h1 = axes;
    set(gca,'FontSize',12)
    hold on
    %plot(X,c2,'ok', 'MarkerSize',5)
    plot(X,c0,'.r', 'MarkerSize',20)
    %plot(X,c1,'.k', 'MarkerSize',20)
    %plot([1:0.01:size(frame,1)/lambda],f([1:0.01:size(frame,1)/lambda]),'-b')
    %plot([1:0.01:size(frame,1)/lambda],g([1:0.01:size(frame,1)/lambda]),'--','Color', [95 95 255]/255)
    if lateralView == false
        plot(X,polyval(f,X),'-k')
        if rootExists == true
            plot([x0(n) x0(n)],[0 0.35], '--', 'Color', [0 0 0] + 0.5) % x = x0
            plot(x0(n),0.35,'d', 'MarkerSize',10, 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
        end
    end
    hold off
    
    % Cosmetic plot stuff.
    if lateralView == false
        xlabel('$i$')
    else
        xlabel('$j$')
    end
    ylabel('Concentrations')
    title(['Concentration of zeros = ' num2str(round(sum(c0)/size(c0,2),2))])
    %legend('$c_{-1}$','$c_{0}$','$c_{+1}$','Location','northwest')
    %legend('Data points','Exponential fit','Linear fit','Location','northwest')
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
    set(gcf,'Units','pixels');
    set(gcf,'Position', [0 0 550 400]*1.5)
    set(gcf,'color','w');
    %set(gca,'Position', [0 0 1 1])
    %pbaspect([1.5 1 1])
    tightfig;
    
    if n == 1
        pause(1);
        if exportGIF == true
            h = gcf;
            filename = [directory '-solventDistribution.gif'];
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',2);
        end
    else
        pause(0.0333); % 0.0167 for 60 FPS, 0.0333 for 30 FPS
        if exportGIF == true
            %set(findall(gca, 'type', 'text'), 'visible', 'off')
            h = gcf;
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',pauseTime);
        end
    end
end
if exportGIF == true
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',5);
end

%%
if timeDep == true && lateralView == false
    frame = importdata([directory '/frame-' num2str(n-skipFrames) '.dat']);
    MCS = MCS(MCS~=0); x0 = x0(x0~=0); %Remove zeros.
    X = 0:0.1:MCS(end);
    f = fit(MCS',x0','power1'); coeffs = coeffvalues(f);
    
    clf;
    % Plotting
    %figure
    h1 = axes;
    set(gca,'FontSize',12)
    hold on
    plot(MCS,x0,'.k', 'MarkerSize',20)
    ws = warning('off','all');  % Turn off warnings.
    plot(X,f(X),'-m')
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
    
    set(gcf,'Units','pixels');
    set(gcf,'Position', [0 0 550 400]*1.5)
    set(gcf,'color','w');
    %set(gca,'Position', [0 0 1 1])
    %pbaspect([1.5 1 1])
    tightfig;
    
    if export == true
        fig = gcf;
        filename = 'timeDependence';
        
        set(fig,'Units','Inches');
        pos = get(fig,'Position');
        set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(fig,filename,'-dpdf','-r0')
    end
end