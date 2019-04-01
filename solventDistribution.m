clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%prefix = '';
prefix = 'automatedRun/256/';
%prefix = 'debug/';
%prefix = 'J_str/';
%prefix = 'PBCvsFBC/';
%prefix = 'solventDistribution/';
%prefix = 'topView/';

folder = 'lambda_4-L_256-J_0.0000_1.0000_2.0000-numIters_2-22-initialDist_60_20_20-FBC';

directory = [prefix folder];

skipFrames = 1;
evapFront = 0.30;
polDeg = 15;
lateralView = false;
timeDep = true; %Shows the time dependence after completion (and exports if export = true).
logLog = true;

F = 'png'; %pdf or png!
export = false; %Turns on the frame export! For GIF exporting, use exporGIF below. DO NOT USE BOTH!
fontSize = 14; % 14 for 0.5\linewidth; 21 for 0.33\linewidth (for 1:1 scale - try 18 if it's too large)
exportGIF = false;
pauseTime = 0.2; %The time between each frame in the GIF.

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
    set(gca,'FontSize',fontSize)
    hold on
    if timeDep == false
        plot(X,c2,'ok', 'MarkerSize',5)
        plot(X,c0,'.r', 'MarkerSize',20)
        plot(X,c1,'.k', 'MarkerSize',20)
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
        xlabel('$j$')
    else
        xlabel('$i$')
    end
    if export ~= true || sum(F == 'png') == 3
        title(['Concentration of zeros = ' num2str(round(sum(c0)/size(c0,2),2)) '; MCS = ' num2str(round(MCSTemp,0))])
    end
    if timeDep == false
        ylabel('Concentration')
        %legend('$c_{0}$','Location','northwest')
        legend('$c_{-1}$','$c_{0}$','$c_{+1}$','Location','northwest')
    else
        ylabel('Concentration')
        legend('$c_{0}$','Polynomial fit','Location','northwest')
    end
    box on
    grid on
    
    xlim([1, size(frame,1)/lambda]);
    %ylim([0 1.0-kLim*n]);
    ylim([0 1.0]);
    xticks([0:size(c0,2)/8:size(c0,2)])
    yticks([0:1/5:1.0])
    
    % yticks([0, 23, 143])
    %yticklabels({'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'})
    
    %     %Rotate ylabel, taking into account its size/centre relation.
    %     ylh = get(gca,'ylabel');
    %     gyl = get(ylh);
    %     ylp = get(ylh, 'Position');
    %     set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
    %tightfig;
    
    set(gca,'FontSize',fontSize)
    if export ~= true || sum(F == 'png') == 3
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
                if timeDep == true
                    filename = sprintf([directory '_MCS_' num2str(round(MCSTemp,0)) '_c0_0%d-solventDistribution-timeDep.' F],str2num(strrep(num2str(round(sum(c0)/size(c0,2),2)),'.','')));
                else
                    filename = sprintf([directory '_MCS_' num2str(round(MCSTemp,0)) '_c0_0%d-solventDistribution.' F],str2num(strrep(num2str(round(sum(c0)/size(c0,2),2)),'.','')));
                end
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
        if timeDep == true
            filename = [directory '-solventDistribution-timeDep.gif'];
        else
            filename = [directory '-solventDistribution.gif'];
        end
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
    cutoff = 15;
    
    frame = importdata([directory '/frame-' num2str(n-skipFrames) '.dat']);
    MCS = MCS(MCS~=0); x0 = x0(x0~=0); %Remove zeros.
    X = 0:0.1:MCS(end);
    
    if logLog == false
        fun = fit(MCS',x0','power1'); coeffs = coeffvalues(fun);
        
        clf;
        % Plotting
        %figure
        h1 = axes;
        set(gca,'FontSize',fontSize)
        hold on
        plot(MCS,x0,'.k', 'MarkerSize',20)
        ws = warning('off','all');  % Turn off warnings.
        plot(X,fun(X),'-m')
        warning(ws)  % Turn them back on.
        hold off
        
        % Cosmetic plot stuff.
        sf = ['Power law fit: $y = ' num2str(round(coeffs(1),2)) '\cdot x^{' num2str(round(coeffs(2),2)) '}$'];
        
        xlabel('MCS')
        ylabel('$j_{\textnormal{demix}}$')
        %title(['Concentration of zeros = ' num2str(round(sum(c0)/size(c0,2),2))])
        legend('Data points',sf,'Location','southeast')
        %legend('Data points','Exponential fit','Linear fit','Location','northwest')
        box on
        grid on
        
        xlim([0, max(MCS)]);
        ylim([1, size(frame,1)/lambda]);
        yticks([0:size(c0,2)/4:size(c0,2)])
        xticks([0:max(MCS)/4:max(MCS)])
    else
        lMCS = log(MCS');
        lx0 = log(x0');
        fun = fit(lMCS(1:end-cutoff),lx0(1:end-cutoff),'poly1'); coeffs = coeffvalues(fun);
        if cutoff > 0
            fun2 = fit(lMCS(end-cutoff+1:end),lx0(end-cutoff+1:end),'poly1'); coeffs2 = coeffvalues(fun2);
            X2 = [min(lMCS(end-cutoff:end)):0.01:max(lMCS(end-cutoff:end))]';
        end
        X = [min(lMCS(1:end-cutoff)):0.01:max(lMCS(1:end-cutoff))]';
        X = [min(lMCS):0.01:max(lMCS)]';
        
        clf;
        % Plotting
        %figure
        h1 = axes;
        set(gca,'FontSize',fontSize)
        hold on
        if cutoff > 0
            plot(lMCS(1:end-cutoff),lx0(1:end-cutoff),'.k', 'MarkerSize',20)
            plot(lMCS(end-cutoff:end),lx0(end-cutoff:end),'.', 'MarkerSize',20, 'Color', [0 0 0] + 0.70)
        else
            plot(lMCS,lx0,'.k', 'MarkerSize',20)
        end
        ws = warning('off','all');  % Turn off warnings.
        plot(X,fun(X),'-m')
        if cutoff > 0
            plot(X2,fun2(X2),'-b')
        end
        warning(ws)  % Turn them back on.
        hold off
        
        % Cosmetic plot stuff.
        sf = ['Fit: $y = ' num2str(round(exp(coeffs(2)),2)) '\cdot x^{' num2str(round(coeffs(1),2)) '}$'];
        
        if cutoff > 0
           sf2 = ['Fit: $y = ' num2str(round(exp(coeffs2(2)),2)) '\cdot x^{' num2str(round(coeffs2(1),2)) '}$']; 
        end
        
        xlabel('MCS $[\ln]$')
        ylabel('$j_{\textnormal{demix}}$ $[\ln]$')
        %title(['Concentration of zeros = ' num2str(round(sum(c0)/size(c0,2),2))])
        if cutoff > 0
            legend('Data points','Excluded data points', sf, sf2, 'Location','northwest')
        else
            legend('Data points',sf,'Location','northwest')
        end
        %legend('Data points','Exponential fit','Linear fit','Location','northwest')
        box on
        grid on
        
        xlim([min(lMCS), max(lMCS)]);
        ylim([min(lx0), max(lx0)]);
        yticks([min(lMCS):(max(lMCS)-min(lMCS))/4:max(lMCS)])
        xticks([min(lx0):(max(lx0)-min(lx0))/4:max(lx0)])
        yticklabels({})
        xticklabels({})
    end
    
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
        filename = sprintf([directory '-timeDep.' F]);
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
