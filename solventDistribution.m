clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%prefix = 'squareScaling/';
prefix = 'PBCvsFBC/';
%prefix = 'evapRand/';
%prefix = '';
if isempty(prefix) == true
    directory = [prefix 'frames'];
else
    directory = [prefix 's_lambda_4-L_256_frames_v3-c_1-J_squareScaled-numIters_2-24-FBC'];
end
lambda = 4;
skipFrames = 1;

exportGIF = true;
pauseTime = 0.2;

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
    for i = 1:lambda:size(frame,1)
        cTemp = 0;
        for x1 = i:1:i+lambda-1
            cTemp = cTemp + size(frame,1) - nnz(frame(x1,:));
        end
        cTemp = cTemp/(lambda^2 * size(frame,1)/lambda);
        c0(k) = cTemp;
        k = k+1;
    end
    
    f = fit([1:1:size(frame,1)/lambda]',c0','exp1');
    g = fit([1:1:size(frame,1)/lambda]',c0','poly1');
    
    clf;
    % Plotting
    %figure
    h1 = axes;
    set(gca,'FontSize',12)
    hold on
    plot([1:1:size(frame,1)/lambda]',c0,'.k')
    plot([1:0.01:size(frame,1)/lambda],f([1:0.01:size(frame,1)/lambda]),'-m')
    plot([1:0.01:size(frame,1)/lambda],g([1:0.01:size(frame,1)/lambda]),'--','Color', [95 95 255]/255)
    hold off
    
    % Cosmetic plot stuff.
    xlabel('$i$')
    ylabel('$c_0$')
    title(['Concentration of zeros = ' num2str(round(sum(c0)/size(c0,2),2))])
    %legend('Data points','Exponential fit','Location','northwest')
    legend('Data points','Exponential fit','Linear fit','Location','northwest')
    box on
    grid on
    
    xlim([1, size(frame,1)/lambda]);
    ylim([0 1.0-kLim*n]);
    xticks([0:size(c0,2)/8:size(c0,2)])
    %xticklabels(split(num2str(log2(unique([sort(x)])))))
    % yticks([0, 23, 143])
    %yticklabels({'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8'})
    
%     %Rotate ylabel, taking into account its size/centre relation.
%     ylh = get(gca,'ylabel');
%     gyl = get(ylh);
%     ylp = get(ylh, 'Position');
%     set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
    %tightfig;
    set(gcf,'Position', [0 0 650 600])
    
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