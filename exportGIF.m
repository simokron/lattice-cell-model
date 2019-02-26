clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

directory = 'se_lambda_32-L_1024_frames_v3-c_1-J_squareScaled-numIters_2-22-FBC_evapRand';
lambda = 1;
numIters = 2^22;
skipFrames = 2;
pauseTime = 0.05;

export = true; %Turns on the export!
height = 825;
width = 800;

sequence = true; %true for whole sequence
once = true; %false for currently running simulations

a = dir([directory '/*.dat']);
b = numel(a);

fprintf(['numFrames = ' num2str(b) '\n'])

go = true; tempPause = true;
while go
    a = dir([directory '/*.dat']);
    b = numel(a);
    
    clc;
    fprintf(['numFrames = ' num2str(b) '\n'])
    
    if sequence == true
        lowLim = 1;
    else
        lowLim = b;
    end
    
    for n = lowLim:skipFrames:b
        frame = importdata([directory '/frame-' num2str(n) '.dat']);
        
        c0 = 1 - nnz(frame)/numel(frame);
        MCS = numIters*(n-1)/(size(frame,1)*size(frame,2));
        
        mymap = [1 1 0
            1 0 0
            0 0 1];
        
        if n == 1
            colormap(mymap)
            imagesc(frame);
            title(['Concentration of zeros = ' num2str(round(c0,2)) '; MCS = ' num2str(MCS)])
            if lambda ~= 1
                xticks([0:lambda:size(frame,1)]+0.5)
                yticks([0:lambda:size(frame,2)]+0.5)
                grid on
            end
            ax = gca;
            ax.GridAlpha = 0.8;
            ax.LineWidth = 1.0;
            ax.YAxis.Visible = 'off';
            ax.XAxis.Visible = 'off';
            set(gcf,'Position', [1920-width/2 720-height/2 width height])
            set(gca,'Position', [0 0 1 width/height])
            %set(gca,'visible','off');
            %set(findall(gca, 'type', 'text'), 'visible', 'on')
            if export == true
                pause(1)
                h = gcf;
                filename = [directory '.gif'];
                % Capture the plot as an image
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % Write to the GIF File
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',2);
                %               fig = gcf;
                %               filename = sprintf([num2str(lambda) '_' num2str(size(frame,1)) '_beta_' strrep(num2str(beta),'.','') '_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                %               fig.PaperUnits = 'points';
                %               fig.PaperPosition = [0 0 300 300];
                %               print(filename,['-d' f]);
                %               set(findall(gca, 'type', 'text'), 'visible', 'on')
            end
            pause(1);
        else
            colormap(mymap)
            imagesc(frame);
            title(['Concentration of zeros = ' num2str(round(c0,2)) '; MCS = ' num2str(MCS)])
            if lambda ~= 1
                xticks([0:lambda:size(frame,1)]+0.5)
                yticks([0:lambda:size(frame,2)]+0.5)
                grid on
            end
            ax = gca;
            ax.GridAlpha = 0.8;
            ax.LineWidth = 1.0;
            ax.YAxis.Visible = 'off';
            ax.XAxis.Visible = 'off';
            if export == true %&& n == 900
                %set(findall(gca, 'type', 'text'), 'visible', 'off')
                h = gcf;
                % Capture the plot as an image
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % Write to the GIF File
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',pauseTime);
            end
            %pause(1)
            pause(0.0333); % 0.0167 for 60 FPS, 0.0333 for 30 FPS
            %pause(0.0167);
        end
    end
    sequence = false;

    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',5);
    
    if once == true
        break;
    end
    
    numPause = 0;
    while tempPause
        pause(5)
        a = dir([directory '/*.dat']);
        if b < numel(a)
            break;
        end
        numPause = numPause + 1;
        if numPause > 2
            go = false;
            break;
        end
    end
end