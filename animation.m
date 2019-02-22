clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

directory = 'frames';
lambda = 1;
beta = 0.6;
numIters = 2^22;

f = 'png'; %pdf or png!
export = false; %Turns on the frame export! For GIF exporting, see separate script.
current = 1.0; %relates to exporting - largest concentration to save (usually best keps at 1.0)
height = 825;
width = 800;

sequence = true; %true for whole sequence
once = false; %false for currently running simulations

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
    
    for n = lowLim:1:b
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
            set(gcf,'Position', [0 0 width height])
            set(gca,'Position', [0 0 1 width/height])
            %set(gca,'visible','off');
            %set(findall(gca, 'type', 'text'), 'visible', 'on')
            if export == true
                pause(1)
                %set(findall(gca, 'type', 'text'), 'visible', 'off')
                fig = gcf;
                filename = sprintf([num2str(lambda) '_' num2str(size(frame,1)) '_beta_' strrep(num2str(beta),'.','') '_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                fig.PaperUnits = 'points';
                fig.PaperPosition = [0 0 300 300];
                print(filename,['-d' f]);
                set(findall(gca, 'type', 'text'), 'visible', 'on')
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
                %             k = 0.11
                for k = 1:9
                    if round(c0,2) == k/10 && k/10 < current
                        current = k/10;
                        %set(findall(gca, 'type', 'text'), 'visible', 'off')
                        fig = gcf;
                        filename = sprintf([num2str(lambda) '_' num2str(size(frame,1)) '_beta_' strrep(num2str(beta),'.','') '_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                        fig.PaperUnits = 'points';
                        fig.PaperPosition = [0 0 300 300];
                        print(filename,['-d' f]);
                        set(findall(gca, 'type', 'text'), 'visible', 'on')
                    end
                end
            end
            %pause(1)
            pause(0.0333); % 0.0167 for 60 FPS, 0.0333 for 30 FPS
            %pause(0.0167);
        end
    end
    sequence = false;
    
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
        if numPause > 20
            go = false;
            break;
        end
    end
end