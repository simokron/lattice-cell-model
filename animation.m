clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

lambda = 1;
directory = ['frames' '/'];
%numIters = 1*10^6;
numIters = 2^19;
f = 'pdf'; %pdf or png!
s = 1; %9 for saving!

current = 0.4;
a = dir([directory '*.dat']);
b = numel(a)

%pause(1)
for n = 1:b
    frame = importdata([directory 'frame-' num2str(n) '.dat']);
    
    c0 = 1 - nnz(frame)/numel(frame);
    
    if n == 1
        mymap = [1 1 0
            1 0 0
            0 0 1];
        colormap(mymap)
        imagesc(frame);
        title(['Concentration of zeros = ' num2str(round(c0,2))])
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
        %set(gca,'visible','off');
        %set(findall(gca, 'type', 'text'), 'visible', 'on')
        if s == 9
            MCS = numIters*(n-1)/(size(frame,1)*size(frame,2));
            pause(1)
            set(findall(gca, 'type', 'text'), 'visible', 'off')
            fig = gcf;
            filename = sprintf([num2str(lambda) '_' num2str(size(frame,1)) '_export_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
            fig.PaperUnits = 'points';
            fig.PaperPosition = [0 0 300 300];
            print(filename,['-d' f]);
            set(findall(gca, 'type', 'text'), 'visible', 'on')
        end
        pause(1);
    else
        imagesc(frame);
        title(['Concentration of zeros = ' num2str(round(c0,2))])
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
        if s == 9 %&& n == 900
%             k = 0.11
           for k = 1:8
               if round(c0,2) == k/10 && k/10 < current
                    current = k/10;
                    MCS = numIters*(n-1)/(size(frame,1)*size(frame,2));
                    set(findall(gca, 'type', 'text'), 'visible', 'off')
                    fig = gcf;
                    filename = sprintf([num2str(lambda) '_' num2str(size(frame,1)) '_export_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                    fig.PaperUnits = 'points';
                    fig.PaperPosition = [0 0 300 300];
                    print(filename,['-d' f]);
                    set(findall(gca, 'type', 'text'), 'visible', 'on')
               end
           end
        end
        %pause(0.1)
        pause(0.0333); % 0.0167 for 60 FPS, 0.0333 for 30 FPS
        %pause(0.0167);
    end
end