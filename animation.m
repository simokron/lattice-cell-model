clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

directory = 'frames/';
a = dir([directory '*.dat']);
b = numel(a)
lambda = 4;
s = 1;

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
        %tightfig;
        if s == 9
            pause(1)
            set(findall(gca, 'type', 'text'), 'visible', 'off')
            fig = gcf;
            filename = sprintf('export_c0_0%d.pdf',str2num(strrep(num2str(round(c0,2)),'.','')));
            fig.PaperUnits = 'points';
            fig.PaperPosition = [0 0 300 300];
            print(filename,'-dpdf');
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
        %         tightfig;
        if s == 9
            if mod(round(c0,2),0.6) == 0 || mod(round(c0,2),0.4) == 0 || mod(round(c0,2),0.2) == 0 || mod(round(c0,2),0.1) == 0
                set(findall(gca, 'type', 'text'), 'visible', 'off')
                fig = gcf;
                filename = sprintf('export_c0_0%d.pdf',str2num(strrep(num2str(round(c0,2)),'.','')));
                fig.PaperUnits = 'points';
                fig.PaperPosition = [0 0 300 300];
                print(filename,'-dpdf');
                set(findall(gca, 'type', 'text'), 'visible', 'on')
            end
        end
        %pause(0.1)
        pause(0.0333); % 0.0167 for 60 FPS, 0.0333 for 30 FPS
        %pause(0.0167);
    end
end