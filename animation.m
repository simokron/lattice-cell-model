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
%prefix = '';
if isempty(prefix) == true
    directory = [prefix 'frames'];
else
    directory = [prefix 's_lambda_4-L_256_frames_v3-c_1-J_squareScaled-numIters_2-24-FBC'];
end
lambda = 4;
numIters = 2^28;

cellVisualisation = false;
linInt = false; m = 1.5;
gridOn = false;

f = 'png'; %pdf or png!
export = false; %Turns on the frame export! For GIF exporting, use exporGIF below. DO NOT USE BOTH!
exportGIF = false;

if exportGIF == true
    skipFrames = 1;
else
    skipFrames = 1;
end
pauseTime = 0.1;

current = 1.0; %relates to exporting - largest concentration to save (usually best keps at 1.0)
height = 838;
width = 800;

sequence = false; %true for whole sequence
once = false; %false for currently running simulations

a = dir([directory '/*.dat']);
b = numel(a);

if b == 0
    fprintf('Empty directory...\n')
    fprintf('Aborting!\n')
    return;
end

fprintf(['numFrames = ' num2str(b) '\n'])

if export == true
    clf;
end

go = true; tempPause = true;
while go
    a = dir([directory '/*.dat']);
    b = numel(a);
    
    clc;
    fprintf(['numFrames = ' num2str(b) '\n'])
    if exportGIF == true
        fprintf('Exporting frames as GIF...\n')
    end
    
    if sequence == true !|| export == true || exportGIF == true
        lowLim = 1;
    else
        lowLim = b;
    end
    
    for n = lowLim:skipFrames:b
        frame = importdata([directory '/frame-' num2str(n) '.dat']);
        
        if cellVisualisation == true
            L = size(frame,1);
            for i = 1:lambda:L
                for j = 1:lambda:L
                    numDown = 0;
                    numZero = 0;
                    numUp = 0;
                    for x2 = i:i+lambda-1
                        for x1 = j:j+lambda-1
                            if frame(x2, x1) == -1
                                numDown = numDown + 1;
                            elseif frame(x2, x1) == 0
                                numZero = numZero + 1;
                            elseif frame(x2, x1) == 1
                                numUp = numUp + 1;
                            end
                        end
                    end
                    cDown_cell((i+lambda-1)/lambda,(j+lambda-1)/lambda) = numDown/lambda^2;
                    cZero_cell((i+lambda-1)/lambda,(j+lambda-1)/lambda) = numZero/lambda^2;
                    cUp_cell((i+lambda-1)/lambda,(j+lambda-1)/lambda) = numUp/lambda^2;
                end
            end
            cDiff = cDown_cell - cUp_cell; %White <---> more down!
            
            map = gray(4096);
            minv = min(cZero_cell(:));
            maxv = max(cZero_cell(:));
            ncol = size(map,1);
            s = round(1+(ncol-1)*(cZero_cell-minv)/(maxv-minv));
            
            mapRed = [0:4096]'./4096;
            mapGreen = [zeros(size(0:4096))]';
            mapBlue = [zeros(size(0:4096))]';
            
            mymap = [mapRed mapGreen mapBlue];
            
            rgb_image = ind2rgb(s,mymap);
            
            map = gray(4096);
            minv = min(cDiff(:));
            maxv = max(cDiff(:));
            ncol = size(map,1);
            s = round(1+(ncol-1)*(cDiff-minv)/(maxv-minv));
            
            old = ind2rgb(s,map);
        end
        
        c0 = 1 - nnz(frame)/numel(frame);
        MCS = numIters*(n-1)/(size(frame,1)*size(frame,2));
        
        if cellVisualisation ~= true
            mymap = [1 1 0
                1 0 0
                0 0 1];
        end
        
        if n == 1
            if cellVisualisation == true
                im = imfuse(old, rgb_image,'blend');
                
                Idouble = im2double(im);
                avg = mean2(Idouble);
                sigma = std2(Idouble);
                if avg-m*sigma > 0
                    if avg+m*sigma < 1
                        im = imadjust(im,[avg-m*sigma avg+m*sigma],[]);
                    else
                        im = imadjust(im,[avg-m*sigma 1],[]);
                    end
                else
                    if avg+m*sigma < 1
                        im = imadjust(im,[0 avg+m*sigma],[]);
                    else
                        im = imadjust(im,[0 1],[]);
                    end
                end
                
                if linInt == true
                    im = im2double(im);
                    imInt = interp3(im,2,'linear');
                    
                    for i = 1:size(imInt,1)
                        for j = 1:size(imInt,2)
                            h = 0;
                            for k = 1:size(imInt,3)/3:size(imInt,3)-size(imInt,3)/3+1
                                imTest(i,j,k+h*(1-size(imInt,3)/3)) = sum(imInt(i,j,k:k+size(imInt,3)/3-1))/(size(imInt,3)/3);
                                h = h + 1;
                            end
                        end
                    end
                    
                    HSV = rgb2hsv(imTest);
                    
                    % "20% more" saturation:
                    HSV(:, :, 2) = HSV(:, :, 2) * 1.25;
                    % or add:
                    % HSV(:, :, 2) = HSV(:, :, 2) + 0.2;
                    HSV(HSV > 1) = 1;  % Limit values
                    imTest = hsv2rgb(HSV);
                    
                    imagesc(imTest);
                    sF = size(imTest,1)/size(im,1);
                    
                else
                    imagesc(im);
                    %imagesc(rgb_image); %DEBUGGING
                    sF = 1;
                    
                end
            else
                colormap(mymap)
                imagesc(frame);
            end
            
            set(gca,'FontSize',14)
            title(['Concentration of zeros = ' num2str(round(c0,2)) '; MCS = ' num2str(MCS)])
            if gridOn == true
                if cellVisualisation == true
                    xticks([0:1:size(cDown_cell,1)]*sF+0.5)
                    yticks([0:1:size(cDown_cell,2)]*sF+0.5)
                else
                    xticks([0:lambda:size(frame,1)]+0.5)
                    yticks([0:lambda:size(frame,2)]+0.5)
                end
                grid on
            end
            ax = gca;
            ax.GridAlpha = 0.8;
            ax.LineWidth = 1.0;
            ax.YAxis.Visible = 'off';
            ax.XAxis.Visible = 'off';
            if export == false
                set(gcf,'Position', [0 0 width height])
                set(gca,'Position', [0 0 1 width/height])
            end
            %set(gca,'visible','off');
            %set(findall(gca, 'type', 'text'), 'visible', 'on')
            if export == true
                pause(1)
                if f == 'png'
                    set(gca,'FontSize',10)
                end
                if f ~= 'png'
                    set(findall(gca, 'type', 'text'), 'visible', 'off')
                end
                %set(gca,'Position', [0 0 1 1])
                fig = gcf;
                %filename = sprintf([num2str(lambda) '_' num2str(size(frame,1)) '_beta_' strrep(num2str(beta),'.','') '_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                filename = sprintf([directory '_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                fig.PaperUnits = 'points';
                fig.PaperPosition = [0 0 300 300];
                print(filename,['-d' f]);
                set(findall(gca, 'type', 'text'), 'visible', 'on')
                %set(gca,'Position', [0 0 1 width/height])
            elseif exportGIF == true
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
            end
            pause(1);
        else
            if cellVisualisation == true
                im = imfuse(old, rgb_image,'blend');
                
                Idouble = im2double(im);
                avg = mean2(Idouble);
                sigma = std2(Idouble);
                if avg-m*sigma > 0
                    if avg+m*sigma < 1
                        im = imadjust(im,[avg-m*sigma avg+m*sigma],[]);
                    else
                        im = imadjust(im,[avg-m*sigma 1],[]);
                    end
                else
                    if avg+m*sigma < 1
                        im = imadjust(im,[0 avg+m*sigma],[]);
                    else
                        im = imadjust(im,[0 1],[]);
                    end
                end
                
                if linInt == true
                    im = im2double(im);
                    imInt = interp3(im,2,'linear');
                    
                    for i = 1:size(imInt,1)
                        for j = 1:size(imInt,2)
                            h = 0;
                            for k = 1:size(imInt,3)/3:size(imInt,3)-size(imInt,3)/3+1
                                imTest(i,j,k+h*(1-size(imInt,3)/3)) = sum(imInt(i,j,k:k+size(imInt,3)/3-1))/(size(imInt,3)/3);
                                h = h + 1;
                            end
                        end
                    end
                    
                    HSV = rgb2hsv(imTest);
                    
                    % "20% more" saturation:
                    HSV(:, :, 2) = HSV(:, :, 2) * 1.25;
                    % or add:
                    % HSV(:, :, 2) = HSV(:, :, 2) + 0.2;
                    HSV(HSV > 1) = 1;  % Limit values
                    imTest = hsv2rgb(HSV);
                    
                    imagesc(imTest)
                    sF = size(imTest,1)/size(im,1);
                    
                else
                    imagesc(im);
                    %imagesc(rgb_image); %DEBUGGING
                    sF = 1;
                end
            else
                colormap(mymap);
                imagesc(frame);
            end
            
            set(gca,'FontSize',14)
            title(['Concentration of zeros = ' num2str(round(c0,2)) '; MCS = ' num2str(MCS)])
            if gridOn == true
                if cellVisualisation == true
                    xticks([0:1:size(cDown_cell,1)]*sF+0.5)
                    yticks([0:1:size(cDown_cell,2)]*sF+0.5)
                else
                    xticks([0:lambda:size(frame,1)]+0.5)
                    yticks([0:lambda:size(frame,2)]+0.5)
                end
                grid on
            end
            ax = gca;
            ax.GridAlpha = 0.8;
            ax.LineWidth = 1.0;
            ax.YAxis.Visible = 'off';
            ax.XAxis.Visible = 'off';
            if export == false
                set(gcf,'Position', [0 0 width height])
                set(gca,'Position', [0 0 1 width/height])
            end
            if export == true %&& n == 900
                for k = 1:9
                    %k = 0.1;
                    if round(c0,2) == k/10 && k/10 < current
                        current = k/10;
                        if f == 'png'
                            set(gca,'FontSize',10)
                        end
                        if f ~= 'png'
                            set(findall(gca, 'type', 'text'), 'visible', 'off')
                        end
                        %set(gca,'Position', [0 0 1 1])
                        fig = gcf;
                        %filename = sprintf([num2str(lambda) '_' num2str(size(frame,1)) '_beta_' strrep(num2str(beta),'.','') '_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                        filename = sprintf([directory '_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                        fig.PaperUnits = 'points';
                        fig.PaperPosition = [0 0 300 300];
                        print(filename,['-d' f]);
                        set(findall(gca, 'type', 'text'), 'visible', 'on')
                        %set(gca,'Position', [0 0 1 width/height])
                    end
                end
            elseif exportGIF == true
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
    
    if once == true || export == true
        break;
    elseif exportGIF == true
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',5);
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
        if numPause > 60
            go = false;
            break;
        end
    end
end