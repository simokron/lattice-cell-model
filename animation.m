clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

prefix = 'squareScaling/';
%prefix = 'PBCvsFBC/';
%prefix = '';
if isempty(prefix) == true
    directory = [prefix 'frames'];
else
    directory = [prefix 's_lambda_4-L_128_frames_v3-c_1-J_ORIGINALsquareScaled-numIters_2-22-PBC'];
end
lambda = 4;
numIters = 2^22;

cellVisualisation = false;
linInt = true; m = 1.5; %m controls the contrast. Decrease to 1.0 for maximum contrast. 1.5 is a decent value.
gridOn = true; %will be disabled if linInt = true.

FourierTransform = false; %disables gridOn an shows the fft image.
FTMap = parula(4096);

f = 'pdf'; %pdf or png!
export = true; %Turns on the frame export! For GIF exporting, use exporGIF below. DO NOT USE BOTH!
gcaOnly = false;
exportGIF = false;
pauseTime = 0.1; %The time between each frame in the GIF.

sequence = true; %true for whole sequence (always true for exports).
once = true; %false for currently running simulations.

if exportGIF == true
    skipFrames = 1;
else
    skipFrames = 1;
end

%The dimensions in pixels of the png/GIF. Note that the height should be increased to account for the title text.
height = 838;
width = 800;
current = 1.0; %relates to exporting - largest concentration to save (usually best keps at 1.0)

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
    
    if sequence == true || export == true || exportGIF == true
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
                
                % "25% more" saturation:
                HSV(:, :, 2) = HSV(:, :, 2) * 1.25;
                % or add:
                % HSV(:, :, 2) = HSV(:, :, 2) + 0.2;
                HSV(HSV > 1) = 1;  % Limit values
                imTest = hsv2rgb(HSV);
                
                imagesc(imTest);
                gridOn = false;
                
                if FourierTransform == true
                    colormap(FTMap);
                    Y = fftshift(fft2(imTest(:,[1:size(imTest,1)])));
                    imagesc(abs(Y));
                else
                    imagesc(imTest);
                end
                
            else
                if FourierTransform == true
                    colormap(FTMap);
                    Y = fftshift(fft2(im(:,[1:size(im,1)])));
                    imagesc(abs(Y));
                    gridOn = false;
                else
                    imagesc(im);
                    %imagesc(rgb_image); %DEBUGGING
                    sF = 1;
                end
                
            end
        else
            if FourierTransform == true
                colormap(FTMap);
                Y = fftshift(fft2(frame(:,[1:size(frame,1)])));
                imagesc(abs(Y));
            else
                colormap(mymap)
                imagesc(frame);
            end
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
        set(gcf,'Units','pixels');
        set(gcf,'Position', [0 0 width height])
        set(gca,'Position', [0 0 1 width/height])
        %set(gca,'visible','off');
        %set(findall(gca, 'type', 'text'), 'visible', 'on')
        if n == 1
            pause(1)
        end
        if export == true
            if gcaOnly == true
                set(gcf,'Position', [0 0 800 800])
                set(gca,'Position', [0 0 1 1])
            elseif sum(f == 'pdf') == 3
                set(gcf,'Position', [0 0 300 300])
                if gridOn == true
                    set(gca,'Position', [0.003 0.005 0.99 0.99])
                else
                    set(gca,'Position', [0 0 1 1])
                end
            end
            for k = 1:9
                if c0 <= 0.1
                    k = 0.1;
                end
                if round(c0,2) == k/10 && k/10 < current
                    current = k/10;
                    fig = gcf;
                    filename = sprintf([directory '_MCS_' num2str(round(MCS,0)) '_c0_0%d.' f],str2num(strrep(num2str(round(c0,2)),'.','')));
                    if sum(f == 'png') == 3
                        frame = getframe(fig);
                        im = frame2im(frame);
                        [imind,cm] = rgb2ind(im,256);
                        imwrite(imind,cm,filename,f);
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
            filename = [directory '.gif'];
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
        pause(0.0333); % 0.0167 for 60 FPS, 0.0333 for 30 FPS
        %pause(0.0167);
    end
    sequence = false;
    
    if exportGIF == true
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',5);
        break;
    end
    
    if once == true || export == true
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