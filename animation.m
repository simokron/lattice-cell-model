clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%prefix = '';
%prefix = 'automatedRun/256/';
%prefix = 'debug/';
%prefix = 'J_str/';
%prefix = 'PBCvsFBC/';
%prefix = 'solventDistribution/';
prefix = 'topView/';

folder = 'lambda_4-L_512-J_0.0000_1.0000_2.0000-numIters_2-26-initialDist_60_35_5-topView';

directory = [prefix folder];

cellVisualisation = true; cD = 16;
linInt = false; %Applies linear interpolation to the frames.
gridOn = false; %will be disabled if linInt = true.

FourierTransform = false; %disables gridOn an shows the fft image.
FTMap = parula(2^cD);

f = 'pdf'; %pdf or png!
export = false; %Turns on the frame export! For GIF exporting, use exporGIF below. DO NOT USE BOTH!
gcaOnly = false;
exportGIF = false;
pauseTime = 0.1; %The time between each frame in the GIF.

sequence = true; %true for whole sequence (always true for exports).
once = false; %false for currently running simulations.

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

ipos = strfind(directory,'lambda_') + strlength("lambda_");
iposLim = strfind(directory,'-L_') - 1;
lambda = str2num(directory(ipos:iposLim));

ipos = strfind(directory,'numIters_2-') + strlength("numIters_2-");
iposLim = strfind(directory,'-initialDist_') - 1;
exponent = str2num(directory(ipos:iposLim));
numIters = 2^exponent;

if b == 0
    fprintf('Empty directory...\n')
    fprintf('Aborting!\n')
    return;
end

T = struct2table(a);
sortedT = sortrows(T, 'date');
sortedA = table2struct(sortedT);

cT = datetime - datetime(getfield(sortedA(b),'date'));
timeSinceLF = seconds(cT);

if timeSinceLF > 60*15
    once = true;
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
        fprintf('Exporting frames as gif...\n')
    elseif export == true
        fprintf(['Exporting frames as ' f '...\n'])
    end
    
    if sequence == true || export == true || exportGIF == true
        lowLim = 1;
    else
        lowLim = b;
    end
    
    for n = lowLim:skipFrames:b
        frame = importdata([directory '/frame-' num2str(n) '.dat']);
        c0 = 1 - nnz(frame)/numel(frame);
        MCS = numIters*(n-1)/(size(frame,1)*size(frame,2));
        
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
            %cDiff = cDiff - cZero_cell;
            
            mapRed = [0:2^cD-1]'./(2^cD-1);
            mapGreen = [zeros(size(0:2^cD-1))]';
            mapBlue = [zeros(size(0:2^cD-1))]';
            solvMap = [mapRed mapGreen mapBlue];
            
            mapRed = [0:2^cD-1]'./(2^cD-1);
            mapGreen = [0:2^cD-1]'./(2^cD-1);
            mapBlue = [0:2^cD-1]'./(2^cD-1);
            diffMap = [mapRed mapGreen mapBlue];
            
            mapRed = [2^cD-1:-1:0]'./(2^cD-1);
            mapGreen = [2^cD-1:-1:0]'./(2^cD-1);
            mapBlue = [2^cD-1:-1:0]'./(2^cD-1);
            upMap = [mapRed mapGreen mapBlue];
            
            minv = 0;
            maxv = 1;
            ncol = size(solvMap,1);
            solv = round(1+(ncol-1)*(cZero_cell-minv)/(maxv-minv));
            solventRGB = ind2rgb(solv,solvMap);
            
            down = round(1+(ncol-1)*(cDown_cell-minv)/(maxv-minv));
            downRGB = ind2rgb(down,diffMap);
            
            up = round(1+(ncol-1)*(cUp_cell-minv)/(maxv-minv));
            upRGB = ind2rgb(up,upMap);
            
            minv = -1;
            maxv = 1;
            diff = round(1+(ncol-1)*(cDiff-minv)/(maxv-minv));
            diffRGB = ind2rgb(diff,diffMap);
            
        end
        
        if cellVisualisation ~= true
            mymap = [1 1 0
                1 0 0
                0 0 1];
        end
        
        if cellVisualisation == true
            im = imadd(solventRGB,diffRGB);
            
            im = im2double(im);
            imSolv = im2double(solventRGB);
            
            HSVim = rgb2hsv(im);
            HSVSolv = rgb2hsv(imSolv);
            
            % Correct the saturation of the merged images.
            HSVimSat = HSVim(:, :, 2);
            HSVSolvSat = HSVSolv(:, :, 3);
            %HSVimSat(HSVimSat <= 0) = eps;
            %HSVimSat = HSVimSat .* (HSVSolv(:,:,3)./HSVimSat);
            %HSVimSat = HSVimSat + 0.9*cZero_cell*HSVSolv(:,:,3);
            %HSVimSat = HSVimSat + 0.5*HSVSolv(:,:,3);
            %HSVimSat = HSVimSat + 0.25*HSVSolv(:,:,3);
            %HSVimSat = HSVSolv(:,:,3);
            %HSVimSat = HSVimSat;
            
            HSVimSat = HSVSolvSat;
            HSVim(:, :, 2) = HSVimSat;
            im = hsv2rgb(HSVim);
            
            if linInt == true
                im = im2double(im); %First convert the image into double
                F = griddedInterpolant(im); %Then create 'a gridded interpolant object for the image'
                F.Method = 'linear'; % Select interpolation method (spline is the best one)
                
                [sx,sy,sz] = size(im); % Record the sizes of the image data
                xq = (1:0.05:sx)'; % Make the grid finer for the "x and y" data
                yq = (1:0.05:sy)';
                zq = (1:sz)'; % Preserve the colour data
                imTest = F({xq,yq,zq}); % Apply the interpolation to the data
                
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
                    %imagesc(downRGB);
                    %imagesc(upRGB);
                    %imagesc(diffRGB); %BLACK AND WHITE
                    %imagesc(solventRGB); %RED
                    
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
        %pause(0.5)
        %pause(0.1)
        pause(0.0333); % 0.0167 for 60 FPS, 0.0333 for 30 FPS
        %pause(0.0167);
        %sum(frame(:))/(size(frame,1)^2)
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
        if numPause > 180
            go = false;
            break;
        end
    end
end