clear; % Clears old variables.
clc; % Clears command window.
clf; % Clears figures.
%close all; % Closes any open windows.

%% LaTeX stuff.
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%prefix = '';
%prefix = 'automatedRun/512/';
%prefix = 'debug/';
%prefix = 'J_str/';
%prefix = 'PBCvsFBC/';
%prefix = 'solventDistribution/';
prefix = 'topView/';

%folder = 'lambda_4-L_256-J_0.0000_1.0000_0.0000-numIters_2-22-initialDist_80_10_10-FBC';

cellVisualisation = false; cD = 16; %cD is the colour-depth (8 for 8 bit, 12 for 12 bit etc).
linInt = false; mag = 20; %Applies linear interpolation to the frames; mag is the magnification (e.g. 20 times).
gridOn = false; %will be disabled if linInt = true.
skipFrames = 10;

FourierTransform = true; %disables gridOn an shows the fft image.
radialDist = true; criticalRegion = false; critUp = 8; critLow = 0;
%FTMap = parula(2^cD);
%FTMap = pink(2^cD);
%FTMap = jet(2^cD);

  %cMap = parula(2^cD);
  %cMap = pink(2^cD);
  cMap = jet(2^cD);
  dataMax = 2^4;
  dataMin = 2^1;
  centerPoint = 1;
  scalingIntensity = 5;
    y = 1:length(cMap); 
  y = y - (centerPoint-dataMin)*length(y)/(dataMax-dataMin);
  y = scalingIntensity * y/max(abs(y));
   y = sign(y).* exp(abs(y));
  y = y - min(y); y = y*511/max(y)+1; 
  FTMap = interp1(y, cMap, 1:512);

f = 'pdf'; %pdf or png!
export = false; %Turns on the frame export! For GIF exporting, use exporGIF below. DO NOT USE BOTH!
gcaOnly = false;
exportGIF = false;
pauseTime = 0.1; %The time between each frame in the GIF.

sequence = true; %true for whole sequence (always true for exports).
once = false; %false for currently running simulations.

if exist('folder') == 0
    % Get a list of all files and folders in this folder.
    files = dir(prefix);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories and remove '.' and '..'.
    subFolders = files(dirFlags);
    for k = 1 : length(subFolders)
        x(k) = sum(subFolders(k).name ~= '.') ~= 0;
    end
    subFolders = subFolders(x~=0);
    % Determine maxmum lenght.
    leng = [];
    for k = 1 : length(subFolders)
        leng = [leng size(subFolders(k).name,2)];
    end
    maxLeng = max(leng);
    % Sort by date modified.
    x = [1:length(subFolders)];
    [sortedDates order] = sort([subFolders(x).datenum],'Descend');
    % Print folder names to command window.
    for k = 1 : length(subFolders)
        fprintf('Folder #%d = %s%s', k, subFolders(order(k)).name, blanks(maxLeng-leng(order(k))));
        fprintf(['\tModified = ', char(datetime(sortedDates(k),'ConvertFrom','datenum','Format','dd/MM'' ''HH'':''mm')),'\n'])
    end
    prompt='\nPlease select a folder...\n';
    x = input(prompt);
    folder = subFolders(order(x)).name;
    clc;
end
directory = [prefix folder];

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

if export == true && exportGIF == true
    fprintf('Pick either "export" or "exportGIF", you pillock!\n')
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

if criticalRegion == true
    load('x0.mat','x0Exp');
    x0 = x0Exp;
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
        if(size(frame,1) ~= size(frame,2))
            break;
        end
        
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
            
            HSVimSat = HSVSolvSat;
            HSVim(:, :, 2) = HSVimSat;
            im = hsv2rgb(HSVim);
            
            if linInt == true
                im = im2double(im); %First convert the image into double
                F = griddedInterpolant(im); %Then create 'a gridded interpolant object for the image'
                F.Method = 'linear'; % Select interpolation method (spline is the best one)
                
                [sx,sy,sz] = size(im); % Record the sizes of the image data
                xq = (1:1/mag:sx)'; % Make the grid finer for the "x and y" data
                yq = (1:1/mag:sy)';
                zq = (1:sz)'; % Preserve the colour data
                imTest = F({xq,yq,zq}); % Apply the interpolation to the data
                
                gridOn = false;
                imagesc(imTest);
            else
                if FourierTransform == true
                    colormap(FTMap);
                    Y = fftshift(fft2(im([1:size(im,1)],:,1))) + fftshift(fft2(im([1:size(im,1)],:,2))) + fftshift(fft2(im([1:size(im,1)],:,3)));
                    if criticalRegion == true
                        if x0(n) > 0
                            temp = sort(unique([floor(x0(n))-critLow:floor(x0(n)) floor(x0(n)):floor(x0(n))+critUp])); temp = temp(temp>0); temp = temp(temp<size(im,1));
                            YTemp = fftshift(fft2(im([temp],:,1))) + fftshift(fft2(im([temp],:,2))) + fftshift(fft2(im([temp],:,3)));
                            %YTemp = fft2(im([temp],:,1)) + fft2(im([temp],:,2)) + fft2(im([temp],:,3));
                            YTemp = abs(YTemp);
                        end
                    else
                        Y = fftshift(fft2(im([1:size(im,1)],:,1))) + fftshift(fft2(im([1:size(im,1)],:,2))) + fftshift(fft2(im([1:size(im,1)],:,3)));
                        YTemp = abs(Y);
                    end
                    gridOn = false;
                    if radialDist == true && exist('YTemp') == 1
                        if criticalRegion ~= true
                            r = [0:1:size(YTemp,1)/2];
                            rAvg = radialAverage(YTemp, size(YTemp,1)/2+1, size(YTemp,1)/2+1, r);
                            %[pks,locs] = findpeaks(rAvg(1:8),r(1:8)); % Find peaks within radius 8.
                            [pks,locs] = findpeaks(rAvg(1:end),r(1:end)); % Find peaks.
                        else
                            [maxValue, linearIndexesOfMaxes] = max(YTemp(:));
                            [rowsOfMaxes colsOfMaxes] = find(YTemp == maxValue);
                            xData = YTemp(rowsOfMaxes,sort(unique([colsOfMaxes:colsOfMaxes+12])));
                            x = [1:size(xData,2)]-1;
                            [pks,locs] = findpeaks(xData,x); % Find peaks along critical x-axis.
                        end
                        if isempty(locs) ~= 1 && x0(n) ~= 0% Save largest peak
                            [maxValue, linearIndexesOfMaxes] = max(pks);
                            rowsOfMaxes = find(pks == maxValue);
                            nSave(n) = n;
                            locsSave(n) = locs(rowsOfMaxes);
                        end
                        clf;
                        hold on
                        if criticalRegion ~= true
                            plot(r,rAvg,'.-k', 'MarkerSize',20);
                        else
                            plot(x,xData,'.-k', 'MarkerSize',20);
                        end
                        xPrime = sort([1:(L/lambda)/8:L/lambda]-1);
                        xPrime = xPrime(xPrime ~= 0);
                        xticks(xPrime);
                        if size(pks,2) ~= 0 && n > 1
                            %                             if size(pks,2) >= 2
                            %                                 plot([locs(1) locs(1)],[0 pks(1)], '--', 'Color', [0 0 0] + 0.5) % x = peak 1
                            %                                 plot([locs(2) locs(2)],[0 pks(2)], '--', 'Color', [0 0 0] + 0.5) % x = peak 2
                            %                                 plot(locs(1:2),pks(1:2),'.m', 'MarkerSize',20);
                            %                                 xticks(sort(unique([xticks locs(1:2)])))
                            %                             else
                            %                                 plot([locs(1) locs(1)],[0 pks(1)], '--', 'Color', [0 0 0] + 0.5) % x = peak 1
                            %                                 plot(locs(1),pks(1),'.m', 'MarkerSize',20);
                            %                                 xticks(sort(unique([xticks locs(1)])))
                            %                             end
                            
                            [maxValue, linearIndexesOfMaxes] = max(pks);
                            rowsOfMaxes = find(pks == maxValue);
                            plot([locs(rowsOfMaxes) locs(rowsOfMaxes)],[0 maxValue], '--', 'Color', [0 0 0] + 0.5) % dominatin peak
                            plot(locs(rowsOfMaxes),maxValue,'.m', 'MarkerSize',20);
                            xticks(sort(unique([xticks locs(rowsOfMaxes)])))
                        end
                        hold off
                    elseif exist('YTemp') == 1
                        imagesc(YTemp);
                    end
                else
                    imagesc(im);
                    if criticalRegion == true && x0(n) ~= 0
                        hold on;
                        %h = plot([0,L], [floor(x0(n)),floor(x0(n))], '-k', 'LineWidth',4); h.Color(4)=0.5;
                        h = rectangle('Position',[0 floor(x0(n))-critLow-1/2 L/lambda+1 critLow+critUp+1], 'FaceColor', 'b'); h.FaceColor(4)=0.3; h.EdgeColor(4)=0.0;
                        if export == true && sum(f == 'pdf') == 3
                            h = rectangle('Position',[0 floor(x0(n))-1/2 L/lambda+1 1], 'FaceColor', 'b'); h.FaceColor(4)=0.305; h.EdgeColor(4)=0.0;
                        else
                            h = rectangle('Position',[0 floor(x0(n))-1/2 L/lambda+1 1], 'FaceColor', 'b'); h.FaceColor(4)=0.3; h.EdgeColor(4)=0.0;
                        end
                        hold off
                    end
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
                Y = fftshift(fft2(frame([1:size(frame,1)],:)));
                YTemp = abs(Y);
%                                 if radialDist == true
%                                     r = [1:1:size(YTemp,2)/2];
%                                     rAvg = radialAverage(YTemp, size(YTemp,2)/2, size(YTemp,2)/2, r);
%                                     [pks,locs] = findpeaks(rAvg(1:8),[1:8]);
%                                     clf;
%                                     hold on
%                                     plot(r,rAvg,'.-k', 'MarkerSize',20);
%                                     if size(pks,2) ~= 0
%                                         plot([locs(1) locs(1)],[0 pks(1)], '--', 'Color', [0 0 0] + 0.5) % x = peak
%                                         plot(locs(1),pks(1),'.m', 'MarkerSize',20);
%                                     end
%                                     hold off
%                                 else
%                                     imagesc(YTemp);
%                                     gridOn = false;
%                                 end
                imagesc(YTemp);
                gridOn = false;
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
        if FourierTransform == true && radialDist == true && exist('r') == 1
            xlim([min(r), max(r)]);
            ylim([0, 2000]*(4/lambda));
            %ylim([0, 5000]);
            yticks([])
            %             xticks(sort([2 1:(L/lambda)/8:L/lambda]-1))
            %             if size(pks,2) ~= 0 && n > 1
            %                 xticks(sort(unique([xticks locs(1)])))
            %             end
            
            if export == true
                set(gca,'FontSize',18) % FOR EXPORTING!!
                title('')
                xlim([min(r), max(r)/2]);
            end
            
            set(gcf,'Units','pixels');
            set(gcf,'Position', [0 0 550 400]*1.5)
            set(gcf,'color','w');
            xlabel('$r$')
            ylabel('Counts (arb.\ units)')
            tightfig;
        elseif FourierTransform == true && criticalRegion == true
            set(gcf,'Units','pixels');
            if n == 1
                if radialDist == false
                    ax = gca;
                    ax.YAxis.Visible = 'off';
                    ax.XAxis.Visible = 'off';
                    set(gcf,'Position', [0 0 width width*((critUp + critLow + 1)/size(Y,1))])
                    set(gca,'Position', [0 0 1 1])
                else
                    yticks([])
                    set(gcf,'Position', [0 0 550 400]*1.5)
                    set(gcf,'color','w');
                    xlabel('$x$')
                    ylabel('Counts (arb.\ untits)')
                    tightfig;
                end
            end
            if exist('Y') == 1 && exist('YTemp') == 1
                if radialDist == false
                    ax = gca;
                    ax.YAxis.Visible = 'off';
                    ax.XAxis.Visible = 'off';
                    set(gcf,'Position', [0 0 width width*(size(YTemp,1)/size(Y,1))])
                    set(gca,'Position', [0 0 1 1])
                else
                    xlim([min(x), max(x)]);
                    ylim([0, 300]);
                    yticks([])
                    set(gcf,'Units','pixels');
                    set(gcf,'Position', [0 0 550 400]*1.5)
                    set(gcf,'color','w');
                    xlabel('$x$')
                    ylabel('Counts (arb.\ untits)')
                    tightfig;
                end
            end
        else
            ax = gca;
            ax.GridAlpha = 0.8;
            ax.LineWidth = 1.0;
            ax.YAxis.Visible = 'off';
            ax.XAxis.Visible = 'off';
            set(gcf,'Units','pixels');
            set(gcf,'Position', [0 0 width height])
            set(gca,'Position', [0 0 1 width/height])
        end
        if n == 1
            pause(1)
        end
        if export == true
            if gcaOnly == true
                set(gcf,'Position', [0 0 800 800])
                set(gca,'Position', [0 0 1 1])
            elseif sum(f == 'pdf') == 3
                if FourierTransform == true
                    if criticalRegion == true && radialDist == true
                        set(gcf,'Units','pixels');
                        set(gcf,'Position', [0 0 550 400])
                        set(gcf,'color','w');
                        tightfig;
                    elseif criticalRegion == true && radialDist == false
                        set(gca,'Position', [0 0 1 1])
                    else
                        set(gcf,'Position', [0 0 300 300])
                        set(gca,'Position', [0 0 1 1])
                    end
                else
                    set(gcf,'Position', [0 0 300 300])
                    if gridOn == true
                        set(gca,'Position', [0.003 0.005 0.99 0.99])
                    else
                        set(gca,'Position', [0 0 1 1])
                    end
                end
            end
            for k = 1:9
                if c0 <= 0.1
                    k = 0.1;
                end
                if round(c0,2) == k/10 && k/10 < current || n == b
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

%% Plotting
if FourierTransform == true && radialDist == true
    clf;
    nSave = nSave(nSave ~= 0); % Strip zeros.
    locsSave = locsSave(nSave ~= 0); % Strip zeros.
    
    plot(nSave,locsSave,'.k', 'MarkerSize',20)
    
    % Cosmetic plot stuff.
    set(gca,'FontSize',14)
    xlabel('MCS (arb.\ units)')
    if criticalRegion ~= true
        ylabel('$r$')
    else
        ylabel('$x$')
    end
    box off
    
    xlim([min(nSave), max(nSave)]);
    ylim([1, 8]);
    yticks([1:8])
    xticks([])
    %     yticklabels({})
    %     xticklabels({})
    
    set(gcf,'Units','pixels');
    set(gcf,'Position', [0 0 550 400])
    set(gcf,'color','w');
    tightfig;
    
    if export == true
        fig = gcf;
        filename = sprintf([directory '-domainCoarsening.' f]);
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