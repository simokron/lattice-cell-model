% This function takes a directory and stores it in 'temp'. It then lists all of the subdirectories and sorts them by the last modified date, prompting the user to pick one.
function folder = listDirs(temp)
    % Get a list of all files and folders in this folder.
    files = dir(temp);
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
    x = input('\nPlease select a folder...\n');
    folder = subFolders(order(x)).name;
    clc;
end