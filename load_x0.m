% This function loads the j_demix vector from the .mat file.
function [x0, MCS] = load_x0(directory,m)
    x0Name = directory;
    fileExistsx0 = exist(x0Name, 'file');
    if fileExistsx0 ~= 2
        fprintf('\n')
        error('The .mat file does not exists!\nPlease re-run solventDistributions with approporiate settings.',class(x0Name))
    else
        if m == 2
            load(x0Name,'x0Exp');
            x0 = x0Exp;
            load(x0Name,'MCSExp');
            MCS = MCSExp;
        else
            fprintf('\nLoading the j_demix vector...')
            load(x0Name,'x0Exp');
            x0 = x0Exp;
            fprintf('\nDone!\n')
        end
    end
end