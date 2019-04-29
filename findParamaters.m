% This function takes a directory and extracts the cell size and the number of iterations per frame (assuming the correct format).
function [lambda, L, numIters] = findParamaters(directory,s)
    if s ~= 1
        fprintf('Extracting the paramaters from the directory name...\n')
    end

    ipos = strfind(directory,'lambda_') + strlength("lambda_");
    iposLim = strfind(directory,'-L_') - 1;
    lambda = str2num(directory(ipos:iposLim));

    ipos = strfind(directory,'L_') + strlength("L_");
    iposLim = strfind(directory,'-J_') - 1;
    L = str2num(directory(ipos:iposLim));

    ipos = strfind(directory,'numIters_2-') + strlength("numIters_2-");
    iposLim = strfind(directory,'-initialDist_') - 1;
    exponent = str2num(directory(ipos:iposLim));
    numIters = 2^exponent;
    
    if s ~= 1
        fprintf(['Done! Found lambda = ' num2str(lambda) ', L = ' num2str(L) ' and numIters = 2^' num2str(log2(numIters)) '.\n'])
    end
end