% Originally written by Hugo Trentesaux here: https://se.mathworks.com/matlabcentral/answers/251714-radial-average-of-2d-matrix#answer_359975

function profile = radialAverage(IMG, cx, cy, w)
    % computes the radial average of the image IMG around the cx,cy point
    % w is the vector of radii starting from zero
    [a,b] = size(IMG);
    [X, Y] = meshgrid( (1:a)-cx, (1:b)-cy);
    R = sqrt(X.^2 + Y.^2);
    profile = [];
    for i = w % radius of the circle
        mask = (i-1<R & R<i+1); % smooth 1 px around the radius
        values = (1-abs(R(mask)-i)) .* double(IMG(mask)); % smooth based on distance to ring
        % values = IMG(mask); % without smooth
        profile(end+1) = mean( values(:) );
    end
end