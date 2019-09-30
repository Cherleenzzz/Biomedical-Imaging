function a = myBSpline(b,c)
switch b
    case 0
        a = ((1-c).^3)/6;
    case 1
        a = (3*c.^3 - 6*c.^2 + 4)/6;
    case 2
        a = (-3*c.^3 + 3*c.^2 + 3*c + 1)/6;
    case 3
        a = (c.^3)/6;
    otherwise
        error('First argument to b-spline must be 0, 1, 2, or 3');
end
end