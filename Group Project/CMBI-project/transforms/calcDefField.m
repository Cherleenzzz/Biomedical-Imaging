function defField = calcDefField(CPG,cxs,cys,czs,rxs,rys,rzs)
%calcDefField function to calculate a deformation field from a B-spline
%control point grid
%
%INPUTS: CPG: a 4D array containing the B-spline control point grid
%           cxs: the x coords of the B-spline control point grid
%           cys: the y coords of the B-spline control point grid
%           czs: the z coords of the B-spline control point grid
%           rxs: the x coords of the reference image (i.e. the points where
%           the deformation field will be calculated)
%           rys: the y coords of the reference image
%           rzs: the z coords of the reference image
%OUTPUTS: defField: a 4D array containing the deformation field

%tol = tolerance on rounding error
tol = 0.01;

%find first cp to effect each x, y, z coord (i, j, k in Rueckert99) and
%distance to each cp (u, v, w in Rueckert99)
i = interp1(cxs,1:length(cxs),rxs);
u = i - floor(i);
i = floor(i) - 1;
j = interp1(cys,1:length(cys),rys);
v = j - floor(j);
j = floor(j) - 1;
%check if 3D
if length(czs) > 1
    k = interp1(czs,1:length(czs),rzs);
    w = k - floor(k);
    k = floor(k) - 1;
end

%check for rounding error
i(1-u < tol) = i(1-u < tol) + 1;
u(1-u < tol) = 0;
j(1-v < tol) = j(1-v < tol) + 1;
v(1-v < tol) = 0;
if length(czs) > 1
    k(1-w < tol) = k(1-w < tol) + 1;
    k(k<1)=1;
    w(1-w < tol) = 0;
end

i(i<1)=1;
j(j<1)=1;

%find value of 4 b-spline basis functions for each x, y, z coord
bspline0X = myBSpline(0,u);
bspline1X = myBSpline(1,u);
bspline2X = myBSpline(2,u);
bspline3X = myBSpline(3,u);
bspline0Y = myBSpline(0,v);
bspline1Y = myBSpline(1,v);
bspline2Y = myBSpline(2,v);
bspline3Y = myBSpline(3,v);
if length(czs) > 1
    bspline0Z = myBSpline(0,w);
    bspline1Z = myBSpline(1,w);
    bspline2Z = myBSpline(2,w);
    bspline3Z = myBSpline(3,w);
end



if length(czs) > 1
    A = zeros(length(rxs),size(CPG,2),size(CPG,3),3);
else
    A = zeros(length(rxs),size(CPG,2),size(CPG,3),2);
end
%sum up 4 cp displacements for each x coord
for x = 1:length(rxs)
    A(x,:,:,:) = bspline0X(x).*CPG(i(x),:,:,:) + bspline1X(x).*CPG(i(x)+1,:,:,:) + bspline2X(x).*CPG(i(x)+2,:,:,:) + bspline3X(x).*CPG(i(x)+3,:,:,:);
end

if length(czs) > 1
    B = zeros(length(rxs),length(rys),size(CPG,3),3);
else
    B = zeros(length(rxs),length(rys),size(CPG,3),2);
end
%sum up 4 cp displacements for each y coord
j(j>size(A,2)-3) = size(A,2)-3;
for y = 1:length(rys)
    B(:,y,:,:) = bspline0Y(y).*A(:,j(y),:,:) + bspline1Y(y).*A(:,j(y)+1,:,:) + bspline2Y(y).*A(:,j(y)+2,:,:) + bspline3Y(y).*A(:,j(y)+3,:,:);
end
clear('A');

if length(czs) > 1
    defField = zeros(length(rxs),length(rys),length(rzs),3);
    
    %sum up 4 cp displacements for each z coord
    for z = 1:length(rzs)
        defField(:,:,z,:) = bspline0Z(z).*B(:,:,k(z),:) + bspline1Z(z).*B(:,:,k(z)+1,:) + bspline2Z(z).*B(:,:,k(z)+2,:) + bspline3Z(z).*B(:,:,k(z)+3,:);
    end
    clear('B')
else
    defField = B;
end