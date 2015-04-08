function data = affine_precompute(I, J)

% Create the output structure
data=struct();

% Compute the px array
data.px = cell(3,1);  

% Generate the coordinate grid
[data.px{1} data.px{2} data.px{3}] = ...
    ndgrid(1:size(I,1), 1:size(I,2), 1:size(I,3));

% Generate the gradient of the moving image
data.gradJ = cell(3,1);
[data.gradJ{2} data.gradJ{1} data.gradJ{3}] = gradient(J);

