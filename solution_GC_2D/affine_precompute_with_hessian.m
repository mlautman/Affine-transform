function data = affine_precompute_with_hessian(I, J)

% Create the output structure
data=struct();

% Compute the px array
data.px = cell(3,1);  

% Generate the coordinate grid
[data.px{1} data.px{2}] = ...
    ndgrid(1:size(I,1), 1:size(I,2));

% Stick ones in the fourth slot of px
data.px{3} = ones(size(data.px{1}));

% Generate the gradient of the moving image
data.gradJ = cell(2,1);
[data.gradJ{2} data.gradJ{1}] = gradient(J);

% Generate the Hessian of the moving image
data.hessJ = cell(2,2);
[data.hessJ{1,2} data.hessJ{1,1}] = gradient(data.gradJ{1});
[data.hessJ{2,2} data.hessJ{2,1}] = gradient(data.gradJ{2});
