function resampled = myTransformImage(fixed, moving, A, b, method)

% Handle optional parameter
if nargin < 5
    method = 'linear*';
end

% Get the X, Y and Z coordinates of all voxels in the fixed image
[px py] = ndgrid(1:size(fixed,1), 1:size(fixed,2));

% Apply the affine transform
qx = A(1,1) * px(:) + A(1,2) * py(:) + b(1);
qy = A(2,1) * px(:) + A(2,2) * py(:) + b(2);

% Clear some memory
clear('px','py');

% Interpolate the moving image at the new coordinates
samples = interpn(moving, qx, qy, method, 0);

% Clear some memory
clear('qx','qy');

% Reorganize the samples into a volume of the same size as the fixed image
resampled = reshape(samples, size(fixed));





