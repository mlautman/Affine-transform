function [f g] = affine_objective(I, J, data, x)

% Get the matrix from x
A = reshape(x(1:9),3,3);
b = x(10:12);

% Coordinate arrays that we will reuse
px = data.px;
qx = cell(3,1); 

% Apply the affine transform to the coordinates
for i = 1:3
    qx{i} = b(i);
    for j = 1:3
        qx{i} = qx{i} + A(i,j) * px{j};
    end
end

% Get the resampled image
K = interpn(J, qx{1}, qx{2}, qx{3}, '*linear', 0);

% Compute the difference image
diff = (K - I);

% Add up the voxels for total mean squared difference
f = sum(sum(sum(diff.^2)));

% Compute the gradient
if nargout > 1
   
    % Compute the gradient with respect to A and b
    gA = zeros(3); gb = zeros(3,1);

    % Interpolate the gradient
    for i = 1:3
        gradK_i = interpn(data.gradJ{i}, qx{1}, qx{2}, qx{3}, '*linear', 0);
        for j = 1:3
            gA(i,j) = 2.0 * sum(sum(sum(diff .* gradK_i .* px{j})));
        end
        gb(i) = 2.0 * sum(sum(sum(diff .* gradK_i)));
    end
    
    g = [gA(:); gb(:)];
    
end

