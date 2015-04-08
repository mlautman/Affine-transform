function [f g H] = affine_masked_objective_with_hessian(I, J, M, data, x)

% Get the matrix from x
A = reshape(x(1:4),2,2);
b = x(5:6);

% Coordinate arrays that we will reuse
qx = cell(2,1); 

% We only care about voxels in the masked region
imask = M > 0;

% Get only the masked coordinates
px = cell(3,1);
for i = 1:3
    px{i} = data.px{i}(imask);
end

% Apply the affine transform to the coordinates
for i = 1:2
    qx{i} = b(i);
    for j = 1:2
        qx{i} = qx{i} + A(i,j) * px{j};
    end
end

% Get the resampled image
K = interpn(J, qx{1}, qx{2}, '*linear', 0);

% Compute the difference image
diff = (K - I(imask));

% Add up the voxels for total mean squared difference
diff2 = diff.^2;
f = sum(diff2(:));

% Compute the gradient
if nargout > 1
   
    idx = [1 2 1 2 1 2];
    qdx = [1 1 2 2 3 3];

    gradK = cell(2,1);
    for i = 1:2
        % Interpolate the gradient image and scale by the x coordinate

        gradK{i} = interpn(data.gradJ{i}, qx{1}, qx{2}, '*linear', 0);
        
    end
    
    % Compute the gradient of the objective
    g = zeros(6,1);
    for i = 1:6
        Gmap = diff .* gradK{idx(i)} .* px{qdx(i)};
        g(i) = 2.0 * sum(Gmap(:));
    end
    
end

% Compute the gradient
if nargout > 2
    
    % Interpolate the Hessian image at Ax + b
    hessK = cell(2,2);
    for i = 1:2
        for j = 1:2
            hessK{i,j} = interpn(data.hessJ{i,j}, qx{1}, qx{2}, '*linear', 0);
        end
    end 
    
    
    
    H = zeros(6);
    for i = 1:6
        Ti = gradK{idx(i)} .* px{qdx(i)};
        for j = 1:6
            Tj = gradK{idx(j)} .* px{qdx(j)};
            Hij = hessK{idx(i),idx(j)} .* px{qdx(i)} .* px{qdx(j)};
            Hmap = Ti .* Tj + diff .* Hij;
            H(i,j) = 2 * sum(Hmap(:));
            % H(i,j) = 2.0 * mean(mean(mean(Ti .* Tj + diff .* Hij)));
        end
    end
            
end