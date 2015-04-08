function filtered = myMeanLPF(image, radius)

% Compute the kernel
w = 2 * radius + 1;
K = ones(w,w) / w^2;

% Compute the result
filtered = imfilter2d(image, K);
