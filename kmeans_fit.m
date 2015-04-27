function [idx, c] = kmeans_fit(data_vec, k, scaling)
    X =  bsxfun(@times, data_vec.X, scaling);
    start = bsxfun(@times, data_vec.C_ave, scaling);
    [idx, c] = kmeans(X, k, 'Start', start);
end