function [idx, c] = kmeans_unsupervised_fit(data_vec, k, scaling)
    X = bsxfun(@times, data_vec.X, scaling);
    [idx, c] = kmeans(X, k);
end