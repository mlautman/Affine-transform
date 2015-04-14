function [idx, c] = kmeans_fit(data_vec, k)
    [idx, c] = kmeans(data_vec.X, k, 'Start', data_vec.C_ave);
end