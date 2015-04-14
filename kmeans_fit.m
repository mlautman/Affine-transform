function [idx, c] = kmeans_fit(data_vec)
    [idx, c] = kmeans(data_vec.X, 4, 'Start', data_vec.C_ave);
end