function plot_kmeans(data_vec, idx, mask)
    figure
    hold on 
    x_v = length(mask) + 2;
    y_v = x_v - 1;
    i = 1; scatter(...
        data_vec.X(find(idx==i),x_v), 240-data_vec.X(find(idx==i),y_v), 'r')
    i = 2; scatter(...
        data_vec.X(find(idx==i),x_v), 240-data_vec.X(find(idx==i),y_v), 'b')
    i = 3; scatter(...
        data_vec.X(find(idx==i),x_v), 240-data_vec.X(find(idx==i),y_v), 'g')
    i = 4; scatter(...
        data_vec.X(find(idx==i),x_v), 240-data_vec.X(find(idx==i),y_v), 'k')
end