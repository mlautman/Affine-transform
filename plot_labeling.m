function plot_labeling(data_vec, idx)
    figure
    hold on 
    x_v = size(data_vec.X,2);
    y_v = x_v - 1;
    max_h = max(data_vec.X(:,y_v));
    i = 1; 
    scatter(...
        data_vec.X(find(idx==i),x_v), data_vec.X(find(idx==i),y_v), 'r')
    i = 2; 
    scatter(...
        data_vec.X(find(idx==i),x_v), data_vec.X(find(idx==i),y_v), 'b')
    if length(unique(idx))>2
        i = 3; 
        scatter(...
            data_vec.X(find(idx==i),x_v), data_vec.X(find(idx==i),y_v), 'g')
    end
    if length(unique(idx))>3
        i = 4; 
        scatter(...
            data_vec.X(find(idx==i),x_v), data_vec.X(find(idx==i),y_v), 'k')
    end
end