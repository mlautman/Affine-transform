function img = gen_img_from_labels(data_vec, idx)
    x_v = size(data_vec.X,2);
    y_v = x_v - 1;
    max_h = max(data_vec.X(:,y_v));
    max_w = max(data_vec.X(:,x_v));
    img = zeros(max_h, max_w);
    for i = 1:length(unique(idx))
        pts = find(idx==i);
        ind_hit = sub2ind([max_h, max_w], ...
            data_vec.X(pts, y_v), ...
            data_vec.X(pts, x_v));
        img(ind_hit)=i;
    end
end