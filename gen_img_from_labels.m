function img = gen_img_from_labels(data_vec, idx)
    x_v = size(data_vec.X,2);
    y_v = x_v - 1;
    max_h = max(data_vec.X(:,y_v));
    max_w = max(data_vec.X(:,x_v));
    img = zeros(max_h, max_w);
    for i = 1:length(unique(idx))
        pts = find(idx==i);
        for k = 1:length(pts)
            img(data_vec.X(pts(k),y_v), data_vec.X(pts(k),x_v))=i;
        end
    end
end