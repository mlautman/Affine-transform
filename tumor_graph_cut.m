function im_new = tumor_graph_cut(data_vec, idx)
    labels = unique(idx);
    im = gen_img_from_labels(data_vec, idx);
    im_new = zeros(size(im));
    for i = 1:length(labels)-1
        e = length(labels) - i;
        l = labels(1:e);
        [~, im_cut] = graph_cut(data_vec,idx, l);
        im_new(find(im_cut~=1)) = l(e);
    end
end