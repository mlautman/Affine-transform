function plot_labeling(data_vec, idx)
    figure
    im = gen_img_from_labels(data_vec, idx);
    im = im/max(max(im));
    imagesc(im);
end