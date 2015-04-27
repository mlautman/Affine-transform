%%
clc; clear all; close all
%%
mask = 1:4;
weight = 1e2;

f_weights = [.1,.1, .1,.1,1,1];

%%
dirs = dir('./data/brain-tumor/train/')';
dirs = dirs(4:length(dirs));
%%
for d = dirs;
    dname = strcat('./data/brain-tumor/train',d.name,'/');
    %%
    disp(dname)
    data = load_tumor_file(d.name);
    data_vec = vectorize_tumor_data(data, mask);
    %%
    idx = GMM_fit(data_vec, f_weights);
    im = gen_img_from_labels(data_vec, idx);
    im_be = im;
    im_be(im==4)=0;
    im_be(im==3)=0;
    [avedist1, hausdist1, dice1] = seg_eval(data.seg, im_be);
    %% 
    im_new = tumor_graph_cut(data_vec, idx, weight);
    im_new_be = im_new;
    im_new_be(im_new==4)=0;
    im_new_be(im_new==3)=0;
    [avedist, hausdist, dice] = seg_eval(data.seg, im_new_be);
%%
    disp([avedist1, avedist, hausdist1, hausdist, dice1,dice]);


    %%
    figure(6); 
    subplot(2,2,1);
    imshow((data.seed)/max(max(data.seed)));
    subplot(2,2,2);
    imshow((data_vec.seg)/max(max(data_vec.seg)));
    subplot(2,2,3);
    imshow((im_be)/max(max(im_be)));
    subplot(2,2,4);
    imshow((im_new_be)/max(max(im_new_be)));
    pause(.1)
end

