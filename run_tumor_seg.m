clear all; 
close all;

%%
mask = 1:4;
weight = 5e1;

f_weights = [.2,.1, .2,.1,1,1];
%%
train_fid = {...
    'sub001', 'sub004', 'sub008', 'sub012', 'sub015', 'sub002', 'sub005',...
    'sub009', 'sub013', 'sub003', 'sub006', 'sub011', 'sub014'};
test_fid = {...
    'sub016', 'sub020', 'sub023', 'sub026', 'sub029', 'sub017', 'sub021',...
    'sub024', 'sub027', 'sub018', 'sub022', 'sub025', 'sub028'};

data_train = cell(length(train_fid),1);
data_vec_train = cell(length(train_fid),1);
for i = 1:length(train_fid)
    fid = train_fid{i};
    data_train{i} = load_tumor_file(fid);
    data_vec_train{i} = vectorize_tumor_data(data_train{i}, mask);
end

%%
data_test = cell(length(test_fid),1);
data_vec_test = cell(length(test_fid),1);
%% test
for i = 1:length(test_fid)
    fid = test_fid{i};
    data_test{i} = load_tumor_file(fid);
    data_vec_test{i} = vectorize_tumor_data(data_test{i}, mask);
    data_vec_test{i}.idx = GMM_fit(data_vec_test{i}, f_weights);
    im = gen_img_from_labels(data_vec_test{i}, data_vec_test{i}.idx);
    im_be = im;
    im_be(im==4)=0;
    im_be(im==3)=0;
    data_test{i}.seg = im_be;
end
%% 
for i =1:length(test_fid)

    im_new = tumor_graph_cut(data_vec_test{i}, data_vec_test{i}.idx, weight);
    im_new_be = im_new;
    im_new_be(im_new==4)=0;
    im_new_be(im_new==3)=0;
    data_test{i}.seg = im_new_be;
end
%% save test results 

if isunix()
    dname = './data/brain-tumor/test/';
    sep = '/';
else
    dname = '.\data\brain-tumor\test\';
    sep = '\';
end
for i =1:length(test_fid)
    i
    imshow(data_test{i}.seg./max(max(data_test{i}.seg)));
    pause(3)
    fpath = strcat(dname,'results',sep,data_test{i}.fid,'_seg.nii');
    myWriteNifti(fpath,double(data_test{i}.seg),data_test{i}.spacing);
end
