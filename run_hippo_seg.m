%%
clear all; 
close all;
clc;
%%
train_fid = {...
   'sub002','sub012','sub003','sub014','sub004','sub015','sub006',...
    'sub016','sub007','sub017','sub008','sub018','sub009','sub019',... 
    'sub011','sub020'};
test_fid = {...
    'sub021', 'sub030', 'sub022', 'sub031', 'sub023', 'sub032', 'sub024',...
    'sub034', 'sub025', 'sub036', 'sub026', 'sub037', 'sub028', 'sub038',...
    'sub029', 'sub040',};

data_train = cell(length(train_fid),1);
for i = 1:length(train_fid)
    fid = train_fid{i};
    data_train{i} = load_hippo_file(fid);
end

data_test = cell(length(test_fid),1);
for i = 1:length(test_fid)
    fid = test_fid{i};
    data_test{i} = load_hippo_file(fid);
end

%%
% Loop over each test image and create a mask for that image by registering
% each training image to that test image and then using the test image's
% mask to add to the initial segmentation location.
sigma = 1;

for i = 1:length(test_fid)
    tic
    for j = 1:length(train_fid)
        tst = data_test{i}.img;
        tst_sp = data_test{i}.spacing(2:3);
        trn = data_train{j}.img;
        trn_sp = data_train{j}.spacing(2:3);
        trn_seg = data_train{j}.seg;
        % mask
        Ms = myGaussianLPF(trn_seg > 0, 4) > 0.001;
        % lpf
        tst = myGaussianLPF(tst, sigma);
        trn = myGaussianLPF(trn, sigma);
        data = affine_precompute_with_hessian(tst, trn);
        % Define the objective function and initial solution
        obj=@(x)(affine_masked_objective_with_hessian(tst,trn,Ms,data,x));

        % plot option for optimization
        plt=@(x,ov,st)(myView(myTransformImage(tst,trn,reshape(x(1:4),2,2),x(5:6)),trn_sp));
        % set option for optimization
        options = optimset(...
            'GradObj','on','Hessian','on','Display','iter',...
            'MaxIter',120);%, 'OutputFcn',@(x,ov,st)(myoptimplot(tst,trn,x,trn_sp)));

        A0 = eye(2); b0 = zeros(2,1); x0 = [A0(:); b0(:)];
        % Run optimization
        [t, xopt,fopt,eflag,output] = evalc('fminunc(obj, x0, options)');
        A = reshape(xopt(1:4), [2,2]);
        b = xopt(5:6);
        mask2 = myTransformImage(tst, trn_seg, A, b);
        %         subplot(2,2,1);  imagesc(trn)
        %         subplot(2,2,2);  imagesc(tst .* (mask2==0) );
        %         subplot(2,2,3);  imagesc(tst)
        data_test{i}.seg = data_test{i}.seg + mask2;
%         pause(1)
    end
    i
    bu = data_test{i}.seg;
    data_test{i}.seg = data_test{i}.seg - min(min(data_test{i}.seg));
    data_test{i}.seg = data_test{i}.seg./max(max(data_test{i}.seg));
    data_test{i}.seg = data_test{i}.seg > .5;

    figure(2)
    subplot(1,3,1); imagesc(data_test{i}.img)
    subplot(1,3,2); imagesc(data_test{i}.seg)
    max_v = max(max(data_test{i}.img));
    subplot(1,3,3); imagesc(data_test{i}.img + max_v.*data_test{i}.seg)
    pause(.1)
    toc
end

%% save test results 

if isunix()
    dname = './data/mri-hippocampus/test/';
else
    dname = '.\data\mri-hippocampus\test\';
end
for i = 1:length(test_fid)
    fpath = strcat(dname,data_test{i}.fid,'_seg.nii');
    myWriteNifti(fpath,double(data_test{i}.seg),data_test{i}.spacing);
end
