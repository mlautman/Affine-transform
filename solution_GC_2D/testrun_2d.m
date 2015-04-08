%% Test with registration between two datasets

clear
clc

% Use practice assignment stuff
addpath /Users/long/Box' Sync'/Study' Abroad'/PhD' Study'/PhD/TA/Spring2015/Grandchallenge2/solution_GC;
addpath /Users/long/Box' Sync'/Study' Abroad'/PhD' Study'/PhD/TA/Spring2015/Grandchallenge2/solution_GC/NIFTI_20110921;
addpath /Users/long/Box' Sync'/Study' Abroad'/PhD' Study'/PhD/TA/Spring2015/Grandchallenge2/solution_GC/solution_MATLAB_1;
addpath /Users/long/Box' Sync'/Study' Abroad'/PhD' Study'/PhD/TA/Spring2015/Grandchallenge2/2D_SELECTED/train;
cd /Users/long/Box' Sync'/Study' Abroad'/PhD' Study'/PhD/TA/Spring2015/Grandchallenge2/

% Read a pair of nifti files
[I1 sp1] = myReadNifti('sub006_mri.nii');
[I2 sp2] = myReadNifti('sub009_mri.nii');
I1 = squeeze(I1);
I2 = squeeze(I2);
sp1 = sp1(2:3);
sp2 = sp2(2:3);

% Load a segmentation image
[S sp1] = myReadNifti('sub006_seg.nii');
S = squeeze(S);
sp1 = sp1(2:3);

%% Histogram match I2 to I1
I2_hm = histogram_match(I1, I2);

%% Smooth the moving image and also compute its gradient

% Set sigma for the experiment
sigma = 1.0;

% Smooth I and J
Is = myGaussianLPF(I1, sigma);
Js = myGaussianLPF(I2_hm, sigma);

% Set up the problem
data = affine_precompute_with_hessian(Is,Js);


%% Test the Gradient and Hessian computation
eps = 0.001;
A0 = eye(2); b0 = zeros(2,1); x0 = [A0(:); b0(:)];
[f_center, g, H] = affine_objective_with_hessian(Is,Js,data,x0);
H_num = zeros(6,6);
g_num = zeros(6,1);

for i = 1:6
    x1 = x0; x1(i) = x0(i) + eps;
    x2 = x0; x2(i) = x0(i) - eps;
    [f1 g1] = affine_objective_with_hessian(Is,Js,data,x1);
    [f2 g2] = affine_objective_with_hessian(Is,Js,data,x2);
    
    g_num(i) = (f1 - f2) / (2 * eps);
    for j = 1:6
        H_num(i,j) = (g1(j) - g2(j)) / (2 * eps);
    end
end

clf;
colormap jet;
subplot(2,2,1); imagesc(H_num); axis image;
subplot(2,2,2); imagesc(H); axis image;
subplot(2,2,3); imagesc(g_num'); axis image;
subplot(2,2,4); imagesc(g'); axis image;


fprintf('Max relative error in Gradient: %f\n', ...
    max(abs(g(:) - g_num(:)) ./ (abs(g(:)))));

fprintf('Max relative error in Hessian: %f\n', ...
    max(abs(H(:) - H_num(:)) ./ (abs(H(:)))));


%% Run the optimization

% Define the objective function and initial solution
obj=@(x)(affine_objective_with_hessian(Is,Js,data,x));

% Set options for optimization
options = optimset(...
    'GradObj','on','Hessian','on','Display','iter',...
    'MaxIter',20, 'OutputFcn',@(x,ov,st)(myoptimplot(Is,Js,x,sp1)));

% Run optimization
tic
[xopt,fopt,eflag,output] = fminunc(obj, x0, options);
toc

%% Do the same, but using the multi-res strategy and histo matching

% Histogram match
[A b] = affine_multires(I1,I2_hm,sp1,[20 10 5],1.0);

% Look at the registration result
myViewInSNAP(I1, sp1)
myViewInSNAP(myTransformImage(I1, I2_hm, A, b), sp1)

%% Now set up a mask

% Expand the segmentation to create a mask. We do this by using mean 
% filtering and thresholding
Ms = myGaussianLPF(S > 0, 4) > 0.001;
myViewInSNAP(Ms .* I1,sp1);

%% Test the derivative computation for masked objective function
eps = 0.001;
A0 = eye(2); b0 = zeros(2,1); x0 = [A0(:); b0(:)];
[f_center g H] = affine_masked_objective_with_hessian(Is,Js,Ms,data,x0);
H_num = zeros(6,6);
g_num = zeros(6,1);

for i = 1:6
    x1 = x0; x1(i) = x0(i) + eps;
    x2 = x0; x2(i) = x0(i) - eps;
    [f1 g1] = affine_masked_objective_with_hessian(Is,Js,Ms,data,x1);
    [f2 g2] = affine_masked_objective_with_hessian(Is,Js,Ms,data,x2);
    
    g_num(i) = (f1 - f2) / (2 * eps);
    for j = 1:6
        H_num(i,j) = (g1(j) - g2(j)) / (2 * eps);
    end
end

clf;
colormap jet;
subplot(2,2,1); imagesc(H_num); axis image;
subplot(2,2,2); imagesc(H); axis image;
subplot(2,2,3); imagesc(g_num'); axis image;
subplot(2,2,4); imagesc(g'); axis image;

fprintf('Max relative error in Gradient: %f\n', ...
    max(abs(g(:) - g_num(:)) ./ (abs(g(:)))));

fprintf('Max relative error in Hessian: %f\n', ...
    max(abs(H(:) - H_num(:)) ./ (abs(H(:)))));


%% Run the optimization with mask

% Define the objective function and initial solution
obj=@(x)(affine_masked_objective_with_hessian(Is,Js,Ms,data,x));

% Set options for optimization
plt=@(x,ov,st)(myView(myTransformImage(Is,Js,reshape(x(1:4),2,2),x(5:6)),sp1));
options = optimset(...
    'GradObj','on','Hessian','on','Display','iter',...
    'MaxIter',20, 'OutputFcn',@(x,ov,st)(myoptimplot(Is,Js,x,sp1)));

% Run optimization
tic
[xopt,fopt,eflag,output] = fminunc(obj, x0, options);
toc

%% Run the multi-res strategy with a mask
[A b] = affine_multires_masked(I1,I2_hm,Ms,sp1,[20 20 10],1.0);

% Look at the registration result
myViewInSNAP(I1, sp1)
myViewInSNAP(myTransformImage(I1, I2_hm, A, b), sp1)
