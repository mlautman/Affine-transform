function [A b] = affine_multires_masked(I, J, M, sp, iter, sigma)
% affine_multires : perform multi-resolution affine registration
% usage:
%       [A b] = affine_multires(I, J, M, sp, iter, sigma)
% parameters:
%       I       Fixed image
%       J       Moving image
%       M       Binary mask image (in the space of image I)
%       sp      Spacing, used only for visualization
%       iter    Vector specifying # iterations at each resolution level
%       sigma   Gaussian smoothing sigma at highest resolution

% Optional parameters
if nargin < 5
    sigma = 1;
end

if nargin < 4
    iter = [20 10];
end

% Initialize the parameters to identity transform
x0 = zeros(6,1); x0(1:4) = eye(2);

% Number of resolutions
nres = length(iter);

% Main loop
for i = 1:nres
    
    % Smooth the image with appropriate sigma
    sub = 2^(nres-i);
    s = sigma * sub;
    
    if iter(i) > 0
    
        % Smooth I and J with increased sigma
        Ism = myGaussianLPF(I, s);
        Jsm = myGaussianLPF(J, s);
        Msm = myGaussianLPF(M, s);

        % Don`t use mask when subsample factor is bigger than 4
        if sub >=4
            Msm = ones(size(Ism));
        end
        
        % Resample to the lower resolution
        [nx ny nz] = size(Ism);
        if i < nres

            % Create a sampling grid
            [rx ry] = ndgrid(...
                linspace(1,nx,round(nx / sub)),...
                linspace(1,ny,round(ny / sub)));

            % Resample to lower resolution
            Isub = interpn(Ism,rx,ry,'*linear',0);
            Jsub = interpn(Jsm,rx,ry,'*linear',0);
            Msub = interpn(Msm,rx,ry,'*linear',0) > 0.5;

        else
            Isub = Ism;
            Jsub = Jsm;
            Msub = M;
        end

        fprintf('Registration at level %d, image size [%d %d %d]\n', i, size(Isub));

        % Set up the optimization
        data = affine_precompute_with_hessian(Isub, Jsub);

        % Define the objective function and initial solution
        obj=@(x)(affine_masked_objective_with_hessian(Isub,Jsub,Msub,data,x));

        % Set options for optimization
        %options = optimset(...
        %    'GradObj','on','Display','iter',...
        %    'MaxIter',iter(i), 'OutputFcn',@(x,ov,st)(myoptimplot(Isub,Jsub,x,sp)));

        options = optimset(...
            'GradObj','on','Hessian','on','Display','iter','MaxIter',iter(i),...
            'OutputFcn',@(x0,ov,st)(myoptimplot(Isub,Jsub,x0,sp)));
  
        %options = optimset(...
        %    'GradObj','on','Hessian','on','Display','iter','MaxIter',iter(i));
        
        % Run optimization
        tic
        xopt = fminunc(obj, x0, options);
        toc
        
        %myoptimplot(Isub, Jsub, xopt, sp);
        
    end
    
    % If not running at final resolution, we need to modify the affine 
    % parameters because the coordinate grid on which the image is defined
    % is going to change. In fact, all we have to change is to scale the
    % translation by the factor of two
    if i < nres
        xopt(5:6) = 2 * xopt(5:6);
        x0 = xopt;
    end
    
end

A = reshape(xopt(1:4),2,2);
b = xopt(5:6);

