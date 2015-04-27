function X = my_icm2(im)

for i=1:3
    mu(i) = mean(im(find(im==i)));
%     sd(i) = std(im(find(im==i)));
end
%% K-Means Clustering

% For k-means we need to define initial means
% mu = [50 100 150];

% Iterate 
for j = 1:12
    
    % For each intensity value, compute square distance to each mean
    dst = zeros(size(im,1), size(im,2), 3);
    for i = 1:3
           dst(:,:,i) = (im - mu(i)).^2;
    end
    % Assign each non-background pixel to the closest cluster
    [dmin seg] = min(dst, [], 3);
    seg(find(im == 0)) = 0;
    % Compute the total variance
    var = sum(dmin(find(im > 0)));
    % Compute the new means
    for i = 1:3
        mu(i) = mean(im(find(seg==i)));
    end
    
end
%% Mixture Modeling

% Initialize alpha, mu and sigma for each class
al = [0.33 0.33 0.34];
mu = [1, 2, 3];
sg = [10 10 10];
%%%%%%%%%%%%%% Used mu from initial GMM input for better results
% Get the array of all non-zero voxels
nz = find(im > 0);
X = im(nz);

% Initialize the Pik array
Pik = zeros(length(X), 3);
Qik = Pik;

for it = 1:20
    % Perform E-step
    for k = 1:3
        Qik(:,k) = al(k) * normpdf(X, mu(k), sg(k)); 
    end
    for k = 1:3
        Pik(:,k) = Qik(:,k) ./ sum(Qik,2);
    end
    % Perform the M-step
    for k = 1:3
        al(k) = mean(Pik(:,k));
        mu(k) = sum(Pik(:,k) .* X) / sum(Pik(:,k));
        sg(k) = sqrt(...
            sum(Pik(:,k) .* (X - mu(k)) .* (X - mu(k))) / ...
            sum(Pik(:,k)));
    end
end


%% MRF: compute the initial segmentation from the GMM
% Change made: Initial input is based on GMM
X0 = im;
%% MRF: iterative conditional modes

% For each non-zero vertex, we need its neighborhood. We construct an 
% nz by 4 array that holds the four neighbors of each non-zero vertex
% or zero if the neighbor is not a non-zero vertex.
Idx=zeros(size(im));
Idx(:)=1:length(im(:));
IdxU=circshift(Idx,[1 0]);
IdxD=circshift(Idx,[-1 0]);
IdxR=circshift(Idx,[0 1]);
IdxL=circshift(Idx,[0 -1]);
N=[IdxU(nz) IdxD(nz) IdxR(nz) IdxL(nz)];
% For each non-zero node, number of non-zero neighbors
Nnnz=sum(im(N(:,:)) > 0,2);
% Get the -log posterior probability maps
LPik = -log(Pik);
% Set the lambda parameter
lambda = 0.5;

%% Perform brute-force iteration
X = X0;
for iter=1:20
    
    % Compute random permutation of the indices
    seq = randperm(length(nz));    
    % Keep track of number of flipped labels
    nflip = 0;    
    % Big loop (basic ICM is not parallelizable)
    for i = seq        
        % Get neighborhood information
        x_nbr = X(N(i,:));
        n_nbr = Nnnz(i);
        i_ctr = nz(i);        
        % For each possible label, compute likelihood
        log_lhd = LPik(i,:)';        
        % For each possible label, compute #nbrs with same label
        n_match = sum((ones(3,1) * x_nbr - [1:3]' * ones(1,4)) == 0,2);
        log_prior = lambda * (n_nbr - n_match);        
        % Get the conditional MAP estimate
        [dummy, map] = min(log_lhd + log_prior);        
        % Replace the current label with the map estimate
        if X(i_ctr) ~= map
            X(i_ctr) = map;
            nflip = nflip + 1;
        end
    end 
%     fprintf('Iter %02d: Flipped %d labels\n', iter, nflip);    
    if(nflip == 0)
        break
    end        
end

end

