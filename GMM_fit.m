function idx = GMM_fit(data_vec, f_weights)
%     model = fitgmdist(data_vec.X,4,'Start', data_vec.Y);
    X = bsxfun(@times, data_vec.X, f_weights);  
    mu = bsxfun(@times, data_vec.C_ave, f_weights);  
    sd = bsxfun(@times, data_vec.C_std, f_weights);  
    
    D = size(sd,2);
    C = size(sd,1);

    % For each voxel, compute the probability assigned by each pdf
    pl=zeros(size(X,1), C);
    for i=1:C
        m = mu(i,:);
        s = diag(sd(i,:));
        s_inv = s\eye(D);
        coeff = -1/2 * log(det(s)) - D/2 * log(2*pi);
        for j = 1:size(data_vec.X,1)
            delta = X(j,:) - m;
            pl(j,i) = coeff - 1/2 * delta * s_inv * delta';
        end
    end
    
    [m,idx] = max(pl,[],2);
    
    
end