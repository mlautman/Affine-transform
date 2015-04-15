function [flow,labels] = graph_cut(data_vec, idx, label)
    weight = 10;
    [C,D] = size(data_vec.C_ave);

    [height, width, ~] = size(data_vec.img);
    N = height*width;
    E = edges4connected(height,width);
    V = abs(data_vec.img(E(:,1)) - data_vec.img(E(:,2)) ) + eps;
    A = sparse(E(:,1), E(:,2), V, N, N, 4*N);
    
    
%     C1 = [round(mean(data_vec.X(find(idx==1),D))), round(mean(data_vec.X(find(idx==1),D-1)))];
%     C2 = [round(mean(data_vec.X(find(idx==2),D))), round(mean(data_vec.X(find(idx==2),D-1)))];
%     C3 = [round(mean(data_vec.X(find(idx==3),D))), round(mean(data_vec.X(find(idx==3),D-1)))];
% 
%     C1_ind = sub2ind(size(data_vec.img),C1(1),C1(2));
%     C2_ind = sub2ind(size(data_vec.img),C2(1),C2(2));
%     C3_ind = sub2ind(size(data_vec.img),C3(1),C3(2));


    ind_hit = sub2ind([height, width], ...
        data_vec.X(find(idx==label),D), ...
        data_vec.X(find(idx==label),D-1));
    ind_miss = sub2ind([height, width], ...
        data_vec.X(find(idx~=label),D), ...
        data_vec.X(find(idx~=label),D-1));
    
    con1 = [...                            
        ind_hit;...  % conn to src1 tumor
        ind_miss;...  % conn to src2 puss
    ]';
    
    con2 = [...
        1*ones(length(ind_hit),1); ... % src1 tumor
        2*ones(length(ind_miss),1); ... % src2 fluid
    ];

    T = sparse(...
        con1,...
        con2, ...  
        ones(length(idx),1)*weight ...
    );

    [flow,labels] = maxflow(A,T);
    labels = reshape(labels,[height width]);
    imagesc(labels); title('labels');
end