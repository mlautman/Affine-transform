function [flow,labels] = graph_cut(data_vec, idx, labels)
    weight = 100;
    [C,D] = size(data_vec.C_ave);

    [height, width, depth] = size(data_vec.img);
    N = height*width;
    E = edges4connected(height,width);
    V = zeros(size(E,1),1);
    for i =1:depth
        im_i = squeeze(data_vec.img(:,:,1));
        V = V + abs(im_i(E(:,1)) - im_i(E(:,2))) + eps;
    end
    A = sparse(E(:,1), E(:,2), V, N, N, 4*N);
    
    
    ind_hit = sub2ind([height, width], ...
        data_vec.X(find(ismember(idx,labels)),D-1), ...
        data_vec.X(find(ismember(idx,labels)),D));
    ind_miss = sub2ind([height, width], ...
        data_vec.X(find(~ismember(idx,labels)),D-1), ...
        data_vec.X(find(~ismember(idx,labels)),D));
    
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
%     imagesc(labels); title('labels');
end