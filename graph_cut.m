function [flow,labels] = graph_cut(data_vec, idx)
    
    [C,D] = size(data_vec.C_ave);

    [height, width] = size(data_vec.img);
    N = height*width;
    E = edges4connected(height,width);
    V = abs(data_vec.img(E(:,1)) - data_vec.img(E(:,2)) ) + eps;
    A = sparse(E(:,1),E(:,2),V,N,N,4*N);
    
    C1 = [round(mean(data_vec.X(find(idx==1),5))), round(mean(data_vec.X(find(idx==1),4)))];
    C2 = [round(mean(data_vec.X(find(idx==2),5))), round(mean(data_vec.X(find(idx==2),4)))];
    C3 = [round(mean(data_vec.X(find(idx==3),5))), round(mean(data_vec.X(find(idx==3),4)))];

    C1_ind = sub2ind(size(data_vec.img),C1(1),C1(2));
    C2_ind = sub2ind(size(data_vec.img),C2(1),C2(2));
    C3_ind = sub2ind(size(data_vec.img),C3(1),C3(2));

    
    con1 = [...                            
            find(data_vec.seed==1);...  % conn to src1 tumor
            find(data_vec.seed==2);...  % conn to src2 puss
        ]';
    
    con2 = [...
            C1_ind*ones(length(find(data_vec.seed==1)),1); ... % src1 tumor
            C2_ind*ones(length(find(data_vec.seed==2)),1); ... % src2 fluid
        ];
    T = sparse(con1,con2, ...  
        ones(...
            length(find(data_vec.seed==1))+length(find(data_vec.seed==2)),1 ...
        )*9e9 ...
    );

    [flow,labels] = maxflow(A,T);
    labels = reshape(labels,[height width]);
    imagesc(labels); title('labels');
end