function data_vec = vectorize_hippo_data(data)
    
    % hippo seed values
    [j_h,i_h] = find(data.seg == 1);    
    h = zeros(length(j_h), 3);
    for i = 1:length(j_h)
        h(i,:) = [...
            data.img(j_h(i), i_h(i)),...
            j_h(i), ...
            i_h(i)...
        ];
    end
    
    % brain seed values
    [j_b,i_b] = find(data.seg == 0 & data.img~=0);
    b = zeros(length(j_b), 3);
    for i = 1:length(j_b)
        b(i,:) = [...
            data.img(j_b(i), i_b(i)),...
            j_b(i), ...
            i_b(i)...
        ];
    end
    
    % unlabled values
    [j_r,i_r] = find(data.img == 0);
    r = zeros(length(j_r), 3);
    for i = 1:length(j_r)
        r(i,:) = [...
            data.img(j_r(i), i_r(i)),...
            j_r(i), ...
            i_r(i)...
        ];
    end
    
    data_vec=struct();
    data_vec.X = [h; b; r];    
    data_vec.Y = [...
		1 * ones(size(h, 1), 1); ...
		2 * ones(size(b, 1), 1); ...
		3 * ones(size(r, 1), 1); ...
	];
    mean_r = [zeros(1, size(r,2)-2), mean(r(:,size(r, 2)-1)), mean(r(:,size(r, 2)))];
    std_r = [ones(1, size(r,2)-2), std(r(:,size(r, 2)-1)), std(r(:,size(r, 2)))];

    data_vec.C_ave = [ mean(h,1); mean(b,1); mean_r];
    data_vec.C_std = [ std(h,1); std(b,1); std_r];
    data_vec.img = data.img;
    data_vec.seg = data.seg;
    data_vec.spacing = data.spacing;
    data_vec.ext = data.ext;
end
    
    