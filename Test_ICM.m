hippo=[2,3,4,6,7,8,9,11,12,14,15,16,17,18,19,20];

GMM_ave=[];
GMM_haus=[];
GMM_dice=[];

for k=1:length(hippo)
    if hippo(k) <10
        dname = ['sub00' num2str(hippo(k))];
    else
        dname = ['sub0' num2str(hippo(k))];
    end
    display(dname)
    data = load_hippo_file(dname);
    data_vec = vectorize_hippo_data(data);

    % % GMM 
    idx = GMM_fit(data_vec);
    im = gen_img_from_labels(data_vec, idx);
    im(im>1)=0
    X = my_icm(im);
    [avedist, hausdist, dice] = seg_eval(data.seg, X);
    GMM_ave=[GMM_ave, avedist(:,2)];
    GMM_haus=[GMM_haus, hausdist(:,2)];
    GMM_dice=[GMM_dice, dice(:,2)];

end

GMM_hippo = [mean(GMM_dice(1,:)); ...
    mean(GMM_ave(1,:)); ...
    mean(GMM_haus(1,:))];

tumor=[1,2,3,4,5,6,8,9,11,12,13,14,15];

GMM_ave=[];
GMM_haus=[];
GMM_dice=[];

for k=1:length(tumor)
    if tumor(k) <10
        dname = ['.\data\brain-tumor\sub00' num2str(tumor(k)) '\'];
    else
        dname = ['.\data\brain-tumor\sub0' num2str(tumor(k)) '\'];
    end
    display(dname)
    data = load_tumor_file(dname);
    mask=1:4;
    data_vec = vectorize_tumor_data(data, mask);

    % % GMM 
    idx = GMM_fit(data_vec);
    im = gen_img_from_labels(data_vec, idx);
    im(im>2) = 0;
    X = my_icm(im);
    [avedist, hausdist, dice] = seg_eval(data.seg, X);
    GMM_ave=[GMM_ave, avedist(:,2)];
    GMM_haus=[GMM_haus, hausdist(:,2)];
    GMM_dice=[GMM_dice, dice(:,2)];

end

GMM_tumor = [mean(GMM_dice(1,2:end)), mean(GMM_dice(2,2:end)); ...
    mean(GMM_ave(1,2:end)), mean(GMM_ave(2,2:end)); ...
    mean(GMM_haus(1,2:end)), mean(GMM_haus(2,2:end)))];

