function data = load_hippo_file(fid)
% eg: load_hippo_file('sub012')

ext= {...
    '_mri.nii', ...
    '_seg.nii'...
};
data.fid = fid;

if str2num(fid(4:length(fid)))<=20
    data.type = 'train';
else 
    data.type = 'test';
end

if isunix()
    dname = strcat('./data/mri-hippocampus/',data.type,'/');
else
    dname = strcat('.\data\mri-hippocampus\',data.type,'\');
end

% init struct 
fpath = char(strcat(dname, fid, ext{1}));
[img, spacing] = myReadNifti(fpath);
[~, n, d] = size(img);
data.ext=ext;
data.img = zeros(n,d);
data.seg = zeros(n,d);
data.spacing = spacing;

fpath = char(strcat(dname, fid, ext{1}));
[data.img(:,:), ~] = myReadNifti(fpath);  

if isequal(data.type, 'train')
    fpath = char(strcat(dname, fid, ext{2}));
    [data.seg(:,:), ~] = myReadNifti(fpath);  
else
    data.seg(:,:) = zeros(size(data.img));  
end

end