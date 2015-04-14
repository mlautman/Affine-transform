function data = load_hippo_file(fid)
ext= {...
    '_mri.nii', ...
    '_seg.nii'...
};

if isunix()
    dname = './data/mri-hippocampus/';
else
    dname = '.\data\mri-hippocampus\';
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

fpath = char(strcat(dname, fid, ext{2}));
[data.seg(:,:), ~] = myReadNifti(fpath);  

end