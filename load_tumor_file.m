function data = load_tumor_file(dname)
% eg: 
% UNIX: data = load_tumor_file('./data/brain-tumor/sub001/')
% WinDOZE: data = load_tumor_file('.\data\brain-tumor\sub001\')
ext= {...
    '_flair.nii', ...
    '_t1.nii', ...
    '_t1c.nii', ...
    '_t2.nii', ...
    '_seed.nii',...
    '_seg.nii'...
};



if isunix()
    pth = strsplit(dname,'/');
    if dname(length(dname)) == '/'
        base = pth(length(pth)-1);
    else
        base = pth(length(pth));
    end
else
    pth = strsplit(dname,'\');
    if dname(length(dname)) == '\'
        base = pth(length(pth)-1);
    else
        base = pth(length(pth));
    end
end

% init struct 
fname = char(strcat(dname, base, ext{1}));
[img, spacing] = myReadNifti(fname);
[n,d] = size(img);

data.ext=ext;
data.img = zeros(n, d, length(ext)-2);
data.seed = zeros(n, d);
data.seg = zeros(n, d);
data.spacing = spacing;

for i = 1:length(ext)-2
    fname = char(strcat(dname, base, ext{i}));
    [data.img(:, :, i), ~] = myReadNifti(fname);  
end

% seed
fname = char(strcat(dname, base, ext{i+1}));
[data.seed(:, :), ~] = myReadNifti(fname); 

% seg
fname = char(strcat(dname, base, ext{i+2}));
[data.seg(:, :), ~] = myReadNifti(fname);  

end