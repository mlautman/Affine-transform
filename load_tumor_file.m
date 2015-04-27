function data = load_tumor_file(fid)
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

data.fid = fid;
if str2num(fid(4:length(fid)))<=15
    subdir = 'train';
    data.type = 'train';
else 
    subdir = 'test';
    data.type = 'test';
end


if isunix()
    dname = strcat('./data/brain-tumor/',subdir,'/', fid, '/');
else
    dname = strcat('.\data\brain-tumor\',subdir,'\', fid, '/');
end

% init struct 
fname = char(strcat(dname, fid, ext{1}));
[img, spacing] = myReadNifti(fname);
[n,d] = size(img);

data.ext=ext;
data.img = zeros(n, d, length(ext)-2);
data.seed = zeros(n, d);
data.seg = zeros(n, d);
data.spacing = spacing;

for i = 1:length(ext)-2
    fname = char(strcat(dname, fid, ext{i}));
    [data.img(:, :, i), ~] = myReadNifti(fname);  
end

% seed
fname = char(strcat(dname, fid, ext{i+1}));
[data.seed(:, :), ~] = myReadNifti(fname); 

% if isequal(data.type, 'train')
    fname = char(strcat(dname, fid, ext{i+2}));
    [data.seg(:, :), ~] = myReadNifti(fname);  
% end

end