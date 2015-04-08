function myViewInSNAP(image,spacing)

% Generate temporary filename
fn=sprintf('%s.nii', tempname());

[nx,ny] = size(image);
image3D = zeros(1,nx,ny);
image3D(1,:,:) = image;
spacing = [1 spacing];

% Save image to temp image
myWriteNifti(fn, double(image3D), spacing);

% SNAP command
snapcmd='/Applications/ITK-SNAP.app/Contents/MacOS/ITK-SNAP';

% Command to execute
cmd=sprintf('%s -g %s&', snapcmd, fn);

% Run command
system(cmd);
