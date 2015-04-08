function stop = myoptimplot(I, J, x, sp)

% Get the matrix from x
A = reshape(x(1:4),2,2);
b = x(5:6);

K = myTransformImage(I, J, A, b);

imagesc(I - K);
colormap('jet');
getframe();

stop = false;