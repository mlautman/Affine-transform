function plot_tumor_seg(data)
% close all;
% plot
figure(1)
for i = 1:length(data.ext)-2
    subplot(3,2,i);
    imshow(...
        reshape(data.img(:,:,i), [240 240]) / ...
        max(max(max(data.img(:,:,i)))));    
    title(data.ext{i}(2:find(data.ext{i}=='.')-1))
end   

i = i+1;
subplot(3,2,i);
imshow(...
    reshape(data.seed(:,:), [240 240]) / ...
    max(max(max(data.seed(:,:)))));    
title(data.ext{i}(2:find(data.ext{i}=='.')-1))

i = i+1;
subplot(3,2,i);
imshow(...
    reshape(data.seg(:,:), [240 240]) / ...
    max(max(max(data.seg(:,:)))));    
title(data.ext{i}(2:find(data.ext{i}=='.')-1))

% figure(2)
% imshow(...
%     reshape(data.img(:,:,1), [240 240]) / ...
%     max(max(max(data.img(:,:,1)))));
% hold on
% [j,i] = find(data.seg == 1);
% scatter(i,j, 'r.')
% [j,i] = find(data.seg == 2);
% scatter(i,j, 'g.')
% [j,i] = find(data.seg == 3);
% scatter(i,j, 'b.')
% 
% [j,i] = find(data.seed == 1);
% scatter(i,j, 'm.')
% [j,i] = find(data.seed == 2);
% scatter(i,j, 'y')
% [j,i] = find(data.seed == 3);
% scatter(i,j, 'b')

end



