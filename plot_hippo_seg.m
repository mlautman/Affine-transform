function plot_hippo_seg(data)

% plot
figure(1)
subplot(1,2,1);
imshow(data.img / max(max(data.img)));    
title(data.ext{1}(2:find(data.ext{1}=='.')-1))

subplot(1,2,2);
imagesc(data.seg / max(max(data.seg)));    
title(data.ext{2}(2:find(data.ext{2}=='.')-1))

figure(2)
imagesc(data.img/max(max(data.img(:,:,1))));
hold on
[j,i] = find(data.seg == 1);
scatter(i,j, 'r*')

end



