function [avedist, hausdist, dice] = seg_eval(gtr, seg)

avedist=[];
hausdist=[];
dice = [];
% %  Need to loop over number of regions in img
for k = 1:length(unique(gtr))-1
    data1 = zeros(size(gtr));
    data1(find(gtr == k)) = 1;
    
    data2 = zeros(size(seg));
    data2(find(seg == k)) = 1;

    GSplus = data1+data2;
    GSintersect = find(GSplus == 2);
    dice =[dice; k, 2*length(GSintersect)/(sum(sum((data1))) + sum(sum(data2)))];
       
    data1 = edge(data1);
    data2 = edge(data2);
    
% %     Distance from set 1 to set 2
    dist12a = bwdist(data1);
    [i,j] = find(data2==1);
    dist12=[];
    for ijk = 1:length(i)
        dist12(ijk) = dist12a(i(ijk), j(ijk));
    end
    
% %     Distance from set 2 to set 1
    dist21a = bwdist(data2);
    [i,j] = find(data1==1);
    dist21=[];
    for ijk = 1:length(i)
        dist21(ijk) = dist21a(i(ijk), j(ijk));
    end
    
    avedist = [avedist; k, mean([mean(dist12), mean(dist21)])];
    hausdist = [hausdist; k, mean([max(dist12), max(dist21)])];
    
end

end
