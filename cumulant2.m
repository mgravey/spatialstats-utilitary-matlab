function [result] = cumulant2(image,maxTemp)

image=image-mean(image(:));

listOfLagVector={};

% create all lag vector
dims=[];
for j=1:size(maxTemp,1)
    dims=[dims norm(maxTemp(j,:),inf)+1];
    listOfLagVector{j}=[0 0];
    for i=1:norm(maxTemp(j,:),inf)
      listOfLagVector{j}=[listOfLagVector{j}; -floor(maxTemp(j,:)./norm(maxTemp(j,:),inf)*i)];
   end
end

% create all shifted map
for j=1:size(maxTemp,1)
    imageStack{j}=nan([size(image) size(listOfLagVector{j},1)]);
    for i=1:size(listOfLagVector{j})
        imageStack{j}(:,:,i)=circshift(image,listOfLagVector{j}(i,:));
        imageStack{j}(end-listOfLagVector{j}(i,1):end,:,i)=nan;
        imageStack{j}(:,end-listOfLagVector{j}(i,2):end,i)=nan;
    end
end

result=nan(dims);
size(result)


if(size(maxTemp,1)==3)
    [v1,v2,v3]=ind2sub(size(result),1:numel(result));
    v=([v1; v2; v3])';
    for l=1:numel(result)
       result(l)=nanmean(reshape(prod(cat(3,imageStack{1}(:,:,v(l,1)),imageStack{2}(:,:,v(l,2)),imageStack{3}(:,:,v(l,3))),3),1,[]));
    end
end

if(size(maxTemp,1)==4)
    [v1,v2,v3,v4]=ind2sub(size(result),1:numel(result));
    v=([v1; v2; v3; v4])';
    for l=1:numel(result)
       result(l)=nanmean(reshape(prod(cat(3,imageStack{1}(:,:,v(l,1)),imageStack{2}(:,:,v(l,2)),imageStack{3}(:,:,v(l,3)),imageStack{4}(:,:,v(l,4))),3),1,[]));
       -nanmean(reshape(prod(cat(3,imageStack{1}(:,:,v(l,1)),imageStack{2}(:,:,v(l,2))),3),1,[]))*nanmean(reshape(prod(cat(3,imageStack{3}(:,:,v(l,3)),imageStack{4}(:,:,v(l,4))),3),1,[]));
       -nanmean(reshape(prod(cat(3,imageStack{1}(:,:,v(l,1)),imageStack{3}(:,:,v(l,3))),3),1,[]))*nanmean(reshape(prod(cat(3,imageStack{2}(:,:,v(l,2)),imageStack{4}(:,:,v(l,4))),3),1,[]));
       -nanmean(reshape(prod(cat(3,imageStack{1}(:,:,v(l,1)),imageStack{4}(:,:,v(l,4))),3),1,[]))*nanmean(reshape(prod(cat(3,imageStack{3}(:,:,v(l,3)),imageStack{2}(:,:,v(l,2))),3),1,[]));
    end
end

result=squeeze(result);

end

