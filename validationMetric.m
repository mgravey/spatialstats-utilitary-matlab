function [output] = validationMetric(imlist,thresold,connectivity, percentile)
%VALIDATIONMETRIC Summary of this function goes here
%   Detailed explanation goes here
if(nargin<2)
    thresold=100;
end
if(nargin<3)
    connectivity=8;
end
if(nargin<4)
    percentile=[0];
end
if(~iscell(imlist))
   imlist={imlist}; 
end
maxValue=-inf;
minValue=inf;

imSize=[];

for i=1:length(imlist)
    locMax=max(imlist{i}(:));
    locMin=min(imlist{i}(:));
    if(locMax>maxValue)maxValue=locMax; end
    if(locMin<minValue)minValue=locMin; end
    imSize=cat(1,imSize,size(imlist{i}));
end
if(length(thresold)==1 && round(thresold)==thresold)
    thresold=linspace(minValue,maxValue,thresold);
end

eulerArray=nan(length(thresold),length(imlist));
eulerArrayC=nan(length(thresold),length(imlist));
connectivityArray=nan(length(thresold),length(imlist));
connectivityArrayC=nan(length(thresold),length(imlist));

for id=1:length(imlist)*length(thresold)
    [j,i]=ind2sub([length(thresold),length(imlist)],id);
    binIm=imlist{i} < thresold(j);
    eulerArray(id)=bweuler(binIm);
    geobodies = bwlabel(binIm, connectivity);
    numGeobodies = max(geobodies(:));
    connectivityArray(id) = sum(binIm(:)) / max(numGeobodies,1) / numel(binIm);
    binIm=~binIm;
    eulerArrayC(id)=bweuler(binIm);
    geobodies = bwlabel(binIm, connectivity);
    numGeobodies = max(geobodies(:));
    connectivityArrayC(id) = sum(binIm(:)) / max(numGeobodies,1) / numel(binIm);
end

output.eulerArray=eulerArray;
output.eulerArrayC=eulerArrayC;
output.connectivityArray=connectivityArray;
output.connectivityArrayC=connectivityArrayC;

% variogram

maxSize=max(imSize,[],1);
output.variogramX=nan(maxSize(2),length(imlist));
output.variogramY=nan(maxSize(1),length(imlist));
output.variogramOmni=nan(max(maxSize),length(imlist));
output.variogramPercentiles=nan(max(maxSize),length(imlist),length(percentile));

for i=1:length(imlist)
   [vario,amount]=variogram(imlist{i});
   output.variogramX(1:ceil(size(vario,2)/2),i)=vario(ceil(end/2),ceil(end/2):end);
   output.variogramY(1:ceil(size(vario,2)/2),i)=vario(ceil(end/2):end,ceil(end/2));
   distImage=bwdist(padarray(1,maxSize-1,'both'));
   for l=1:max(maxSize)
       mask=and(distImage<=l,distImage>(l-1));
       output.variogramOmni(l,i)=sum(vario(mask).*amount(mask))/sum(amount(mask));
       output.variogramPercentiles(l,i,:)=prctile(vario(mask),percentile);
   end
end

output.percentiles=percentile;
output.thresold=thresold;

end