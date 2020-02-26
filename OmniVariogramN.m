function [variogramOmni] = OmniVariogramN(image)
%OMNIVARIOGRAMN Summary of this function goes here
%   Detailed explanation goes here

[vario, amount]=variogramN(image);
maxSize=max(size(image),[],1);
distImage=bwdist(padarray(1,maxSize-1,'both'));
variogramOmni=nan(1,max(maxSize));

for l=1:max(maxSize)
   mask=and(distImage<=l,distImage>(l-1));
   variogramOmni(l)=sum(vario(mask).*amount(mask))/sum(amount(mask));
end


end

