function [result, ind, arch] = IsCloseToSample(xOpt,xSamp,arch)
[disClose, ind] = min(pdist2(xOpt,xSamp),[],2);%Closest to candidate
sampleTemp = xSamp;
sampleTemp(ind,:) = [];
[clsDis, ~] = min(pdist2(xSamp(ind,:),sampleTemp),[],2);% Closest to closest
radius = clsDis/2;
diff = disClose - radius;
if (diff < 0) %|| (isMem == 1)
    result = 1;
else
    arch = [arch; xSamp(ind,:)]; 
    result = 0;
end
end%function