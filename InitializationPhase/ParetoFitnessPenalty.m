function pf = ParetoFitnessPenalty(fSample)

numSam = size(fSample,1);

pf = zeros(numSam,1);

k = linspace(1,numSam,numSam);

fMin = min(fSample,[],1);

fMax = max(fSample,[],1);

for sam = 1 : numSam
    tempK = k;
    tempK(sam) = [];
    other = tempK;
    fSam = fSample(sam,:);
    fSamNorm = (fSam - fMin) ./ (fMax - fMin);
    fOther = fSample(other,:);
    fMinTemp = ones(size(fOther,1),1) * fMin;
    fMaxTemp = ones(size(fOther,1),1) * fMax;
    fOtherNorm = (fOther - fMinTemp) ./ (fMaxTemp - fMinTemp);
    fSamNormTemp = ones(size(fOther,1),1) * fSamNorm;
    fDiffTemp = fSamNormTemp - fOtherNorm;
    minAll = min(fDiffTemp,[],2);
    pf(sam) = 1 - max(minAll);
end%sam
maxPf = max(pf);
pf(pf < 1) = (1-pf(pf < 1))+2*maxPf;
end%fun