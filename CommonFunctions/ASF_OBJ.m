function [ASF_Val, fSample]=ASF_OBJ(x,Ref,DataInfo,Tag)
fSample=P_objective('value',DataInfo.Prob,DataInfo.NumObj,x);

if Tag==0%%ASF evaluation for a ref and a set of x    
    TempRef = repmat(Ref,size(x,1),1);
    TempW=repmat(DataInfo.ASFWeight,size(x,1),1);
    ASF_Val = max(TempW .* (fSample - TempRef),[],2);
else %%ASF evaluation for samples in obj and a set of ref
    fSampleTemp=repmat(fSample,size(Ref,1),1);
    TempW=repmat(DataInfo.ASFWeight,size(Ref,1),1);
    ASF_Val = max(TempW .* (fSampleTemp - Ref),[],2);
end

end