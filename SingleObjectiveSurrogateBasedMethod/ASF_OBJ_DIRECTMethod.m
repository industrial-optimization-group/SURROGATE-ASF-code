function [ASF_Val, fSample]=ASF_OBJ_DIRECTMethod(x,DataInfo)
x=x';
fSample=P_objective('value',DataInfo.Prob,DataInfo.NumObj,x);
%ASF evaluation for a ref and a set of x
TempRef = repmat(DataInfo.RefForASF,size(x,1),1);
TempW=repmat(DataInfo.ASFWeight,size(x,1),1);
ASF_Val = max(TempW .* (fSample - TempRef),[],2);

end