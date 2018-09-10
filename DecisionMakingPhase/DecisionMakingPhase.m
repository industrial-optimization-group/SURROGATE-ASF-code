function DataInfo=DecisionMakingPhase(DataInfo)
askDM=1;
DMReferencesHistory=Inf*ones(DataInfo.MaxFunEvlForDM+1,DataInfo.NumObj);
PreferredSolutionInObjectiveSpaceHistory=Inf*ones(DataInfo.MaxFunEvlForDM+1,DataInfo.NumObj);
PreferredSolutionInDecisionSpaceHistory=Inf*ones(DataInfo.MaxFunEvlForDM+1,DataInfo.NumVar);
Count=1;
clc;
disp(['You can only provide at most ' num2str(DataInfo.MaxFunEvlForDM) ' reference points'])
while askDM   
    %%Give ranges for reference points     
    disp('Explorable ranges in the objective space for objective functions: ');
    for Obj=1:DataInfo.NumObj        
        disp([num2str(DataInfo.Ideal(Obj)) ' <= f' num2str(Obj) ' <= ' num2str(DataInfo.Nadir(Obj))])
    end   
    InRef = input(['Input the reference point number ' num2str(Count) ' within the explorable ranges with a dimension 1 x ' num2str(DataInfo.NumObj) ' ']);
    if (length(InRef) ~= DataInfo.NumObj) 
        disp(['The reference point must have a dimension 1 x ' num2str(DataInfo.NumObj) '.'])
    elseif  any(DataInfo.Ideal > InRef) || ...
            any(InRef > DataInfo.Nadir)
        disp('The componet of the reference point must be within the explorable ranges.')
    else
        DMReferencesHistory(Count,:)=InRef;
        ZProj=FindProjectedReference(DMReferencesHistory(Count,:),DataInfo);
        RegionInd=FindSubregion(ZProj,DataInfo);
        [PreferredSolutionInDecisionSpaceHistory(Count,:), PreferredSolutionInObjectiveSpaceHistory(Count,:)]=OptimizeSurrogate(RegionInd,ZProj,DataInfo);
        disp(['Given reference point = [' num2str(DMReferencesHistory(Count,:)) '], and '])
        disp(['corresponding prefered solution in the objecive space = [' num2str(PreferredSolutionInObjectiveSpaceHistory(Count,:)) '].'])        
        if Count >= DataInfo.MaxFunEvlForDM
            disp('Maximum number of function evaluations is reached. Decision making phase is terminated.')        
            Tag=0;
        else
            if Count+1 == DataInfo.MaxFunEvlForDM
                disp('You can only provide one more reference point.')
            end
            Tag=input('Would you like to provide another reference point? 1/0 ');
        end
        if Tag == 0
            askDM=0;
        else            
            Count=Count+1;
        end
    end
end
DMReferencesHistory(Count+1:end,:)=[];
PreferredSolutionInObjectiveSpaceHistory(Count+1:end,:)=[];
PreferredSolutionInDecisionSpaceHistory(Count+1:end,:)=[];
TheMostPreferredSolutionInDecisionSpace= PreferredSolutionInDecisionSpaceHistory(Count,:);
TheMostPreferredSolutionInObjectiveSpace= PreferredSolutionInObjectiveSpaceHistory(Count,:);
DataInfo.DMReferencesHistory=DMReferencesHistory;
DataInfo.PreferredSolutionInObjectiveSpaceHistory=PreferredSolutionInObjectiveSpaceHistory;
DataInfo.PreferredSolutionInDecisionSpaceHistory=PreferredSolutionInDecisionSpaceHistory;
DataInfo.TheMostPreferredSolutionInDecisionSpace=TheMostPreferredSolutionInDecisionSpace;
DataInfo.TheMostPreferredSolutionInObjectiveSpace=TheMostPreferredSolutionInObjectiveSpace;
end

function ZProj=FindProjectedReference(DMref,DataInfo)

[~, RefTemp]=ReferenceVectorGenerator(300,0,DataInfo.NumObj);

RefTemp=MapSamples(RefTemp,[min(DataInfo.ExtremesObj); ... 
    max(DataInfo.ExtremesObj)],[zeros(1,DataInfo.NumObj);ones(1,DataInfo.NumObj)]);

[~,IndProj]=min(sqrt(sum(((repmat(DMref,size(RefTemp,1),1)-RefTemp) ./ ... 
    repmat(DataInfo.ASFWeight,size(RefTemp,1),1)).^2,2)));
ZProj=RefTemp(IndProj,:);
end

function RegionInd=FindSubregion(ZProj,DataInfo)
NormalizedZProj=(ZProj-DataInfo.Ideal) ./ (DataInfo.Nadir-DataInfo.Ideal);
[DistVal,IndRef]=sort(pdist2(NormalizedZProj, ... 
    DataInfo.SubRegionNumberingReferencePoints),'ascend');
if min(DistVal) > 0.0001
    IndRef=sort(IndRef(1:size(DataInfo.ExtremesObj,1)));
    RegionInd=find(ismember(DataInfo.HyperBoxesSubRegionExtremesIndices(:, ... 
        1:size(DataInfo.ExtremesObj,1)),IndRef,'rows')==1);
else
    PredRefInd=find(ismember(DataInfo.PredeterminedReferencePoints,ZProj,'rows')==1);
    [TempInd,~]=find(DataInfo.HyperBoxesSubRegionExtremesIndices==PredRefInd);
    RegionInd=TempInd(1);
end
end

function [xOpt, fOpt]=OptimizeSurrogate(RegionInd,ZProj,DataInfo)
opts.es = 1e-4;
opts.maxevals = 500;
opts.maxits = 500;
opts.maxdeep = 1000;
opts.testflag = 0;
opts.showits = 0;
Problem.f ='SurrogatePrediction';
Struc.Model = DataInfo.SurrogateInfo{RegionInd}.Surrogate;
Struc.SurrogateType=DataInfo.SurrogateType;
Struc.ZRef = ZProj;
bound = [DataInfo.SurrogateInfo{RegionInd}.lb' DataInfo.SurrogateInfo{RegionInd}.ub'];
[~, xOpt] = Direct(Problem,bound,opts,Struc);
xOpt=xOpt';
fOpt = P_objective('value',DataInfo.Prob,DataInfo.NumObj,xOpt);

end