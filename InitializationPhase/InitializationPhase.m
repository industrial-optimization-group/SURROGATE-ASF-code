function DataInfo=InitializationPhase(DataInfo)
DataInfo=FindInitialExtremes(DataInfo);
DataInfo=FindPredeterminedReferencePointsAndNonDominatedSolutions(DataInfo);
DataInfo=FindHyperBoxesSubResionsExtremesIndices(DataInfo);
DataInfo=BuildSurrogates(DataInfo);
delete('Results.Mat');
end

function DataInfo=FindInitialExtremes(DataInfo)
ExtTemp=zeros(DataInfo.NumObj,DataInfo.NumObj);
for i =1:DataInfo.NumObj
    ExtTemp(i,i)=DataInfo.InitialRanges(i,2);
    IndTemp=1:DataInfo.NumObj;
    IndTemp(i)=[];
    ExtTemp(i,IndTemp)=DataInfo.InitialRanges(IndTemp,1);
end
DataInfo.InitialExtremes=ExtTemp;
end


function DataInfo=GenerateUniformReferencePoints(DataInfo)
if strcmp(DataInfo.RefTypeGeneration,'Predetermined')
    [~,DataInfo.SubRegionNumberingReferencePoints]= ... 
        ReferenceVectorGenerator(DataInfo.NumDivPreRef,0,DataInfo.NumObj);      
    DataInfo.PredeterminedReferencePoints=MapSamples( ...
        DataInfo.SubRegionNumberingReferencePoints,[DataInfo.Ranges(:,1) DataInfo.Ranges(:,2)]', ...
        [zeros(1,DataInfo.NumObj); ones(1,DataInfo.NumObj)]);
elseif strcomp(DataInfo.RefTypeGeneration,'ReferenceSamples')
    ...
end
end

function DataInfo=FindPredeterminedReferencePointsAndNonDominatedSolutions(DataInfo)
t1=cputime;
%%%Frist: adjusting extermes in the objective space
ExtremTempObj=zeros(DataInfo.NumObj,DataInfo.NumObj);
ExtremTempDec=zeros(DataInfo.NumObj,DataInfo.NumVar);
for Ext=1:DataInfo.NumObj
    DataInfo.RefForASF=DataInfo.InitialExtremes(Ext,:);
    [ExtremTempDec(Ext,:),ExtremTempObj(Ext,:)]= ... 
        SolveComputationallyExpensiveSingleObjective(DataInfo,'ASF');
end
if size(unique(ExtremTempObj,'rows'),1)~=size(ExtremTempObj,1)
    TempObj=floor(rand(size(ExtremTempObj))* 10^(4))/10^(8);
    ExtremTempObj = ExtremTempObj+TempObj;
    TempDec=floor(rand(size(ExtremTempDec))* 10^(4))/10^(8);
    ExtremTempDec=ExtremTempDec+TempDec;    
end
DataInfo.ExtremesObj=ExtremTempObj;
DataInfo.ExtremesDec=ExtremTempDec;
DataInfo.Ideal=min(ExtremTempObj);
DataInfo.Nadir=max(ExtremTempObj);
DataInfo.Ranges=[DataInfo.Ideal;DataInfo.Nadir]';
DataInfo.ASFWeight = 1 ./ (DataInfo.Nadir-DataInfo.Ideal);
%%% Generating predetermined reference points
DataInfo.RefTypeGeneration='Predetermined';
DataInfo=GenerateUniformReferencePoints(DataInfo);

%%% Finding nondominated solutions corresponding to
%%% predetermined reference points
NondominatedSolObjTemp=zeros(size(DataInfo.PredeterminedReferencePoints));
NondominatedSolDecTemp=zeros(size(DataInfo.PredeterminedReferencePoints,1) ...
    ,DataInfo.NumVar);
IndExt=zeros(DataInfo.NumObj,1);
for Ext=1:DataInfo.NumObj
    [~,IndExt(Ext)]=min(pdist2(DataInfo.ExtremesObj(Ext,:), ... 
        DataInfo.PredeterminedReferencePoints));
end

NondominatedSolObjTemp(IndExt,:)=DataInfo.PredeterminedReferencePoints(IndExt,:);
NondominatedSolDecTemp(IndExt,:)=DataInfo.ExtremesDec;
for PreDetRef=1:size(DataInfo.PredeterminedReferencePoints,1)
    if ~ismember(PreDetRef,IndExt)
        DataInfo.RefForASF=DataInfo.PredeterminedReferencePoints(PreDetRef,:);
        [NondominatedSolDecTemp(PreDetRef,:),NondominatedSolObjTemp(PreDetRef,:)]= ...
            SolveComputationallyExpensiveSingleObjective(DataInfo,'ASF');
    end
end
DataInfo.NondominatedSolDec=NondominatedSolDecTemp;
DataInfo.NondominatedSolObj=NondominatedSolObjTemp;
t2=cputime-t1;
disp(t2)

end



function [xBest,fBest]=SolveComputationallyExpensiveSingleObjective(DataInfo,ProblemName)
if strcmp(ProblemName,'ASF')
    data_file = 'datainput_SURROGATEASF'; %name of DataInfo file
else
    data_file = 'datainput_EachObjective';
end
Solver=DataInfo.SingleSurrogateBasedSolver;
DataInfo.Tag = 0;
save DataInfo DataInfo
if strcmp(Solver,'Mat')
    maxeval = DataInfo.MaxFunEvlForEachPreRef; %maximum number of allowed function evaluations
    surogate_model = 'RBFcub'; %selected response surface
    sampling_technique = 'CANDloc'; %global randomized sampling strategy
    initial_design = 'LHS'; %symmetric Latin hypercube design as initial design
    number_startpoints = []; %default number of points for initial experimental design
    starting_point = []; %no user-specified points to be added to the initial design
    NumberNewSamples = 1; %1 new point is selected in each iteration    
    [xBest,~] = MATSuMoTo(data_file,maxeval,surogate_model,sampling_technique,...
        initial_design,number_startpoints,starting_point,NumberNewSamples);
elseif strcmp(Solver,'StochasticRBF')    
    Solutions=StochasticRBF(data_file,DataInfo.MaxFunEvlForEachPreRef,1,0,1);
    xBest=Solutions.BestPoints;
elseif strcmp(Solver,'DIRECT')
    opts.es = 1e-4;
    opts.maxevals = DataInfo.MaxFunEvlForEachPreRef;
    opts.maxits = 500;
    opts.maxdeep = 1000;
    opts.testflag = 0;
    opts.showits = 0;
    Problem.f ='ASF_OBJ_DIRECTMethod';
    bound = [DataInfo.LowerBounds' DataInfo.UpperBounds'];
    [~, xBest] = Direct(Problem,bound,opts,DataInfo);
    xBest=xBest';
end
fBest=P_objective('value',DataInfo.Prob,DataInfo.NumObj,xBest);
end

function DataInfo=FindHyperBoxesSubResionsExtremesIndices(DataInfo)
%%%Form subregions
HyperBoxesSubRegionExtremesIndicesTemp=[];
LevelsTemp=sort(unique(DataInfo.SubRegionNumberingReferencePoints(:,DataInfo.NumObj)),'descend');

for l=1:length(LevelsTemp)
    RefTempInd=find((DataInfo.SubRegionNumberingReferencePoints(:,DataInfo.NumObj)==LevelsTemp(l))==1);
    RefCanDonw=find((DataInfo.SubRegionNumberingReferencePoints(:,DataInfo.NumObj) < LevelsTemp(l))==1)';
    RefCanUp=find((DataInfo.SubRegionNumberingReferencePoints(:,DataInfo.NumObj) > LevelsTemp(l))==1)';
    for RefInd=1:length(RefTempInd)
        if length(RefCanDonw) > 1 %%Finding down extrems
            DistTempDown=pdist2(DataInfo.SubRegionNumberingReferencePoints(RefTempInd(RefInd),:), ...
                DataInfo.SubRegionNumberingReferencePoints(RefCanDonw,:));
            [~, DistTempIndDown]=sort(DistTempDown);
            HyperBoxesSubRegionExtremesIndicesTempDown=sort([RefTempInd(RefInd) ...
                RefCanDonw(DistTempIndDown(1:size(DataInfo.ExtremesObj,1)-1))]);
            HyperBoxesSubRegionExtremesIndicesTemp= ... 
                [HyperBoxesSubRegionExtremesIndicesTemp;HyperBoxesSubRegionExtremesIndicesTempDown];
        end
        if length(RefCanUp) > 1
            DistTempUp=pdist2(DataInfo.SubRegionNumberingReferencePoints(RefTempInd(RefInd),:), ...
                DataInfo.SubRegionNumberingReferencePoints(RefCanUp,:));
            [~, DistTempIndUp]=sort(DistTempUp);
            RefUpTempInd=RefCanUp(DistTempIndUp(1:size(DataInfo.ExtremesObj,1)-1));
            VectT=DataInfo.SubRegionNumberingReferencePoints(RefUpTempInd,:)- ...
                repmat(DataInfo.SubRegionNumberingReferencePoints(RefTempInd(RefInd),:), ...
                length(RefUpTempInd),1);
            if std(sqrt(sum(VectT.^2,2))) < 0.0001
                HyperBoxesSubRegionExtremesIndicesTempUp=sort([RefTempInd(RefInd) ...
                    RefCanUp(DistTempIndUp(1:size(DataInfo.ExtremesObj,1)-1))]);%%1 means extreme is down
                HyperBoxesSubRegionExtremesIndicesTemp=[HyperBoxesSubRegionExtremesIndicesTemp; ... 
                    HyperBoxesSubRegionExtremesIndicesTempUp];
            end
        end
    end
end

HyperBoxesSubRegionExtremesIndicesTemp=unique(HyperBoxesSubRegionExtremesIndicesTemp,'rows');
DataInfo.HyperBoxesSubRegionExtremesIndices=HyperBoxesSubRegionExtremesIndicesTemp;
end

function DataInfo=BuildSurrogates(DataInfo)
NumHypBox=size(DataInfo.HyperBoxesSubRegionExtremesIndices,1);
SurrogateInfo{NumHypBox}=[];
for Region=1:NumHypBox
    %%Generating reference sample points
    if DataInfo.NumObj==3
        ExtInd=DataInfo.HyperBoxesSubRegionExtremesIndices(Region,:);
        ProjMat=zeros(size(DataInfo.ExtremesObj,2),DataInfo.NumObj);
        I=eye(size(DataInfo.ExtremesObj,2),DataInfo.NumObj);
        I(end,:)=[];
        ProjMat(size(DataInfo.ExtremesObj,2),:)=DataInfo.PredeterminedReferencePoints( ...
            ExtInd(end),:);
        IndTemp=ExtInd(ExtInd~=ExtInd(end));
        [~,Ind_]=min(sqrt(sum((repmat(I(1,:),size(DataInfo.ExtremesObj,2)-1,1)- ...
            DataInfo.PredeterminedReferencePoints(IndTemp,:)).^2,2)));
        ProjMat(1,:)=DataInfo.PredeterminedReferencePoints(IndTemp(Ind_),:);
        IndTemp(Ind_)=[];
        ProjMat(2,:)=DataInfo.PredeterminedReferencePoints(IndTemp,:);
    elseif DataInfo.NumObj==2
        ExtInd=DataInfo.HyperBoxesSubRegionExtremesIndices(Region,:);
        ProjMat=DataInfo.PredeterminedReferencePoints(ExtInd,:);
    end
    [~,RefSamTemp]=ReferenceVectorGenerator(DataInfo.NumDivRefSam,0,DataInfo.NumObj);
    RefSamTemp=RefSamTemp*ProjMat;
    lb=min(DataInfo.NondominatedSolDec(ExtInd,:));
    ub=max(DataInfo.NondominatedSolDec(ExtInd,:));
    t=ub-lb<0.001;
    lb(t)=lb(t)-0.001;
    ub(t)=ub(t)+0.001;
    t= lb < DataInfo.LowerBounds;
    lb(t) = DataInfo.LowerBounds(t);
    t= lb > DataInfo.UpperBounds;
    ub(t) = DataInfo.UpperBounds(t);
    xSample=lhsdesign(2*(DataInfo.NumVar+1),DataInfo.NumVar,'criterion','correlation');
    xSample=MapSamples(xSample,[lb;ub],[zeros(1,DataInfo.NumVar);ones(1,DataInfo.NumVar)]);
    %%Building surrogates
    [TrainingPoints, ASFTrainingVal, fSample]=FormTrainingPoints(xSample, ... 
        RefSamTemp,DataInfo);
    SurrogateTemp=TrainSurrogate(TrainingPoints,ASFTrainingVal,DataInfo);
    [xSample,fSample,TrainingPoints,ASFTrainingVal,Surrogate]= ... 
        UpdateSurrogate(xSample,fSample,RefSamTemp,TrainingPoints, ... 
        ASFTrainingVal,SurrogateTemp,lb,ub,DataInfo);
    SurrogateInfo{Region}.lb=lb;
    SurrogateInfo{Region}.ub=ub;
    SurrogateInfo{Region}.xSample=xSample;
    SurrogateInfo{Region}.fSample=fSample;
    SurrogateInfo{Region}.RefSamples=RefSamTemp;
    SurrogateInfo{Region}.TrainingPoints=TrainingPoints;
    SurrogateInfo{Region}.ASFTraningVal=ASFTrainingVal;
    SurrogateInfo{Region}.Surrogate=Surrogate;
    SurrogateInfo{Region}.ExtremsIndices=ExtInd(1:DataInfo.NumObj);    
end
DataInfo.SurrogateInfo=SurrogateInfo;
end

function [TrainingPoints, ASFTrainingVal, fSample]=FormTrainingPoints(xSample,RefSam,DataInfo)
TrainingPoints=[];
ASFTrainingVal=[];
fSample=zeros(size(xSample,1),DataInfo.NumObj);
for SamInd=1:size(xSample,1);
    [ASFValTemp, fSample(SamInd,:)]=ASF_OBJ(xSample(SamInd,:),RefSam,DataInfo,1);
    TrainingPoints=[TrainingPoints;[repmat(xSample(SamInd,:),size(RefSam,1),1) RefSam]];
    ASFTrainingVal=[ASFTrainingVal;ASFValTemp];
end
end

function Model=TrainSurrogate(InputValue,OutputValue,DataInfo)

switch DataInfo.SurrogateType
    case 'RBF'
        Model = rbfcreate(InputValue',OutputValue','RBFFunction','cubic');
    case 'SVR'
        Model = fitrsvm(InputValue,OutputValue,'KernelFunction', ... 
            'gaussian','KernelScale','auto','Standardize',true);
    case 'KRIGING'
        theta = repmat(10, 1, size(InputValue,2)); lob = repmat(1e-1, ...
            1, size(InputValue,2)); upb = repmat(20, 1, size(InputValue,2));
        [Model, ~] = dacefit(InputValue, OutputValue, @regpoly0, @corrgauss, ... 
            theta, lob, upb);
end
end

function [xSample,fSample,TrainingPoints,ASFTrainingVal,Surrogate]= ...
    UpdateSurrogate(xSample,fSample,RefSam,TrainingPoints,ASFTrainingVal,Surrogate,lb,ub,DataInfo)
AddPenalty=0;
xSamplePF=xSample;
fSamplePF=fSample;
Update=1;
while Update
    ParetoFitnessPenaltyVal=ParetoFitnessPenalty(fSamplePF);
    if AddPenalty
        ParetoFitnessPenaltyVal(end) = max(ParetoFitnessPenaltyVal)+ParetoFitnessPenaltyVal(end);
        AddPenalty=0;
    end
    ParetoFitnessSurrogate=TrainSurrogate(xSamplePF,ParetoFitnessPenaltyVal,DataInfo);
    opts.es = 1e-4;
    opts.maxevals = 500;
    opts.maxits = 500;
    opts.maxdeep = 1000;
    opts.testflag = 0;
    opts.showits = 0;
    Problem.f ='PredictParetoFitnessPenalty';
    Struc.Model = ParetoFitnessSurrogate;
    Struc.SurrogateType=DataInfo.SurrogateType;
    bound = [lb' ub'];
    [~, xNewCandidate] = Direct(Problem,bound,opts,Struc);
    xNewCandidate = xNewCandidate';
    [NewCandidateTrainingPoints, NewCandidateASFTrainingVal, fNewCandidate] ... 
        =FormTrainingPoints(xNewCandidate,RefSam,DataInfo);
    PredictedAsfNewCandidate=PredictSurrogate(NewCandidateTrainingPoints,Surrogate,DataInfo);
    R2 = 1-(sum((NewCandidateASFTrainingVal - PredictedAsfNewCandidate).^2) ...
        / (sum((NewCandidateASFTrainingVal-mean(NewCandidateASFTrainingVal)).^2)));
    RMSE = sqrt(sum((PredictedAsfNewCandidate - NewCandidateASFTrainingVal).^2) ... 
        /size(NewCandidateASFTrainingVal,1));
    if  size(fSample,1) >= DataInfo.MaxFunEvlForEachHyp || (R2 > 0.95) || (RMSE < 0.001) 
        if ~ismember(xNewCandidate,xSample)
            xSample = [xSample; xNewCandidate];
            fSample = [fSample; fNewCandidate];
            TrainingPoints = [TrainingPoints; NewCandidateTrainingPoints];
            ASFTrainingVal = [ASFTrainingVal; NewCandidateASFTrainingVal];
            Surrogate=TrainSurrogate(TrainingPoints,ASFTrainingVal,DataInfo);
        else
            Surrogate=TrainSurrogate(TrainingPoints,ASFTrainingVal,DataInfo);
        end
        Update=0;
    else
        xSample = [xSample; xNewCandidate];
        fSample = [fSample; fNewCandidate];
        
        TrainingPoints = [TrainingPoints; NewCandidateTrainingPoints];
        ASFTrainingVal = [ASFTrainingVal; NewCandidateASFTrainingVal];
        Surrogate=TrainSurrogate(TrainingPoints,ASFTrainingVal,DataInfo);
        
        xSamplePF = [xSamplePF; xNewCandidate];
        fSamplePF = [fSamplePF; fNewCandidate];
    end%
end
end