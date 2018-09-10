function Data = datainput_EachObjective
load DataInfo
Data.xlow=DataInfo.LowerBounds; %lower variable bounds
Data.xup=DataInfo.UpperBounds; %upper variable bounds
Data.dim=DataInfo.NumVar; %problem dimension
Data.integer=[]; %indices of integer variables
Data.continuous=(1:DataInfo.NumVar); %indices of continuous variables
Obj=DataInfo.Obj;
Data.objfunction=@(x)ObjectiveFunction(x,Obj,DataInfo); %handle to objective function
end%function datainput_SURROGATEASF

function y=ObjectiveFunction(x,Obj,DataInfo)
Output=P_objective('value',DataInfo.Prob,DataInfo.NumObj,x);
y=Output(Obj);
end