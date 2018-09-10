function Data = datainput_SURROGATEASF
load DataInfo
Data.xlow=DataInfo.LowerBounds; %lower variable bounds
Data.xup=DataInfo.UpperBounds; %upper variable bounds
Data.dim=DataInfo.NumVar; %problem dimension
Data.integer=[]; %indices of integer variables
Data.continuous=(1:DataInfo.NumVar); %indices of continuous variables
Ref = DataInfo.RefForASF; 
Data.objfunction=@(x)ASF_OBJ(x,Ref,DataInfo); %handle to objective function
end%function datainput_SURROGATEASF