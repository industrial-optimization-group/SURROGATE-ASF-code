function Data = datainput_SURROGATEASF
load DataInfo
Data.xlow=DataInfo.LowerBounds; %lower variable bounds
Data.xup=DataInfo.UpperBounds; %upper variable bounds
Data.dim=DataInfo.NumVar; %problem dimension
Data.integer=[]; %indices of integer variables
Data.continuous=(1:DataInfo.NumVar); %indices of continuous variables
Ref = DataInfo.RefForASF; 
Tag=DataInfo.Tag;
Data.objfunction=@(x)ASF_OBJ(x,Ref,DataInfo,Tag); %handle to objective function
end%function datainput_SURROGATEASF