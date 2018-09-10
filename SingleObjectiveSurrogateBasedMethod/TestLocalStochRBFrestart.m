function Solution=TestLocalStochRBFrestart(Data,maxeval,Ntrials,...
    NumberNewSamples)
%----------------********************************--------------------------
% Copyright (C) 2013 Cornell University
% This file is part of the program StochasticRBF.m
%
%    StochasticRBF.m is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    StochasticRBF.m is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with StochasticRBF.m.  If not, see <http://www.gnu.org/licenses/>.
%----------------********************************--------------------------
%
%
%
%----------------*****  Contact Information *****--------------------------
%   Primary Contact (Implementation Questions, Bug Reports, etc.):
%   Juliane Mueller: juliane.mueller2901@gmail.com
%       
%   Secondary Contact:
%       Christine A. Shoemaker: cas12@cornell.edu
%
%----------------********************************--------------------------
%--------------------------------------------------------------------------
% TestLocalStochRBFrestart calls the surrogate model optimization algorithm
% Ntrials times, and collects the results of each trial in the
% sruct-variable "Solution"
%
% Input:
% Data: struct-variable with problem information (variable bounds,
%       objective/simulation function handle, etc.)
% maxeval: maximum number of function evaluations for each trial
% Ntrials: number of trials for running algorithm
% NumberNewSamples: number of points where objective/simulation function
%                   is evaluated in every iteration of the algorithm; if
%                   NumberNewSamples > 1, evaluations are in parallel, in
%                   which case Matlab Parallel Computing Toolbox is
%                   required.
%
% Output:
% Solution: struct-variable containing all the results of every run of the
% algorithm
%--------------------------------------------------------------------------

warning off;

% initializations; struct "Solution" contains results of all trials
Solution.BestPoints = zeros(Ntrials,Data.dim); %best point in each trial
Solution.BestValues = zeros(Ntrials,1); %best function value in each trial
Solution.NumFuncEval = zeros(Ntrials,1); %number of fevals in each trial
Solution.AvgFuncEvalTime = zeros(Ntrials,1); %average times for function evaluations
Solution.FuncVal=zeros(maxeval,Ntrials); %array with all object function values of each trial
Solution.DMatrix=zeros(maxeval,Data.dim,Ntrials); %array with all sample points of each trial
Solution.NumberOfRestarts= zeros(Ntrials,1); %number of algorithm restarts in each trial (=number local opt found)

%loop over all algorithm trials
for jj = 1:Ntrials %for each trial   
    rand('state',jj); %set random number seed according to jj'th trial
    randn('state',jj);
    
    %call the surrogate optimization function
    Data_temp=LocalStochRBFrestart(Data,maxeval,NumberNewSamples);
        
    %gather results in "Solution" struct-variable
    Solution.BestValues(jj)=Data_temp.Fbest;
    Solution.BestPoints(jj,:)=Data_temp.xbest;
    Solution.NumFuncEval(jj,1)=Data_temp.NumberFevals;
    Solution.AvgFuncEvalTime(jj)=mean(Data_temp.fevaltime);
    Solution.FuncVal(:,jj) = Data_temp.Y;
    Solution.DMatrix(:,:,jj) =Data_temp.S;
    Solution.NumberOfRestarts(jj,1)=Data_temp.NumberOfRestarts;
end

end %function