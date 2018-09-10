function Data = LocalStochRBFrestart(Data,maxeval,NumberNewSamples)  
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
% LocalStochRBFrestart is the surrogate model algorithm that uses an RBF
% interpolant to approximate the simulation model and uses this information
% to decide where to sample in every iteration.
% The method may stop at a local optimum, and in case a local optimum has
% been found in fewer than maxeval function evaluations, the algorithm
% starts from scratch until totally maxeval function evaluations of the 
% current trial are used up.
%
% Input: 
% Data: struct-variable with problem information (variable bounds,
%       objective/simulation function handle, etc.)
% maxeval: maximum number of function evaluations for each trial
% NumberNewSamples: number of points where objective/simulation function
%                   is evaluated in every iteration of the algorithm; if
%                   NumberNewSamples > 1, evaluations are in parallel, in
%                   which case Matlab Parallel Computing Toolbox is
%                   required.
%
% Output:
% Data: updated struct-variable containing all the results of the current
%       trial
%--------------------------------------------------------------------------

m = 2*(Data.dim+1);  % number of points in initial experimental design
%initialize arrays for collecting results of current trial
% number of restarts (when algorithm restarts within one trial after
% encountering local optimum)
numstart = 0;
Y_all=[]; %collect all objective function values of the current trial here
S_all=[]; %collect all sample points of the current trial here
value = inf; %best objective function value found so far in the current trial 
numevals = 0; %number of function evaluations done so far
Fevaltime_all=[]; %collect all objective function evaluation times of the current trial here

while numevals < maxeval %do until max. number of allowed f-evals reached
    numstart = numstart + 1; %increment number of algorithm restarts
    
    % create initial experimental design by symmetric Latin hypercube
    % sampling
    %for cubic and thin-plate spline RBF: rank_P must be Data.dim+1
    rank_P = 0;     
    %regenerate initial experimental design until matrix rank is
    %dimension+1
    while rank_P ~= Data.dim+1 
        Data.S = SLHDstandard(Data.dim,m);
        P = [ones(m,1),Data.S]; %matrix augmented with vector of ones for computing RBF model parameters
        rank_P = rank(P);
    end
    
    
    % for the current number of starts, run local optimization
    Data= LocalStochRBFstop(Data,maxeval-numevals,NumberNewSamples);
    
    %update best solution found if current solution is better than best
    %point found so far
    if Data.Fbest < value 
        solution = Data.xbest; %best point 
        value = Data.Fbest; %best function value
    end 
    
    %after LocalStochRBFstop.m stops (in local min, update vectors) 
    Fevaltime_all=[Fevaltime_all;Data.fevaltime];   %collect function evaluation times  
    Y_all=[Y_all;Data.Y]; %collect function values
    S_all=[S_all;Data.S]; %collect sample sites
    numevals = numevals + Data.NumberFevals; %update number of function evaluations done so far  
end

%update Data struct with the results of the current trial
Data.S=S_all; %sample sites
Data.Y=Y_all; %function values
Data.fevaltime=Fevaltime_all; %function evaluation times
Data.xbest=solution; %best point found in this trial
Data.Fbest=value; %best function value found in this trial
Data.NumberFevals=numevals; %number of function evaluations done in this trial
Data.NumberOfRestarts=numstart; %number of restarts of local search

end%function