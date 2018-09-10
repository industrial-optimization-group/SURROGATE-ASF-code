function Solution=StochasticRBF(data_file,maxeval, Ntrials,PlotResult,NumberNewSamples)
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
%--------------------------------------------------------------------------
% StochasticRBF - Unconstrained optimization using surrogate models.
% StochasticRBF attempts to solve problems of the form:
%         min  F(X)   subject to   LB <= X <= UB
% where 
% F(X) is a computationally expensive objective function whose analytical
% description is not available (black-box simulation model)
% LB and UB are the lower and upper variable bounds, i.e. only
% box-constrained problems can be solved.
% 
%StochasticRBF.m is an optimization algorithm that uses a cubic RBF
%surrogate model. The implementation is based on the two published papers:
%a) "A Stochastic Radial Basis Function Method for the Global Optimization 
%    of Expensive Functions" by R.G. Regis and C.A. Shoemaker, 2007, 
%    INFORMS Journal on Computing, vol. 19, pp. 497-509
%b) "Parallel Stochastic Global Optimization Using Radial Basis Functions" 
%    by R.G. Regis and C.A. Shoemaker, 2009, INFORMS Journal on Computing,
%    vol. 21, pp. 411-426
%These papers and this code should be referenced whenever they are used to 
%generate results for the user's own research. 
%
%Contributers to this matlab implementation are J. Mueller, R.G. Regis, and 
%C.A.Shoemaker. See also http://people.sju.edu/~rregis/pages/software.html for
%a basic implementation of the algorithm.
%
%----------------*****  Code Author Information *****----------------------
%   Code Author (Implementation Questions, Bug Reports, etc.): 
%       Juliane Mueller: juliane.mueller2901@gmail.com
%   Co-Authors:
%       Rommel G. Regis: rregis@sju.edu
%       Christine A. Shoemaker: cas12@cornell.edu
%**************************************************************************
%   Please refer with all questions, comments, bug reports, etc. to
%   juliane.mueller2901@gmail.com
%----------------********************************--------------------------
%
%
% Input:
% data_file: an m-file that contains all problem information. See for
%            example the file "datainput_hartman3.m"
% maxeval: the maximal number of allowed function evaluations in each trial
% Ntrials: the number of trials (e.g. if Ntrials=3, the algorithm runs
%          three times)
% PlotResult: 0 - if no plot of the solution is desired
%             1 - plots average objective function value over Ntrials vs 
%                 number of function evaluations
% NumberNewSamples: the number of function evaluation points to be selected
%                   in every iteration. If NumberNewSamples > 1, function
%                   evaluations are done in parallel (requires matlab
%                   parallel computing toolbox)
% 
% The algorithm saves the results in the file Results.mat to the current
% working directory.
%
% For further information, see the manual.
%
%--------------------------------------------------------------------------


%% ------------------------------------------------------------------------
close all; %close all curently open figures
%% ----------------------- START INPUT CHECK-------------------------------
if nargin < 1 ||  ~ischar(data_file) %check if input data file present
    error(['You have to supply a file name with your data. \n' ...
    'See example files and tutorial for information how to define problems.'])
end
Data=feval(data_file); % load problem data

%check if dimension is defined in input data file
if ~isfield(Data,'dim')
    error('You must provide the problem dimension.')
elseif abs(Data.dim-round(Data.dim))>0 || Data.dim <= 0 
    error('Dimension must be positive integer.')
end

%check if all upper and lower bounds are given in input data file
if ~isfield(Data,'xlow') || ~isfield(Data,'xup') || ...
        length(Data.xlow) ~= Data.dim || length(Data.xup) ~= Data.dim
    error('Vector length of lower and upper bounds must equal problem dimension')
end

%check if lower bounds < upper bounds
if any(Data.xlow >= Data.xup)
    error('Lower bounds have to be lower than upper bounds.')
end

%check if maximum number of allowed function evaluations is given
if nargin < 2  || isempty(maxeval) 
    fprintf(['No maximal number of allowed function evaluations given. \n '...
            'I use default value maxeval = 20 * dimension.'])
    maxeval= 20*Data.dim;
elseif  maxeval<0 || abs(maxeval-round(maxeval))>0   
    error('Maximal number of allowed function evaluations must be positive integer.\n');
end

%check if number of trials is given
if nargin < 3 || isempty(Ntrials) 
    fprintf(['No maximal number of trials given. \n '...
            'I use default value NumberOfTrials=1.'])
    Ntrials= 1;
elseif  Ntrials<0  || abs(Ntrials-round(Ntrials))>0     
    error('Maximal number of trials must be positive integer.\n');        
end
 
%check if user indicated if s/he wants to see results in plot form
if nargin < 4 || isempty(PlotResult)
    fprintf(['No indication if result plot wanted. \n '...
            'I use default value PlotResult=1.'])
    PlotResult= 1;
elseif abs(PlotResult)>0
    PlotResult=1;
end

%check if user provided desred number of new sample sites in each run
if nargin < 5 || isempty(NumberNewSamples) 
    fprintf(['No number of desired new sample sites given. \n '...
            'I use default value NumberNewSamples=1.'])
    NumberNewSamples= 1;
elseif  NumberNewSamples<0  || abs(NumberNewSamples-round(NumberNewSamples))>0     
    error('Number of new sample sites must be positive integer.\n');        
end

if NumberNewSamples > 1 %requires parallel computing toolbox
    v = ver;
    if ~any(strcmp('Parallel Computing Toolbox', {v.Name}))
        error(['You need Matlab Parallel Computing Toolbox if you choose \n '...
            'NumberNewSamples > 1'])
    end
end

 
%% ----------------------end input check-----------------------------------

%% ------------------------ Optimization  ---------------------------------
%parameter settings
Data.Ncand=500*Data.dim; % Traditional LMSRBF - number of candidate points
Data.phifunction='cubic'; %radial basis function type
Data.polynomial='linear'; %polynomial tail for cubic RBF

%optimization procedure
Solution=TestLocalStochRBFrestart(Data,maxeval,...
    Ntrials,NumberNewSamples);
%% -------------------- End Optimization  ---------------------------------

%% --------------------- Plot results -------------------------------------
if PlotResult %in case plotting of results is wanted
    %initialize vector of best function values found so far
    Y_cur_best=zeros(maxeval,Ntrials); 
    for ii = 1:Ntrials %go through all trials
        Y_cur=Solution.FuncVal(:,ii); % function values of current trial (trial ii)
        Y_cur_best(1,ii)=Y_cur(1); %first best function value is first function value computed
        for jj=2:maxeval %go through remaining function values
            if Y_cur(jj) < Y_cur_best(jj-1,ii) %in case f-value in iteration jj is better than current best in iteration jj-1
                Y_cur_best(jj,ii)=Y_cur(jj); % set current best value at iteration jj to this value
            else %otherwise current best in iteration jj is same as in iteration jj-1
                Y_cur_best(jj,ii)=Y_cur_best(jj-1,ii);
            end
        end
    end
    %compute means over matrix of current best values (Y_cur_best has dimension 
    % maxeval x Ntrials)
    Ymean=mean(Y_cur_best, 2); 
    Yplot=zeros(maxeval,1); %initialize vector for plotting results
    %sort results according to best point found till iteration jj
    Yplot(1)=Ymean(1);
    for jj=2:maxeval
        if Ymean(jj) < Yplot(jj-1)
            Yplot(jj)=Ymean(jj);
        else
            Yplot(jj)=Yplot(jj-1);
        end
    end
      
    figure %make a new figure
    plot(1:maxeval, Yplot,'b') %plot the line of best f-values found so far
    xlabel('Number Of Function Evaluations') 
    ylabel(sprintf('Average Best Objective Function Value In %d Trials',Ntrials));

end
%% ----------------- End Plot results -------------------------------------

%save the result to mat-file into current working directory
save Results.mat Solution

end %function