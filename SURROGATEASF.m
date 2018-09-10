function SURROGATEASF()
%----------------********************************--------------------------
% Copyright (C) 2017 By Mohammad Tabatabaei 
%
% This file is part of the program SURROGATEASF.m
%
%    SURROGATEASF.m is a free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    SURROGATEASF.m is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with SURROGATEASF.m.  If not, see <http://www.gnu.org/licenses/>.
%----------------********************************--------------------------
%
%--------------------------------------------------------------------------
% SURROGATE-ASF: An interactive surrogate-based method for solving 
% computationally expensive %multiobjective optimization problems of the form 
%         min {f_1(x),…f_k(x)}    subject to   LB <= x <= UB
% where 
% f_i, i=1,..k (k=2,3), are computationally expensive objective function 
% whose closed-form formulas are not available (i.e., black-box functions). 
%
% LB and UB are the lower and upper bounds of the decision variables. 
% SURROAGTE-ASF aims at finding the most preferred solution for the 
% decision maker (user) within a  limited number of function evaluations. 
% SURROGATE-ASF consists of two phases: i.e., initialization and 
% decision making phases. In the initialization phase, SURROGATE-ASF employs 
% a surrogate-based single  objective optimization method to get some 
% information regarding the location of the solutions in the 
% decision/objective space. Then the original problem is decomposed into 
% a finite number of single-objective surrogate problems. In the decision 
% making phase, the interaction with the decision maker is conducted to find 
% the most preferred solution. The decision maker expresses his/her 
% preferences in the form of a reference point. 
% Please refer with all questions, comments, bug reports, etc. to
% tabatabaei62@yahoo.com
%----------------********************************--------------------------
%
%
% Input:
% The multiobjective optimization problem is defined in “P_objective.m”. 
% P_objective.m is a part of the implementation of RVEA method developed in 
%
% R. Cheng, Y. Jin, M. Olhofer and B. Sendhoff, A Reference Vector Guided 
% Evolutionary Algorithm for Many-objective Optimization, IEEE Transactions
% on Evolutionary Computation, 2016.
%
% In this implementation of SURROGATE-ASF, two surrogate-based methods 
% developed for computationally expensive single-objective optimization 
% problems are incorporated as listed below:
%
% MATSuMoTo developed in J. Muller and C. A. Shoemaker. Influence of ensemble 
% surrogate models and sampling strategy on the solution quality of algorithms 
% for computationally expensive black-box global optimization problems. 
% Journal of Global Optimization, 60(2):123-144, 2014.
%
% and
%
% StochasticRBF developed in R.G. Regis and C.A. Shoemaker, A Stochastic 
% Radial Basis Function Method for the Global Optimization of Expensive 
% Functions, INFORMS Journal on Computing, vol. 19, pp. 497-509, 2007 
%
% and 
%
% R.G. Regis and C.A. Shoemaker, Parallel Stochastic Global Optimization 
% Using Radial Basis Functions, INFORMS Journal on Computing, vol. 21, pp. 
% 411-426, 2009.
%
% In order to solve single-objective surrogate problems, the DIRECT method 
% developed in 
%
% D. R. Jones, C. D. Perttunen, and B. E. Stuckman. Lipschitzian optimization 
% without the Lipschitz constant. Journal of Optimization Theory and Applications,
% 79(1):157-181, 1993.” is incorporated. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Feel free to replace these methods with other methods.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%How to use this code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To use this code, run "SURROGATE-ASF.m". It starts to ask the name of   %
% the problem defined in "P_objective.m", the number of objective         %
% functions, number of decision variables, the lower and upper bounds of  %
% the decision variables, and the maximum number of function evaluations. %
% Then, the method starts to solve the problem.                           %
% In order to provide preferred ranges for objective functions, a matrix k%
% by 2 is assigned to DataInfo.InitialRanges. Each row of this matrix     %
% contains the minimum and maximum of the preferred range for the         %
% corresponsing objective function. In this code, the defult preferred    %
% range for each objective function is [-50, 50].                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             When using SURROGATE-ASF, please cite the following paper:  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Mohammad Tabatabaei, Markus Hartikainen, Karthik Sindhya, Jussi Hakanen, 
% Kaisa Miettinen, An Interactive Surrogate-based Method for Computationally 
% Expensive Multiobjective Optimization, submitted.

clc;
close all;
warning off;
delete('Results.Mat');
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
DataInfo=InputData();
switch DataInfo.NumObj
    case 2
        DataInfo.NumDivRefSam=10;
        DataInfo.NumDivPreRef=3;
        DataInfo.InitialRanges=[-50 50;-50 50];
        FunEvl = floor(linprog([-4 -3],[4 3],DataInfo.MaxFunEvlBudget-10,[],[], ...
            [2*(DataInfo.NumVar+1) 2*(DataInfo.NumVar+1)], ...
            [DataInfo.MaxFunEvlBudget-10 DataInfo.MaxFunEvlBudget-10]));
        DataInfo.MaxFunEvlForEachPreRef=FunEvl(1);
        DataInfo.MaxFunEvlForEachHyp=FunEvl(2);
        DataInfo.MaxFunEvlForDM = DataInfo.MaxFunEvlBudget - sum([4 3] .* FunEvl');
    case 3
        DataInfo.NumDivRefSam=5;
        DataInfo.NumDivPreRef=4;
        DataInfo.InitialRanges=[-50 50;-50 50;-50 50];
        FunEvl = floor(linprog([-6 -4],[6 4],DataInfo.MaxFunEvlBudget-10,[],[], ...
            [2*(DataInfo.NumVar+1) 2*(DataInfo.NumVar+1)], ...
            [DataInfo.MaxFunEvlBudget-10 DataInfo.MaxFunEvlBudget-10]));
        DataInfo.MaxFunEvlForEachPreRef=FunEvl(1);
        DataInfo.MaxFunEvlForEachHyp=FunEvl(2);
        DataInfo.MaxFunEvlForDM = DataInfo.MaxFunEvlBudget - sum([6 4] .* FunEvl');
end

DataInfo.ASFWeight=1./(DataInfo.InitialRanges(:,2)-DataInfo.InitialRanges(:,1))';
DataInfo.SingleSurrogateBasedSolver='Mat';%'Mat','DIRECT','StochasticRBF'
DataInfo.SurrogateType='RBF';%'RBF';'KRIGING'
DataInfo=InitializationPhase(DataInfo); 
DataInfo=DecisionMakingPhase(DataInfo);
clc;
disp('The most preferred solution in the objective space is') 
disp(['[' num2str(DataInfo.TheMostPreferredSolutionInObjectiveSpace) '].'])
disp(' ')
disp('The most preferred solution in the decision space is') 
disp(['[' num2str(DataInfo.TheMostPreferredSolutionInDecisionSpace) '].'])
save DataInfo DataInfo
end

function DataInfo=InputData()

DataInfo.Prob=input('Input the problem name as an string ');

ask=1;
while ask
   DataInfo.NumObj=input('Input objective numbers (greater than 1 and less than or equal 3) ');
   if DataInfo.NumObj < 2 ||  DataInfo.NumObj > 3
       disp('The number of objectives must be 2 or 3.')
   else 
       ask = 0;
   end       
end
ask=1;
while ask
   DataInfo.NumVar=input('Input decision variables numbers (greater than 0) ');
   if ~isscalar(DataInfo.NumVar) || DataInfo.NumVar < 1
       disp('The number of decision variables must be greater than 0.')
   else 
       ask = 0;
   end       
end
ask=1;
while ask
    DataInfo.LowerBounds=input(['Input a vector of lower bound with a dimension 1 x ' num2str(DataInfo.NumVar) ' ']);
    if length(DataInfo.LowerBounds)==DataInfo.NumVar
        ask = 0;
    else
        disp(['The vector must have a dimension 1 x ' num2str(DataInfo.NumVar)]'.')
    end
end

ask=1;
while ask
    DataInfo.UpperBounds=input(['Input a vector of upper bounds with a dimension 1 x ' num2str(DataInfo.NumVar) ' ']);
    if length(DataInfo.UpperBounds)~=DataInfo.NumVar 
        disp(['The vector must have a dimension 1 x ' num2str(DataInfo.NumVar)]'.')
    elseif any(DataInfo.LowerBounds >= DataInfo.UpperBounds)
        disp('Upper bound must be greater than or equal lower bound.')
    else
        ask=0;
    end
end
ask=1;
while ask
    DataInfo.MaxFunEvlBudget = input('Input the maximum number of function evaluations ');  
    if DataInfo.MaxFunEvlBudget > 0
        ask = 0;
    end
end
end