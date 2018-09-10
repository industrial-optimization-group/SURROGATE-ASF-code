function [RBFVALUE, NORMVALUE] = ComputeRBF(CandPoint,Data)
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
% ComputeRBF predicts the objective function values of the candidate points
% and also returns the distance of each candidate point to all already
% sampled points
%
% Input: 
% CandPoint: (Ncand x dimension) matrix with candidate points for next
%            expensive function evaluation
% Data: struct-variable with all problem information
%
% Output:
% RBFVALUE: objective function value predicted by RBF model
% NORMVALUE: matrix with distances of all candidate points to already
%            sampled points
%--------------------------------------------------------------------------

numpoints = size(CandPoint,1);%determine number of candidate points

%compute pairwise distances between candidates and already sampled points
NORMVALUE=pdist2(CandPoint,Data.S(1:Data.m,:))'; 

%compute radial basis function value for distances
U_Y = phi(NORMVALUE,Data.phifunction);

%determine the polynomial tail (depending on rbf model)
switch Data.polynomial
case 'none'
    PolyPart = zeros(numpoints,1);
case 'constant'
    PolyPart = ctail*ones(numpoints,1);
case 'linear'
    PolyPart = [ones(numpoints,1),CandPoint]*Data.ctail;
case 'quadratic'
    temp = [ones(numpoints,1),CandPoint,zeros(numpoints,(Data.dim*(Data.dim+1))/2)];
    columnpos = Data.dim+1;
    for i = 1:Data.dim
        for j = i:Data.dim
            columnpos = columnpos + 1;
            temp(:,columnpos) = CandPoint(:,i).*CandPoint(:,j);
        end
    end
    PolyPart = temp*Data.ctail;
otherwise
    disp('Error: Invalid polynomial tail.');
    return;
end

%predict objective function values at candidate points
RBFVALUE = (U_Y')*Data.lambda + PolyPart;

end %function