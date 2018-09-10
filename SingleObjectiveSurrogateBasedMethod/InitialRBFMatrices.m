function [PHI,phi0,P,pdim] = InitialRBFMatrices(maxevals,Data,PairwiseDistance)
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
% set up matrices for computing parameters of RBF model based on points in 
% initial experimental design
%
% Input:
% maxevals: maximal number of allowed function evaluations
% Data: struct-variable with all problem information such as sampled points
% PairwiseDistance: pairwise distances between points in initial
%                   experimental design
%
% Output:
% PHI: matrix containing pairwise distances of all points to each other,
%      will be updated in following iterations
% phi0: PHI-value of two equal points (depends on RBF model!)
% P: sample site matrix, needed for determining parameters of polynomial
%    tail
% pdim: dimension of P-matrix (number of columns)
%--------------------------------------------------------------------------


% initial PHI matrix - radial basis function value for pairwise distances
PHI = zeros(maxevals);
switch Data.phifunction
case 'linear'
    PairwiseDistance= PairwiseDistance;
case 'cubic'
    PairwiseDistance= PairwiseDistance.^3;
case 'thinplate'
    PairwiseDistance=PairwiseDistance.^2.*log(PairwiseDistance+realmin);
    PairwiseDistance(logical(eye(size(PairwiseDistance)))) = 0;
end

PHI(1:Data.m,1:Data.m)=PairwiseDistance;
phi0 = phi(0,Data.phifunction); %phi-value where distance of 2 points =0 (diagonal entries)

% initial P matrix
switch Data.polynomial
case 'none'
    pdim = 0;
    P = [];
case 'constant'
    pdim = 1;
    P = ones(maxevals,1);
case 'linear'
    pdim = Data.dim + 1;
    P = [ones(maxevals,1),Data.S];
case 'quadratic'
    pdim = ((Data.dim+1)*(Data.dim+2))/2;
    P = [ones(maxevals,1),Data.S,zeros(maxevals,(Data.dim*(Data.dim+1))/2)];
    columnpos = Data.dim+1;
    for ii = 1:Data.dim
        for jj = ii:Data.dim
            columnpos = columnpos + 1;
            P(1:Data.m,columnpos) = Data.S(1:Data.m,ii).*Data.S(1:Data.m,jj);
        end
    end
otherwise
    disp('Error: Invalid polynomial tail.');
    return;
end

end %function