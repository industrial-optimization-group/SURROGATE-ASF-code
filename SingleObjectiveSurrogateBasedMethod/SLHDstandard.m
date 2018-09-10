function InitialPoints = SLHDstandard(d,m)
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
% SLHD creates a symmetric latin hypercube design. d is the dimension of 
% the input and m is the number of initial points to be selected.
%
% InitialPoints is a matrix with the initial sample points, scaled to the
% unit hypercube
%--------------------------------------------------------------------------

delta = (1/m)*ones(1,d);

X = zeros(m,d);
for j = 1:d
    for i = 1:m
        X(i,j) = ((2*i-1)/2)*delta(j);
    end
end

P = zeros(m,d);
P(:,1) = (1:m)';
if (mod(m,2) == 0)
   k = m/2;
else
   k = (m-1)/2;
   P(k+1,:) = (k+1)*ones(1,d);
end

for j = 2:d
   P(1:k,j) = randperm(k)';
   for i = 1:k
      if (rand(1) <= 0.5)
         P(m+1-i,j) = m+1-P(i,j);
      else
         P(m+1-i,j) = P(i,j);
         P(i,j) = m+1-P(i,j);
      end
   end
end

InitialPoints = zeros(m,d);
for j = 1:d
    for i = 1:m
        InitialPoints(i,j) = X(P(i,j),j);
    end
end