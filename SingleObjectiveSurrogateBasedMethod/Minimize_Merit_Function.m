function [xselected,normval]=Minimize_Merit_Function(Data,CandPoint,...
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
% Minimize_Merit_Function computes the distance and response surface
% criteria for every candidate point. The values are scaled to [0,1], and
% the candidate with the best weighted score of both criteria becomes the
% new sample point. If there are more than one new sample point to be
% selected, the distances of the candidate points to the previously
% selected candidate point have to be taken into account.
%
% Input:
% Data: struct-variable with all problem information
% CandPoint: matrix with candidate points
% NumberNewSamples: number of points to be selected for next costly evaluation 
%
% Output:
% xselected: matrix with all seleced points for next evaluation
% normval: cell array with distances to previously evaluated points and
% other selected candidate points, needed for updating PHI matrix later
%--------------------------------------------------------------------------    
    
% predict objective function value of the candidate points, and compute
% their distances to all already sampled points
[CandValue, NormValue] = ComputeRBF(CandPoint,Data);

%scale predicted objective function values to [0,1]
MinCandValue = min(CandValue); %find min of predIcted objective function value
MaxCandValue = max(CandValue); %find maximum of predicted objective function value
if MinCandValue == MaxCandValue  %compute scaled objective function value scores
    ScaledCandValue=ones(length(CandValue),1);
else
    ScaledCandValue = (CandValue-MinCandValue)/(MaxCandValue-MinCandValue);
end


if NumberNewSamples == 1 %only one point is selected
    valueweight=0.95; %weight for response surface criterion 
    %scale distances to already evaluated points to [0,1]
    CandMinDist = (min(NormValue,[],1))'; %minimum distance of every candidate point to already sampled points  
    MaxCandMinDist = max(CandMinDist); %maximum of distances
    MinCandMinDist = min(CandMinDist); %minimum of distances
    if  MaxCandMinDist ==MinCandMinDist  %compute distance criterion scores
        ScaledCandMinDist =ones(length(CandMinDist),1);
    else
        ScaledCandMinDist = (MaxCandMinDist-CandMinDist)/(MaxCandMinDist-MinCandMinDist);
    end

    %compute weighted score for all candidates
    CandTotalValue = valueweight*ScaledCandValue + (1 - valueweight)*ScaledCandMinDist;

    %assign bad scores to candidate points that are too close to already sampled
    %points
    CandTotalValue(CandMinDist < Data.tolerance) = Inf; 

    %find candidate with best score -> becomes new sample point
    [MinCandTotalValue,selindex] = min(CandTotalValue);
    xselected(1,:) = CandPoint(selindex,:);  %new sample point
    normval{1} = NormValue(:,selindex)'; %new sample point's distance to already evaluated points

else  %more than one new sample point wanted
    wp_id=0;
    if NumberNewSamples < 4
        if NumberNewSamples ==2
            weightpattern=[0.5,0.95];
        elseif NumberNewSamples == 3
            weightpattern=[0.3,0.5,0.95];
        end
    else
        weightpattern=[0.3,0.5,0.8,0.95];
    end
    for ii =1:NumberNewSamples
        wp_id=wp_id+1;
        if wp_id > length(weightpattern)
            wp_id=1;
        end
        valueweight=weightpattern(wp_id);
        if ii == 1 %select first candidate point
            CandMinDist = (min(NormValue,[],1))'; %minimum distance of every candidate point to already sampled points  
            MaxCandMinDist = max(CandMinDist); %maximum of distances
            MinCandMinDist = min(CandMinDist); %minimum of distances
            if  MaxCandMinDist ==MinCandMinDist  %compute distance criterion scores
                ScaledCandMinDist =ones(length(CandMinDist),1);
            else
                ScaledCandMinDist = (MaxCandMinDist-CandMinDist)/(MaxCandMinDist-MinCandMinDist);
            end

            %compute weighted score for all candidates
            CandTotalValue = valueweight*ScaledCandValue + (1 - valueweight)*ScaledCandMinDist;

            %assign bad scores to candidate points that are too close to already sampled
            %points
            CandTotalValue(CandMinDist < Data.tolerance) = Inf; 

            %find candidate with best score -> becomes new sample point
            [MinCandTotalValue,selindex] = min(CandTotalValue);
            xselected(1,:) = CandPoint(selindex,:);  %new sample point
            normval{1} = NormValue(:,selindex)'; %new sample point's distance to already evaluated points
        else    
            %compute distance of all candidate points to the previously selected
            %candidate point
            NormValueP=sqrt(sum((repmat(xselected(ii-1,:),size(CandPoint,1),1)-CandPoint).^2,2));
            NormValue=[NormValue;NormValueP']; %augment distance matrix

            %re-scale distance values to [0,1]
            CandMinDist = (min(NormValue,[],1))'; %minimum distance of every candidate point to already sampled points  
            MaxCandMinDist = max(CandMinDist); %maximum of distances
            MinCandMinDist = min(CandMinDist); %minimum of distances
            if  MaxCandMinDist ==MinCandMinDist  %compute distance criterion scores
                ScaledCandMinDist =ones(length(CandMinDist),1);
            else
                ScaledCandMinDist = (MaxCandMinDist-CandMinDist)/(MaxCandMinDist-MinCandMinDist);
            end

            %compute weighted score for all candidates
            CandTotalValue = valueweight*ScaledCandValue + (1 - valueweight)*ScaledCandMinDist;    
            %assign bad values to points that are too close to already
            %evaluated/chosen points
            CandTotalValue(CandMinDist < Data.tolerance)=Inf;

            [MinCandTotalValue,selindex] = min(CandTotalValue); %find best candidate
            xselected(ii,:) = CandPoint(selindex,:);   %select best candidate
            normval{ii} = NormValue(:,selindex)';  %distance of candidate to all already evaluated and selected points       
        end
    end
end


ii=1; %delete identical selected sample points
while ii <= size(xselected,1)-1 
    jj=ii+1;
    while jj <= size(xselected,1)
        if sqrt(sum((xselected(ii,:)-xselected(jj,:)).^2)) < Data.tolerance
            xselected(jj,:)=[];
        else
            jj=jj+1;
        end
    end
    ii=ii+1;
end

if isempty(xselected) %regenerate potential sample points if no new point selected
    while 1
        xselected = Data.xlow+(Data.xup-Data.xlow).*rand(1,Data.dim);
        dist_new=min(pdist2(xselected,Data.S));
        if dist_new>Data.tolerance
            break
        end
    end
end    
end%function