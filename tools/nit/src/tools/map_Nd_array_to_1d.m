function R_1d = map_Nd_array_to_1d(R_Nd)
%%% *function R_1d = map_Nd_array_to_1d(R_Nd)*
%%%
%%% ### Description
%%% This function maps a N-dimensional array defined by `R_Nd` to 1D vector.
%%% The function is useful when analysing multi-dimensional activity (either
%%% in time on a single unit or across units) when interested in understanding
%%% if specific activation patterns result in different information levels.
%%%
%%% The multi-dimensional input is mapped on a single dimensional vector of
%%% the same size as the possible combinations of all responses in the
%%% input multi-dimensional activity.
%%%
%%% As an example, the 2-dimensional input `R_nd` below:
%%%
%%%     R_nd = [ 1 2 1; 2 3 1]
%%%
%%% Is mapped on the rows of the matrix of all possible 2-dimensional responses:
%%%
%%%     R_nd = [ 1 1; 2 1; 1 2; 2 2; 1 3; 2 3]
%%%
%%% The output `R_1d` is then:
%%%
%%% 	R_nd = [3 6 1]
%%%
%%% ### Inputs:
%%% - *R_Nd*: *nDims X nTrials* input response matrix
%%%
%%% ### Outputs:
%%% - *R_1D*: output 1D response based on the N-dimensional space defined by `R_Nd`
%%%
%
%  This source code is part of:
%  NIT - Neuroscience Information Toolbox
%  Copyright (C) 2020  Roberto Maffulli, Miguel Angel Casal Santiago
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.

% define variability range of response for each dimension
Ndims = length(R_Nd(:,1));
for d = 1:Ndims
    resps{d} = unique(R_Nd(d,:));
end

% 1d input is returned identical
if Ndims == 1 
    R_1d = R_Nd ;
else
    ii = 1:Ndims;
    % build grid based on responses at each dimension
    [resps_grid{ii}] = ndgrid(resps{ii}) ;
    % concatenate grid in single matrix
    resps_grid = reshape(cat(Ndims+1,resps_grid{:}), [], Ndims);
    % R_1d is mapping responses to indices of resps_grid
    R_1d = zeros(1,length(R_Nd(1,:)));
    for i = 1:length(R_Nd(1,:))
        [~,R_1d(i)] = ismember(R_Nd(:,i)',resps_grid,'rows');
    end
end
end

