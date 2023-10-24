function [p_src, p_crs, n_S, n_R, n_C, n_singleton_dims] = build_p(S, R, C)
%%% *function [p_src, p_crs, n_S, n_R, n_C, n_singleton_dims] = build_p(S, R, C)*
%%%
%%% ### Description
%%% The function estimates the probability distribution p(s,r,c) from the 3D histogram of the input (S, R, C) occurrences.
%%%
%%% ### Inputs:
%%% - *S*: must be a one dimensional array of *n_trials X 1* elements representing the discrete value of the stimulus presented in each trial.
%%% - *R*: must be a one dimensional array of *n_trials X 1* elements representing the response at each trial, here we suppose the response has been binned already.
%%% - *C*: must be a one dimensional array of *n_trials X 1* elements representing the discrete value of the choice made by the subject in each trial.
%%%
%%% ### Outputs:
%%% - *p_src*: joint probability p(s,r,c).
%%% - *p_crs*: joint probability p(c,r,s).
%%% - *n_S*: number of stimuli.
%%% - *n_R*: number of responses.
%%% - *n_C*: number of choices.
%%% - *n_singleton_dims*: number of singleton dimensions.
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
%

R = reshape(R,length(R),1);
R_discrete_values=unique(R);
S_values=unique(S);
C_values=unique(C);

N_trials=length(R(:,1));
n_R=numel(R_discrete_values);
n_S=numel(S_values);
n_C=numel(C_values);

R_stored=false(N_trials,n_R);
S_stored=false(N_trials,n_S);

p_src=zeros(n_S,n_R,n_C);
n_singleton_dims = 3-length(size(p_src))+sum(size(p_src)==1);

for cc=1:n_C
    for rr=1:n_R
        if cc == 1
            % store for performance
            R_stored(:,rr)=(R==R_discrete_values(rr));
        end
        for ss=1:n_S
            if rr == 1
                % store for performance
                S_stored(:,ss)=(S==S_values(ss));
            end
            p_src(ss,rr,cc)=sum(R_stored(:,rr).*S_stored(:,ss).*(C==C_values(cc)));
        end
    end
end

p_src=p_src/sum(p_src(:));
p_crs=permute(p_src,[3 2 1]);

end