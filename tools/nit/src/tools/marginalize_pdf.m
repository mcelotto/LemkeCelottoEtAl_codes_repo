function [marginalized_pdf] = ...
    marginalize_pdf(pdf,survivingdims,iscont,varargin)
%%% *[marginalized_pdf] = marginalize_pdf(pdf,survivingdims)*
%%% 
%%% ### Description
%%% This function reduces the number of variables in a multivariate pdf to the ones specified in survivingdims. E.g. for a 3-variate input pdf with survivingdims = [1,2] the resulting reduced_pdf will be the joint probability (1,2).
%%%
%%% ### Inputs:
%%% - *pdf*: multi-variate input pdf, specified as a multi-dimensional array such as each component of `size(pdf)` is equal to the length of the corresponding random variable.
%%% - *survivingdims*: vector (or single value) of variables that will be used in the resulting reduced_pdf. If a single value the function returns the corresponding marginal pdf of such variable.
%%% - *iscont*: vector of boolean values of same length of `length(size(pdf))` (number of random variables in original pdf). True if i-th variable is continuous, false if it is discrete.
%%% - *varargin*: vectors of $$x_i$$ used to sample the input pdf. Used for the integration of the pdf over the reduced variables. There should be a $$x_i$$ vector for each of the dimensions of the input joint pdf.
%%%
%%% ### Outputs:
%%% - *marginalized_pdf*: output joint (or marginal if survivingdims is not a vector) probability between variables in survivingdims.

assert(length(size(survivingdims)) <= 2, ...
    'survivingdims can only be a single value or a 1-D vector.');

assert(length(unique(survivingdims)) == length(survivingdims), ...
    'The variable indices in survivingdims cannot be repeated');

assert(length(varargin) == length(size(pdf)),...
    'The last n arguments should be equal to number of variables in pdf');

for i = 1:length(varargin)
    pdf_size = size(pdf);
    assert(pdf_size(i) == length(varargin{i}),...
        'Size of varargin(%d) does not match pdf size.',...
        'Expected number of elements: %d',length(varargin{i}),pdf_size(i));
    assert(isvector(varargin{i}),...
        'Each varargin should be a one dimensional vector ');
end

vars = 1:ndims(pdf);
for i=1:length(survivingdims)
    d = survivingdims(i);
    d_index = find(vars==d);
    vars([i,d_index]) = vars([d_index,i]);
    iscont([i,d_index]) = iscont([d_index,i]);
%     dx([i,d_index]) = dx([d_index,i]);
end

marginalized_pdf = permute(pdf,vars);
permuted_x = cell(1,length(varargin));
for i = 1:length(vars)
    permuted_x{i} = varargin{vars(i)};
end
for i=ndims(pdf):-1:length(survivingdims)+1
    if iscont(i)
        marginalized_pdf = trapz(permuted_x{i},marginalized_pdf,i);
    else
        marginalized_pdf = sum(marginalized_pdf,i);
    end
end
end

