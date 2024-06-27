function params = addNetworkModelParameter(params, setup, varargin)
% Simple utility function for adding a new graph parameter to a list of
% parameters.
% See  NetworkModelParameter.m for more information about the parameters class.
%  
% Modified from the ModelParameter equivalent, addParameter.m

p = NetworkModelParameter(setup, varargin{:});
if isempty(params)
    params = {p};
elseif iscell(params)
    params = [params, {p}];
elseif isa(params, class(p))
    params = [{params}, {p}];
else
    error('Unknown format of input ''params''');
end
end
