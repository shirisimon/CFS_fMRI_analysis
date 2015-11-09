function [bv, vmult] = tal2bv(tal, offset, res)
% tal2bv  - converting TAL system coordinates to BV system
%
% FORMAT:       [bv, vmult] = tal2bv(tal [, offset [, res]])
%
% Input fields:
%
%       tal         TAL compliant coordinates
%       offset      1x3 coordinate offset (default: [-128, -128, -128])
%       res         1x3 resolution (default: [1, 1, 1])
%
% Output fields:
%
%       bv          BV system coordinates
%       vmult       multiplication matrix
 
% Version:  v0.5c
% Build:    6120415
% Date:     Dec-04 2006, 3:15 PM CET
% Author:   Jochen Weber, Brain Innovation, B.V., Maastricht, NL
% URL/Info: http://wiki.brainvoyager.com/BVQXtools
 
% argument check
if nargin < 1 || ...
   ~isa(tal, 'double') || ...
   ~any(size(tal) == 3)
    error( ...
        'BVQXtools:BadArgument', ...
        'Bad or missing first argument.' ...
    );
end
if size(tal, 2) == 3 && ...
    size(tal, 1) ~= 3
    tal = tal';
end
if nargin < 3 || ...
   ~isa(res, 'double') || ...
   (numel(res) ~= 1 && numel(res) == 3) || ...
    any(isinf(res) | isnan(res))
    res = [1, 1, 1];
elseif numel(res) == 1
    res = [res, res, res];
else
    res = res(:)';
end
if nargin < 2 || ...
   ~isa(offset, 'double') || ...
    numel(offset) ~= 3 || ...
    any(isinf(offset) | isnan(offset))
    offset = [-128, -128, -128] .* res;
else
    offset = offset(:);
end
 
% make 4x4 compliant and prepare mult matrix
tal(4, :) = 1;
vmult = [ ...
            0 -1/res(3)         0  -offset(1)/res(3); ...
            0         0 -1/res(2)  -offset(2)/res(2); ...
    -1/res(1)         0         0  -offset(3)/res(1); ...
            0         0         0          1];
 
% multiply coordinates
bv = vmult * tal;
bv(4, :) = [];