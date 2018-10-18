function c = struct2cellarg(s)
% function c = struct2cellarg(s)
%
% Converts struct into paired cell array {'field1', value1, 'field2', value2, ...}.
% that can be passed as varargin into mat2mod (which uses MIRT toolbox varargin parser).
%
% Input:
%   s     struct 
%           s.field1 = value1;
%           s.field2 = value2; etc
%
% Output:
%   c     paired cell array
%
% Basically does the opposite of built-in struct() function, which converts
% paired cell array {'field1', 'value1', 'field2', 'value2', ...} into a struct.
%
% Example usage:
%  >> lims.MaxSlew = 200;
%  >> lims.SlewUnit = 'T/m/s';
%  >> mat2mod('rf', rf, 'gx', gx, 'system', struct2cellarg(lims));
%
% $Id: struct2cellarg.m,v 1.2 2018/10/18 12:35:17 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/lib/v2/struct2cellarg.m,v $

fields = fieldnames(s);
vals = struct2cell(s);
c = [];
for ii = 1:length(fields)
	c = [c fields(ii) vals(ii)];
end
