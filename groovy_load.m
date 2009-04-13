function [names, values] = groovy_load(filename, sep, missing)
% Load file with first text value and subsequent numeric values
% FORMAT [names, values] = groovy_load(filename, sep, missing)
%
% Inputs
% filename       - text file to load
% sep            - optional separator character (tab default)
% missing        - optional value to indicate missing input [NaN]
%
% Outputs
% names          - cell array of names loaded
% values         - matrix of numeric values, MxN, where M is the number
%                  of rows in the file, and N is the maximum number of
%                  numeric values loaded (can differ across lines)
%
% If there are different numbers of values per line, absent values for 
% columns where are there are less than N values will be filled with 
% values of 'missing' input

if nargin < 2
  sep = [];
end
if isempty(sep)
  sep = sprintf('\t');
end
if nargin < 3
  missing = [];
end
if isempty(missing)
  missing = NaN;
end
filename = deblank(filename);
fid = fopen(filename, 'rt');
if fid == -1
  error(['Cannot open ' filename]);
end
names = {};
values = [];

loop_i = 1;
while 1
  tline = fgetl(fid);
  if ~ischar(tline)
    break
  end
  fields = sf_split(tline, sep);
  if length(fields) == 1 % no sep
    continue
  end
  names{loop_i} = sscanf(fields{1}, '%s');
  numbers = [];
  for nno = 2:length(fields)
    numbers = [numbers sscanf(fields{nno}, '%f')];
  end
  column_diff = length(numbers)-size(values, 2);
  if column_diff < 0
    numbers = [numbers repmat(missing, 1, -column_diff)];
  elseif column_diff > 0
    values = [values repmat(missing, size(values,1), column_diff)];
  end
  values = [values; numbers];
  loop_i = loop_i + 1;
end
fclose(fid);
return

function fields = sf_split(str, sep)
% Splits line into files at sep
sep_inds = strfind(str, sep);
starts = [1 sep_inds+1];
stops = [sep_inds-1 length(str)];
fields = {};
for fno = 1:length(sep_inds)+1
  fields{fno} = str(starts(fno):stops(fno));
end
return
