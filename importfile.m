
%% Importing data.txt file
function d = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  D = IMPORTFILE(FILENAME) reads data from text file FILENAME for the
%  default selection.  Returns the data as column vectors.
%
%  D = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  d = importfile("/Users/amrinaferdous/Desktop/Tools/Git/Frontier/data/data.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 04-Oct-2019 00:17:48

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 1
    dataLines = [1, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "d";
opts.VariableTypes = "double";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
d = tbl.d;
end