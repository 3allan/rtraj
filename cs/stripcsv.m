function newcsv = stripcsv(csv, cols, firstRow)
% 
% Matt Werner (m.werner@vt.edu) - Dec 9, 2020
% 
% Strip away all raw data from a .csv file except for those data existing
% in the columns specified by cols. The indicated columns from the .csv
% file will be placed in order of appearance into the new .csv file and 
% saved, ensuring that pre-existing files in the current working directory
% are not overwritten.
% 
%    Inputs:
% 
%               csv - Name of the .csv file from which to strip away
%                     columns of data. If not in the current working
%                     directory, CSV must lead to the path, either
%                     absolutely (/home/...) or relatively (./../...).
%                     Size: 1-by-1 (string)
%                     Units: ?
% 
%              cols - Indicated columns of data from the original .csv file
%                     to keep and place into the new .csv file in order of
%                     appearance.
%                     Size: 1-by-m (scalar)
%                     Units: - (N/A)
% 
%          firstRow - Indication as to which row in the csv file should be
%                     read first for each and every column. This quantity
%                     is helpful in the event that the csv file contains
%                     headers (text). In this case, the firstRow should be
%                     set to the row number at which the data actually
%                     begins to be listed. Otherwise, initial data may be
%                     simply skipped over.
% 

% Leave checks to Matlab

% Extract the .csv file's data
csvtmp = csvread(csv, firstRow, 0);

% Export only the requested columns.
newcsv = csvtmp(:, cols);