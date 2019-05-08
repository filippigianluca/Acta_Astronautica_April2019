function [data,f] = readFileImportdata(name_file)

% Function to open a file that contains, together with the data, comments
% and other strings. The datas are loaded in the field .data while the
% comments in the field .textdata.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the file

f=0;
if exist(name_file, 'file')
 dati = importdata(name_file);
 data = dati.data;
 f = 1;
 disp(['Reading from file ' name_file ' successfull']);
else
  disp(['File ' name_file ' not found']);       
end

    
end
