function [num, text, raw,fileS, scales] = getFiles(PathName) 

revFile = 'reverseFrames.xlsx'; % file for coordinates correction of corrupted vids
[num, text, raw] = xlsread([PathName,revFile]);

scalesFile = fopen(strcat(PathName,'vidScales.txt')); % file of video scaling
C_text = textscan(scalesFile,'%s',2,'Delimiter',',');
C_data = textscan(scalesFile,'%s %n', 'Delimiter',',');
fclose(scalesFile);
fileS = C_data{1};
scales = C_data{2};

end