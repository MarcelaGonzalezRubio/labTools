function [header,outmat] = JSONtxt2cell(filename)
%JSONtxt2cell uses json_parse to open a txt file and parse it into a cell
%array
%   ...
fid = fopen(filename);
bigstring = fread(fid);
fclose(fid);
% smallstring = char(bigstring(1:150)');
bigstring = char(bigstring');
%figure out how many "[" and "]" there are
[IB] = ismember(bigstring,'[');
[IC] = ismember(bigstring,']');
startindex = find(IB);
stopindex = find(IC);
outcell = cell(length(startindex),1);

%the header
header = JSON.parse(bigstring(startindex(1):stopindex(1)));

for z = 2:length(startindex)
    temp = bigstring(startindex(z):stopindex(z));
%     g = JSON.parse(temp);
    [ID] = ismember(temp,',');%find out how many items there are
    commaindex = find(ID);
    g{1} = str2double(temp(2:commaindex(1)-1));
    for zz = 2:length(commaindex)
       g{zz} = str2double(temp(commaindex(zz-1)+1:commaindex(zz)-1));
    end
    g{end+1} = str2double(temp(commaindex(end)+1:end-1));
    outcell{z-1} = g;
    clear g
end

outcell = outcell(~cellfun('isempty',outcell));
outmat = zeros(length(outcell),length(outcell{1}));

for z = 1:length(outcell)
   outmat(z,:) = cell2mat(outcell{z}); 
end
clear outcell

%save('Pyton.mat','header','outmat') 
end
