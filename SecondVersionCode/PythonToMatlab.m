[filename] = uigetfile('.txt','MultiSelect','on');

files=what('./'); 
fileList=files.mat;
subject=fileList{1};
subject={subject(1:end-4)};

for i=1:size(filename,2)
    [header,outmat] = JSONtxt2cell(filename{i});
    save(['PytonData' num2str(i) '.mat'],'header','outmat');
    clear header
    clear outmat
end


results=SyncPython_forStatic(subject{1});
