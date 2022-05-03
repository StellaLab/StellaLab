clear all; close all; clc
%% Part 1 Read in and process initial mass and time data
% Pre.txt is the raw pre-reading data and post.txt for the post-reading
% data. These will need to line up with an additional conditions file to
% append data appropriately for analysis

fds = fileDatastore('*pre.txt', 'ReadFcn', @importdata); 
fds2 = fileDatastore('*post.txt', 'ReadFcn', @importdata);
fds3 = fileDatastore('*cdns.xlsx','ReadFcn', @importdata);
fds4 = fileDatastore('Graph_112_ATP_Other.xlsx','ReadFcn', @importdata);
fullFilePre = fds.Files;
fullFilePost = fds2.Files;
fullFileCdns = fds3.Files;
fullFileGraph = fds4.Files;

numGraph = length(fullFileGraph);                                          %This is the number graph excel files. 
numFiles = length(fullFilePre);                                            %This is the number of post and pre data

datetoday = datestr(datetime('today'));                                    % We can attach this varible to any output name so that it has the date you ran this analyses on it.

% This loop currently does not work. In theory it would go through all the
% excel files and graph them all. Currently there is a bug that it can only
% do it once. Once this is fixed though, the code will generate all the
% analyses needed according to this loop.


for i = 1:numGraph
graph = read(fds4);                                                        % Read the graph excel file
cdnsName = graph.textdata(:,1);                                                 % The text data of this file are the conditons we need to find/graph
cdnsnumber = length(cdnsName);                                             % Find the length of the conditions for this graph
indicator = graph.textdata(:,2);                                               % Indicator will help us find what pre/post data files contain the right data we need
indicator = indicator(~cellfun('isempty',indicator))
cell_dat_split = cellfun(@(x)regexp(x,',','split'),indicator,'UniformOutput',0)
basilIndicator = graph.data(:,1);                                          % basilIndicator will tell us what conditions need to be subtracted from each other for the basal calculations

% This section is to find the name of this graph, i.e. Graph_20.xlsx will
% take Graph_20 to name output files with the same prefix Graph_20
graphTitle = fds4.Files{i};
nameGraph = strsplit(graphTitle, '\');
nameGraph  = string(nameGraph{7});
nameGraph = strsplit(nameGraph, '.xlsx');
title_name = nameGraph(1);

%This inital loop is to create text files that are associated with each
%individual experiment. In the text file there will be data in x columns
%for x conditions prompted. These files will go through the second part of
%this program with a series of analysis.


for k = 1:numFiles

    pre = read(fds);
    post = read(fds2);
    cdns = read(fds3);
    
    % Reading in the data
    x = pre.data;
    y = post.data;
    [rownum,colnum]=size(y);                                                    % Call the number of columns & rows
    
    %With how the reader sets up the columns. We want to take data from
    %array 3 to end as the numbers we'll be doing elemetary calculations
    %on.
    f= x(:,3:end);
    g= y(:,3:end);
    for i=1:colnum
          Mn = mean(f, 1);                                                  % Loop through to take the mean of the first four rows
          B = (g - Mn)./Mn;                                                 % Take A from columns 3 on to subtract Mn and divide.
    end
    
    %We read the time from array 4 to end - 2 and indicate the format so it
    %can read it as a datetime
    T = post.textdata(4:(end-2),1);
    Time = datetime(T,'InputFormat','mm:ss');
    
    %A small glitch not sure how to get around in a more efficent way.
    %Basically sometimes the last file has a time that is LESS than the max
    %time. So I won't be able to create PrismTable because of this, it
    %needs the max time so that I can create a buffered table etc. SO
    %what I did to fix this was grab a specific file for the time to
    %reference later.
    if k == 23
        Time2 = datetime(T,'InputFormat','mm:ss');
    end
    
    %Change the \ to / if running on a mac. This is to grab the appropriate
    %name as x so that I can name the output text file with the prefix of
    %the date
    x = fds.Files{k};
    name = strsplit(x, '\');
    name  = string(name{7});
    name = strsplit(name, '_');
    x = name(1);

     
%% Plot Prompt 
    B( :, all( isnan( B ), 1 ) ) = [];                                      %If there are nan values replace them with a blank
    C = reshape(B,rownum,3,[]);                                             %Reshape to take a mean between every 3rd column, this is how the reader
    D = mean(C,2);                                                          %Organizes readings. Seems to conventional.
    E = reshape(D,rownum,[]);
    
    [C,ia,idx] = unique(cdns(:,1),'stable');                                %Cool part! If there are repeating conditions we can combine them
    idxLength = max(size(idx));                                             %We use these to find the indexes of where the repeating condtions are
    
    for m = 1:rownum
    tableEE(:,m) = accumarray(idx,E(m,:),[],@mean);                         %Creating a new table where the repeating conditons are averaged and representative in a new
                                                                            %tableEE stored as E
    end
    E = tableEE';
    %% if cdns repeats take corresponding values in table E and average, this is new table E
    
    tableE = array2table(E);                                                %Now we take array to table and find the variable names
    vars = 1:width(tableE);                                                 %Renaming the columns in the table with the variable names(This is created by appending
    tableE = renamevars(tableE,vars,C');                                    %The conditions that the user inputed as cdns                     
    Date = table(Time);                                                     %Earlier we found Time as a datetime and we will attach this to the table
    PrismTable = [Date,tableE];                                             %Now we have a table fit to input data into PRISM
    clear tableEE

%     writetable(PrismTable, strcat(x,'_','PrismData.xlsx'));


       status = mkdir(title_name);                                          %Find the name so that later we can just pull all these into a new folder with
       status2 = fullfile(title_name);                                      %The title name that user was prompted for
        %% 
       D=zeros(max(size(T)),cdnsnumber);                                    %Create empty array to pad any data that is less than the max length
       D2=zeros(max(size(T)),1);                                            %We want to do this so we can output the data all as the same length otherwise
                                                                            %Matlab will freak out
      
      % As we loop throguh all the post and pre data, in one file you will
      % Go through all the conditions and see if the post/pre data has
      % those conditons, if yes store them in that appropriate column 
      % that lines up with the condition, if it doesn't have data fill that
      % column with zeroes. Output is a data excel file that will be used
      % for later calculations.                                                                     
       for o = 1:cdnsnumber
       index = find(strcmp(cdns,cdnsName(o)));                              %This is to go through and place zeroes where there is no data
       
       if index >= 1
           cdnsnewName{o}= cdnsName(o);
           cdnsnewName = [cdnsnewName{:}];
           T2 = table2array(PrismTable(:,cdnsnewName(:)))
           clearvars cdnsnewName;
           D(:,o) = T2;
       end 
       if isempty(index); 
           D(:,o) = D2
           continue
       end 
       end
       
        
        writematrix(D, strcat(x,'_',title_name,'data.txt'));
 
end


%% Final calculations after parsing through all the data
        folder = pwd;
        txtFiles = dir(fullfile(folder,'*data.txt'));                      %load in/save text files saved from previous loop into struct
        numFiles;                                                          %number of files loaded

clearvars myData;

for m = 1:numFiles                                                         %myData holds all the imported data
     tx = txtFiles(m).name;
     myData{m} = [importdata(tx)];
end

% cnm = cdnsnumber - 4;                                                      % At the moment,cnm represents the amount of conditions (only in the plateBASE code)  
                                                                           % that need to be present in an experiment file to use it for the correct calculations


%% For loop with if statements that will store data from all experiments that pass through the loop
countCondTable = zeros(cdnsnumber,1);                                      % Make a table that can store the count of the conditions, while going through this loop 
                                                                           % it will store the numbers in this table. Good to referece to know if the code is parsing through data correctly.
% indicatorList = [1:max(indicator)];                                        % Indicator variable will have numbers to group conditions together, conditions that need to be subtracted from each 
                                                                           % each other. We need to find the individual numbers so i.e. 1,2,3,4 when in reality the indicator list will have 1,1,2,2, etc
                                                                           % This will make a for loop below function more efficiently.
for j = 1:numFiles
    myDatas = cell2mat(myData(j));
%     Ddatas{j} = zeros(181,cdnsnumber);                                                                       
                                                                           %This if statement is to only take the data if it has present, cnm number of nonzero columns, equivilent to say cnm number of conditions data.
            for i = 1:length(cell_dat_split)
                idnx = str2double(cell_dat_split{i,1})
                lengthIdnx = length(idnx)
                if sum(all(myDatas(:,idnx) ~= 0)) == lengthIdnx            % If the sum of myDatas at this finalIdnx == 4 then proceed
                    
                    Dsize = 181-max(size(myDatas));                         %So if the data has less than cnm columns == 0, as in it has data for DMSO and another conditons
                    for i = 1:lengthIdnx
                        indexx = idnx(i);
                        Ddatas{j}(:,indexx) = padarray(myDatas(:,indexx),[Dsize 0],0,'post');      %Take it and pad the array for averaging later. But if it only has DMSO as in ==0 is greater than cnm
                                                                         %Skip. We don't want to average controls in an experiment that didn't run the conditions we're interested in.
                                                                           %All data that passed the above conditons will be stored in Ddatas
                    end

                end   
                
            end
          
        
            
%          if sum(all(myDatas == 0)) >= cnm
%             continue
%          end
            
        
        

end

for i = 1:length(Ddatas)
    if isempty(Ddatas{i}) == 0
                datassize = size(Ddatas{i},2)
                if datassize < cdnsnumber            
                    Ddatas{i} = padarray(Ddatas{i},[0,cdnsnumber-datassize],0,'post')
                end
    end
end
%% 
countCondTable = zeros(cdnsnumber,1);
for i = 1:length(Ddatas)
    if isempty(Ddatas{i})== 0
                for j = 1:cdnsnumber
                        if sum(Ddatas{i}(:,j) ~= 0);                             %Loop to count the conditions, this only works if it passed through the previous conditions.
                            countCondTable(j)= countCondTable(j)+ 1;
                        end
                end
    end
end
%% Find Basal

% This loop will go through all the Ddatas to calculate the basal. 

for i = 1:length(Ddatas)
    cdnsData = Ddatas(i)
    if isempty(cdnsData{:}) == 0                                           % Make sure the data is not empty, if it is it won't pass the if statement    
    cdnsData = cell2mat(cdnsData);                                         % Make the data into numbers, so that we can do calculations
    for j = 1:cdnsnumber
        idnx = basilIndicator(j);                                          % Same as before, find the index so that you can index what columns to subtract from each other for basal.
        if idnx == 1                                                       % We have two conditions here, the first is if the index is 1. 
            basil = cdnsData(:,j);                                         % basil will be the cdnsData at j in the loop, and Cdn will be the idnx (indexed) data for what to subtract basil from.
            Cdn = cdnsData(:,idnx);

            t = basil < 0 ;                                                % For where the data is negative, do the follow calculation and store it in the variable basil
            basil(t) = basil(t) - Cdn(t);

            t2 = basil > 0;                                                % For where the data is positive, ""
            basil(t2) = basil(t2) - Cdn(t2);

            cdnsNewData(:,j) = basil;                                      %cdnsNewData will store these new calculations and these will be used for the else statement below to find the basal for index outside of 1
            cdnsNewData(:,1) = 0;                                          % DMSO column will always be zero, unfortunately, the way we do if >0 or < 0 subtact/add it makes the DMSO column not zero.
        else 
            basil = cdnsData(:,j);                                         % Everything is the same as the statment above except now you will subtact from the indexed cdnsNewData stored in Cdn.
            Cdn = cdnsNewData(:,idnx);

            t = basil < 0 ;
            basil(t) = basil(t) - Cdn(t);

            t2 = basil > 0;
            basil(t2) = basil(t2) - Cdn(t2);

            cdnsNewData(:,j) = basil;
        end
    end
    
    masterBasil{:,i} = cdnsNewData;                                        % This variable is where we store all our basal caculations after going through the whole length of Ddatas.                            
    end
end

%AVG of Basal
myData = cell2mat(masterBasil);
avgDatas = reshape(myData,181,cdnsnumber,[]);
avgData = sum(avgDatas,3) ./ sum(avgDatas~=0,3);
basalDataTable = cell2table(num2cell((avgData)));


%AVG of all data
myData2 = cell2mat(Ddatas);
avgDatas2 = reshape(myData2,181,cdnsnumber,[]);
avgData2 = sum(avgDatas2,3) ./ sum(avgDatas2~=0,3);
avgDataTable = cell2table(num2cell((avgData2)));
avgDataTable.Properties.VariableNames = cdnsName;

%%find 
n = array2table(countCondTable);

%find SEM
SD = std(avgDatas,[],3);
sdlength = max(size(SD));
SEM = SD./sqrt(sdlength);
SEMtable = cell2table(num2cell(SEM));

%% Final Output Table with BASIL and SEM columns
%Name Columns of Table
for t = 1:cdnsnumber
    nName(t,:) = [strcat("N count",num2str(t))];
    semName(t,:) = [strcat("SEM", num2str(t))];
    basilName(t,:) = [strcat("Basil", num2str(t))];

end

%Assign names to table
% n.Properties.VariableNames = nName;
SEMtable.Properties.VariableNames = semName;
basalDataTable.Properties.VariableNames = basilName;

%Create a final data table with the columns you'd like for each cndn
for t = 1:cdnsnumber
    dataTable{t} =[avgDataTable(:,t) basalDataTable(:,t) SEMtable(:,t)];
end


datasTable = cat(2, dataTable{:});                                         %This will just make it so the datatable is now a table and not made up of cells.

%% DOSE RESPONSE , output a kinetic table that will be up to the time point specified by Nephi. Input it in where you see kineticmin.
%Dose Response AVG Table
drTime = table(Time2);
drTable = [drTime,datasTable];

kineticMin = ismember(Time2, datetime('00:30:00'));
kineticMin2 = find(kineticMin, 1);

kineticNew = drTable(1:kineticMin2,:);
writetable(kineticNew, strcat(title_name,'_','kinetic_', datetoday,'.xlsx')); %Output the kinetic table

%For mean of response and mean of SD new table
% drTable2 = [drTime,avgDataTable];

%% MOA
% To create a x-y minutes MOA table. Store index of where x minute is and 
% where y minute is and this will be used later to index for MOA tables 
% with different time points.

%% 4-5 MIN
fourMin = ismember(Time2, datetime('00:04:00'));
fourMin2 = find(fourMin, 1);
fiveMin = ismember(Time2, datetime('00:05:00'));
fiveMin2 = find(fiveMin, 1);
moaTable = drTable(fourMin2:fiveMin2,:);                                    %Don't really need these lines of code, similarly as below. Just in case you'd like to reference to make sure its including all data that needs to be averaged

%% 25-30 MIN
twentyfiveMin = ismember(Time2, datetime('00:25:00'));
twentyfiveMin2 = find(twentyfiveMin, 1);
thirtyMin = ismember(Time2, datetime('00:30:00'));
thirtyMin2 = find(thirtyMin, 1);
moaTable = drTable(twentyfiveMin2:thirtyMin2,:);
%% 6-9 MIN
sixMin = ismember(Time2, datetime('00:06:00'));
sixMin2 = find(sixMin, 1);
nineMin = ismember(Time2, datetime('00:09:00'));
nineMin2 = find(nineMin, 1);
moaNineTen = drTable(sixMin2:nineMin2,:);

%% 5-8 MIN
fiveMin = ismember(Time2, datetime('00:05:00'));
fiveMin2 = find(fiveMin, 1);
eightMin = ismember(Time2, datetime('00:08:00'));
eightMin2 = find(eightMin, 1);
moaTwentyNineThirty = drTable(fiveMin2:eightMin2,:);

%% Minute BASIL & MEAN Table
% Now take the indexes stored in the above variables, to index the
% basalDataTable we have. We take the mean. Index the SEM table similarly
% and take the mean, include the count that we stored in n as a table. Name
% columns as you'd like and then output the table as an excel.
%% 4-5
avgResponse = table(mean(basalDataTable{fourMin2:fiveMin2,:},1)');
avgResponse.Properties.VariableNames = {'Mean Basil'};

meanSEM = table(mean(SEMtable{fourMin2:fiveMin2,:},1)');
meanSEM.Properties.VariableNames = {'Mean(SEM)'};

n.Properties.VariableNames = {'Count'};

cdnsMOA = table(cdnsName);
cdnsMOA.Properties.VariableNames = {'Conditions'};

MOATable = [cdnsMOA avgResponse meanSEM n];
writetable(MOATable, strcat(title_name,'_','4-5-Min_', datetoday,'.xlsx'));

%% 25-30
avgResponse = table(mean(basalDataTable{twentyfiveMin2:thirtyMin2,:},1)');
avgResponse.Properties.VariableNames = {'Mean Basil'};

meanSEM = table(mean(SEMtable{twentyfiveMin2:thirtyMin2,:},1)');
meanSEM.Properties.VariableNames = {'Mean(SEM)'};

n.Properties.VariableNames = {'Count'};

cdnsMOA = table(cdnsName);
cdnsMOA.Properties.VariableNames = {'Conditions'};

MOATable = [cdnsMOA avgResponse meanSEM n];
 writetable(MOATable, strcat(title_name,'_','25-30-Min_', datetoday,'.xlsx'));
%% 4-5
% avgResponse = table(mean(basilTable{fourMin2:fiveMin2,:},1)');
% avgResponse.Properties.VariableNames = {'Mean Basil'};
% 
% meanSEM = table(mean(SEMtable{fourMin2:fiveMin2,:},1)');
% meanSEM.Properties.VariableNames = {'Mean(SEM)'};
% 
% n.Properties.VariableNames = {'Count'};
% 
% cdnsMOA = table(cdnsName);
% cdnsMOA.Properties.VariableNames = {'Conditions'};
% 
% MOATable = [cdnsMOA avgResponse meanSEM n];
% writetable(MOATable, strcat(title_name,'_','4-5-Min.xlsx'));

% %% 5-8
% avgResponse = table(mean(basilTable{fiveMin2:eightMin2,:},1)');
% avgResponse.Properties.VariableNames = {'Mean Basil'};
% 
% meanSEM = table(mean(SEMtable{fiveMin2:eightMin2,:},1)');
% meanSEM.Properties.VariableNames = {'Mean(SEM)'};
% 
% n.Properties.VariableNames = {'Count'};
% 
% cdnsMOA = table(cdnsName);
% cdnsMOA.Properties.VariableNames = {'Conditions'};
% 
% MOATable = [cdnsMOA avgResponse meanSEM n];
% writetable(MOATable, strcat(title_name,'_','5-8-Min.xlsx'));

%% Excel for slope
oneMin = ismember(Time2, datetime('00:01:00'));
oneMin2 = find(oneMin, 1);
threeMin = ismember(Time2, datetime('00:03:00'));

slope = table((avgDataTable{oneMin,:}- avgDataTable{threeMin,:})');
slope.Properties.VariableNames = {'Slope'};
slopeTable = [cdnsMOA slope]
writetable(slopeTable, strcat(title_name,'_','1-3-Min_Slope_', datetoday,'.xlsx'));
%% PLOT averages

%OG Avg Table
% dataTable = horzcat(dataTable{:});
% writetable(dataTable, strcat(title_name,'Avg_Data.xlsx'));


        fig = figure(1);
        lineProps.col{1}=[1 0 0];
        plot(Time2, avgData);


        hold on;
        title(strcat('Avg',title_name));
        legend(cdnsName,'Location','northeast');
        xlabel('Time (MM:SS)');
        hold off 
%         saveas(fig, strcat(title_name,'_',datetoday, '.pdf'));
%% Move all files containing title name to title name folder
title = convertStringsToChars(title_name)        
movefile(['*',title,'*'],status2);


end
