%% Replicating the ANALYSIS FOR THE DESI Scores

% script that read data from CSV file and:
%1. compute some descriptive statistics separately for each indicator and
%2. Standardize the data with min-max for aggregation purposes
%3. Replciate the scores results from DESI 2022 (with updated data)
%4. Calculate the scores for the new proposed adjustments
%5. Sensitivity analysis for the sub-dimensional weigthts for the major
%   scenario
%5. Plot country scores for National and CB

clear all
clc
opengl software %this command solves the problem with the video card that can cause figures to be black (all black)

% Read dataset with variable names in first row
% year in 1st column (equal to year 7 in all rows)
% country label in 2nd column
% indicator names from column 3rd till the end
dataset=readtable('DESI_Y7_2.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_d=dataset.Properties.VariableNames;
dim=size(dataset);   
n_vars_d=dim(2)-2;   % minus year and country columns
n_countries=dim(1)-1;  % Excluding EU
country_names=dataset(1:n_countries,2);

country_names=dataset(1:n_countries,2);

% save EU line and delete it from the table:
EU = dataset(dim(1),:);
data=table2array(dataset(1:n_countries,3:end)); % excluding EU from the calculations

% Add new indicators
dataset_ni=readtable('new_indicators_egov_2023.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_ni=dataset_ni.Properties.VariableNames;

dim=size(dataset_ni);
n_vars_ni=dim(2)-1;
n_countries_ni=dim(1); 

data_ni = table2array(dataset_ni(1:n_countries_ni,2:9));  % data_ni containes the new indicators

% Table containing both indicators form DESI 2022 and new indicators
data_tot= [data data_ni];

% Merging toghether the variable names from 
% the first and the second table
var_names = [var_names_d var_names_ni(2:9)]; 
n_vars = n_vars_d + n_vars_ni;   % updating the number of total variable 33 + 
                                 % the 8 new indicators

%The function grpstats computes different descriptive stats on sets of indicators by groups (for instance by group of countries)
%The instruction [], tells Matlab not to cluster data into any group
%grpstats treats NaNs and NaTs as missing values, and removes them
[average,sigma,skew]=grpstats(data_tot,[],{'mean','std','skewness'});
CV=sigma./average;                           
[max_data,ID_max]=max(data_tot);
country_max=dataset(ID_max,2);
country_max=table2cell(country_max);
country_max=country_max'; 
[min_data,ID_min]=min(data_tot);
country_min=dataset(ID_min,2);
country_min=table2cell(country_min);
country_min=country_min';

% Load min-max.xlsx which contains:
% min-max value for each indicator, to min max standardize the data
% weights for each indicator
dataset2 = readtable('min_max_2023.xlsx');
min_val = dataset2([1:33 37:44],"Min");
max_val = dataset2([1:33 37:44],"Max");
weights = dataset2([1:33 37:44],"Weight");  % original weight for each indicator 

% Load weights of each subdimension: contained in Weigths_Subdim.xlsx
dataset3 = readtable('Weigths_Subdim.xlsx');
weights_subdim = table2array(dataset3(:,"Weights"));

fid = fopen('output_all_2023.txt','w'); %open text file to put results
for i= 1:n_vars
% display results variable by variable
        fprintf(fid,'%s%s\n','results for: ',var_names{1,i+2});
        string=strcat({'These are results for indicator:  '},var_names(1,i+2));
        disp(string);
        string=strcat({'Mean value:  '},num2str(average(1,i)));
        fprintf(fid,'%s\n','mean value:');
        fprintf(fid,'%12.4f \n',average(1,i));
        disp(string);
        string=strcat({'Standard deviation (unbiased):  '},num2str(sigma(1,i)));
        fprintf(fid,'%s\n','standard deviation:');
        fprintf(fid,'%12.4f \n',sigma(1,i));
        disp(string);
        string=strcat({'Coefficient of variation:  '},num2str(CV(1,i)));
        fprintf(fid,'%s\n','coefficient of variation:');
        fprintf(fid,'%12.4f \n',CV(1,i));
        disp(string);
        string=strcat({'Skewness:  '},num2str(skew(1,i)));
        fprintf(fid,'%s\n','Skewness:');
        fprintf(fid,'%12.4f \n',skew(1,i));
        disp(string);
        string=strcat({'Maximum value:  '},num2str(max_data(1,i)));
        fprintf(fid,'%s\n','maximum');
        fprintf(fid,'%12.4f \n',max_data(1,i));
        disp(string);
        string=strcat({'Country corresponding to maximum value:  '},country_max(1,i));
        fprintf(fid,'%s\n','Country corresponding to maximum value:');
        fprintf(fid,'%s\n',country_max{1,i});
        disp(string);
        string=strcat({'Minimum value:  '},num2str(min_data(1,i)));
        fprintf(fid,'%s\n','minimum:');
        fprintf(fid,'%12.4f \n',min_data(1,i));
        disp(string);
        string=strcat({'Country corresponding to minimum value:  '},country_min(1,i));
        fprintf(fid,'%s\n','Country corresponding to minimum value:');
        fprintf(fid,'%s\n',country_min{1,i});
        disp(string);
        
        min_vec = table2array(min_val(i,1));
        max_vec = table2array(max_val(i,1));

        data_std_minmax(:,i) = (data_tot(:,i)-min_vec)./(max_vec-min_vec);
        
end %end of cycle on variables
fclose(fid);

%save matrix of standardized indicators: minmax 
save DESI_std_minmax_DPS_2023.txt data_std_minmax -ascii


%% Replicated analysis:
% Human capital-----------------------------------------
t = table2array(weights);
weights_HC1 = t(1:3,:)/sum(t(1:3,:));
weights_HC2 = t(4:7,:)/sum(t(4:7,:));

sub_HC1 = data_std_minmax(:,1:3)*weights_HC1;
sub_HC2 = data_std_minmax(:,4:7)*weights_HC2;

sub_HC = [sub_HC1 sub_HC2];

HC_Index = 100*sub_HC*weights_subdim(1:2);

% Connectivity-------------------------------------------
weights_C1 = t(8:10,:)/sum(t(8:10,:));
weights_C2 = t(11:13,:)/sum(t(11:13,:));
weights_C3 = t(14:16,:)/sum(t(14:16,:));
weights_C4 = t(17,:)/sum(t(17,:));

sub_C1 = data_std_minmax(:,8:10)*weights_C1;
sub_C2 = data_std_minmax(:,11:13)*weights_C2;
sub_C3 = data_std_minmax(:,14:16)*weights_C3;
sub_C4 = data_std_minmax(:,17)*weights_C4;

sub_C = [sub_C1 sub_C2 sub_C3 sub_C4];

C_Index = 100*sub_C*weights_subdim(3:6); 

% Integration of Digital Technology-----------------------
weights_I1 = t(18,:)/sum(t(18,:));
weights_I2 = t(19:25,:)/sum(t(19:25,:));
weights_I3 = t(26:28,:)/sum(t(26:28,:));

sub_I1 = data_std_minmax(:,18)*weights_I1;
sub_I2 = data_std_minmax(:,19:25)*weights_I2;
sub_I3 = data_std_minmax(:,26:28)*weights_I3;

sub_I = [sub_I1 sub_I2 sub_I3];

I_Index = 100*[sub_I1 sub_I2 sub_I1]*weights_subdim(7:9); 
 

% Digital public services--------------------------------
weights_DPS = t(29:33,:)/sum(t(29:33,:));
sub_DPS = data_std_minmax(:,29:33)*weights_DPS;

% Measuring internal consistency in the only dimension 
% of Digital public services:
InternalDPS_within = cronbach_alpha(data_std_minmax(:,29:33));  %0.8272 

DPS_Index = 100*weights_subdim(10)*sub_DPS;

% DESI final:

DESI = 0.25*(HC_Index+C_Index+I_Index+DPS_Index);

%% Score calculations based on the minor adjustment:

temp_DPS_min = data_std_minmax(:,[29:32 36:38]);     % (data_std_minmax contains the minmax values for the original 33 indicators
                                                     %  + minmax values for the newly introduced 8 indicators)
                                                     %  Selecting: 29:32) 4a1 new, 4a2, 4a3, 4a4 +
                                                     %             36:38) Mobile fr, User support National,
                                                     %                    Transparency
% Need to substitute User support National , with User support aggregated
% and standardize User support aggregated 
user_support = 0.5* data_ni(:,4) + 0.5*data_ni(:,8); % Obtained as 0.5 national user support and 0.5 cb user support
min_user = 0.75*min(user_support);
max_user = 1.25*max(user_support);
user_support_minmax = (user_support - min_user)/(max_user - min_user);
temp_DPS_m = [temp_DPS_min(:,[1:5]) user_support_minmax temp_DPS_min(:,7)]; 


% Doing the same with the weights
tt = table2array(weights([29:32 36:38],1));        % taking only the weights for the
                                                    % selected indicators + new indicators                                                 
                                                   

weights_DPS_new = tt/sum(tt);  
DPS_Index_new = 100*temp_DPS_m*weights_DPS_new;

DESI_new_DPS = 0.25*(HC_Index+C_Index+I_Index+DPS_Index_new);  % total DESI scores, changing DPS, 
                                                               % keeping the other dimensions fixed


total = [DPS_Index DPS_Index_new];                                                             
%save matrix of standardized indicators: minmax 
writematrix(total,'DPS_scores_minor_2.xls');

%% Scores excluding Open Data only (no new indicators)
temp_DPS_min = data_std_minmax(:,[29:32]);    

tt = table2array(weights([29:32],1)); 

weights_DPS_new = tt/sum(tt);  
DPS_Index_new_noOD = 100*temp_DPS_min*weights_DPS_new;

OD =  data_std_minmax(:,33);

corr(OD,DPS_Index_new_noOD);

%% Score calculations based on the MINOR adjustment DIVIDING 
%  Into Digital Public services dimension and other indicators

temp_DPS_min = data_std_minmax(:,[29:32 36:38]);     % (data_std_minmax contains the minmax values for the original 33 indicators
                                                     %  + minmax values for the newly introduced 8 indicators)
                                                     %  Selecting: 29:32) 4a1 new, 4a2, 4a3, 4a4 +
                                                     %             36:38) Mobile fr, User support National,
                                                     %                    Transparency
% Need to substitute User support National , with User support aggregated
% and standardize User support aggregated 
user_support = 0.5* data_ni(:,4) + 0.5*data_ni(:,8); % Obtained as 0.5 national user support and 0.5 cb user support
min_user = 0.75*min(user_support);
max_user = 1.25*max(user_support);
user_support_minmax = (user_support - min_user)/(max_user - min_user);
temp_DPS_m = [temp_DPS_min(:,[1:5]) user_support_minmax temp_DPS_min(:,7)]; 


tt = table2array(weights([29:32 36:38],1));         % taking only the weights for the
                                                    % selected indicators + new indicators 

g1_H = [3 4];          % Digital public services sub-dimension
                                                     
g2_H = [1:2 5:7];         % Other indicators

weights_DPS1 = t(g1_H,:)/sum(t(g1_H,:));  
weights_DPS2 = t(g2_H,:)/sum(t(g2_H,:));  
sub_DPS1 = temp_DPS_m(:,g1_H)*weights_DPS1;
sub_DPS2 = temp_DPS_m(:,g2_H)*weights_DPS2;

% Measuring internal consistency in DPS:
corr_DPS_sub = corr(sub_DPS1, sub_DPS2);  

DPS_Index_new_sub = 100*(2/3*sub_DPS1 + 1/3*sub_DPS2); % 2/3 assigned to DPS sub, 1/3 assigned to the rest of indicators

DESI_new_DPS_sub = 0.25*(HC_Index+C_Index+I_Index+DPS_Index_new_sub);  % total DESI scores, changing DPS, 
                                                                       % keeping the other dimensions fixed                                                


total = [DPS_Index DPS_Index_new DPS_Index_new_sub];                                                             
%save matrix of standardized indicators: minmax 
writematrix(total,'DPS_scores_minor_3.xls');


%% Score calculations based on the MAJOR adjustment:
% DPS) Digital Public Services: extract the retained indicators from DPS (29:30)+
%      append the new introduced 6 indicators (34:41)

temp_DPS = data_std_minmax(:,[29:30 34:41]);     % (data_std_minmax contains the minmax values for the original 33 indicators
                                                 %  + minmax values for the newly introduced 8 indicators)
                                                 % Selecting: 29:30) 4a1 new, 4a2,
                                                 %            34:37) Online av Citiz, Online av Bus, Mob friend, User supp Nat
                                                 %            38:41) Transparency, CB OA cit, CB OA bus, CB User supp

% Doing the same with the weights
tt = table2array(weights([29:30 34:41],1));         % taking only the weights for the
                                                    % selected indicators + new indicators
                                                    % t =  repelem(1,8), EQUAL WEIGTHING
                                                    % (weights contains the weights of all the original indicators
                                                    %  
                                                    %  + the new 6 indicators of DPS)
g1_H = [1:7];          % NCB indicators wich positions are
                                                     
g2_H = [8:10];         % CB indicators

weights_DPS1_new = t(g1_H,:)/sum(t(g1_H,:));  
weights_DPS2_new = t(g2_H,:)/sum(t(g2_H,:));  
sub_DPS1_new = temp_DPS(:,g1_H)*weights_DPS1_new;
sub_DPS2_new = temp_DPS(:,g2_H)*weights_DPS2_new;

% Measuring internal consistency in DPS:
Corr_DPS_new = corr(sub_DPS1_new, sub_DPS2_new);  % 

DPS_Index_new_05 = 50*(sub_DPS1_new + sub_DPS2_new); % 0.5 on each subdim and 100 to change scale

DESI_new_DPS_05 = 0.25*(HC_Index+C_Index+I_Index+DPS_Index_new_05);  % total DESI scores, changing DPS, 
                                                                     % keeping the other dimensions fixed

% Changing weights:
% 60% NCB, 40% CB
DPS_Index_new_06 = 60*sub_DPS1_new  + 40*sub_DPS2_new; 

DESI_new_DPS_06 = 0.25*(HC_Index+C_Index+I_Index+DPS_Index_new_06);  


% 70% NCB, 30% CB
DPS_Index_new_07 = 70*sub_DPS1_new  + 30*sub_DPS2_new; 

DESI_new_DPS_07 = 0.25*(HC_Index+C_Index+I_Index+DPS_Index_new_07);  


% 80% NCB, 20% CB
DPS_Index_new_08 = 80*sub_DPS1_new  + 20*sub_DPS2_new; 

DESI_new_DPS_08 = 0.25*(HC_Index+C_Index+I_Index+DPS_Index_new_08);    

%% Compare: Index for DPS at 2 levels, EQUAL WEIGHTING of SUBDIMENSIONS
%           - DIMENSION: Digital Public Services
%           - OVERALL: Changing DPS dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for DPS: both on the subdimension and on the total
% score

% DIMENSION LEVEL: equal weigting for NCB and CB---------------------------
% Score differences for DPS
diff_DPS_dim = DPS_Index_new_05 - DPS_Index;

% ranking new DPS index, and original DPS index
ranking_DPS_new = get_rank(DPS_Index_new_05)';
ranking_DPS = get_rank(DPS_Index)';

% Difference in rankings for DPS
diff_rank_DPS = ranking_DPS - ranking_DPS_new;

% Final results:
DPS_results_dim = [country_names array2table(DPS_Index) array2table(DPS_Index_new_05) array2table(diff_DPS_dim) array2table(ranking_DPS) array2table(ranking_DPS_new) array2table(diff_rank_DPS)];

% Plotting the dimension level results
% 1. Plot of differences for Digital Public services only,
%    compating the reference indicator for DPS with the new one
figure
labels = table2array(country_names);
scatter(DPS_Index, DPS_Index_new_05,25,'b','o')
axis([20 100 20 100])
hline = refline([1 0]);
set(hline,'Color','r')
text(DPS_Index, DPS_Index_new_05, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Digital Public Services: DPS original scores from DESI 2022 vs New scores')

% Correlation between new old scores
corr(DPS_Index, DPS_Index_new_05)

% 2. Plotting the rankings differences for DPS only
% ranking DPS rankings differences
[value_DPS_rank,p_diff_DPS_rank] = sort(diff_rank_DPS,'asc');
ranking_countries_DPS_rank = country_names(p_diff_DPS_rank,:);

labels = table2array(ranking_countries_DPS_rank(:,"Country"));
figure
barh(value_DPS_rank)
axis([-8 8 0 28])
dx = (value_DPS_rank)>=0;
text(value_DPS_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_DPS_rank) < 0;
text(value_DPS_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Digital Public Services- Ranking differences for DPS dimension')


%% Plottin results obtained with different weights
% Scores:
ts0 = timeseries(DPS_Index);
ts1 = timeseries(DPS_Index_new_05);
ts2 = timeseries(DPS_Index_new_06);
ts3 = timeseries(DPS_Index_new_07);
ts4 = timeseries(DPS_Index_new_08);

labels=table2array(country_names(:,1));
figure
plot(ts0, '-*')
hold on
plot(ts1, '-o')
hold on
plot(ts2)
hold on
plot(ts3)
hold on
plot(ts4)
axis([0 26 20 100])
xlabel('Country')
ylabel('DPS score value')
legend('DESI 2022','equal weighting','60% NCB, 40% CB','70% NCB, 30% CB','80% NCB, 20% CB')
text(0:26, DPS_Index_new_05, labels,'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
title('Effect of different weighting set for CB and NCB sub-dimensions')
set(gca,'xtick',[])

scores_tot = [DPS_Index_new_05 DPS_Index_new_06 DPS_Index_new_07 DPS_Index_new_08 DPS_Index];
writematrix(scores_tot,'DPS_scores_Y7.xls');

%% Ranking differences:
% Score differences for DPS:

% -- Equal weighting
diff_DPS_dim_05 = DPS_Index_new_05 - DPS_Index;

% ranking new DPS index, and original DPS index
ranking_DPS_new_05 = get_rank(DPS_Index_new_05)';
ranking_DPS_05 = get_rank(DPS_Index)';

% Difference in rankings for DPS
diff_rank_DPS_05 = ranking_DPS_05 - ranking_DPS_new_05;

% -- 60, 40:
% Score differences for DPS
diff_DPS_dim_06 = DPS_Index_new_06 - DPS_Index;

% ranking new DPS index, and original DPS index
ranking_DPS_new_06 = get_rank(DPS_Index_new_06)';
ranking_DPS_06 = get_rank(DPS_Index)';

% Difference in rankings for DPS
diff_rank_DPS_06 = ranking_DPS_06 - ranking_DPS_new_06;

% -- 70, 30:
% Score differences for DPS
diff_DPS_dim_07 = DPS_Index_new_07 - DPS_Index;

% ranking new DPS index, and original DPS index
ranking_DPS_new_07 = get_rank(DPS_Index_new_07)';
ranking_DPS_07 = get_rank(DPS_Index)';

% Difference in rankings for DPS
diff_rank_DPS_07 = ranking_DPS_07 - ranking_DPS_new_07;

% -- 80, 20:
% Score differences for DPS
diff_DPS_dim_08 = DPS_Index_new_08 - DPS_Index;

% ranking new DPS index, and original DPS index
ranking_DPS_new_08 = get_rank(DPS_Index_new_08)';
ranking_DPS_08 = get_rank(DPS_Index)';

% Difference in rankings for DPS
diff_rank_DPS_08 = ranking_DPS_08 - ranking_DPS_new_08;

%% Comparison of OA_city, OA_bus and User_support, divided into national and CB
%  For the 3 countries with the highest change in scores and the 3
%  countries with the lowest change in score


% Plotting the dimension level differences score differences for DPS
diff_DPS = DPS_Index_new_05 - DPS_Index;
[value_DPS,p_diff_DPS] = sort(diff_DPS,'asc');
ranking_countries_DPS_rank = country_names(p_diff_DPS,:);

labels = table2array(ranking_countries_DPS_rank(:,"Country"));
figure
barh(value_DPS)
axis([-15 15 0 28])
dx = (value_DPS)>=0;
text(value_DPS(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_DPS) < 0;
text(value_DPS(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Digital Public Services- Score differences for DPS dimension')


% Plotting OA_cit
%table2array(country_names)' temp_DPS(:,[3 8])

% Malta
figure
MT = [temp_DPS(20,[3 8]); temp_DPS(20,[4 9]); temp_DPS(20,[6 10])];
bar(MT)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Malta');
set(gca,'xtick',[1:3],'xticklabel',["Public services for citizens"; "Public services for businesses"; "User support"]);

% Luxembourg
figure
LU = [temp_DPS(18,[3 8]); temp_DPS(18,[4 9]); temp_DPS(18,[6 10])];
bar(LU)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Luxembourg');
set(gca,'xtick',[1:3],'xticklabel',["Public services for citizens"; "Public services for businesses"; "User support"]);

% Belgium
figure
BE = [temp_DPS(2,[3 8]); temp_DPS(2,[4 9]); temp_DPS(2,[6 10])];
bar(BE)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Belgium');
set(gca,'xtick',[1:3],'xticklabel',["Public services for citizens"; "Public services for businesses"; "User support"]);


% Poland
figure
PL = [temp_DPS(22,[3 8]); temp_DPS(22,[4 9]); temp_DPS(22,[6 10])];
bar(PL)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Poland');
set(gca,'xtick',[1:3],'xticklabel',["Public services for citizens"; "Public services for businesses"; "User support"]);

% EL
figure
EL = [temp_DPS(9,[3 8]); temp_DPS(9,[4 9]); temp_DPS(9,[6 10])];
bar(EL)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Greece');
set(gca,'xtick',[1:3],'xticklabel',["Public services for citizens"; "Public services for businesses"; "User support"]);

% FR
figure
FR = [temp_DPS(12,[3 8]); temp_DPS(12,[4 9]); temp_DPS(12,[6 10])];
bar(FR)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for France');
set(gca,'xtick',[1:3],'xticklabel',["Public services for citizens"; "Public services for businesses"; "User support"]);



