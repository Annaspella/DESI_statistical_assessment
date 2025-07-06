%% Replicating the ANALYSIS FOR THE DESI Scores

%script that read data from XLS file and:
%1. compute some descriptive statistics separately for each indicator and
%2. draw histograms for each indicator
%3. trasform indicators

clear all
clc
opengl software %this command solves the problem with the video card that can cause figures to be black (all black)

% read dataset with variable names in first row
% country in 1st column
% region in 2nd column
% region ID in 3rd column
% indicator names from column 4th till the end
dataset=readtable('DESI_Y6.csv', ReadVariableNames=true, VariableNamingRule='preserve');
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
dataset_ni=readtable('new_indicators_egov.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_ni=dataset_ni.Properties.VariableNames;


dim=size(dataset_ni);
n_vars_ni=dim(2)-1;
n_countries_ni=dim(1); 

data_ni=table2array(dataset_ni(1:n_countries_ni,2:9));

% Table containing both original and new indicators
data_tot= [data data_ni];

% Merging toghether the variable names from 
% the first and the second table
var_names = [var_names_d var_names_ni(2:9)]; 
n_vars = n_vars_d + n_vars_ni;  % updating the number of total variable 33 + 
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
dataset2 = readtable('min_max.xlsx');
min_val = dataset2([1:33 37:44],"Min");
max_val = dataset2([1:33 37:44],"Max");
weights = dataset2([1:33 37:44],"Weight");  % original weight for each indicator 

% Load weights of each subdimension: contained in Weigths_Subdim.xlsx
dataset3 = readtable('Weigths_Subdim.xlsx');
weights_subdim = table2array(dataset3(:,"Weights"));


% computation of percentage of missing data  
perc_miss=ones(1,n_vars)*(-1);

for i=1:n_vars
    temp=find(isnan(data_tot(:,i)));
    temp1=size(temp);
    number_nan=temp1(1);
    perc_miss(1,i)=number_nan/n_countries*100;
end

fid = fopen('output_all.txt','w'); %open text file to put results
for i= 1:n_vars
% display results variable by variable
        fprintf(fid,'%s%s\n','results for: ',var_names{1,i+2});
        string=strcat({'These are results for indicator:  '},var_names(1,i+2));
        disp(string);
        string=strcat({'Percentage of missing values:  '},num2str(perc_miss(1,i)));
        fprintf(fid,'%s\n','percentage of missing values:');
        fprintf(fid,'%12.4f \n',perc_miss(1,i));
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
save DESI_std_minmax_DPS.txt data_std_minmax -ascii


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

%% Score calculations based on the proposed adjustment:
% DPS) Digital Public Services: extract the retained indicators from DPS +
%      append the new introduced 6 indicators

temp_DPS = data_std_minmax(:,[29:30 34:41]);     % (data_std_minmax contains the minmax values for the original indicators
                                                 %  + minmax values for the newly introduced 8 indicators)

% Doing the same with the weights
tt = table2array(weights([29:30 34:41],1));         % taking only the weights for the
                                                    % selected indicators + new indicators
                                                    % t =  repelem(1,8), EQUAL WEIGTHING
                                                    % (weights contains the weights of all the original indicators
                                                    %  
                                                    %  + the new 6 indicators of DPS)
g1_H = [1:7];          % NCB indicators wich positions are
                                                     
g2_H = [8:10];          % CB indicators

weights_DPS1_new = t(g1_H,:)/sum(t(g1_H,:));  
weights_DPS2_new = t(g2_H,:)/sum(t(g2_H,:));  
sub_DPS1_new = temp_DPS(:,g1_H)*weights_DPS1_new;
sub_DPS2_new = temp_DPS(:,g2_H)*weights_DPS2_new;

% Measuring internal consistency in DPS:
InternalDPS_across_new = cronbach_alpha([sub_DPS1_new sub_DPS2_new]);  % 0.7

DPS_Index_new_05 = 50*(sub_DPS1_new  + sub_DPS2_new); % 0.5 on each subdim and 100 to change scale

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
title('Digital Public Services- Ranking differences for Human Capital dimension')



% DESI LEVEL-----------------------------------------------------
% Score differences for DESI
diff_DPS_overall =  DESI_new_DPS_05 - DESI;

% ranking new DESI index with DPS changed, and DESI original
ranking_DPS_overall = get_rank(DESI_new_DPS_05)';
ranking_DESI = get_rank(DESI)';

% Difference in rankings for DESI_DPS
diff_rank_DESI_DPS = ranking_DESI-ranking_DPS_overall;

% Final results
DPS_results_overall = [country_names array2table(DESI) array2table(DESI_new_DPS_05) array2table(diff_DPS_overall) array2table(ranking_DESI) array2table(ranking_DPS_overall) array2table(diff_rank_DESI_DPS)];


% Plotting the results at DESI LEVEL---------------------------------------------------
% 1. Plot of score differences between DESI and DESI with human capital new subdimensions, 
%    keeping the other dimensions fixed
figure
labels = table2array(country_names);
scatter(DESI, DESI_new_DPS_05, 25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(DESI, DESI_new_DPS_05, labels,'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Human Capital: Original scores of DESI 2022 vs New scores')

% 2. Plotting the rankings differences for DESI and DESI with DPS changed
%    ranking DESI rankings differences (plot purpose)
[value_DESI_rank,p_diff_DESI_rank] = sort(diff_rank_DESI_DPS,'asc');
ranking_countries_DESI_rank = country_names(p_diff_DESI_rank,:);
labels = table2array(ranking_countries_DESI_rank(:,"Country"));
figure
barh(value_DESI_rank)
axis([-6 6 0 28])
dx = (value_DESI_rank)>=0;
text(value_DESI_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_DESI_rank) < 0;
text(value_DESI_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Human Capital- Ranking differences for DESI')


%% Plottin results obttained with different weights

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
writematrix(scores_tot,'DPS_scores.xls');

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

% Plotting the dimension level results
% 1. Plot of differences for human capital only,
%    comparing the reference indicator for DPS with the new one
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
title('Human Capital- Ranking differences for Human Capital dimension')


% Plotting OA_cit
table2array(country_names)' temp_DPS(:,[3 8])

% Malta
figure
MT = [temp_DPS(20,[3 8]); temp_DPS(20,[4 9]); temp_DPS(20,[6 10])];
bar(MT)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Malta');
set(gca,'xtick',[1:3],'xticklabel',["Online availability citizens"; "Online availability businesses"; "User support"]);

% Luxembourg
figure
LU = [temp_DPS(18,[3 8]); temp_DPS(18,[4 9]); temp_DPS(18,[6 10])];
bar(LU)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Luxembourg');
set(gca,'xtick',[1:3],'xticklabel',["Online availability citizens"; "Online availability businesses"; "User support"]);

% Belgium
figure
BE = [temp_DPS(2,[3 8]); temp_DPS(2,[4 9]); temp_DPS(2,[6 10])];
bar(BE)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Belgium');
set(gca,'xtick',[1:3],'xticklabel',["Online availability citizens"; "Online availability businesses"; "User support"]);


% Poland
figure
PL = [temp_DPS(22,[3 8]); temp_DPS(22,[4 9]); temp_DPS(22,[6 10])];
bar(PL)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Poland');
set(gca,'xtick',[1:3],'xticklabel',["Online availability citizens"; "Online availability businesses"; "User support"]);

% Serbia
figure
SE = [temp_DPS(25,[3 8]); temp_DPS(25,[4 9]); temp_DPS(25,[6 10])];
bar(SE)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Serbia');
set(gca,'xtick',[1:3],'xticklabel',["Online availability citizens"; "Online availability businesses"; "User support"]);

% Lituania
figure
LT = [temp_DPS(17,[3 8]); temp_DPS(17,[4 9]); temp_DPS(17,[6 10])];
bar(LT)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Lituania');
set(gca,'xtick',[1:3],'xticklabel',["Online availability citizens"; "Online availability businesses"; "User support"]);

% Romania:
figure
RO = [temp_DPS(24,[3 8]); temp_DPS(24,[4 9]); temp_DPS(24,[6 10])];
bar(RO)
axis([0 4 0 1.5]);
ylabel('Score value')
legend('National','Cross Border')
title('National vs CB indicators for Romania');
set(gca,'xtick',[1:3],'xticklabel',["Online availability citizens"; "Online availability businesses"; "User support"]);


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


