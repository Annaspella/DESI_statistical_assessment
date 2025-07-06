%% Replicating the ANALYSIS FOR THE DESI Scores and HUMAN CAPITAL

% script that read data from CSV file and:
%1. compute some descriptive statistics separately for each indicator and
%2. Standardize the data with min-max standardization
%3. Replciate the scores results from DESI 2022
%4. Calculate the scores for the new proposed adjustments
%5. Plot the differences in scores and rankings 

clear all
clc
opengl software %this command solves the problem with the video card that can cause figures to be black (all black)

% Read dataset with variable names in first row
% year in 1st column (equal to year 6 in all rows)
% country label in 2nd column
% indicator names from column 3rd till the end
dataset=readtable('DESI_Y6.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_d=dataset.Properties.VariableNames;
dim=size(dataset);   
n_vars_d=dim(2)-2;    % minus year and country columns
n_countries=dim(1)-1;  % Excluding EU
country_names=dataset(1:n_countries,2);

% save EU line and delete it from the table:
EU = dataset(dim(1),:);
data=table2array(dataset(1:n_countries,3:end)); % excluding EU from the calculations

% Standardizing data
dataset_ni=readtable('new_indicators_HumanCapital.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_ni=dataset_ni.Properties.VariableNames;


new_ind = table2array(dataset_ni(1:n_countries, 2:3));

% Table containing both original and new indicators
data_tot= [data new_ind];

% Merging toghether the variable names from 
% the first and the second table
var_names = [var_names_d var_names_ni(2:3)]; 
n_vars = n_vars_d + 2;  % updating the number of total variable with 
                       % the 2 new indicators

% computation of percentage of missing data  
perc_miss=ones(1,n_vars)*(-1);

for i=1:n_vars
    temp=find(isnan(data_tot(:,i)));
    temp1=size(temp);
    number_nan=temp1(1);
    perc_miss(1,i)=number_nan/n_countries*100;
end

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
min_val = dataset2(:,"Min");
max_val = dataset2(:,"Max");
weights = dataset2(:,"Weight");  % original weight for each indicator 

% Load weights of each subdimension: contained in Weigths_Subdim.xlsx
dataset3 = readtable('Weigths_Subdim.xlsx');
weights_subdim = table2array(dataset3(:,"Weights"));


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

        if i ~= 34   % 34 corresponding to never used internet, opposite direction, need to standardize it differently
            data_std_minmax(:,i) = (data_tot(:,i)-min_vec)./(max_vec-min_vec);
        else
            data_std_minmax(:,i) = 1-((data_tot(:,i)-min_vec)./(max_vec-min_vec));  % attenzione, quantità cambiata, vedere risultati
        end
end %end of cycle on variables
fclose(fid);

%save matrix of standardized indicators: minmax 
save DESI_std_minmax_HC.txt data_std_minmax -ascii



%% Replicated analysis:
% Human capital
t = table2array(weights);
weights_HC1 = t(1:3,:)/sum(t(1:3,:));
weights_HC2 = t(4:7,:)/sum(t(4:7,:));

sub_HC1 = data_std_minmax(:,1:3)*weights_HC1;
sub_HC2 = data_std_minmax(:,4:7)*weights_HC2;

sub_HC = [sub_HC1 sub_HC2];

% Measuring internal consistency in HC:
InternalC_across_sq = cronbach_alpha(sub_HC);   % 0.7377

HC_Index = 100*(weights_subdim(1:2)'*sub_HC')';

% Connectivity
weights_C1 = t(8:10,:)/sum(t(8:10,:));
weights_C2 = t(11:13,:)/sum(t(11:13,:));
weights_C3 = t(14:16,:)/sum(t(14:16,:));
weights_C4 = t(17,:)/sum(t(17,:));

sub_C1 = data_std_minmax(:,8:10)*weights_C1;
sub_C2 = data_std_minmax(:,11:13)*weights_C2;
sub_C3 = data_std_minmax(:,14:16)*weights_C3;
sub_C4 = data_std_minmax(:,17)*weights_C4;

sub_C = [sub_C1 sub_C2 sub_C3 sub_C4];

C_Index = 100*(weights_subdim(3:6)'*sub_C')'; 

% Integration
weights_I1 = t(18,:)/sum(t(18,:));
weights_I2 = t(19:25,:)/sum(t(19:25,:));
weights_I3 = t(26:28,:)/sum(t(26:28,:));

sub_I1 = data_std_minmax(:,18)*weights_I1;
sub_I2 = data_std_minmax(:,19:25)*weights_I2;
sub_I3 = data_std_minmax(:,26:28)*weights_I3;

sub_I = [sub_I1 sub_I2 sub_I3];

I_Index = 100*(weights_subdim(7:9)'*sub_I')'; 
 

% Digital public services
weights_D = t(29:33,:)/sum(t(29:33,:));
sub_D = data_std_minmax(:,29:33)*weights_D;

D_Index = 100*weights_subdim(10)*sub_D;

% DESI final:

DESI = 0.25*(HC_Index+C_Index+I_Index+D_Index);

%% Score calculations based on the new subdimension:
% HC) Human capital: major adjustment-------------------------------------

t = table2array(weights);
g1_H = [1:3 6 34 35];                              % add never_used_intenet indicator 
                                                   % which position in 34
                                                   % and frequency_of_use
                                                   % which position at 35
g2_H = [4 7];
weights_HC1_new = t(g1_H,:)/sum(t(g1_H,:));  
weights_HC2_new = t(g2_H,:)/sum(t(g2_H,:));  
sub_HC1_new = data_std_minmax(:,g1_H)*weights_HC1_new;
sub_HC2_new = data_std_minmax(:,g2_H)*weights_HC2_new;

% Measuring internal consistency in HC:
InternalC_across_new = cronbach_alpha([sub_HC1_new sub_HC2_new]);  % 0.7689

HC_Index_new = 50*(sub_HC1_new  + sub_HC2_new); % 0.5 on each subdim and 100 to change scale

DESI_new_HC = 0.25*(HC_Index_new+C_Index+I_Index+D_Index);  % total DESI scores, changing HC, 
                                                            % keeping the other dimensions fixed
 


              %% Compare: Index for Human capital at 3 levels
%           - SUBDIMENSION: Basic skills and Advanced skilles
%           - DIMENSION: Human capital
%           - OVERALL: Changing Human Capital dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for HC: both on the subdimension and on the total
% score

% SUBDIMENSION LEVEL-----------------------------------------------------
sub_HC10 = 100*sub_HC1;
sub_HC10_new = 100*sub_HC1_new;
diff_HC1_subd = sub_HC10_new - sub_HC10;
sub_HC20 = 100*sub_HC2;
sub_HC20_new = 100*sub_HC2_new;
diff_HC2_subd = sub_HC20_new - sub_HC20;

% Ranking new HC subdimensions index, and original HC index
ranking_HC1_new = get_rank(sub_HC10_new)';
ranking_HC2_new = get_rank(sub_HC20_new)';

ranking_HC1 = get_rank(sub_HC10)';
ranking_HC2 = get_rank(sub_HC20)';

% Difference in rankings for HC1 and HC2 subdimension
diff_rank_HC1 = ranking_HC1 - ranking_HC1_new;
diff_rank_HC2 = ranking_HC2 - ranking_HC2_new;

% Final results
HC1_results_subd = [country_names array2table(sub_HC10) array2table(sub_HC10_new) array2table(diff_HC1_subd) array2table(ranking_HC1) array2table(ranking_HC1_new) array2table(diff_rank_HC1)];
HC2_results_subd = [country_names array2table(sub_HC20) array2table(sub_HC20_new) array2table(diff_HC2_subd) array2table(ranking_HC2) array2table(ranking_HC2_new) array2table(diff_rank_HC2)];

% Plotting the results
% HC1 SUBDIMESION: ----------------------------------------------
% 1. Plot the differences for human capital subdimension 1
%    x axis: score 2022
%    y axis: new scores

%    ranking HC score differences
figure
labels = table2array(country_names);
scatter(sub_HC10, sub_HC10_new, 25,'b','o')
axis([0 90 0 90])
hline = refline([1 0]);
set(hline,'Color','r')
text(sub_HC10, sub_HC10_new, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Human Capital: Basic skills original scores from DESI 2022 vs New scores')

% 2. Plotting the rankings differences for subdimension 1 only
% ranking HC rankings differences
[value_HC_rank,p_diff_HC_rank] = sort(diff_rank_HC1,'asc');
ranking_countries_HC_rank = country_names(p_diff_HC_rank,:);

labels = table2array(ranking_countries_HC_rank(:,"Country"));
figure
barh(value_HC_rank)
axis([-10 10 0 28])
dx = (value_HC_rank)>=0;
text(value_HC_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_HC_rank) < 0;
text(value_HC_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Human Capital: Basic skills subdimension ranking differences')


% HC2 SUBDIMESION: ----------------------------------------------
% 1. Plot the differences for human capital subdimension 1
%    x axis: score 2022
%    y axis: new scores

%    ranking HC score differences
figure
labels = table2array(country_names);
scatter(sub_HC20, sub_HC20_new, 25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(sub_HC20, sub_HC20_new, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Human Capital: Advanced skills original scores from DESI 2022 vs New scores')

% 2. Plotting the rankings differences for subdimension 1 only
% ranking HC rankings differences


% DIMENSION LEVEL-----------------------------
% Score differences for HC
diff_HC_dim = HC_Index_new - HC_Index;

% ranking new HC index, and original HC index
ranking_HC_new = get_rank(HC_Index_new)';
ranking_HC = get_rank(HC_Index)';

% Difference in rankings for HC
diff_rank_HC = ranking_HC - ranking_HC_new;

% Final results:
HC_results_dim = [country_names array2table(HC_Index) array2table(HC_Index_new) array2table(diff_HC_dim) array2table(ranking_HC) array2table(ranking_HC_new) array2table(diff_rank_HC)];

% Plotting the dimension level results
% 1. Plot of differences for human capital only,
%    compating the reference indicator for HC with the new one

figure
labels = table2array(country_names);
scatter(HC_Index, HC_Index_new,25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(HC_Index, HC_Index_new, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Human Capital: Human capital original scores from DESI 2022 vs New scores')

% Correlation between new old scores
corr(HC_Index, HC_Index_new)

% 2. Plotting the rankings differences for HC only
% ranking HC rankings differences
[value_HC_rank,p_diff_HC_rank] = sort(diff_rank_HC,'asc');
ranking_countries_HC_rank = country_names(p_diff_HC_rank,:);

labels = table2array(ranking_countries_HC_rank(:,"Country"));
figure
barh(value_HC_rank)
axis([-8 8 0 28])
dx = (value_HC_rank)>=0;
text(value_HC_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_HC_rank) < 0;
text(value_HC_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Human Capital- Ranking differences for Human Capital dimension')


% DESI LEVEL-----------------------------------------------------
% Score differences for DESI
diff_HC_overall =  DESI_new_HC - DESI;

% ranking new DESI index with HC changed, and DESI original
ranking_HC_overall = get_rank(DESI_new_HC)';
ranking_DESI = get_rank(DESI)';

% Difference in rankings for DESI_HC
diff_rank_DESI_HC = ranking_DESI-ranking_HC_overall;

% Final results
HC_results_overall = [country_names array2table(DESI) array2table(DESI_new_HC) array2table(diff_HC_overall) array2table(ranking_DESI) array2table(ranking_HC_overall) array2table(diff_rank_DESI_HC)];


% Plotting the results at DESI LEVEL---------------------------------------------------
% 1. Plot of score differences between DESI and DESI with human capital new subdimensions, 
%    keeping the other dimensions fixed

figure
labels = table2array(country_names)
scatter(DESI, DESI_new_HC,25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(DESI, DESI_new_HC, labels,'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Human Capital: Original scores of DESI 2022 vs New scores')

% 2. Plotting the rankings differences for DESI and DESI with HC changed
%    ranking DESI rankings differences (plot purpose)
[value_DESI_rank,p_diff_DESI_rank] = sort(diff_rank_DESI_HC,'asc');
ranking_countries_DESI_rank = country_names(p_diff_DESI_rank,:);

labels = table2array(ranking_countries_DESI_rank(:,"Country"));
figure
barh(value_DESI_rank)
axis([-5 5 0 28])
dx = (value_DESI_rank)>=0;
text(value_DESI_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_DESI_rank) < 0;
text(value_DESI_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Human Capital- Ranking differences for DESI')

                                             


