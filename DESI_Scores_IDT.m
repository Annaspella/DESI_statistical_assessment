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
var_names=dataset.Properties.VariableNames;

dim=size(dataset);
n_vars=dim(2) -2;
n_countries=dim(1)-1;  % Excluding EU

country_names=dataset(1:n_countries,2);

%computation of percentage of missing data present 
perc_miss=ones(1,n_vars)*(-1);

%save EU line and delete:
EU = dataset(dim(1),:);

% Change the sign to the new indicator: increasing it's values, the
% subdimension values descrese. Make it consistenz with the other
% indicators in HC

% Scelte per scelta min e max: prendi come min l'utopian 0 e come max
% controlla lo storico degli ultimi 6 anni e confornta con il max del 2017:
% se max 2017 è più alto allora tieni quello, altimenti prendi il max
% storico. In ogni caso ci si aspetta che il max diminuisca nel tempo

% max set equal to 35 = level of the index for Bulgaria in 2015
% ma level in 2021 is equal to 20
% see Never_internet_15_21 table to see the historical values

data=table2array(dataset(1:n_countries,3:end)); % excluding EU from the calculations

for i=1:n_vars
    temp=find(isnan(data(:,i)));
    temp1=size(temp);
    number_nan=temp1(1);
    perc_miss(1,i)=number_nan/n_countries*100;
end


%The function grpstats computes different descriptive stats on sets of indicators by groups (for instance by group of countries)
%The instruction [], tells Matlab not to cluster data into any group
%grpstats treats NaNs and NaTs as missing values, and removes them

[average,sigma,skew]=grpstats(data,[],{'mean','std','skewness'});
CV=sigma./average;                           % ./ element-wise right division
[max_data,ID_max]=max(data);
country_max=dataset(ID_max,2);
country_max=table2cell(country_max);
country_max=country_max'; 
[min_data,ID_min]=min(data);
country_min=dataset(ID_min,2);
country_min=table2cell(country_min);
country_min=country_min';

% Load min-max value for each indicator, to min max standardize the data
dataset2 = readtable('min_max.xlsx');
min_val = dataset2(:,"Min");
max_val = dataset2(:,"Max");
weights = dataset2(:,"Weight");  % weights relative to the indicators linear combination

% Load weights of each subdimension
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

        data_std_minmax(:,i) = (data(:,i)-min_vec)./(max_vec-min_vec);

end %end of cycle on variables
fclose(fid);

%save matrix of standardized indicators: minmax 
save DESI_std_minmax_HC.txt data_std_minmax -ascii



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

% Measuring internal consistency in HC:
InternalC_across_sq = cronbach_alpha(sub_I);   %0.8589

I_Index = 100*[sub_I1 sub_I2 sub_I1]*weights_subdim(7:9); 
 

% Digital public services--------------------------------
weights_D = t(29:33,:)/sum(t(29:33,:));
sub_D = data_std_minmax(:,29:33)*weights_D;

D_Index = 100*weights_subdim(10)*sub_D;

% DESI final:

DESI = 0.25*(HC_Index+C_Index+I_Index+D_Index);


%% Score calculations based on the new subdimension: Original weigthing, only one group
% I) Integration of Digital Technology:------------------------------------
temp_I=data_std_minmax(:,18:28);
temp_I_new=temp_I(:,[1:6 8:11]); % Excluding 3b6

g_I = 1:10;


weights_I_new = t(g_I,:)/sum(t(g_I,:));
sub_I_new = temp_I_new(:,g_I)*weights_I_new;

% Meausring dimension index
I_Index_new = 100*sub_I_new; 

% Measuring overall index
DESI_new_I = 0.25*(HC_Index+C_Index+I_Index_new+D_Index);  % total DESI scores, changing scores for I,
                                                           % keeping the other dimensions fixed

%% Compare: Index for Integration od Digital Technology at 2 levels: Original weighting
%           - DIMENSION: Integration od Digital Technology
%           - OVERALL: Changing Integration od Digital Technology dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for Integration od Digital Technology: both on the 
% subdimension and on the total score

% DIMENSION LEVEL-----------------------------
% Score differences for IDT
diff_I_dim = I_Index_new - I_Index;

% ranking new IDT index, and original HC index
ranking_I_new = get_rank(I_Index_new)';
ranking_I = get_rank(I_Index)';

% Difference in rankings for IDT
diff_rank_I = ranking_I - ranking_I_new;

% Final results:
I_results_dim = [country_names array2table(I_Index) array2table(I_Index_new) array2table(diff_I_dim) array2table(ranking_I) array2table(ranking_I_new) array2table(diff_rank_I)];

% Plotting the dimension level results
% 1. Plot of differences for human capital only,
%    compating the reference indicator for HC with the new one

figure
labels = table2array(country_names);
scatter(I_Index, I_Index_new,25,'b','o')
axis([10 70 10 70])
hline = refline([1 0]);
set(hline,'Color','r')
text(I_Index, I_Index_new, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Integration of Digital Technology: IDT scores from DESI 2022 vs New scores')

corr(I_Index, I_Index_new)

% 2. Plotting the rankings differences for IDT only
% ranking IDT rankings differences ONLY FOR PLOTTING PURPOSE
[value_I_rank,p_diff_I_rank] = sort(diff_rank_I,'asc');
ranking_countries_I_rank = country_names(p_diff_I_rank,:);

labels = table2array(ranking_countries_I_rank(:,"Country"));
figure
barh(value_I_rank)
axis([-8 8 0 28])
dx = (value_I_rank)>=0;
text(value_I_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_I_rank) < 0;
text(value_I_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Integration of Digital Technology - Ranking differences for IDT dimension')


% DESI LEVEL-----------------------------------------------------
% Score differences for DESI
diff_I_overall =  DESI_new_I - DESI;

% ranking new DESI index with IDT changed, and DESI original
ranking_I_overall = get_rank(DESI_new_I)';
ranking_DESI = get_rank(DESI)';

% Difference in rankings for DESI_IDT
diff_rank_DESI_I = ranking_DESI-ranking_I_overall;

% Final results
I_results_overall = [country_names array2table(DESI) array2table(DESI_new_I) array2table(diff_I_overall) array2table(ranking_DESI) array2table(ranking_I_overall) array2table(diff_rank_DESI_I)];


% Plotting the results at DESI LEVEL---------------------------------------------------
% 1. Plot of score differences between DESI and DESI with Integration od Digital Technology new subdimensions, 
%    keeping the other dimensions fixed

figure
labels = table2array(country_names);
scatter(DESI, DESI_new_I,25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(DESI, DESI_new_I, labels,'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Integration of Digital Technology: Original scores of DESI 2022 vs New scores')

% 2. Plotting the rankings differences for DESI and DESI with IDT changed
%    ranking DESI rankings differences (plot purpose)
[value_DESI_rank,p_diff_DESI_rank] = sort(diff_rank_DESI_I,'asc');
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


%% Score calculations based on the new subdimension: New weigthing
% - Weightimg 3b1 and 3b4 as 2, 
% - the rest is 1

% I) Integration of Digital Technology:------------------------------------
t = [2 1 1 1 2 1 1 1 1 1]';      


temp_I=data_std_minmax(:,18:28);
temp_I_new=temp_I(:,[1:6 8:11]);  % Excluding 3b6
g_I = 1:10;


weights_I_new = t(g_I,:)/sum(t(g_I,:));
sub_I_new = temp_I_new(:,g_I)*weights_I_new;

% Meausring dimension index
I_Index_new_2 = 100*sub_I_new; 


% Measuring overall index
DESI_new_I_2 = 0.25*(HC_Index+C_Index+I_Index_new_2+D_Index);  % total DESI scores, changing scores for I,
                                                               % keeping the other dimensions fixed

%% Compare: Index for Integration od Digital Technology at 2 levels: Original weighting
%           - DIMENSION: Integration od Digital Technology
%           - OVERALL: Changing Integration od Digital Technology dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for Integration od Digital Technology: both on the 
% subdimension and on the total score

% DIMENSION LEVEL-----------------------------
% Score differences for IDT
diff_I_dim_2 = I_Index_new_2 - I_Index;

% ranking new IDT index, and original HC index
ranking_I_new_2 = get_rank(I_Index_new_2)';
ranking_I = get_rank(I_Index)';

% Difference in rankings for IDT
diff_rank_I_2 = ranking_I - ranking_I_new_2;

% Final results:
I_results_dim_2 = [country_names array2table(I_Index) array2table(I_Index_new_2) array2table(diff_I_dim_2) array2table(ranking_I) array2table(ranking_I_new_2) array2table(diff_rank_I_2)];

% Plotting the dimension level results
% 1. Plot of differences for human capital only,
%    compating the reference indicator for HC with the new one

figure
labels = table2array(country_names);
scatter(I_Index, I_Index_new_2,25,'b','o')
axis([10 70 10 70])
hline = refline([1 0]);
set(hline,'Color','r')
text(I_Index, I_Index_new_2, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Integration of Digital Technology: Human capital original scores from DESI 2022 vs New scores')

% 2. Plotting the rankings differences for IDT only
% ranking IDT rankings differences ONLY FOR PLOTTING PURPOSE
[value_I_rank_2,p_diff_I_rank_2] = sort(diff_rank_I_2,'asc');
ranking_countries_I_rank_2 = country_names(p_diff_I_rank_2,:);

labels = table2array(ranking_countries_I_rank_2(:,"Country"));
figure
barh(value_I_rank_2)
axis([-8 8 0 28])
dx = (value_I_rank_2)>=0;
text(value_I_rank_2(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_I_rank_2) < 0;
text(value_I_rank_2(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Integration of Digital Technology - Ranking differences for Human Capital dimension')


% DESI LEVEL-----------------------------------------------------
% Score differences for DESI
diff_I_overall_2 =  DESI_new_I_2 - DESI;

% ranking new DESI index with IDT changed, and DESI original
ranking_I_overall_2 = get_rank(DESI_new_I_2)';
ranking_DESI = get_rank(DESI)';

% Difference in rankings for DESI_IDT
diff_rank_DESI_I_2 = ranking_DESI-ranking_I_overall_2;

% Final results
I_results_overall_2 = [country_names array2table(DESI) array2table(DESI_new_I_2) array2table(diff_I_overall_2) array2table(ranking_DESI) array2table(ranking_I_overall_2) array2table(diff_rank_DESI_I_2)];


% Plotting the results at DESI LEVEL---------------------------------------------------
% 1. Plot of score differences between DESI and DESI with Integration od Digital Technology new subdimensions, 
%    keeping the other dimensions fixed

figure
labels = table2array(country_names);
scatter(DESI, DESI_new_I_2,25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(DESI, DESI_new_I_2, labels,'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Integration of Digital Technology: Original scores of DESI 2022 vs New scores')

% 2. Plotting the rankings differences for DESI and DESI with IDT changed
%    ranking DESI rankings differences (plot purpose)
[value_DESI_rank_2,p_diff_DESI_rank_2] = sort(diff_rank_DESI_I_2,'asc');
ranking_countries_DESI_rank_2 = country_names(p_diff_DESI_rank_2,:);

labels = table2array(ranking_countries_DESI_rank_2(:,"Country"));
figure
barh(value_DESI_rank_2)
axis([-5 5 0 28])
dx = (value_DESI_rank_2)>=0;
text(value_DESI_rank_2(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_DESI_rank_2) < 0;
text(value_DESI_rank_2(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Integration of Digital Technology - Ranking differences for DESI')


%% Score calculations based on the new subdimension: Equal weigthing
% (Weighting all the indicators as 1)
% I) Integration of Digital Technology:------------------------------------
temp_I=data_std_minmax(:,18:28);
temp_I_new=temp_I(:,[1:6 8:11]); % Excluding 3b6

g1_I = 1;
g2_I = 2:7;
g3_I = 8:10;

t = repelem(1, 10)';

weights_I1_ew = t(g1_I,:)/sum(t(g1_I,:));
sub_I1_ew = temp_I_new(:,g1_I)*weights_I1_ew;

weights_I2_ew = t(g2_I,:)/sum(t(g2_I,:));
sub_I2_ew = temp_I_new(:,g2_I)*weights_I2_ew;

weights_I3_ew = t(g3_I,:)/sum(t(g3_I,:));
sub_I3_ew = temp_I_new(:,g3_I)*weights_I3_ew;

% Measuring internal consistency in IDT:
InternalC_across_ew = cronbach_alpha([sub_I1_ew sub_I2_ew sub_I3_ew]); 

% Meausring dimension index
I_Index_ew = 100*[sub_I1_ew sub_I2_ew sub_I3_ew]*weights_subdim(7:9); 

% Measuring overall index
DESI_new_I_ew = 0.25*(HC_Index+C_Index+I_Index_ew+D_Index);  % total DESI scores, changing scores for I,
                                                               % keeping the other dimensions fixed

%% Compare: Index for Integration od Digital Technology at 2 levels: Equal weighting
%           - DIMENSION: Integration od Digital Technology
%           - OVERALL: Changing Integration od Digital Technology dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for Integration od Digital Technology: both on the 
% subdimension and on the total score

% DIMENSION LEVEL-----------------------------
% Score differences for IDT
diff_I_dim_ew = I_Index_ew - I_Index;

% ranking new IDT index, and original HC index
ranking_I_ew = get_rank(I_Index_ew)';
ranking_I = get_rank(I_Index)';

% Difference in rankings for IDT
diff_rank_I_ew = ranking_I - ranking_I_ew;

% Final results:
I_results_dim_ew = [country_names array2table(I_Index) array2table(I_Index_ew) array2table(diff_I_dim_ew) array2table(ranking_I) array2table(ranking_I_ew) array2table(diff_rank_I_ew)];

% Plotting the dimension level results
% 1. Plot of differences for human capital only,
%    compating the reference indicator for HC with the new one

figure
labels = table2array(country_names);
scatter(I_Index, I_Index_ew,25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(I_Index, I_Index_ew, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Integration of Digital Technology: Human capital original scores from DESI 2022 vs New scores')

% 2. Plotting the rankings differences for IDT only
% ranking IDT rankings differences ONLY FOR PLOTTING PURPOSE
[value_I_rank_ew,p_diff_I_rank_ew] = sort(diff_rank_I_ew,'asc');
ranking_countries_I_rank_ew = country_names(p_diff_I_rank_ew,:);

labels = table2array(ranking_countries_I_rank_ew(:,"Country"));
figure
barh(value_I_rank_ew)
axis([-8 8 0 28])
dx = (value_I_rank_ew)>=0;
text(value_I_rank_ew(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_I_rank_ew) < 0;
text(value_I_rank_ew(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Integration of Digital Technology - Ranking differences for Human Capital dimension')


% DESI LEVEL-----------------------------------------------------
% Score differences for DESI
diff_I_overall_ew =  DESI_new_I_ew - DESI;

% ranking new DESI index with IDT changed, and DESI original
ranking_I_overall_ew = get_rank(DESI_new_I_ew)';
ranking_DESI = get_rank(DESI)';

% Difference in rankings for DESI_IDT
diff_rank_DESI_I_ew = ranking_DESI-ranking_I_overall_ew;

% Final results
I_results_overall_ew = [country_names array2table(DESI) array2table(DESI_new_I_ew) array2table(diff_I_overall_ew) array2table(ranking_DESI) array2table(ranking_I_overall_ew) array2table(diff_rank_DESI_I_ew)];


% Plotting the results at DESI LEVEL---------------------------------------------------
% 1. Plot of score differences between DESI and DESI with Integration od Digital Technology new subdimensions, 
%    keeping the other dimensions fixed

figure
labels = table2array(country_names);
scatter(DESI, DESI_new_I_ew,25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(DESI, DESI_new_I_ew, labels,'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Integration of Digital Technology: Original scores of DESI 2022 vs New scores')

% 2. Plotting the rankings differences for DESI and DESI with IDT changed
%    ranking DESI rankings differences (plot purpose)
[value_DESI_rank_ew,p_diff_DESI_rank_ew] = sort(diff_rank_DESI_I_ew,'asc');
ranking_countries_DESI_rank_ew = country_names(p_diff_DESI_rank_ew,:);

labels = table2array(ranking_countries_DESI_rank_ew(:,"Country"));
figure
barh(value_DESI_rank_ew)
axis([-5 5 0 28])
dx = (value_DESI_rank_ew)>=0;
text(value_DESI_rank_ew(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_DESI_rank_ew) < 0;
text(value_DESI_rank_ew(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Human Capital- Ranking differences for DESI')


