%% Replicating the ANALYSIS FOR THE DESI Scores

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
n_vars_d=dim(2)-2;     % minus year and country columns
n_countries=dim(1)-1;  % Excluding EU
country_names=dataset(1:n_countries,2);

% save EU line and delete it from the table:
EU = dataset(dim(1),:);
data=table2array(dataset(1:n_countries,3:end)); % excluding EU from the calculations

%------------------------------
% Add new indicator: 5G stations
% Choices for min e max: min taken equal to 0 (present in our data)
% while max taken equal to 2 times the max in our data
% Standardizing data
dataset_ni=readtable('new_indicators_connectivity.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_ni=dataset_ni.Properties.VariableNames;

dim=size(dataset_ni);
n_vars_ni=dim(2)-1;
n_countries_ni=dim(1); 

data_ni=table2array(dataset_ni(1:n_countries_ni,2));   % data_ni containes the new indicators containing NaN values
                                                       % and raw scores
                          
%---------------------------
% Table containing original indicators plus 5G stations
data_tot= [data data_ni];

% Merging toghether the variable names from 
% the first and the second table
var_names = [var_names_d var_names_ni(2)];
n_vars = n_vars_d + n_vars_ni;

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
CV=sigma./average;                           % ./ element-wise right division
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
% Selecting all 33 indicators + 5G stations [1:33 36]
dataset2 = readtable('min_max.xlsx');
min_val = dataset2([1:33 36],"Min");
max_val = dataset2([1:33 36],"Max");
weights = dataset2([1:33 36],"Weight");   % original weight for each indicator

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

        data_std_minmax(:,i) = (data_tot(:,i)-min_vec)./(max_vec-min_vec);    

end %end of cycle on variables
fclose(fid);

%save matrix of standardized indicators: minmax 
save DESI_std_minmax_C.txt data_std_minmax -ascii



%% Replicated analysis:
% Human capital----------------------------------------
t = table2array(weights);
weights_HC1 = t(1:3,:)/sum(t(1:3,:));
weights_HC2 = t(4:7,:)/sum(t(4:7,:));

sub_HC1 = data_std_minmax(:,1:3)*weights_HC1;
sub_HC2 = data_std_minmax(:,4:7)*weights_HC2;

sub_HC = [sub_HC1 sub_HC2];

HC_Index = 100*sub_HC*weights_subdim(1:2);

% Connectivity-----------------------------------------
weights_C1 = t(8:10,:)/sum(t(8:10,:));
weights_C2 = t(11:13,:)/sum(t(11:13,:));
weights_C3 = t(14:16,:)/sum(t(14:16,:));
weights_C4 = t(17,:)/sum(t(17,:));

sub_C1 = data_std_minmax(:,8:10)*weights_C1;
sub_C2 = data_std_minmax(:,11:13)*weights_C2;
sub_C3 = data_std_minmax(:,14:16)*weights_C3;
sub_C4 = data_std_minmax(:,17)*weights_C4;

sub_C = [sub_C1 sub_C2 sub_C3 sub_C4];

% Measuring internal consistency in Connectivity:
InternalC_across_sq = cronbach_alpha(sub_C);   % -0.3265

C_Index = 100*sub_C*weights_subdim(3:6); 


% Integration of Digital Technology-------------------
weights_I1 = t(18,:)/sum(t(18,:));
weights_I2 = t(19:25,:)/sum(t(19:25,:));
weights_I3 = t(26:28,:)/sum(t(26:28,:));

sub_I1 = data_std_minmax(:,18)*weights_I1;
sub_I2 = data_std_minmax(:,19:25)*weights_I2;
sub_I3 = data_std_minmax(:,26:28)*weights_I3;

sub_I = [sub_I1 sub_I2 sub_I3];

I_Index = 100*[sub_I1 sub_I2 sub_I1]*weights_subdim(7:9); 
 

% Digital public services------------------------------
weights_D = t(29:33,:)/sum(t(29:33,:));
sub_D = data_std_minmax(:,29:33)*weights_D;

D_Index = 100*weights_subdim(10)*sub_D;

% DESI final:

DESI = 0.25*(HC_Index+C_Index+I_Index+D_Index);


%% 1) Score calculations based on the new subdimension: MINOR with ORIGINAL WEIGHTs
% C) Connectivity: minor adjustment-------------------------------------
temp_C = data_std_minmax(:,[8:17 34]);     % temp_C: extract indicators from connectivity + 
                                           % 5G stations indicator which column index is 34
                                           % (data_std_minmax contains the minmax values for the original indicators
                                           %  + minmax values for the newly introduced indicator 5G stations)

temp_C_new=temp_C(:,[2 5:9 11]);           % temp_C_new: retains only the a priori 
                                           % selected indicators + 5G stations

% Doing the same with the weights
tt = table2array(weights([8:17 34],1));        % taking only the weights for the
t =  tt([2 5:9 11]);                           % selected indicators + 5G stations
                                               % t =  repelem(1,7), EQUAL WEIGTHING
                                               % (weights contains the
                                               % weights of all the original indicators
                                               % in Connectivity
                                               % + the new indicator for Connecti
                                               
% Defining the two new sub-groups:
g1_H = [1:3];                            
g2_H = [4:7];
weights_C1_new = t(g1_H,:)/sum(t(g1_H,:));  
weights_C2_new = t(g2_H,:)/sum(t(g2_H,:));  
sub_C1_new = temp_C_new(:,g1_H)*weights_C1_new;
sub_C2_new = get_score(temp_C_new(:,g2_H),weights_C2_new);

% Measuring internal consistency in HC:
InternalC_across_new_min = cronbach_alpha([sub_C1_new sub_C2_new]);  % -0.9141
corr_min = corr(sub_C1_new, sub_C2_new);

C_Index_new = 50*(sub_C1_new + sub_C2_new); % 0.5 on each subdim and 100 to change scale

DESI_new_C = 0.25*(HC_Index + C_Index_new + I_Index + D_Index);  % total DESI scores, changing Connect, 
                                                                 % keeping the other dimensions fixed


%% 1.1) Compare: Index for Connectivity at 2 levels: MINOR with ORIGINAL WEIGHTS
%           - DIMENSION: Connectivity
%           - OVERALL: Changing Connectivity dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for Connectivity: both on the subdimension and on the total
% score

% DIMENSION LEVEL-----------------------------
% Score differences for Connectivity
diff_C_dim = C_Index_new - C_Index;

% ranking new HC index and original Connectivity index
ranking_C_new = get_rank(C_Index_new)';
ranking_C = get_rank(C_Index)';

% Difference in rankings for Connectivity
diff_rank_C = ranking_C - ranking_C_new;

% Final results:
C_results_dim = [country_names array2table(C_Index) array2table(C_Index_new) array2table(diff_C_dim) array2table(ranking_C) array2table(ranking_C_new) array2table(diff_rank_C)];

% Plotting the dimension level results
% 1. Plot of differences for human capital only,
%    compating the reference indicator for C with the new one
figure
labels = table2array(country_names);
scatter(C_Index, C_Index_new,25,'b','o')
axis([30 90 30 90])
hline = refline([1 0]);
set(hline,'Color','r')
text(C_Index, C_Index_new, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Connectivity: Connectivity original scores from DESI 2022 vs New scores')

% 2. Plotting the rankings differences for C only
% ranking C rankings differences
[value_C_rank,p_diff_C_rank] = sort(diff_rank_C,'asc');
ranking_countries_C_rank = country_names(p_diff_C_rank,:);

labels = table2array(ranking_countries_C_rank(:,"Country"));
figure
barh(value_C_rank)
axis([-10 10 0 28])
dx = (value_C_rank)>=0;
text(value_C_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_C_rank) < 0;
text(value_C_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Connectivity - Ranking differences for Human Capital dimension')


% DESI LEVEL-----------------------------------------------------
% Score differences for DESI
diff_C_overall =  DESI_new_C - DESI;

% ranking new DESI index with C changed, and DESI original
ranking_C_overall = get_rank(DESI_new_C)';
ranking_DESI = get_rank(DESI)';

% Difference in rankings between DESI and DESI_C
diff_rank_DESI_C = ranking_DESI-ranking_C_overall;

% Final results
C_results_overall = [country_names array2table(DESI) array2table(DESI_new_C) array2table(diff_C_overall) array2table(ranking_DESI) array2table(ranking_C_overall) array2table(diff_rank_DESI_C)];


% Plotting the results at DESI LEVEL---------------------------------------------------
% 1. Plot of score differences between DESI and DESI with connectivity new subdimensions, 
%    keeping the other dimensions fixed
figure
labels = table2array(country_names);
scatter(DESI, DESI_new_C,25,'b','o')
axis([20 80 20 80])
hline = refline([1 0]);
set(hline,'Color','r')
text(DESI, DESI_new_C, labels,'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Connectivity: Original scores of DESI 2022 vs New scores')

% 2. Plotting the rankings differences for DESI and DESI with C changed
%    ranking DESI rankings differences (only for plot purpose)
[value_DESI_rank,p_diff_DESI_rank] = sort(diff_rank_DESI_C,'asc');
ranking_countries_DESI_rank = country_names(p_diff_DESI_rank,:);

labels = table2array(ranking_countries_DESI_rank(:,"Country"));
figure
barh(value_DESI_rank)
axis([-5 5 0 28])
dx = (value_DESI_rank)>=0;
text(value_DESI_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_DESI_rank) < 0;
text(value_DESI_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Connectivity - Ranking differences for DESI')



%% ADDING NEWDATA
% Substituting 2a2,2a3 data and 2d1
dataset_nd=readtable('new_indicators_connectivity_replace.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_nd=dataset_nd.Properties.VariableNames([1 5:7]); 

dim=size(dataset_nd);
n_countries_nd=dim(1);  

data_nd=table2array(dataset_nd(1:n_countries_nd,5:7));

% Standardizing min-max and substituitng to old values:

data_std_minmax_new = data_std_minmax;
j=1;

for i= [9 10 17]
        min_vec = table2array(min_val(i,1));
        max_vec = table2array(max_val(i,1));

        data_std_minmax_new(:,i) = (data_nd(:,j)-min_vec)./(max_vec-min_vec);  
        j = j+1;
end %end of cycle on variables


%% 2) Score calculations based on the new subdimension: NEWDATA, NO PRICE INDEX
% (no price index)
% C) Connectivity: -------------------------------------------------------
temp_C_1 = data_std_minmax_new(:,[8:17 34]);     % temp_C: extract indicators from connectivity + 
                                           % 5G stations indicator which column index is 34
                                           % In total it contains 11
                                           % indicators                                   

temp_C_new_1=temp_C_1(:,[2 5:9 11]);           % temp_C_new: retains only the
                                               % selected indicators: 2a2, 2b2, 2b3
                                               % 2c1, 2c2, 2c3 + 5G stations



% Doing the same with the weights
tt = table2array(weights([8:17 34],1));        % taking only the weights for the
t =  tt([2 5:9 11]);                           % selected indicators + 5G stations

% Defining the two new sub-groups:
g1_C = [1:3];                            
g2_C = [4:7];
weights_C1_new = t(g1_C,:)/sum(t(g1_C,:));  
weights_C2_new = t(g2_C,:)/sum(t(g2_C,:));  
sub_C1_new = temp_C_new_1(:,g1_C)*weights_C1_new;
sub_C2_new = get_score(temp_C_new_1(:,g2_C),weights_C2_new);

% Measuring internal consistency in HC:
InternalC_across_new_min = cronbach_alpha([sub_C1_new sub_C2_new]);  % -0.9098
corr_min = corr(sub_C1_new, sub_C2_new);

C_Index_new = 50*(sub_C1_new + sub_C2_new); % 0.5 on each subdim and 100 to change scale

DESI_new_C = 0.25*(HC_Index + C_Index_new + I_Index + D_Index);  % total DESI scores, changing Connect, 
                                                                 % keeping the other dimensions fixed
    

%% 2.1) Compare: Index for Connectivity at 2 levels: NEWDATA, NO PRICE INDEX
%           - DIMENSION: Connectivity 
%           - OVERALL: Changing Connectivity dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for Connectivity: both on the subdimension and on 
% the total score

% DIMENSION LEVEL-----------------------------
% Score differences for Connectivity
diff_C_dim = C_Index_new - C_Index;

% ranking new HC index and original Connectivity index
ranking_C_new = get_rank(C_Index_new)';
ranking_C = get_rank(C_Index)';

% Difference in rankings for Connectivity
diff_rank_C = ranking_C - ranking_C_new;

% Final results:
C_results_dim = [country_names array2table(C_Index) array2table(C_Index_new) array2table(diff_C_dim) array2table(ranking_C) array2table(ranking_C_new) array2table(diff_rank_C)];


% Plotting the DIMENSION level results:
% 1. Plot of differences for Conn only,
%    compating the reference indicator for C with the new one
figure
labels = table2array(country_names);
scatter(C_Index, C_Index_new,25,'b','o')
axis([30 90 30 90])
hline = refline([1 0]);
set(hline,'Color','r')
text(C_Index, C_Index_new, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Connectivity: Connectivity original scores from DESI 2022 vs New scores')

% Correlation between new old scores
corr(C_Index, C_Index_new)

% 2. Plotting the rankings differences for C only
% ranking C rankings differences
[value_C_rank,p_diff_C_rank] = sort(diff_rank_C,'asc');
ranking_countries_C_rank = country_names(p_diff_C_rank,:);

labels = table2array(ranking_countries_C_rank(:,"Country"));
figure
barh(value_C_rank)
axis([-15 15 0 28])
dx = (value_C_rank)>=0;
text(value_C_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_C_rank) < 0;
text(value_C_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Connectivity - Ranking differences for Connectivity dimension')



%% 3) Score calculations based on the new subdimension: NEW DATA, INCLUDING PRICE INDEX in FIXED
% C) Connectivity: minor adjustment-------------------------------------
temp_C_2 = data_std_minmax_new(:,[8:17 34]);

temp_C_new_2 = temp_C_2(:,[2 5:11]);           % temp_C_new: retains only the a priori 
                                               % selected indicators + 5G stations

% Doing the same with the weights
tt = table2array(weights([8:17 34],1));      % taking only the weights for the
t =  tt([2 5:11]);                           % selected indicators + 5G stations

% Defining the two new sub-groups:
g1_C = [1:3 7];     % inclusion of d21 (7) in the fixed subdimension                          
g2_C = [4 5 6 8];

weights_C1_new = t(g1_C,:)/sum(t(g1_C,:));  
weights_C2_new = t(g2_C,:)/sum(t(g2_C,:));

sub_C1_new = temp_C_new_2(:,g1_C)*weights_C1_new;
sub_C2_new = get_score(temp_C_new_2(:,g2_C),weights_C2_new);


% Measuring internal consistency in HC:
InternalC_across_new_min = cronbach_alpha([sub_C1_new sub_C2_new]);  %-0.8108
corr_min = corr([sub_C1_new sub_C2_new]);
writematrix(corr_min, "corr.xls");

C_Index_new = 50*sub_C1_new + 50*sub_C2_new; 

DESI_new_C = 0.25*(HC_Index + C_Index_new + I_Index + D_Index);  % total DESI scores, changing Connectivity, 
                                                                 % keeping the other dimensions fixed
    

%% 3.1) Compare: Index for Connectivity at 2 levels: NEW DATA, INCLUDING PRICE INDEX in FIXED
%           - DIMENSION: Connectivity 
%           - OVERALL: Changing Connectivity dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for Connectivity: both on the subdimension and on the total
% score

% DIMENSION LEVEL-----------------------------
% Score differences for Connectivity
diff_C_dim = C_Index_new - C_Index;

% ranking new HC index and original Connectivity index
ranking_C_new = get_rank(C_Index_new)';
ranking_C = get_rank(C_Index)';

% Difference in rankings for Connectivity
diff_rank_C = ranking_C - ranking_C_new;

% Final results:
C_results_dim = [country_names array2table(C_Index) array2table(C_Index_new) array2table(diff_C_dim) array2table(ranking_C) array2table(ranking_C_new) array2table(diff_rank_C)];

% Plotting the dimension level results
% 1. Plot of differences for human capital only,
%    compating the reference indicator for C with the new one
figure
labels = table2array(country_names);
scatter(C_Index, C_Index_new,25,'b','o')
axis([30 90 30 90])
hline = refline([1 0]);
set(hline,'Color','r')
text(C_Index, C_Index_new, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Connectivity: Connectivity original scores from DESI 2022 vs New scores')

% Correlation between new old scores
corr(C_Index, C_Index_new)


% 2. Plotting the rankings differences for C only
% ranking C rankings differences
[value_C_rank,p_diff_C_rank] = sort(diff_rank_C,'asc');
ranking_countries_C_rank = country_names(p_diff_C_rank,:);

labels = table2array(ranking_countries_C_rank(:,"Country"));
figure
barh(value_C_rank)
axis([-15 15 0 28])
dx = (value_C_rank)>=0;
text(value_C_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_C_rank) < 0;
text(value_C_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Connectivity - Ranking differences for Connectivity dimension')


%% 4) Score calculations based on the new subdimension: NEWDATA + PRICE INDEX as a separated dimension
% C) Connectivity: minor adjustment-------------------------------------
temp_C_2 = data_std_minmax_new(:,[8:17 34]);

temp_C_new_2 = temp_C_2(:,[2 5:11]);           % temp_C_new: retains only the a priori 
                                               % selected indicators + 5G stations

% Doing the same with the weights
tt = table2array(weights([8:17 34],1));      % taking only the weights for the
t =  tt([2 5:11]);                           % selected indicators + 5G stations

% Defining the two new sub-groups:
g1_C = [1:3];                            
g2_C = [4 5 6 8];
g3_C = 7;
weights_C1_new = t(g1_C,:)/sum(t(g1_C,:));  
weights_C2_new = t(g2_C,:)/sum(t(g2_C,:));

sub_C1_new = temp_C_new_2(:,g1_C)*weights_C1_new;
sub_C2_new = get_score(temp_C_new_2(:,g2_C),weights_C2_new);
sub_C3_new = temp_C_new_2(:,g3_C);  % price index indicator


% Measuring internal consistency in HC:
InternalC_across_new_min = cronbach_alpha([sub_C1_new sub_C2_new sub_C3_new]);  % 0.0463
corr_min = corr([sub_C1_new sub_C2_new sub_C3_new]);
writematrix(corr_min, "corr.xls");

C_Index_new = 50*sub_C1_new + 40*sub_C2_new + 10*sub_C3_new; 

DESI_new_C = 0.25*(HC_Index + C_Index_new + I_Index + D_Index);  % total DESI scores, changing Connectivity, 
                                                                 % keeping the other dimensions fixed
    

%% 4.1) Compare: Index for Connectivity at 2 levels: NEWDATA + PRICE INDEX as a separated dimension
%           - DIMENSION: Connectivity 
%           - OVERALL: Changing Connectivity dimension inside DESI, the
%                      rest stays fixed
% Analysis on the results for Connectivity: both on the subdimension and on the total
% score

% DIMENSION LEVEL-----------------------------
% Score differences for Connectivity
diff_C_dim = C_Index_new - C_Index;

% ranking new HC index and original Connectivity index
ranking_C_new = get_rank(C_Index_new)';
ranking_C = get_rank(C_Index)';

% Difference in rankings for Connectivity
diff_rank_C = ranking_C - ranking_C_new;

% Final results:
C_results_dim = [country_names array2table(C_Index) array2table(C_Index_new) array2table(diff_C_dim) array2table(ranking_C) array2table(ranking_C_new) array2table(diff_rank_C)];

% Plotting the dimension level results
% 1. Plot of differences for human capital only,
%    compating the reference indicator for C with the new one
figure
labels = table2array(country_names);
scatter(C_Index, C_Index_new,25,'b','o')
axis([30 90 30 90])
hline = refline([1 0]);
set(hline,'Color','r')
text(C_Index, C_Index_new, labels, 'FontSize',8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('New scores')
xlabel('Original scores')
title('Connectivity: Connectivity original scores from DESI 2022 vs New scores')

% 2. Plotting the rankings differences for C only
% ranking C rankings differences
[value_C_rank,p_diff_C_rank] = sort(diff_rank_C,'asc');
ranking_countries_C_rank = country_names(p_diff_C_rank,:);

labels = table2array(ranking_countries_C_rank(:,"Country"));
figure
barh(value_C_rank)
axis([-15 15 0 28])
dx = (value_C_rank)>=0;
text(value_C_rank(dx),(27-sum(dx)+1):27 , labels(dx),'VerticalAlignment','HorizontalAlignment','left')
sx = (value_C_rank) < 0;
text(value_C_rank(sx),1:(27-sum(dx)) , labels(sx),'VerticalAlignment','HorizontalAlignment','right')
title('Connectivity - Ranking differences for Connectivity dimension')


