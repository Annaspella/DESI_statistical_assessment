% Cleaning and standardization of the data of CONNECTIVITY dimension

% script that read data from CSV files and:
% 1. compute some descriptive statistics separately for each indicator and
%    set-up a text file with results
% 2. draw histograms for each indicator
% 3. trasform indicators
% 4. Perform PCA over subsets of indicators

clear all
clc
opengl software %this command solves the problem with the video card that can cause figures to be black (all black)

% Read dataset with variable names in first row
% year in 1st column (equal to year 6 in all rows)
% country label in 2nd column
% indicator names from column 3rd on
dataset=readtable('DESI_Y6.csv', ReadVariableNames=true, VariableNamingRule='preserve'); % file DESI_Y6, contains data with already imputed values
var_names=dataset.Properties.VariableNames;

dim=size(dataset);
n_vars=dim(2)-2;       % Minus two columns corresponding to year and country
n_countries=dim(1)-1;  % Not considering EU

% Save EU line and delete from table data
EU = dataset(dim(1),:);

% data: table containing the data with imputed values for DESI year 6 (2021)
data=table2array(dataset(1:n_countries,3:end)); % excluding EU from the analysis 
                                                % (EU last row, taking the first 27)

% Computation of percentage of missing data:
perc_miss=ones(1,n_vars)*(-1);
                                                
% Reading data with missing value to calculate the percentage 
% of missing values
% data2: table containing the raw data for DESI year 6 (2021)
dataset2=readtable('DESI_missing.csv', ReadVariableNames=true, VariableNamingRule='preserve'); % file DESI_missing, contains raw data with missing values
data2=table2array(dataset2(1:n_countries,3:end)); % excluding EU from the analysis

for i=1:n_vars
    temp=find(isnan(data2(:,i)));
    temp1=size(temp);
    number_nan=temp1(1);
    perc_miss(1,i)=number_nan/n_countries*100;
end

% Standardizing data without missing values: 
% Z scores normalization for PCA 
data_std_Z=zscore(data);

%The function grpstats computes different descriptive stats on sets of indicators by groups (for instance by group of countries)
%The instruction [], tells Matlab not to cluster data into any group
%grpstats treats NaNs and NaTs as missing values, and removes them

[average,sigma,skew]=grpstats(data,[],{'mean','std','skewness'});
CV=sigma./average;                          
[max_data,ID_max]=max(data);
country_max=dataset(ID_max,2);
country_max=table2cell(country_max);
country_max=country_max'; 
[min_data,ID_min]=min(data);
country_min=dataset(ID_min,2);
country_min=table2cell(country_min);
country_min=country_min';

% Load original weights for each indicator
dataset3 = readtable('min_max.xlsx');
weights = dataset3(:,"Weight"); 

% Load original weights for each sub-dimension
dataset4 = readtable('Weigths_Subdim.xlsx');
weights_subdim = table2array(dataset4(:,"Weights")); 

% save matrix of standardized indicators: Zscore
save DESI_std.txt data_std_Z -ascii

%% Connectivity:
fid = fopen('output_connectivity.txt','w'); %open text file to put results
for i=8:17   % 8:17, indexing the 10 indicators of Connectivity
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

        % Plotting histograms
         figure
            hist(data(:,i),50)
            string=strcat({'Indicator: '},var_names{1,i+2},{' '},{'skewness='},num2str(skew(1,i)));
            title(string)
            
end %end of cycle on variables

fclose(fid);

%----------------------------------------------
% Status quo
%----------------------------------------------
%  Classical PCA over all indicators in Connectivity
temp_C= data_std_Z(:,8:17);      % temp_C: contains the standardized data for Connectivity
i_rows=[];
[i_rows,tt]=find(isnan(temp_C)); %find rows WITH missing values
temp_C(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_C,p_C]=corrcoef(temp_C);
[loadings_C,scores_C,eigenvalues_C]=pca(temp_C);


% Cronbach alpha:
c_C1 = cronbach_alpha(temp_C);
eigen_C1 = 100*(eigenvalues_C./sum(eigenvalues_C));

% Plot of PC variances
figure
plot(eigenvalues_C,'bo-')
axis([1 10 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 6 7 8 9 10 11])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C./sum(eigenvalues_C),'bo-')
axis([1 10 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 6 7 8 9 10])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')

%------------ Robust analysis-----------------------------
% Make use of functions contained in folder LIBRA_20160628 
% ROBPCA: to avoid the curse of dimensionality, n=28 must 
%         be greater than 5k thus we set k_max equal to 4
k=10;
results_C = robpca(temp_C,'k',10, 'alpha',0.75, 'plots',1);
results_C.L./sum(results_C.L)

% Loadings: k fixed to 2
results_C = robpca(temp_C,'k',2, 'alpha',0.75, 'plots',1);
results_C.P
writematrix(results_C.P, 'loadings_robC.xls');

% eigenvalues 
results_C.L
results_C.L./sum(results_C.L)
writematrix(results_C.L, 'eigen_robC.xls');
writematrix(results_C.L./sum(results_C.L), 'explained_var_robC.xls');

%grafici per ROBPCA: score plot 
results_C_rob = robpca(temp_C, 'k', 2, 'alpha',0.75,'classic', 0);

%grafici per CPCA: score plot CPCA
robpca(temp_C, 'k', 2, 'classic', 1);

%ellipsplot
mcd_C = mcdcov(results_C_rob.T, 'plots', 0);
dist = repelem(1,27);
ellipsplot(mcd_C.center, mcd_C.cov, results_C_rob.T, dist=dist)

mcd_C_class = mcdcov(scores_C(:,1:2), 'plots', 0);
ellipsplot(mcd_C_class.center, mcd_C_class.cov, scores_C(:,1:2), dist=dist)


%----------------------------------------------
% Replicating the Classical PCA over 2a1, 2a2, 2a3 
% separated from 2b1, 2b2, 2b3 and 2c1, 2c2, 2c3 and 2d1
%------------------------
% group 1:  2a1, 2a2, 2a3 
temp_CC=temp_C(:,1:3);
i_rows=[];
[i_rows,tt]=find(isnan(temp_CC)); %find rows WITH missing values
temp_CC(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_C1,p_C1]=corrcoef(temp_CC);
[loadings_C1,scores_C1,eigenvalues_C1]=pca(temp_CC);

% Cronbach alpha:
c_C2 = cronbach_alpha(temp_CC);
eigen_C2 = 100*(eigenvalues_C1./sum(eigenvalues_C1));

% Plot of PC variances
figure
plot(eigenvalues_C1,'bo-')
axis([1 3 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C1./sum(eigenvalues_C1),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')

%-----------------------
% group 2: 2b1, 2b2, 2b3
temp_CC=temp_C(:,4:6);
i_rows=[];
[i_rows,tt]=find(isnan(temp_CC)); %find rows WITH missing values
temp_CC(i_rows,:)=[];             %keep only rows without missing values, 
                           
[r_C2,p_C2]=corrcoef(temp_CC);
[loadings_C2,scores_C2,eigenvalues_C2]=pca(temp_CC);

% Cronbach alpha:
c_C3 = cronbach_alpha(temp_CC);
eigen_C3 = 100*(eigenvalues_C2./sum(eigenvalues_C2));


% Plot of PC variances
figure
plot(eigenvalues_C2,'bo-')
axis([1 3 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C2./sum(eigenvalues_C2),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')


%-------------------------------------------------------
% Replicating the Classical PCA over 2c1, 2c2, 2c3:
temp_CC=temp_C(:,7:9);
i_rows=[];
[i_rows,tt]=find(isnan(temp_CC)); %find rows WITH missing values
temp_CC(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_C3,p_C3]=corrcoef(temp_CC);
[loadings_C3,scores_C3,eigenvalues_C3]=pca(temp_CC);

% Cronbach alpha:
c_C4 = cronbach_alpha(temp_CC);
eigen_C4 = 100*(eigenvalues_C3./sum(eigenvalues_C3));


% Plot of PC variances
figure
plot(eigenvalues_C3,'bo-')
axis([1 3 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C3./sum(eigenvalues_C3),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')


%-----------------------------------------------
% Adjustment - MINOR
%-----------------------------------------------
%  PCA_new: Replicating the Classical analysis without 2a1 2a3 2b1 2d1

temp_C_new=temp_C(:,[2 5:9]); % temp_C contains standardized indicators: 
                              % taking only 2a2, 2b2, 2b3, 2c1, 2c2 and 2c3
                              % for MAJOR adjustment use command
                              % temp_C_new=temp_C(:,[2 5:8]);
i_rows=[];
[i_rows,tt]=find(isnan(temp_C_new)); %find rows WITH missing values
temp_C_new(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_C4,p_C4]=corrcoef(temp_C_new);
[loadings_C4,scores_C4,eigenvalues_C4]=pca(temp_C_new);

% Cronbach alpha:
c_C5 =cronbach_alpha(temp_C_new);
eigen_C5 = 100*(eigenvalues_C4./sum(eigenvalues_C4));
 
% Plot of PC variances
figure
plot(eigenvalues_C4,'bo-')
axis([1 6 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C4./sum(eigenvalues_C4),'bo-')
axis([1 6 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')

%-----------------------
% Including 5G STATIONS
% Upload of the new indicator
dataset_ni=readtable('new_indicators_connectivity.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_ni=dataset_ni.Properties.VariableNames;

dim=size(dataset_ni);
n_vars_ni=dim(2)-1;
n_countries_ni=dim(1);  

data_ni=table2array(dataset_ni(1:n_countries_ni,2:end));

[average_ni,sigma_ni,skew_ni]=grpstats(data_ni,[],{'mean','std','skewness'});
CV_ni=sigma_ni./average_ni;                         
[max_ni,ID_max_ni]=max(data_ni);
country_max_ni=dataset_ni(ID_max_ni,1);
country_max_ni=table2cell(country_max_ni);
country_max_ni=country_max_ni'; 
[min_ni,ID_min_ni]=min(data_ni);
country_min_ni=dataset_ni(ID_min_ni,1);
country_min_ni=table2cell(country_min_ni);
country_min_ni=country_min_ni';

fid = fopen('output_connectivity_ni.txt','w'); %open text file to put results
for i=1
% display results variable by variable
        fprintf(fid,'%s%s\n','results for: ',var_names_ni{1,i+1});
        string=strcat({'These are results for indicator:  '},var_names_ni(1,i+1));
        disp(string);
        string=strcat({'Mean value:  '},num2str(average_ni(1,i)));
        fprintf(fid,'%s\n','mean value:');
        fprintf(fid,'%12.4f \n',average_ni(1,i));
        disp(string);
        string=strcat({'Standard deviation (unbiased):  '},num2str(sigma_ni(1,i)));
        fprintf(fid,'%s\n','standard deviation:');
        fprintf(fid,'%12.4f \n',sigma_ni(1,i));
        disp(string);
        string=strcat({'Coefficient of variation:  '},num2str(CV_ni(1,i)));
        fprintf(fid,'%s\n','coefficient of variation:');
        fprintf(fid,'%12.4f \n',CV_ni(1,i));
        disp(string);
        string=strcat({'Skewness:  '},num2str(skew_ni(1,i)));
        fprintf(fid,'%s\n','Skewness:');
        fprintf(fid,'%12.4f \n',skew_ni(1,i));
        disp(string);
        string=strcat({'Maximum value:  '},num2str(max_ni(1,i)));
        fprintf(fid,'%s\n','maximum');
        fprintf(fid,'%12.4f \n',max_ni(1,i));
        disp(string);
        string=strcat({'Country corresponding to maximum value:  '},country_max_ni(1,i));
        fprintf(fid,'%s\n','Country corresponding to maximum value:');
        fprintf(fid,'%s\n',country_max_ni{1,i});
        disp(string);
        string=strcat({'Minimum value:  '},num2str(min_ni(1,i)));
        fprintf(fid,'%s\n','minimum:');
        fprintf(fid,'%12.4f \n',min_ni(1,i));
        disp(string);
        string=strcat({'Country corresponding to minimum value:  '},country_min_ni(1,i));
        fprintf(fid,'%s\n','Country corresponding to minimum value:');
        fprintf(fid,'%s\n',country_min_ni{1,i});
        disp(string);
        % Plotting histograms
         figure
            hist(data_ni(:,i),50)
            string=strcat({'Indicator: '},var_names_ni{1,i+1},{' '},{'skewness='},num2str(skew_ni(1,i)));
            title(string)
            
end %end of cycle on variables
fclose(fid);

% Be carefull: 5G stations indicator contains NaN values
% Thus, when standardizing the variables inly using the non Nan rows

% Standardizing the new indicators
% Z scores of the new indicator
temp_C_1 = [temp_C_new data_ni];
i_rows=[];
[i_rows,tt]=find(isnan(data_ni)); %find rows WITH missing values
temp_C_1(i_rows,:)=[];            %keep only rows without missing values, 
                           
data_std_Z_ni = [temp_C_1(:,1:6), zscore(temp_C_1(:,7))]; % data_std_Z_ni: 
                                                          % Columns contains standardized
                                                          % indicators 2a2, 2b2, 2b3,
                                                          % 2c1, 2c2 and 2c3 + 5G stations
                                                          % rows are relative to the countries
                                                          % without an NaN values for 5G stations
                                                          % MAJOR: data_std_Z_ni = [temp_C_1(:,1:5), zscore(temp_C_1(:,6))]; 
size(data_std_Z_ni)  % 23 x 7                                                        

%-----------------------------------------------
% 1) PCA over indicators: 2a2, 2b2, 2b3, 2c1, 2c2 and 2c3 and 5G stations (based on 23 countries only)                                                        
[r_C5,p_C5]=corrcoef(data_std_Z_ni);
[loadings_C5,scores_C5,eigenvalues_C5]=pca(data_std_Z_ni);

c_C6 =cronbach_alpha(data_std_Z_ni);
eigen_C6 = 100*(eigenvalues_C5./sum(eigenvalues_C5));

writematrix(loadings_C5, "load.xls");
writematrix(eigenvalues_C5, "eig.xls");
writematrix(eigen_C6, "eig_perc.xls");
writematrix(r_C5, "corrc.xls");

% Plot of PC variances
figure
plot(eigenvalues_C5,'bo-')
axis([1 6 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 ])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C5./sum(eigenvalues_C5),'bo-')
axis([1 6 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')


%-----------------------------------------------
% PCA  grouping 2a2 2b2 2b3 and 2c1 2c2 2c3 5G stations separately
%---------------------
% group 1: 2a2 2b2 2b3                                                  
temp_CC=temp_C_new(:,1:3);  % Using data from temp_C_new which don't contaianes the values for 5G
                            % stations (PCA based on 27 countries)
i_rows=[];
[i_rows,tt]=find(isnan(temp_CC)); %find rows WITH missing values
temp_CC(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_C6,p_C6]=corrcoef(temp_CC);
[loadings_C6,scores_C6,eigenvalues_C6]=pca(temp_CC);  %computation of PCA on the standardized indicators -C


% Cronbach alpha:
c_C7 =cronbach_alpha(temp_CC);
eigen_C7 = 100*(eigenvalues_C6./sum(eigenvalues_C6));


% Plot of PC variances
figure
plot(eigenvalues_C6,'bo-')
axis([1 3 0 3])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C6./sum(eigenvalues_C6),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')


%-------------------------------
% group 2: 2c1, 2c2, 2c3 and 5G stations                                                
temp_CC=data_std_Z_ni(:,4:7); % Using table data_std_Z_ni 
                              % (PCA based on 23 rows only)
% For major adjustment, delete 2c3 
% temp_CC=data_std_Z_ni(:,4 5 7);

[r_C7,p_C7]=corrcoef(temp_CC);
[loadings_C7,scores_C7,eigenvalues_C7]=pca(temp_CC);

% Cronbach alpha:
c_C8 =cronbach_alpha(temp_CC);
eigen_C8 = 100*(eigenvalues_C7./sum(eigenvalues_C7));

% Plot of PC variances
figure
plot(eigenvalues_C7,'bo-')
axis([1 4 0 3])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C7./sum(eigenvalues_C7),'bo-')
axis([1 4 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')

%-----------------------------------------------
% 2) Adjustment - NEW DATA for 2a2, 2a3 (reincluded) and try to reinclude 2d1 (new data) 
%-----------------------------------------------

% PCA_new: Replicating the Classical analysis without 2a1 2b1

temp_C_new=temp_C(:,[2 3 5:10]);  % temp_C contains standardized indicators: 
                                   % taking only 2a2, 2a3, 2b2, 2b3, 2c1, 2c2,
                                   % 2c3 and 2d1
                                   % Need of substituting data for 2a2 and 2a3
                                   % with new versions
i_rows=[];
[i_rows,tt]=find(isnan(temp_C_new)); %find rows WITH missing values
temp_C_new(i_rows,:)=[];             %keep only rows without missing values, 



%% ADDING NEW DATA: 
% 1) substituting 2a2,2a3, 2d1 data: 
%    new_indicators_connectivity_replace.csv contains raw data
%    1st column: country indicator, 
%    5th column: the 2a2 new data, 
%    6th column: the 2a3 new data,
%    7th column: the 2d1 new data
%   The number of countries is in this case 27 (EU not included in the table)

dataset_nd=readtable('new_indicators_connectivity_replace.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_nd=dataset_nd.Properties.VariableNames([1 5:7]); 

dim=size(dataset_nd);
n_countries_nd=dim(1);  

data_nd=table2array(dataset_nd(1:n_countries_nd,5:7));

[average_nd,sigma_nd,skew_nd]=grpstats(data_nd,[],{'mean','std','skewness'});
CV_nd=sigma_nd./average_nd;                         
[max_nd,ID_max_nd]=max(data_nd);
country_max_nd=dataset_nd(ID_max_nd,1);
country_max_nd=table2cell(country_max_nd);
country_max_nd=country_max_nd'; 
[min_nd,ID_min_nd]=min(data_nd);
country_min_nd=dataset_nd(ID_min_nd,1);
country_min_nd=table2cell(country_min_nd);
country_min_nd=country_min_nd';

fid = fopen('output_connectivity_nd.txt','w'); %open text file to put results
for i=1:3
% display results variable by variable
        fprintf(fid,'%s%s\n','results for: ',var_names_nd{1,i+1});
        string=strcat({'These are results for indicator:  '},var_names_nd(1,i+1));
        disp(string);
        string=strcat({'Mean value:  '},num2str(average_nd(1,i)));
        fprintf(fid,'%s\n','mean value:');
        fprintf(fid,'%12.4f \n',average_nd(1,i));
        disp(string);
        string=strcat({'Standard deviation (unbiased):  '},num2str(sigma_nd(1,i)));
        fprintf(fid,'%s\n','standard deviation:');
        fprintf(fid,'%12.4f \n',sigma_nd(1,i));
        disp(string);
        string=strcat({'Coefficient of variation:  '},num2str(CV_nd(1,i)));
        fprintf(fid,'%s\n','coefficient of variation:');
        fprintf(fid,'%12.4f \n',CV_nd(1,i));
        disp(string);
        string=strcat({'Skewness:  '},num2str(skew_nd(1,i)));
        fprintf(fid,'%s\n','Skewness:');
        fprintf(fid,'%12.4f \n',skew_nd(1,i));
        disp(string);
        string=strcat({'Maximum value:  '},num2str(max_nd(1,i)));
        fprintf(fid,'%s\n','maximum');
        fprintf(fid,'%12.4f \n',max_nd(1,i));
        disp(string);
        string=strcat({'Country corresponding to maximum value:  '},country_max_nd(1,i));
        fprintf(fid,'%s\n','Country corresponding to maximum value:');
        fprintf(fid,'%s\n',country_max_nd{1,i});
        disp(string);
        string=strcat({'Minimum value:  '},num2str(min_nd(1,i)));
        fprintf(fid,'%s\n','minimum:');
        fprintf(fid,'%12.4f \n',min_nd(1,i));
        disp(string);
        string=strcat({'Country corresponding to minimum value:  '},country_min_nd(1,i));
        fprintf(fid,'%s\n','Country corresponding to minimum value:');
        fprintf(fid,'%s\n',country_min_nd{1,i});
        disp(string);
        % Plotting histograms
         figure
            hist(data_nd(:,i),50)
            string=strcat({'Indicator: '},var_names_nd{1,i+1},{' '},{'skewness='},num2str(skew_nd(1,i)));
            title(string)
            
end %end of cycle on variables
fclose(fid);

% Substituting with the new data indicators for 2a2, 2a3 and 2d1
% Standardizing with zscore:
temp_C_new(:, [1 2 8]) = zscore(data_nd);

% 2) Add 5G stations:
%    new_indicators_connectivity.csv contains raw data
%    1st column: country indicator, 
%    2nd column: the 5G stations indicator new data, 
%   The number of countries is in this case 27 (EU not included in the table)
dataset_ni=readtable('new_indicators_connectivity.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_ni=dataset_ni.Properties.VariableNames;

dim=size(dataset_ni);
n_vars_ni=dim(2)-1;
n_countries_ni=dim(1);  

data_ni=table2array(dataset_ni(1:n_countries_ni,2:end));

[average_ni,sigma_ni,skew_ni]=grpstats(data_ni,[],{'mean','std','skewness'});
CV_ni=sigma_ni./average_ni;                         
[max_ni,ID_max_ni]=max(data_ni);
country_max_ni=dataset_ni(ID_max_ni,1);
country_max_ni=table2cell(country_max_ni);
country_max_ni=country_max_ni'; 
[min_ni,ID_min_ni]=min(data_ni);
country_min_ni=dataset_ni(ID_min_ni,1);
country_min_ni=table2cell(country_min_ni);
country_min_ni=country_min_ni';

fid = fopen('output_connectivity_ni.txt','w'); %open text file to put results
for i=1
% display results variable by variable
        fprintf(fid,'%s%s\n','results for: ',var_names_ni{1,i+1});
        string=strcat({'These are results for indicator:  '},var_names_ni(1,i+1));
        disp(string);
        string=strcat({'Mean value:  '},num2str(average_ni(1,i)));
        fprintf(fid,'%s\n','mean value:');
        fprintf(fid,'%12.4f \n',average_ni(1,i));
        disp(string);
        string=strcat({'Standard deviation (unbiased):  '},num2str(sigma_ni(1,i)));
        fprintf(fid,'%s\n','standard deviation:');
        fprintf(fid,'%12.4f \n',sigma_ni(1,i));
        disp(string);
        string=strcat({'Coefficient of variation:  '},num2str(CV_ni(1,i)));
        fprintf(fid,'%s\n','coefficient of variation:');
        fprintf(fid,'%12.4f \n',CV_ni(1,i));
        disp(string);
        string=strcat({'Skewness:  '},num2str(skew_ni(1,i)));
        fprintf(fid,'%s\n','Skewness:');
        fprintf(fid,'%12.4f \n',skew_ni(1,i));
        disp(string);
        string=strcat({'Maximum value:  '},num2str(max_ni(1,i)));
        fprintf(fid,'%s\n','maximum');
        fprintf(fid,'%12.4f \n',max_ni(1,i));
        disp(string);
        string=strcat({'Country corresponding to maximum value:  '},country_max_ni(1,i));
        fprintf(fid,'%s\n','Country corresponding to maximum value:');
        fprintf(fid,'%s\n',country_max_ni{1,i});
        disp(string);
        string=strcat({'Minimum value:  '},num2str(min_ni(1,i)));
        fprintf(fid,'%s\n','minimum:');
        fprintf(fid,'%12.4f \n',min_ni(1,i));
        disp(string);
        string=strcat({'Country corresponding to minimum value:  '},country_min_ni(1,i));
        fprintf(fid,'%s\n','Country corresponding to minimum value:');
        fprintf(fid,'%s\n',country_min_ni{1,i});
        disp(string);
        % Plotting histograms
         figure
            hist(data_ni(:,i),50)
            string=strcat({'Indicator: '},var_names_ni{1,i+1},{' '},{'skewness='},num2str(skew_ni(1,i)));
            title(string)
            
end %end of cycle on variables
fclose(fid);

% Be carefull: 5G stations indicator contains NaN values
% Thus, when standardizing the variables inly using the non Nan rows

% Standardizing the new indicators
% Z scores of the new indicator
temp_C_1 = [temp_C_new data_ni];
i_rows=[];
[i_rows,tt]=find(isnan(data_ni)); %find rows WITH missing values
temp_C_1(i_rows,:)=[];            %keep only rows without missing values, 
                           
data_std_Z_ni = [temp_C_1(:,[1 3:8]) zscore(temp_C_1(:,9))]; % temp_C_new: 
                                                          % columns contains indicators 
                                                          % 1] 2a2 new, (if we want to add 2a3 new add index 2) 
                                                          % 3:8] 2b2, 2b3, 2c1, 2c2, 2c3, 2d1 
                                                          % zscore(temp_C_1(:,9))] 5G stations
                                                          % rows are relative to the countries
                                                          % without an NaN values for 5G stations
                                                          
size(data_std_Z_ni)  % 23 x 7        

%% ----------------------------------------------
% PCA over indicators: 2a2 new, 2b2, 2b3, 2c1, 2c2, 2c3, 2d1 and 5G stations (based on 23 countries only)                                                        
[r_C8,p_C8]=corrcoef(data_std_Z_ni);
[loadings_C8,scores_C8,eigenvalues_C8]=pca(data_std_Z_ni);

c_C9 =cronbach_alpha(data_std_Z_ni);
eigen_C9 = 100*(eigenvalues_C8./sum(eigenvalues_C8));

writematrix(loadings_C8,'load_CC.xls');


% Plot of PC variances
figure
plot(eigenvalues_C8,'bo-')
axis([1 6 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 ])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C8./sum(eigenvalues_C8),'bo-')
axis([1 6 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')


%-----------------------------------------------
% PCA  grouping 2a2 new 2b2 2b3  and 2c1 2c2 2c3 5G stations separately
%---------------------
% group 1: 2a2 new, 2b2, 2b3                                                  
temp_CC = temp_C_new(:,[1 3 4]);  % Using data from temp_C_new which don't contain the values for 5G
                                  % stations (PCA based on 27 countries)
i_rows=[];
[i_rows,tt]=find(isnan(temp_CC)); %find rows WITH missing values
temp_CC(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_C9,p_C9]=corrcoef(temp_CC);
[loadings_C9,scores_C9,eigenvalues_C9]=pca(temp_CC);  %computation of PCA on the standardized indicators -C

writematrix(loadings_C9, 'loadings_C9.xls');
writematrix(eigenvalues_C9, 'eigenvalues_C9.xls');

% Cronbach alpha:
c_C10 = cronbach_alpha(temp_CC);
eigen_C10 = 100*(eigenvalues_C9./sum(eigenvalues_C9));


% Plot of PC variances
figure
plot(eigenvalues_C9,'bo-')
axis([1 3 0 3])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C9./sum(eigenvalues_C9),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')

% 3) Adding 2d1---------------
% group 1: 2a2 new, 2b2 2b3 2d1                                                
temp_CC = temp_C_new(:,[1 3 4 8]);  % Using data from temp_C_new which don't contaianes the values for 5G
                                    % stations (PCa based on 27 countries)
i_rows=[];
[i_rows,tt]=find(isnan(temp_CC)); %find rows WITH missing values
temp_CC(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_C9,p_C9]=corrcoef(temp_CC);
[loadings_C9,scores_C9,eigenvalues_C9]=pca(temp_CC);  %computation of PCA on the standardized indicators -C

writematrix(loadings_C9, 'loadings_C9.xls');
writematrix(eigenvalues_C9, 'eigenvalues_C9.xls');

% Cronbach alpha:
c_C10 = cronbach_alpha(temp_CC);
eigen_C10 = 100*(eigenvalues_C9./sum(eigenvalues_C9));


% Plot of PC variances
figure
plot(eigenvalues_C9,'bo-')
axis([1 4 0 4])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C9./sum(eigenvalues_C9),'bo-')
axis([1 4 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')

%-------------------------------
% group 2: 2c1, 2c2, 2c3 and 5G stations                                                
temp_CC=data_std_Z_ni(:,[4:6 8]); % Using table data_std_Z_ni:
                                  % taking indicators 2c1,2c2,2c3 and 5G
                                  % stations
                                  % (PCA based on 23 rows only)
% For major adjustment, delete 2c3 
% temp_CC=data_std_Z_ni(:,4 5 8);

[r_C10,p_C10]=corrcoef(temp_CC);
[loadings_C10,scores_C10,eigenvalues_C10]=pca(temp_CC);

writematrix(loadings_C10, 'loadings_C10.xls');
writematrix(eigenvalues_C10, 'eigenvalues_C10.xls');

% Cronbach alpha:
c_C11 =cronbach_alpha(temp_CC);
eigen_C11 = 100*(eigenvalues_C10./sum(eigenvalues_C10));

% Plot of PC variances
figure
plot(eigenvalues_C10,'bo-')
axis([1 4 0 3])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C10./sum(eigenvalues_C10),'bo-')
axis([1 4 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')

% Adding 2d1---------------
% group 2: 2c1, 2c2, 2c3, 2d1 and 5G stations                                          
temp_CC=  [data_std_Z_ni(:,[4:6 8]) data_std_Z_ni(:,7)] ; % Using table data_std_Z_ni 
                                                          % (PCA based on 23 rows only)
% For major adjustment, delete 2c3 
% temp_CC=data_std_Z_ni(:,4 5 7);

[r_C10,p_C10]=corrcoef(temp_CC);
[loadings_C10,scores_C10,eigenvalues_C10]=pca(temp_CC);

writematrix(loadings_C10, 'loadings_C10.xls');
writematrix(eigenvalues_C10, 'eigenvalues_C10.xls');

% Cronbach alpha:
c_C11 =cronbach_alpha(temp_CC);
eigen_C11 = 100*(eigenvalues_C10./sum(eigenvalues_C10));

% Plot of PC variances
figure
plot(eigenvalues_C10,'bo-')
axis([1 4 0 3])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 1 2 3 4 5])
title('Connectivity-Variances')

% Plot of explained variances
figure
plot(eigenvalues_C10./sum(eigenvalues_C10),'bo-')
axis([1 4 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Connectivity-Explained variances')



