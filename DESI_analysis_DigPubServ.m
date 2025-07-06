%% Cleaning and standardization of the data
clear all
clc
opengl software %this command solves the problem with the video card that can cause figures to be black (all black)

%read dataset with variable names in first row
%country in 1st column
%region in 2nd column
%region ID in 3rd column
%indicator names from column 4th on
dataset=readtable('DESI_Y6.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names=dataset.Properties.VariableNames;

dim=size(dataset);
n_vars=dim(2)-2;
n_countries=dim(1)-1;  % Not considering EU

%computation of percentage of missing data present 
perc_miss=ones(1,n_vars)*(-1);

%save EU line and delete:
EU = dataset(dim(1),:);

data=table2array(dataset(1:n_countries,3:end)); % excluding EU from the analysis

% reading data with missing value to calculate the percentage of missing
% values
dataset2=readtable('DESI_missing.csv', ReadVariableNames=true, VariableNamingRule='preserve');
data2=table2array(dataset2(1:n_countries,3:end)); % excluding EU from the analysis

for i=1:n_vars
    temp=find(isnan(data2(:,i)));
    temp1=size(temp);
    number_nan=temp1(1);
    perc_miss(1,i)=number_nan/n_countries*100;
end

% Standardizing data without missing values
data_std_Z=zscore(data);
data_dst_minmax=data;

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



%% Digital Public Services: getting the summarizing statistics for indicators
%  in DPS
fid = fopen('output_D.txt','w'); %open text file to put results
for i=29:33
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

        % Plotting histograms
         figure
            hist(data(:,i),50)
            string=strcat({'Indicator: '},var_names{1,i+2},{' '},{'skewness='},num2str(skew(1,i)));
            title(string)
            
end %end of cycle on variables
fclose(fid);

% PCA_1:
temp_D=data_std_Z(:,29:33);
i_rows=[];
[i_rows,tt]=find(isnan(temp_D)); %find rows WITH missing values
temp_I(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_D,p_D]=corrcoef(temp_D);
[loadings_D,scores_D,eigenvalues_D]=pca(temp_D);%computation of PCA on the standardized indicators -D

writematrix(r_D,'corr_D.xls')
writematrix(loadings_D(:,1:3),'loadings_D.xls')
writematrix(eigenvalues_D,'eigen_D.xls')
writematrix(100*(eigenvalues_D./sum(eigenvalues_D)),'cumulative_D.xls')

% Cronbach alpha:
c_D1 =cronbach_alpha(temp_D);

% Plot of PC variances
figure
plot(eigenvalues_D,'bo-')
axis([1 5 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 6 7 8])
set(gca,'YTick',[0 1 2 3 4 5])
title('DigitalPublicServices- Variances')

% Plot of explained variances
figure
plot(eigenvalues_D./sum(eigenvalues_D),'bo-')
axis([1 5 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('DigitalPublicServices-Explained variances')


% ROBPCA: to avoid the curse of dimensionality, n=28 must be greater than 5k
%         thus we set k_max equal to 4number of total variables
k=5;
results_D = robpca(temp_D,'k',k, 'alpha',0.75, 'plots',1);

% Loadings:
results_D.P
writematrix(results_D.P(:,1:3), 'loadings_ROB_I.xls')

% eigenvalues 
results_D.L
results_D.L ./sum(results_D.L)

%grafici per ROBPCA: score plot 
results_D_rob = robpca(temp_D, 'k', 2, 'alpha',0.75,'classic', 0);

%grafici per CPCA: score plot CPCA
robpca(temp_D, 'k', 2, 'classic', 1);

%ellipsplot
mcd_D = mcdcov(results_D_rob.T, 'plots', 0);
dist = repelem(1,27);
ellipsplot(mcd_D.center, mcd_D.cov, results_D_rob.T, dist=dist)

mcd_D_class = mcdcov(scores_D(:,1:2), 'plots', 0);
ellipsplot(mcd_D_class.center, mcd_D_class.cov, scores_D(:,1:2), dist=dist)

%-----------------------------------------------
% PCA_2: Replicating the Classical analysis with all indicators except 4a5
temp_DD=temp_D(:,[1:4]); 
i_rows=[];
[i_rows,tt]=find(isnan(temp_DD)); %find rows WITH missing values
temp_DD(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_DD,p_DD]=corrcoef(temp_DD);
[loadings_DD,scores_DD,eigenvalues_DD]=pca(temp_DD);%computation of PCA on the standardized indicators -C

writematrix(r_DD,'corr_DD.xls');
writematrix(loadings_DD(:,1:3),'loadings_DD.xls')

% Cronbach alpha:
c_D2 =cronbach_alpha(temp_DD);

% Plot of PC variances
figure
plot(eigenvalues_DD,'bo-')
axis([1 4 0 4])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 1 2 3 4])
title('Digital Public Services-Variances')

% Plot of explained variances
figure
plot(eigenvalues_DD./sum(eigenvalues_DD),'bo-')
axis([1 4 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Public Services -Explained variances')

%-----------------------
% Including Mobile friendliness, user support,
%           Transparency (aggregate), CB online availabiltiy, CB user support
% Upload of the new indicator
dataset_ni = readtable('new_indicators_egov.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_ni = dataset_ni.Properties.VariableNames;

dim = size(dataset_ni);
n_vars_ni=dim(2)-1;
n_countries_ni = dim(1);
n_vars_ni = length(var_names_ni)-1;

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

fid = fopen('output_Digital_P_S_ni.txt','w'); %open text file to put results
for i=1:n_vars_ni
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

% Standardizing the new indicators
% Z scores of the new indicator
temp_D_1 = [data_ni];
i_rows=[];
[i_rows,tt]=find(isnan(data_ni)); %find rows WITH missing values
temp_D_1(i_rows,:)=[];            %keep only rows without missing values, 
temp_D_1_Z = zscore(temp_D_1);
                           
data_std_Z_ni = [temp_DD(:,[1 2]) temp_D_1_Z];           % merging togheter already std
                                                        % data for 4a1, 4a2, 4a3 and 4a4
                                                        % (excluding 4a3 because
                                                        % 50% online ava 50% CB online ava
                                             % with std data for new indicators: 
                                             % Online availability,
                                             % Mobile friendliness, 
                                             % User support,
                                             % Transparency (aggregate), 
                                             % CB online availabiltiy
                                             % CB user support
                                             % --> Total 10 indicators
size(data_std_Z_ni)  % 27 x 10             

%-----------------------------------------------
% PCA_3: Replicating the Classical analysis with all indicators except 4a5
%        Including new indicators: Mobile friendliness, user support,
%        transparency (aggregate), CB online availabiltiy, CB user support
[r_D3,p_D3]=corrcoef(data_std_Z_ni);
[loadings_D3,scores_D3,eigenvalues_D3]=pca(data_std_Z_ni); % computation of PCA on the standardized indicators 

writematrix(r_D3,'corr_D3.xls');
writematrix(loadings_D3(:,1:3),'loadings_D3.xls')

% Cronbach alpha:
c_D3 =cronbach_alpha(data_std_Z_ni);

% Plot of PC variances
figure
plot(eigenvalues_D3,'bo-')
axis([1 10 0 7])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
set(gca,'YTick',[0 1 2 3 4 5 6 7])
title('Digital Public Services-Variances')

% Plot of explained variances
figure
plot(eigenvalues_D3./sum(eigenvalues_D3),'bo-')
axis([1 10 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Public Services - Explained variances')

% PCA_3.2: 
%        Separate into 2 sub-dimensions
%----------- National sub-dimension -----------
[r_D4, p_D4] = corrcoef(data_std_Z_ni(:,1:7));
[loadings_D4,scores_D4,eigenvalues_D4]=pca(data_std_Z_ni(:,1:7)); % computation of PCA on the standardized indicators -C

writematrix(r_D4,'corr_D4.xls');
writematrix(loadings_D4(:,1:3),'loadings_D4.xls');

% Cronbach alpha:
c_D4 = cronbach_alpha(data_std_Z_ni(:,1:7));

% Plot of PC variances
figure
plot(eigenvalues_D4,'bo-')
axis([1 7 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7])
set(gca,'YTick',[0 1 2 3 4 5 6 ])
title('Digital Public Services-Variances')

% Plot of explained variances
figure
plot(eigenvalues_D4./sum(eigenvalues_D4),'bo-')
axis([1 7 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Public Services - Explained variances')

%------------ CB sub-dimension ----------------
[r_D5,p_D5]=corrcoef(data_std_Z_ni(:,8:10));
[loadings_D5,scores_D5,eigenvalues_D5]=pca(data_std_Z_ni(:,8:10)); % computation of PCA on the standardized indicators -C

writematrix(r_D5,'corr_D5.xls');
writematrix(loadings_D5(:,1:3),'loadings_D5.xls');

% Cronbach alpha:
c_D5 =cronbach_alpha(data_std_Z_ni(:,8:10));

% Plot of PC variances
figure
plot(eigenvalues_D5,'bo-')
axis([1 3 0 3])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 1 2 3 4 5 6])
title('Digital Public Services-Variances')

% Plot of explained variances
figure
plot(eigenvalues_D5./sum(eigenvalues_D5),'bo-')
axis([1 3 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Public Services -Explained variances')


% PCA_4: Replicating the Classical analysis with all indicators in
%          DESI2022 +  Mobile friendliness + User support (aggregated)
%          Transparency (aggregate)
user_support = 0.5* data_ni(:,4)+0.5*data_ni(:,8);
min_user = min(user_support);
max_user = max(user_support);
temp_DPS = [temp_D temp_D_1_Z(:,3) zscore(user_support) temp_D_1_Z(:,5)];
[r_D5]=corrcoef(temp_DPS);
[loadings_D5,scores_D5,eigenvalues_D5]=pca(temp_DPS); % computation of PCA on the standardized indicators 

writematrix(r_D5,'corr_D5.xls');
writematrix(loadings_D5(:,1:3),'loadings_D5.xls');

% Cronbach alpha:
c_D5 = cronbach_alpha(temp_DPS);

% Plot of PC variances
figure
plot(eigenvalues_D5,'bo-')
axis([1 8 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7 8])
set(gca,'YTick',[0 1 2 3 4 5 6 ])
title('Digital Public Services- Variances')

% Plot of explained variances
figure
plot(eigenvalues_D5./sum(eigenvalues_D5),'bo-')
axis([1 8 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7 8])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Public Services - Explained variances')

% PCA_4.2: Removing Open Data
temp_DPS_2 = [temp_D(:,1:4) temp_D_1_Z(:,3) zscore(user_support) temp_D_1_Z(:,5)];
[r_D6,p_D6]=corrcoef(temp_DPS_2);
[loadings_D6,scores_D6,eigenvalues_D6]=pca(temp_DPS_2); % computation of PCA on the standardized indicators -C

writematrix(r_D6,'corr_D6.xls');
writematrix(loadings_D6(:,1:3),'loadings_D6.xls');

% Cronbach alpha:
c_D6 = cronbach_alpha(temp_DPS_2);

% Plot of PC variances
figure
plot(eigenvalues_D6,'bo-')
axis([1 7 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7])
set(gca,'YTick',[0 1 2 3 4 5 6])
title('Digital Public Services-Variances')

% Plot of explained variances
figure
plot(eigenvalues_D6./sum(eigenvalues_D6),'bo-')
axis([1 7 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Public Services -Explained variances')


