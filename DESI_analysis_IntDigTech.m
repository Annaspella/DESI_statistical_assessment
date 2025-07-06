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
weights = dataset2(:,"Weight");  %weights relative to the indicators linear combination

% Load weights of each subdimension
dataset3 = readtable('Weigths_Subdim.xlsx');
weights_subdim = table2array(dataset3(:,"Weights"));


fid = fopen('output_all.txt','w'); %open text file to put results
for i=1:n_vars
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
        string=strcat({'Maximum value:  '},num2str(max(1,i)));
        fprintf(fid,'%s\n','maximum');
        fprintf(fid,'%12.4f \n',max(1,i));
        disp(string);
        string=strcat({'Country corresponding to maximum value:  '},country_max(1,i));
        fprintf(fid,'%s\n','Country corresponding to maximum value:');
        fprintf(fid,'%s\n',country_max{1,i});
        disp(string);
        string=strcat({'Minimum value:  '},num2str(min(1,i)));
        fprintf(fid,'%s\n','minimum:');
        fprintf(fid,'%12.4f \n',min(1,i));
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

%save matrix of standardized indicators: Zscore
save DESI_std.txt data_std_Z -ascii

%% Integration of Digital Technology:
fid = fopen('output_IDT.txt','w'); %open text file to put results
for i=18:28
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
temp_I=data_std_Z(:,18:28);
i_rows=[];
[i_rows,tt]=find(isnan(temp_I)); %find rows WITH missing values
temp_I(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_I,p_I]=corrcoef(temp_I);
[loadings_I,scores_I,eigenvalues_I]=pca(temp_I);%computation of PCA on the standardized indicators -H

writematrix(r_I,'corr_I.xls')
writematrix(loadings_I(:,1:3),'loadings_I.xls')
writematrix(eigenvalues_I,'eigen_I.xls')
writematrix(100*(eigenvalues_I./sum(eigenvalues_I)),'cumulative_I.xls')

% Cronbach alpha:
c_I1 =cronbach_alpha(temp_I);

% Plot of PC variances
figure
plot(eigenvalues_I,'bo-')
axis([1 11 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 6 7 8 9 10 11])
set(gca,'YTick',[0 1 2 3 4 5])
title('IntergrationDigitalTechnology- Variances')

% Plot of explained variances
figure
plot(eigenvalues_I./sum(eigenvalues_I),'bo-')
axis([1 11 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 6 7 8 9 10 11])
set(gca,'YTick',[0 1 2 3 4 5])
title('IntergrationDigitalTechnology - Explained variances')


%------------ Robust analysis-----------------------------
% Make use of functions contained in folder LIBRA_20160628 
% ROBPCA: to avoid the curse of dimensionality, n=28 must be greater than 5k
%         thus we set k_max equal to 4
k=11;
results_I = robpca(temp_I,'k',k, 'alpha',0.75, 'plots',1);

% Loadings:
results_I.P
writematrix(results_I.P(:,1:3), 'loadings_ROB_I.xls')

% eigenvalues 
results_I.L
results_I.L ./sum(results_I.L)

%grafici per ROBPCA: score plot 
results_I_rob = robpca(temp_I, 'k', 2, 'alpha',0.75,'classic', 0);

%grafici per CPCA: score plot CPCA
robpca(temp_I, 'k', 2, 'classic', 1);

%ellipsplot
mcd_I = mcdcov(results_I_rob.T, 'plots', 0);
dist = repelem(1,27);
ellipsplot(mcd_I.center, mcd_I.cov, results_I_rob.T, dist=dist)

mcd_I_class = mcdcov(scores_I(:,1:2), 'plots', 0);
ellipsplot(mcd_I_class.center, mcd_I_class.cov, scores_I(:,1:2), dist=dist)

%-------------------------------------
% Status quo
%-------------------------------------
% PCA: Replicating the Classical PCA over 3b1, 3b2, 3b3, 3b4, 3b5, 3b6, 3b7 
% separated from 3c1, 3c2 and 3c3:
%
% group 1: 3b1, 3b2, 3b3, 3b4, 3b5, 3b6, 3b7

temp_II=temp_I(:,2:8); 
i_rows=[];
[i_rows,tt]=find(isnan(temp_II)); %find rows WITH missing values
temp_II(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_II,p_II]=corrcoef(temp_II);
[loadings_II,scores_II,eigenvalues_II]=pca(temp_II); %computation of PCA on the standardized indicators IDT

writematrix(r_II,'corr_II.xls');
writematrix(loadings_II(:,1:3),'loadings_II.xls')
writematrix(eigenvalues_II,'eigen_II.xls')
writematrix(100*(eigenvalues_II./sum(eigenvalues_II)),'cumulative_II.xls')

% Cronbach alpha:
c_I2 =cronbach_alpha(temp_II);

% Plot of PC variances
figure
plot(eigenvalues_II,'bo-')
axis([1 7 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7])
set(gca,'YTick',[0 1 2 3 4 5])
title('Digital Technologies-Variances')

% Plot of explained variances
figure
plot(eigenvalues_II./sum(eigenvalues_II),'bo-')
axis([1 7 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Technologies -Explained variances')


% group 2: 3c1, 3c2, 3c3
temp_II=temp_I(:,9:11); 
i_rows=[];
[i_rows,tt]=find(isnan(temp_II)); %find rows WITH missing values
temp_II(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_II,p_II]=corrcoef(temp_II);
[loadings_II,scores_II,eigenvalues_II]=pca(temp_II);%computation of PCA on the standardized indicators -C

writematrix(r_II,'corr_II.xls');
writematrix(loadings_II(:,1:3),'loadings_II.xls')
writematrix(eigenvalues_II,'eigen_II.xls')
writematrix(100*(eigenvalues_II./sum(eigenvalues_II)),'cumulative_II.xls')


% Cronbach alpha:
c_I3 =cronbach_alpha(temp_II);

% Plot of PC variances
figure
plot(eigenvalues_II,'bo-')
axis([1 3 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 1 2 3 4 5])
title('Digital Technologies-Variances')

% Plot of explained variances
figure
plot(eigenvalues_II./sum(eigenvalues_II),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Technologies -Explained variances')


%-------------------------------------
% Adjustment 
%-------------------------------------
% PCA_2: Replicating the Classical analysis with all indicators except 3b6:
temp_I_new=temp_I(:,[1:6 8 9 10 11]); 
i_rows=[];
[i_rows,tt]=find(isnan(temp_I_new)); %find rows WITH missing values
temp_I_new(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_II,p_II]=corrcoef(temp_I_new);
[loadings_II,scores_II,eigenvalues_II]=pca(temp_I_new);%computation of PCA on the standardized indicators -C

writematrix(r_II,'corr_II.xls');
writematrix(loadings_II(:,1:3),'loadings_II.xls')
writematrix(eigenvalues_II,'eigen_II.xls')
writematrix(100*(eigenvalues_II./sum(eigenvalues_II)),'cumulative_II.xls')

% Cronbach alpha:
c_I4 =cronbach_alpha(temp_II);

% Plot of PC variances
figure
plot(eigenvalues_II,'bo-')
axis([1 10 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
set(gca,'YTick',[0 1 2 3 4 5])
title('Digital Technologies-Variances')

% Plot of explained variances
figure
plot(eigenvalues_II./sum(eigenvalues_II),'bo-')
axis([1 10 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Technologies -Explained variances')

%-----------------------------------------------
% PCA: Replicating the Classical PCA over 3b1, 3b2, 3b3, 3b4, 3b5, 3b7 separated from 3c1, 3c2 and 3c3:

% group 1: 3a1, 3b1, 3b2, 3b3, 3b4, 3b5, 3b7
temp_II=temp_I_new(:,1:7); 
i_rows=[];
[i_rows,tt]=find(isnan(temp_II)); %find rows WITH missing values
temp_II(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_II,p_II]=corrcoef(temp_II);
[loadings_II,scores_II,eigenvalues_II]=pca(temp_II);%computation of PCA on the standardized indicators -C

writematrix(r_II,'corr_II.xls');
writematrix(loadings_II(:,1:3),'loadings_II.xls')
writematrix(eigenvalues_II,'eigen_II.xls')
writematrix(100*(eigenvalues_II./sum(eigenvalues_II)),'cumulative_II.xls')

% Cronbach alpha:
c_I5 =cronbach_alpha(temp_II);

% Plot of PC variances
figure
plot(eigenvalues_II,'bo-')
axis([1 6 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'YTick',[0 1 2 3 4 5])
title('Digital Technologies-Variances')

% Plot of explained variances
figure
plot(eigenvalues_II./sum(eigenvalues_II),'bo-')
axis([1 6 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Technologies -Explained variances')


% group 2: 3c1, 3c2, 3c3
temp_II=temp_I_new(:,8:10); 
i_rows=[];
[i_rows,tt]=find(isnan(temp_II)); %find rows WITH missing values
temp_II(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_II,p_II]=corrcoef(temp_II);
[loadings_II,scores_II,eigenvalues_II]=pca(temp_II); %computation of PCA on the standardized indicators -C

writematrix(r_II,'corr_II.xls');
writematrix(loadings_II(:,1:3),'loadings_II.xls')
writematrix(eigenvalues_II,'eigen_II.xls')
writematrix(100*(eigenvalues_II./sum(eigenvalues_II)),'cumulative_II.xls')


% Cronbach alpha:
c_I6 =cronbach_alpha(temp_II);

% Plot of PC variances
figure
plot(eigenvalues_II,'bo-')
axis([1 3 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 1 2 3 4 5])
title('Digital Technologies-Variances')

% Plot of explained variances
figure
plot(eigenvalues_II./sum(eigenvalues_II),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Digital Technologies -Explained variances')



