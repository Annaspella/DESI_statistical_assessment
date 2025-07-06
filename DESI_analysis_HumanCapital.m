%% Cleaning and standardization of the data of HUMAN CAPITAL dimension

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

% Reading data with missing value to calculate the percentage of
% missing values
% data2: table containing the raw data for DESI year 6 (2021)
dataset2=readtable('DESI_missing.csv', ReadVariableNames=true, VariableNamingRule='preserve');
data2=table2array(dataset2(1:n_countries,3:end)); % excluding EU from the analysis

for i=1:n_vars
    temp=find(isnan(data2(:,i)));
    temp1=size(temp);
    number_nan=temp1(1);
    perc_miss(1,i)=number_nan/n_countries*100;
end

% Standardizing data without missing values
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

% save matrix of standardized indicators: Zscore
save DESI_std.txt data_std_Z -ascii

%% Human Capital:
fid = fopen('output_humancapital.txt','w'); %open text file to put results
for i=1:7  % 1:7, indexing the 7 indicators of Human Capital
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
% Classical PCA over all indicators in Human Capital
temp_H=data_std_Z(:,1:7);        % temp_H: contains the standardized data for Human Capital
i_rows=[];
[i_rows,tt]=find(isnan(temp_H)); %find rows WITH missing values
temp_H(i_rows,:)=[];             %keep only rows without missing values, 
                               
[r_H,p_H]=corrcoef(temp_H);
[loadings_H,scores_H,eigenvalues_H]=pca(temp_H); 

% Cronbach alpha:
c_H1 =cronbach_alpha(temp_H);
eigen_H1 = 100*(eigenvalues_H./sum(eigenvalues_H));

% Plot of PC variances
figure
plot(eigenvalues_H,'bo-')
axis([1 7 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 6 7 ])
set(gca,'YTick',[0 1 2 3 4 5])
title('Human Capital-Variances')

% Plot of explained variances
figure
plot(eigenvalues_H./sum(eigenvalues_H),'bo-')
axis([1 7 0 1])
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 6 7])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Human Capital-% of variances')


% Plotting  the standardize scores for Female ICT specialists against ICT specialists
% To show how the two are not correlated
country_names=dataset(1:n_countries,2);

figure
labels = table2array(country_names);
scatter(temp_H(:,4), temp_H(:,5),25,'b','o')
axis([-3 3 -3 3])
hline = refline([0 0]);
vline = xline(0);
set(hline,'Color','r')
set(vline,'Color','r')
text(temp_H(:,4), temp_H(:,5), labels, 'FontSize', 8, 'VerticalAlignment','bottom','HorizontalAlignment','left')
ylabel('Std index values for Female ICT specialiss')
xlabel('Std index values for ICT specialists')
title('Human Capital: ICT specialists vs Female ICT specialsits')


%------------ Robust analysis-------------------------------
% Make use of functions contained in folder LIBRA_20160628 
% ROBPCA: to avoid the curse of dimensionality, n=28 must be 
%         greater than 5k thus we set k_max equal to 4
k=4;
results = robpca(temp_H,'k',k, 'alpha',0.75, 'plots',1);

% Loadings:
results.P

% eigenvalues 
results.L
results.L./sum(results.L);

%grafici per ROBPCA: score plot 
results_rob = robpca(temp_H, 'k', 2, 'alpha',0.75,'classic', 0);

%grafici per CPCA: score plot CPCA, select 2 dimensions
robpca(temp_H, 'k', 2, 'classic', 1);

%ellipsplot
mcd = mcdcov(results_rob.T, 'plots', 0);
dist = repelem(1,27);
ellipsplot(mcd.center, mcd.cov, results_rob.T, dist=dist)

mcd_class = mcdcov(scores_H(:,1:2), 'plots', 0);
ellipsplot(mcd_class.center, mcd_class.cov, scores_H(:,1:2), dist=dist)


%----------------------------------------------
% PCA: Replicating the Classical PCA over 1a1, 1a2, 1a3 
% separated from 1b1, 1b2, 1b3, 1b4:
%-----------------------
% group 1: 1a1, 1a2, 1a3 
temp_HH=temp_H(:,[1:3]);
i_rows=[];
[i_rows,tt]=find(isnan(temp_HH)); %find rows WITH missing values
temp_HH(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_H1,p_H1]=corrcoef(temp_HH);
[loadings_H1,scores_H1,eigenvalues_HH]=pca(temp_HH); 

% Cronbach alpha and percentage of variance explained by first component:
c_H2 =cronbach_alpha(temp_HH);
eigen_H2 = 100*(eigenvalues_HH./sum(eigenvalues_HH));

% Plot of PC variances
figure
plot(eigenvalues_HH,'bo-')
axis([1 3 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 1 2 3 4 5])
title('Human Capital-Variances')

% Plot of explained variances
figure
plot(eigenvalues_HH./sum(eigenvalues_HH),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Human Capital-Explained variances')

%-----------------------------
% group 2: 1b1, 1b2, 1b3, 1b4:
temp_HH=temp_H(:,[4 5 6 7]);
i_rows=[];
[i_rows,tt]=find(isnan(temp_HH)); %find rows WITH missing values
temp_HH(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_H2,p_H2]=corrcoef(temp_HH);
[loadings_H2,scores_H2,eigenvalues_H2]=pca(temp_HH); 

% Cronbach alpha and percentage of variance explained by first component:
c_H3 = cronbach_alpha(temp_HH);
eigen_H3 = 100*(eigenvalues_H2./sum(eigenvalues_H2));

% Plot of PC variances
figure
plot(eigenvalues_H2,'bo-')
axis([1 4 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 1 2 3 4 5])
title('Human Capital - Variances')

% Plot of explained variances
figure
plot(eigenvalues_H2./sum(eigenvalues_H2),'bo-')
axis([1 4 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Human Capital - Explained variances')

%-----------------------------------------------
% Minor Adjustment final results
%-----------------------------------------------
% PCA_2: Replicating the Classical PCA over 1a1, 1a2, 1a3, 1b3 
% separated from 1b1, 1b4:
%----------------------------
% group 1: 1a1, 1a2, 1a3, 1b3
temp_HH=temp_H(:,[1:3 6]);
i_rows=[];
[i_rows,tt]=find(isnan(temp_HH)); %find rows WITH missing values
temp_HH(i_rows,:)=[];             %keep only rows without missing values, 
                               
[r_H3,p_H3]=corrcoef(temp_HH);
[loadings_H3,scores_H3,eigenvalues_H3]=pca(temp_HH); 

% Cronbach alpha:
c_H4 =cronbach_alpha(temp_HH);
eigen_H4 =100*(eigenvalues_H3./sum(eigenvalues_H3));

% Plot of PC variances
figure
plot(eigenvalues_H3,'bo-')
axis([1 4 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 1 2 3 4 5])
title('Human Capital-Variances')

% Plot of explained variances
figure
plot(eigenvalues_H3./sum(eigenvalues_H3),'bo-')
axis([1 4 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Human Capital-Explained variances')

%-------------------
% group 2: 1b1, 1b4:
temp_HH=temp_H(:,[4 7]);
i_rows=[];
[i_rows,tt]=find(isnan(temp_HH)); %find rows WITH missing values
temp_HH(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_H4,p_H4]=corrcoef(temp_HH);
[loadings_H4,scores_H4,eigenvalues_H4]=pca(temp_HH); 

writematrix(r_HH,'corr_HH.xls');
writematrix(loadings_HH(:,1:2),'loadings_HH.xls')

% Cronbach alpha:
c_H5 = cronbach_alpha(temp_HH);
eigen_H5 = 100*(eigenvalues_H4./sum(eigenvalues_H4));

% Plot of PC variances
figure
plot(eigenvalues_H4,'bo-')
axis([1 2 0 5])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2])
set(gca,'YTick',[0 1 2 3 4 5])
title('Human Capital - Variances')

% Plot of explained variances
figure
plot(eigenvalues_H4./sum(eigenvalues_H4),'bo-')
axis([1 3 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Human Capital - Explained variances')


%--------------------------------------------------
% Major adjustment : FINAL VERSION
% -------------------------------------------------
% Adding new indicators to Human Capital:
%
% new_indcators_ HumanCapital.csv contains
% 1st row Country labels
% 2nd row Never used internet indicator (Opposite orientation)
% 3rd row Frequency of internet use indicator
dataset_ni=readtable('new_indicators_HumanCapital.csv', ReadVariableNames=true, VariableNamingRule='preserve');
var_names_ni=dataset_ni.Properties.VariableNames;

dim=size(dataset_ni);
n_vars_ni=dim(2)-1;          % Removing country column
n_countries_ni=dim(1)-1;     % Not considering EU
EU_ni = dataset(dim(1),:);
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

fid = fopen('output_humancapital_ni.txt','w'); %open text file to put results
for i=1:2
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

% Standarsizing the new indicators:
% Z scores of the new indicators, kept with their original orientations
data_std_Z_ni = zscore(data_ni);

%---------------------------
% PCA over ALL INDICATORS (except ICT female specialist) + new indicators
temp_HH = [temp_H(:,[1:4 6 7]) data_std_Z_ni];
i_rows=[];
[i_rows,tt]=find(isnan(temp_HH)); %find rows WITH missing values
temp_HH(i_rows,:)=[];             %keep only rows without missing values, 
                               
[r_H5,p_H5]=corrcoef(temp_HH);
[loadings_H5,scores_H5,eigenvalues_H5]=pca(temp_HH); 

% Cronbach alpha:
c_H6 = cronbach_alpha(temp_HH);
eigen_H6 = 100*(eigenvalues_H5./sum(eigenvalues_H5));

% Plot of PC variances
figure
plot(eigenvalues_H5,'bo-')
axis([1 8 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7 8])
set(gca,'YTick',[0 1 2 3 4 5 6 7 8])
title('Human Capital-Variances')

% Plot of explained variances
figure
plot(eigenvalues_H5./sum(eigenvalues_H5),'bo-')
axis([1 8 0 1]);
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6 7 8])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Human Capital-Explained variances')


% Replicating the Classical PCA over 1a1, 1a2, 1a3, 1b3 + new indicators
% separated from 1b1, 1b4:
%----------------------------
% group 1: 1a1, 1a2, 1a3, 1b3 + new indicators
temp_HH= [temp_H(:,[1:3 6]) data_std_Z_ni];
i_rows=[];
[i_rows,tt]=find(isnan(temp_HH)); %find rows WITH missing values
temp_HH(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_H6,p_H6]=corrcoef(temp_HH);
[loadings_H6,scores_H6,eigenvalues_H6]=pca(temp_HH);%computation of PCA on the standardized indicators -C

writematrix(r_HH,'corr_HH.xls');
writematrix(loadings_HH(:,1:3),'loadings_HH.xls')

% Cronbach alpha:
c_H7 = cronbach_alpha(temp_HH);
eigen_H7 = 100*(eigenvalues_H6./sum(eigenvalues_H6));

% Plot of PC variances
figure
plot(eigenvalues_H6,'bo-')
axis([1 6 0 6])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'YTick',[0 1 2 3 4 5])
title('Human Capital-Variances')

% Plot of explained variances
figure
plot(eigenvalues_H6./sum(eigenvalues_H6),'bo-')
axis([1 6 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 3 4 5 6])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Human Capital-Explained variances')


%-------------------
% group 2: 1b1, 1b4
temp_HH= temp_H(:,[4 7]);
i_rows=[];
[i_rows,tt]=find(isnan(temp_HH)); %find rows WITH missing values
temp_HH(i_rows,:)=[];             %keep only rows without missing values, 
                              
[r_H7,p_H7]=corrcoef(temp_HH);
[loadings_H7,scores_H7,eigenvalues_H7]=pca(temp_HH);%computation of PCA on the standardized indicators

writematrix(r_HH,'corr_HH.xls');
writematrix(loadings_HH(:,1:2),'loadings_HH.xls')

% Cronbach's alpha:
c_H8 = cronbach_alpha(temp_HH);
eigen_H8 = 100*(eigenvalues_H7./sum(eigenvalues_H7));

% Plot of PC variances
figure
plot(eigenvalues_H7,'bo-')
axis([1 2 0 2])
hline = refline([0 1]);
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2 ])
set(gca,'YTick',[0 1 2])
title('Human Capital-Variances')

% Plot of explained variances
figure
plot(eigenvalues_H7./sum(eigenvalues_H7),'bo-')
axis([1 2 0 1])
set(hline,'Color','r')
xlabel('dimension')
ylabel('eigenvalue')
set(gca,'XTick',[1 2])
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
title('Human Capital-Explained variances')

