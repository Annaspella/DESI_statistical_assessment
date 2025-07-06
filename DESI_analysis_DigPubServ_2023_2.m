%% Cleaning and standardization of the data
clear all
clc
opengl software %this command solves the problem with the video card that can cause figures to be black (all black)

%read dataset with variable names in first row
%country in 1st column
%region in 2nd column
%region ID in 3rd column
%indicator names from column 4th on
dataset=readtable('DESI_Y7_2.csv', ReadVariableNames=true, VariableNamingRule='preserve'); % DESI_Y7_2 contains data 
                                                                                           % data updated to 2022
                                                                                           % 4a1 substituted by I_GOVANYS 
var_names=dataset.Properties.VariableNames;

dim=size(dataset);
n_vars=dim(2)-2;       % Minus two columns corresponding to year and country 
n_countries=dim(1)-1;  % Not considering EU

% Save EU line and delete from table "data"
EU = dataset(dim(1),:);

% data: table containing the data with imputed values for DESI year 7 (2022)
data=table2array(dataset(1:n_countries,3:end)); % excluding EU from the analysis
                                                % (EU last row, taking the first 27)


% Standardizing data values: 
% Z scores normalization for PCA 
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

% Load NEW weights for each indicator
% load NEW min-max for each indicator, (see "calculations minmax egov 2023data.xlsx")
dataset2 = readtable('min_max_2023.xlsx');
min_val = dataset2(:,"Min");
max_val = dataset2(:,"Max");

% save matrix of standardized indicators: Zscore
save DESI_std.txt data_std_Z -ascii

%% Digital Public Services: getting the summary statistics for indicators in DPS only

fid = fopen('output_DPS_2023_2.txt','w'); % open text file to save results
for i=29:33                               % 29:33 - indexing the 5 indicators in DPS
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
        data_std_minmax(:,i) = (data(:,i)-min_vec)./(max_vec-min_vec);    

        % Plotting histograms
         figure
            hist(data(:,i),50)
            string=strcat({'Indicator: '},var_names{1,i+2},{' '},{'skewness='},num2str(skew(1,i)));
            title(string)
            
end  % end of cycle on variables
fclose(fid);

%----------------------------------------------
% Status quo
%----------------------------------------------
%  Classical PCA over all indicators in DPS

temp_D=data_std_Z(:,29:33);

%i_rows=[];
%[i_rows,tt]=find(isnan(temp_D)); %find rows WITH missing values
%temp_I(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_D,p_D]=corrcoef(temp_D);
[loadings_D,scores_D,eigenvalues_D]=pca(temp_D);%computation of PCA on the standardized indicators -D

writematrix(r_D,'corr_DPS_Y3.xls')
writematrix(loadings_D(:,1:3),'loadings_DPS_Y3.xls')
writematrix(eigenvalues_D,'eigen_DPS_Y3.xls')
writematrix(100*(eigenvalues_D./sum(eigenvalues_D)),'cumulative_DPS_Y3.xls')

% Cronbach alpha:
c_D1 = cronbach_alpha(temp_D);

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
title('DigitalPublicServices - Scree plot')


%-----------------------------------------------
% Replicating the Classical analysis with all indicators except 4a5
%-----------------------------------------------
temp_DD=temp_D(:,[1:4]); 
i_rows=[];
[i_rows,tt]=find(isnan(temp_DD)); %find rows WITH missing values
temp_DD(i_rows,:)=[];             %keep only rows without missing values, 
                               %no necessary because imputation of missing
                               %data was previously done. We keep it as a
                               %control
[r_DD,p_DD]=corrcoef(temp_DD);
[loadings_DD,scores_DD,eigenvalues_DD]=pca(temp_DD);%computation of PCA on the standardized indicators -C

writematrix(loadings_DD(:,1:3),'loadings_DPS_2.xls')

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
title('Digital Public Services - Scree plot')

%-------------------------------------
% Uploading Mobile friendliness, user support, Online availability,
%           Transparency (aggregate), CB online availabiltiy, CB user support
% ------------------------------------
% Upload of the new indicator
dataset_ni = readtable('new_indicators_egov_2023.csv', ReadVariableNames=true, VariableNamingRule='preserve');
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

fid = fopen('output_DPS_ni_2023.txt','w'); %open text file to put results
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
                           
data_std_Z_ni = [temp_DD(:,[1 2]) temp_D_1_Z];          % merging togheter (already std
                                                        % data) for 4a1, 4a2, 4a3 and 4a4
                                                        % (excluding 4a3 and 4a4 because
                                                        % 50% online availability
                                                        % 50% CB online ava)
                                             % with std data for new indicators: 
                                             % Online availability,
                                             % Mobile friendliness, 
                                             % User support,
                                             % Transparency (aggregate), 
                                             % CB online availabiltiy
                                             % CB user support
                                             % --> Total 10 indicators
size(data_std_Z_ni)  % 27 x 10          

%--------------------------------------------
% Minor adjustment:   Replicating the Classical analysis with all indicators in
%                     DESI2022 - Open data +  Mobile friendliness + User support (aggregate)
%                     Transparency (aggregate)
%---------------------------------------------
user_support = 0.5* data_ni(:,4) + 0.5* data_ni(:,8); % Obtained as 0.5 national user support and 0.5 cb user support

temp_DPS = [temp_D(:,1:4) temp_D_1_Z(:,3) zscore(user_support) temp_D_1_Z(:,5)]; % taking standardized values for:
                                                                                 % 4a1, 4a2, 4a3, 4a4, Mobile friendliness 
                                                                                 % User support(aggregate),Transparency(aggregate)
[r_D6,p_D6]=corrcoef(temp_DPS);
[loadings_D6,scores_D6,eigenvalues_D6]=pca(temp_DPS); % computation of PCA on the standardized indicators -C

% 
% -------------------------------------------------
% Major adjustment: Replicating the Classical analysis with all indicators in
%                   DESI2022 - Open data CB 
%                   Including new indicators: Mobile friendliness, User
%                   support, Online availability bus, Online availability cit,
%                   Transparency (aggregate), CB online availabiltiy, CB user support
%---------------------------------------------------
[r_D3,p_D3]=corrcoef(data_std_Z_ni);
[loadings_D3,scores_D3,eigenvalues_D3]=pca(data_std_Z_ni); % computation of PCA on the standardized indicators 

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
title('Digital Public Services - Scree plot')

% ---------------------------------------------
%        Separate into 2 sub-dimensions
%----------- National sub-dimension -----------
[r_D4, p_D4] = corrcoef(data_std_Z_ni(:,1:7));
[loadings_D4,scores_D4,eigenvalues_D4]=pca(data_std_Z_ni(:,1:7)); % computation of PCA on the standardized indicators -C

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
title('Digital Public Services - Scree plot')

%------------ CB sub-dimension ----------------
[r_D5,p_D5]=corrcoef(data_std_Z_ni(:,8:10));
[loadings_D5,scores_D5,eigenvalues_D5]=pca(data_std_Z_ni(:,8:10)); % computation of PCA on the standardized indicators -C

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
title('Digital Public Services - Scree plot')
