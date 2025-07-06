function [a, item_variances, item_number, subject_number] = cronbach_alpha(item_scores_mat)
    % Calculate Cronbach's alpha for item_scores_mat is an SxI numeric matrix of raw scores with a row for each subject and a column for each item.
    % All items should be on the same scale.
    % item_scores_mat cannot contain any missing values (this would make some items less reliable than others).
    % Cronbach's alpha is an index of internal consistency (e.g. for a questionnaire scale)
    % Cronbach's alpha >=0.80 is considered good internal consistency.
    % Note that alpha can sometimes be negative. 
    % This could happen if you have items that should have been reversed scored.
    % E Zakreski 2018
    
    % Check input
    assert(isnumeric(item_scores_mat),...
        'cronbach_alpha:input_not_numeric',...
        'item_scores_mat must be numeric.');
    assert(ismatrix(item_scores_mat),...
        'cronbach_alpha:input_not_matrix',...
        'item_scores_mat must be a numeric matrix with a row for each subject and item for each column.');
    % Warn if any scores are missing 
    if any(any(isnan(item_scores_mat)))
        error('cronbach_alpha:input_has_missing_scores',...
            'item_scores_mat contains one or more missing values (NaN''s).');
    end

    % How many items and subject
    [subject_number, item_number] = size(item_scores_mat);
    
    % Get the variance of each item across subjects
    item_variances = var(item_scores_mat,1,1); %use biased variance (w=1)
    
    % Get each subject's total score
    subject_total_scores = sum(item_scores_mat,2);
    
    % Cronbach's alpha
    a = 1 - sum(item_variances)/var(subject_total_scores,1); %use biased variance (w=1)
    
    a = item_number/(item_number - 1)*a;

end