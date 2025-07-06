%% The function get_score:
%  Receives as an input a matrix of scores (containing NaN values),
%  and a vector of raw weights for the indicators contained in the matrix
%  (number of columns of score = number of elements of weights)

%  Returns a vector containing the scores at the subdimension levele
%  where the NaN values has been ignored and the final sub_score values
%  are calculated on the basis of the non NaN values by row 
% (weights are standardize differently for each row)

function sub_score = get_score(scores, weights)
    [r,~] = size(scores);
    sub_score = [];
    
    for i = 1:r                                          % Looping inside rows
        ind = isnan(scores(i,:));                        % Returns a vector of 0 and 1, 1 if is NaN
       
        weights_new = weights(~ind)/sum(weights(~ind));  % Weights standardization: ignoring the NaN indicators
        scores_new = scores(i,~ind);                     
        sub_score = [sub_score scores_new*weights_new];  % Sub-scores weights calculation: only based on the non NaN indicators
    end  
        sub_score = sub_score';
    return 
end