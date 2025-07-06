%% The function get_rank: 
% Receives as an input a list of scores
% Returns as an output the list of the ranking

function ranking = get_rank(Index)
    [~,ordered_index] = sort(Index,'descend');     % Returns two list: one of sorted elements 
    num = 1:length(ordered_index);                 % and on with the corresponging indexes. We take the second one
    ranking= [];                                    
    for i = 1:length(num)                          % The function returns the vector containeing the 
        j=1;                                       % the list of the ranking 
        while j <= length(num)
            if ordered_index(j)==i
                ranking = [ranking num(j)];
                break
            else 
                j = j+1;
            end
            
        end    
     end  
    return
end
