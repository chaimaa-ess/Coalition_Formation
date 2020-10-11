function [ buyers, sellers ] = split_buyers_sellers( D, C )
%This function splits the list of participants into e subsets: 
% sellers with an excess of energy to sell
% buyers with a defficiency of energy to buy

    sellers=int16.empty(1,0);        
    buyers=int16.empty(1,0);
    for d=1:length(C)
        if D(C(d)) < 0
            sellers=[sellers, C(d)]; %#ok<AGROW>
        elseif D(C(d)) > 0
            buyers=[buyers, C(d)]; %#ok<AGROW>
        end
    end

end

