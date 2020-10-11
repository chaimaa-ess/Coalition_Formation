function [ PrefB, PrefS, rank ] = distance_rank( buyers, sellers, dist, distU )
%This function returns 3 matrices:
% PrefB: PrefB(b,:): lists the sellers given in order of preference of
% buyer b; the nearest seller to b is the most prefered
% PrefS: PrefS(s,:): lists the buyers given in order of preference of
% seller s
% Rank(s,index): returns the rank of b (which is buyers(index)) from s's perspective   

    PrefS=zeros(length(sellers),length(buyers));
        for s=1:length(sellers)
            [~, IndiceS]= sort(dist(sellers(s),buyers()));
            for indexS=1:length(IndiceS)
                PrefS(s, indexS)=buyers(IndiceS(indexS));
            end
        end

        PrefB=zeros(length(buyers), length(sellers));     
        for b=1:length(buyers)
            [~, IndiceB]= sort(dist(buyers(b),sellers()));
            for indexB=1:length(IndiceB)
                s = sellers(IndiceB(indexB)) ;
                if dist(b,s) > distU(b)
                    PrefB(b, indexB) = 0;
                else
                    PrefB(b, indexB)=s;
                end
            end
        end
        
        rank=zeros(length(sellers),length(buyers));
        for s=1:length(sellers)
            for b=1:length(buyers)
                for k=1:length(buyers)
                    if PrefS(s,k)==buyers(b)
                        rank(s,b)=k;
                    end
                end
            end  
        end

end

