%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Cooperative Game Theory                               %
%                  Gale Shapley (Matching)                              %
%               Preference based on distance                            %
%       Preference based on losses shown in comments 42                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%Remark 1 = Loop over all the sellers-buyers is coded here               %
%Remark 2 = Loop over one buyer-seller association is coded in matching1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

function [l,lU,e, iter]  = matching(dist,distU,D,C,N)
l=0;   %losses in the networked µgrid
lU=0;  %losses due to exchange with the main grid 
iter=0;%nbr of iterations til convergence 

%Parameters of the distribution line system   
beta=0.02;                  % transformer loss rate
U0=50*10^3;                 % (kV) voltage of the utility energy 
%U1=22*10^3;                 % (kV) voltage of the local electric network
R= 0.2;                     % (ohm/km) losses in the main grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(C) > 1
    
    %Separate buyers from sellers
    [ buyers, sellers ] = split_buyers_sellers( D, C );    

    DUp=D;               %Updated matrix of energy
    sellersUp=sellers;   %Updated list of sellers
    buyersUp=buyers;     %Updated list of buyers
    loss=zeros(N, N);    %matrix of losses during the exchange
    Exch=zeros(N, N);    %matrix of energy exchanged between µgrids
    
   %%%    Execute the matching until the sum of energy to sell/buy =0       
    while ~(isempty(sellersUp) || isempty(buyersUp))
    
        %%% Matching Game buyer_optimized       
        [ PrefB, ~, rank ] = distance_rank( buyersUp, sellersUp, dist, distU );
        %%% In case rank based on losses
        %[ PrefB, PrefS, rank ] = loss_rank( buyersUp, sellersUp, dist )
        
        buyer_free=zeros(1,length(buyersUp));
        seller_partner=zeros(1,length(sellersUp));
        condition=zeros(1,1);

        while min(condition)==0
            for n=1:length(buyersUp)
                if buyer_free(n)~=0
                    continue
                end
                
                next=find(PrefB(n,:)>0,1); 
                seller=PrefB(n,next);
                PrefB(n,next)=0;
                
                selId=find(sellersUp==seller);
                if seller_partner(selId)==0
                    seller_partner(selId)=buyersUp(n);         
                    buyer_free(n)=1;
                else
                    comp=seller_partner(selId); % the buyer itself and not its index   
                    competitor= find(buyersUp==comp); %index of competitor in buyersUp
                    if rank(selId,n) < rank(selId, competitor)
                        seller_partner(selId)=buyersUp(n);     
                        buyer_free(n)=1;
                        buyer_free(competitor)=0;
                    end
                end
            end
            
            if length(buyersUp) <= length(sellersUp)
                condition=buyer_free;
            else
                condition=seller_partner;
            end
            
        end
 
    %%% Energy Exchange
       [ Exch, sellersUp, buyersUp, DUp, loss ] = Energy_Exchange( DUp, sellersUp, buyersUp, seller_partner, dist, loss, Exch );
        iter=iter+1;
        
    end       
   
    %losses of energy exchange between the subset and the utility
    for pf=1:length(C)          
        lU = lU + (DUp(C(pf))^2)*distU(C(pf))*R/(U0^2) + ( beta*abs(DUp(C(pf))) );       
    end
    
    %add internal losses
     l=sum(sum(loss));
     e=sum(sum(Exch))/sum(abs(D));

else 
    lU = ( D(C)^2*distU(C)*R/(U0^2) ) + ( beta*abs(D(C)) );
    e=0;
end