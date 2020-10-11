function [ Exch, sellersUp, buyersUp, DUp, loss ] = Energy_Exchange( DUp, sellersUp, buyersUp, seller_partner, ....
                                                         dist, loss, Exch )
%Energy_Exchange calculates the energy to be exchanged between each pair
%(b,s)
%  
sellersUp2=sellersUp;
buyersUp2=buyersUp;
    
%Parameters of our distribution line system   
U1=22*10^3;                 % (kV) voltage of the local electric network
R= 0.2; %*10^(-3);          % (ohm/km) losses in the main grid


    for s=1:length(sellersUp)  
            b=seller_partner(s);    % b is the buyer itself
            if b~=0
                Eqn=[(R*dist(sellersUp(s),b)/(U1^2)) -1 DUp(b)];
                r=roots(Eqn);
                if isreal(r)
                    r(r<=0)=nan;
                    Sol=min(r); 
                else %Delta<0
                    Sol=U1^2/(2*R*dist(sellersUp(s),b)); 
                end
                
                E = min (Sol, abs(DUp(sellersUp(s))) );
                loss(b,sellersUp(s)) = E^2*R*dist(sellersUp(s),b)/(U1^2);
                DUp(sellersUp(s)) = DUp(sellersUp(s)) + E; 
                DUp(b)= DUp(b)- ( E - loss(b,sellersUp(s)) ); 

                if DUp(b)==0
                   bId=buyersUp2==b;
                   buyersUp2(bId)=[]; 
                end
                if DUp(sellersUp(s))==0
                   SId= sellersUp2 == sellersUp(s);
                   sellersUp2(SId)=[];
                end            
    
                Exch(sellersUp(s),b)=E;
            end
    end
    sellersUp=sellersUp2;
    buyersUp=buyersUp2;

