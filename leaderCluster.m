%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Cooperative Game Theory                               %
%                 Leader Cluster Algorithm                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=5;                       % length of the square
beta=0.02;                  % transformer loss rate
U0=50*10^3;                 % (kV) voltage of the utility energy 
U1=22*10^3;                 % (kV) voltage of the local electric network
R= 0.2; %*10^(-3);          % (ohm/km) losses in the main grid
%U1/2R the maximum to be sent = 242KWh with d=L=5

% Averaged statistical results for coalitions
Num=9;
Loss=zeros(1,Num);
Exchanged_Energy=zeros(1,Num);
avg_iter=zeros(1,Num);

% Averaged statistical results for grand coalition
LossGC=zeros(1,Num);
Exchanged_Energy_GC=zeros(1,Num);
avg_Exchanged_Energy_GC=zeros(1,Num);

%Loss in the classical network
Loss0=zeros(1,Num);

S=5000; % Nbr of time we run the simulation to get refined statistical results
 
n=1; %index to refere to the elements of the vector N

 for N=[ 5 7 10 13 15 17 20 25 30]          % # of MGs  
    iter=0; 
   for statistic = 1:S

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% MGs Placement %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X= L*rand(1,N);          % x coordonate array of MGs
        Y= L*rand(1,N);
        
    %%%%%%%%%%%%%%%%%%%%%% Distance Matrix   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dist=zeros(N,N);
        distU=zeros(1,N);
        for i=1:N
            for j=1:N
                dist(i,j)= sqrt ( (X(i)-X(j))^2 + (Y(i)-Y(j))^2 );
            end
            distU(i)= sqrt ( (X(i)-L/2)^2 + (Y(i)-L/2)^2 );
        end
        
    %%%%%%%%%%%%%%% Net Demand vector establishment %%%%%%%%%%%%%%%%%%%%%%%
    % Di is a Gaussian random variable with zero mean and a variance that is 
    % uniformly distributed between 10 and 1000 KW
        %sigmamin=3.162; sigmamax=31.62;    
        %sigmamin=1; sigmamax=3.162; %[1 10]
        sigmamin=sqrt(10); sigmamax=sqrt(100); 
        sigma = sigmamin + (sigmamax-sigmamin)*rand;
        D=normrnd(0, sigma, [1,N])*10^3*10^4;
        
         for pf=1:N
             Loss0(n)= Loss0(n) + ( (D(pf)^2)*distU(pf)*R/(U0^2) ) + ( beta*abs(D(pf)) ); %(n)           
         end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Coalition={1};
        Centroid = {[X(1) Y(1)]};
        for elt=2:N
            %dmin=distU(elt);
            value_ref=(D(elt)^2)*distU(elt)*R/(U0^2) + ( beta*abs(D(elt)) );           
            for c=1:length(Coalition)
                C_length=length(Coalition);
                d= sqrt ( (Centroid{c}(1)-X(elt))^2 + (Centroid{c}(2)-Y(elt))^2 );
%                 if d < dmin
%                    ref=c; dmin=d; 
%                 end
                if  value_ref >= (D(elt)^2)*d*R/(U1^2)
                   ref=c;
                   value_ref = (D(elt)^2)*d*R/(U1^2);
                end
               iter=iter+1;
            end
            if value_ref ~= (D(elt)^2)*distU(elt)*R/(U0^2) + ( beta*abs(D(elt)) )  %dmin ~= distU(elt) %
                Coalition{ref}=[Coalition{ref} elt];
                Centroid{ref}(1) = sum(X(Coalition{ref}))/length(Coalition{ref});
                Centroid{ref}(2) = sum(Y(Coalition{ref}))/length(Coalition{ref});
            else
                Coalition{1, length(Coalition)+1} = elt;
                Centroid{1,length(Coalition)}=[X(elt) Y(elt)];
            
            end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% Apply matching on the resulting Coalition %%%%%%%%%%%%
        for c=1:length(Coalition)
            C=Coalition{1,c}; 
            [l, e] = matching(dist, distU,D,C,N);
            Loss(n) = Loss(n) + l; %(n)
            Exchanged_Energy(n) = Exchanged_Energy(n) + e; %(n)
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
    %%%%%%%%%%%%%% Apply matching on GC (grand coalition) %%%%%%%%%%%%%%%%%
        [lGC, eGC] = matching(dist, distU, D, (1:N), N);
        LossGC(n) = LossGC(n) + lGC ; %(n)
        Exchanged_Energy_GC(n) = Exchanged_Energy_GC(n) + eGC; %(n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
   end %for statistic
%     
    Loss(n) = Loss(n)/(N*S);
    Exchanged_Energy(n) = Exchanged_Energy(n)/(S);
    avg_iter(n)=iter/S;
    
    LossGC (n) = LossGC(n)/(N*S);
    Exchanged_Energy_GC(n) = Exchanged_Energy_GC(n)/(N*S);
%     
    Loss0(n)= Loss0(n)/(N*S);
    
   n=n+1; 
end
save('ALC_D_10-8_loss_matchloss_exch.mat') 

xlab=[5 7 10 13 15 17 20 25 30];
figure
grid on
plot(xlab, Loss, 'g');
hold on
plot (xlab, LossGC);
hold on
plot (xlab, Loss0, 'r');
xlabel('N° of MGs'); ylabel('Average Loss per MG (KW)')
legend('New method', 'Grand Coalition', 'Classical Model' )

figure
grid on
plot(xlab, Exchanged_Energy, 'g');
hold on
plot(xlab, Exchanged_Energy_GC);
xlabel('N° of MGs'); ylabel('Average Exchange per MG')
legend('New method', 'Grand Coalition')
% 
figure
grid on
plot(xlab, avg_iter);
xlabel('N° of MGs'); ylabel('Average iter')
