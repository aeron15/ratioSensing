function dy = SimulatingGrowthSequntiolEq(t,y,g1,g2,d1,d2,ts,ti,cost,r)




for i = 1:5
    y(i) = max([0 y(i)]);
end
dy = zeros(5,1);    % a column vector
%% Define states

%%
 % Cell is not switchd (cant grow on 2)
 if y(4)<ts
     dy(5) = 0; % induction timer is off
     % Glc >0
     if y(2)>0
         g = g1;
         u1 = g1*d1;
         u2 = 0;
     end
     % Glc = 0
     if y(2)==0
         g = 0;
         u1 = 0;
         u2 = 0;
     end
 end

%%
% Cell switchd
if y(4)>=ts
    dy(5) = 1; % induction timer is on
    % induction is NOT done
    if y(5) < ti
        % Glc >0
        if y(2)>0
            g = g1/cost;
            u1 = (g1/cost)*d1*r;
            u2 = 0;
        end
        % Glc = 0
        if y(2)==0
            g = 0;
            u1 = 0;
            u2 = 0;
        end
    end
    
    % induction is DONE
    if y(5) >= ti
        % Glc > N1*d
        if y(2)>y(1)*d1
            g = g1/cost;
            u1 = (g1/cost)*d1*r;
            u2 = 0;
        end
        % Glc = < N1*d
        if y(2) <= y(1)*d1
            % Gal > 0
            if y(3) > y(1)*d1
                g = g1/cost * y(2)/(y(1)*d1) + g2*(1 -y(2)/(y(1)*d1)) ;
                u1 =  g1 * y(2)/(y(1)*d1)*d1 ;
                u2 =  g2*(1 -y(2)/(y(1)*d1))*d1 ;
            end
            % Gal = 0
            if y(3) <= y(1)*d1
                g =g1/cost * y(2)/(y(1)*d1) + g2*y(3)/(y(1)*d1);
                u1 =  g1 * y(2)/(y(1)*d1)*d1  ;
                u2 =  g2 * y(3)/(y(1)*d1)*d1  ;
            end
        end
    end
    
    
    
    
    
    
    
    
end
%%
                 
    dy(1) =  y(1)*g; % N
    dy(2) = -y(1)*u1; % Glc
    dy(3) = -y(1)*u2; % GAl
    dy(4) = 1; % time while y(5) is time from induction
    


