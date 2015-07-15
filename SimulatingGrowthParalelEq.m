function dy = SimulatingGrowthSequntiolEq(t,y,g1,g2,d1,d2,ts,ti,cost,r)




for i = 1:5
    y(i) = max([0 y(i)]);
end
dy = zeros(5,1);    % a column vector


%% Creat a step fucntion growth rate
% if y(2) < 2^-6   % glc
%     g1 = 0;
% end
% if y(3) < 2^-4  % gal
%     g2 = 0;
% end


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

        g = (y(2)*g1 + y(3)*g2)/(y(2)+y(3));
        u1 = d1*r*(y(2)*g1)/(y(2)+y(3));
        u2 = d2*r*(y(3)*g2)/(y(2)+y(3));
    end
    

end
%%
                 
    dy(1) =  y(1)*g; % N
    dy(2) = -y(1)*u1*g; % Glc
    dy(3) = -y(1)*u2*g; % GAl
    dy(4) = 1; % time while y(5) is time from induction
    


