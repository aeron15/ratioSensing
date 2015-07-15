function dy = SimulatingGrowthEq(t,y,ra1,ra2,rh1,rh2,tha,thh,ti,d,cost)




for i = 1:4
    y(i) = max([0 y(i)]);
end
dy = zeros(6,1);    % a column vector
%% Define states
% A
if y(3)>tha*y(4)
    ra = ra1;
    d1a = d;
    d2a = 0;
end
if y(3)<=tha*y(4)
    dy(5) = 1;
    
        if y(3)>0
            ra = ra1*cost;
            d1a = d;
            d2a = 0;
        end
        
        if y(3)<=0
            
            if y(5)<ti
                ra = 0;
                d1a = 0;
                d2a = 0;
            end
            
            if y(5)>=ti
                
                if y(4)>0
                    ra = ra2;
                    d1a = 0;
                    d2a = d;
                end
                if y(4)<=0
                    ra=0;
                d1a = 0;
                d2a = 0;
                end
            end
        end
end

     

% H
if y(3)>thh*y(4)
    rh = rh1;
    d1h = d;
    d2h = 0;
end
if y(3)<=thh*y(4)
    dy(6) = 1;
    
        if y(3)>0
            rh = rh1*cost;
            d1h = d;
            d2h = 0;
        end
        
        if y(3)<=0
            
            if y(6)<ti
                rh = 0;
                d1h = 0;
                d2h = 0;
            end
            
            if y(6)>=ti
                
                if y(4)>0
                    rh = rh2;
                    d1h = 0;
                    d2h = d;
                end
                if y(4)<=0
                    rh=0;
                d1h = 0;
                d2h = 0;
                end
            end
        end
end

                       
           
    dy(1) = y(1)*ra;
    dy(2) = y(2)*rh;
    dy(3) = -y(1)*d1a - y(2)*d1h;
    dy(4) = -y(1)*d2a - y(2)*d2h;
    
    


