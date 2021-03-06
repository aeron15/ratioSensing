%% Go over mulitplate raw data and produce a 24*24 matrix

function [E_area,E_prec,E_mean] = Plates2matMch(plates,data,plates_hists,d,map,th_const,off_peak)



for p = 1:length(plates)
    names = fieldnames( plates_hists.(plates{p}));
    xx = struct2cell(data{p});
    for i = 1:8
        for j =1:12
            n = i+d(p,1);
            m = j+d(p,2);
            
        if p==1;
            y_off_mCh = xx{1}.mCh_yfp.f;
            x_off_mCh = xx{1}.mCh_yfp.xi;

        end
        
            if~isempty(find(strcmp(names,map{i,j})))
                
                k = find(strcmp(names,map{i,j}));
                
                if length(xx{k}.mCh_yfp.f)>1
                    [prec_const,perc_adaptive,perc_area_temp] = AnalyzeHist_YS(xx{k}.mCh_yfp.f,xx{k}.mCh_yfp.xi,th_const,off_peak,y_off_mCh,x_off_mCh);
                    
                    E_area(n,m) = perc_area_temp;
                    E_prec(n,m) = prec_const;           
                    E_mean(n,m) = xx{k}.mCh_yfp.mean;

                else
                    E_area(n,m) = nan;
                    E_prec(n,m) = xx{k}.mCh_yfp.perc_ind;
                    E_mean(n,m) = xx{k}.mCh_yfp.mean;
                end
                
                
            else
                E_area(n,m)  = nan;
                E_prec(n,m) = nan;
                E_mean(n,m) = nan;

            end
        end

    end
end

