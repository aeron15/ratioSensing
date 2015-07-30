%% Go over mulitplate raw data and produce a 24*24 matrix

function [E_area,E_prec,E_mean] = Plates2mat2(plates,data,plates_hists,d,map,th_const,off_peak)



for p = 1:length(plates)
    names = fieldnames( plates_hists.(plates{p}));
    n=0;
    m=0;
    xx = struct2cell(data{p});
    for i = 1:8
        for j =1:12
            n = i+d(p,1);
            m = j+d(p,2);
            
        if p==1;
            y_off_bfp = xx{1}.bfp_yfp.f;
            x_off_bfp = xx{1}.bfp_yfp.xi;

        end
        
            if~isempty(find(strcmp(names,map{i,j})))
                
                k = find(strcmp(names,map{i,j}));
                
                if length(xx{k}.bfp_yfp.f)>1
                    [prec_const,perc_adaptive,perc_area_temp] = AnalyzeHist_YS(xx{k}.bfp_yfp.f,xx{k}.bfp_yfp.xi,th_const,off_peak,y_off_bfp,x_off_bfp);
                    
                    E_area{p}(n,m) = perc_area_temp;
                    E_prec{p}(n,m) = prec_const;           
                    E_mean{p}(n,m) = xx{k}.bfp_yfp.mean;

                else
                    E_area{p}(n,m) = nan;
                    E_prec{p}(n,m) = xx{k}.bfp_yfp.perc_ind;
                    E_mean{p}(n,m) = xx{k}.bfp_yfp.mean;
                end
                
                
            else
                E_area{p}(n,m)  = nan;
                E_prec{p}(n,m) = nan;
                E_mean{p}(n,m) = nan;

            end
        end

    end
end

