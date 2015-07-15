function [x,y] = GetHistogram(plates,data,plates_hists,plate_name,well_name)


    p = find(strcmp(plates,plate_name));
    names = fieldnames( plates_hists.(plates{p}));
    k = find(strcmp(names,well_name));
    xx = struct2cell(data{p});

    y = xx{k}.bfp_yfp.f;
    x = xx{k}.bfp_yfp.xi;
