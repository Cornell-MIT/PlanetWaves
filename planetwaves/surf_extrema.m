function surf_extrema(my_array,variable_name,model,extrema_type)
% make a surf plot of the max and min values of a variable of interest

    %my_array = squeeze(my_array(round(model.LonDim/2),round(model.LatDim/2),:,:));

    reshaped_array = reshape(my_array,model.LonDim*model.LatDim,model.Fdim*model.Dirdim);
    
    if strcmp(extrema_type,'tallest')
        val_ar = max(reshaped_array,[],2);
    elseif strcmp(extrema_type,'smallest')
        val_ar = min(reshaped_array,[],2);
    else
        error('SURF_EXTREMA: extrema type not specified')
    end

    valgrid = reshape(val_ar,model.LonDim,model.LatDim);
    

    figure;
    surf(valgrid,'FaceColor','interp')
    colorbar
    title([variable_name,' ',extrema_type,' at center of grid'])
    view(2)
    hold on

end