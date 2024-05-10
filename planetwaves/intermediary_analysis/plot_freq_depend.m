function plot_freq_depend(my_array,variable_name,depth,frequencies,model)
% plot the frequency dependence of a variable in the deepest and shallowest section of the grid

    [deep_spot_depth, deep_spot_li] = max(depth(:));
    [deep_ri,deep_ci] = ind2sub(size(depth),deep_spot_li);

    [shallow_spot_depth, shallow_spot_li] = min(depth(:));
    [shallow_ri,shallow_ci] = ind2sub(size(depth),shallow_spot_li);



    figure;
    plot(frequencies,squeeze(my_array(deep_ri,deep_ci,:,round(model.Dirdim/2))),'LineWidth',3)
    hold on
    plot(frequencies,squeeze(my_array(shallow_ri,shallow_ci,:,round(model.Dirdim/2))),'LineWidth',3)
    xline(frequencies(model.cutoff_freq),'-','cutoff frequency')
    legend({strcat('Deep=',num2str(deep_spot_depth)), strcat('Shallow=',num2str(shallow_spot_depth))},'Location', 'Best')
    xlabel('frequency [Hz]')
    ylabel(variable_name)
    set(gca, 'XScale', 'log')
    hold off;
    

end