function plot_prolMean(prol_mean,app)
    % Plot proliferantion mean
    max_val=max(max(prol_mean));
    min_val=min(min(prol_mean));
    clim(app,[min_val,max_val]);
    imagesc(prol_mean,'Interpolation','bilinear','parent',app);
end