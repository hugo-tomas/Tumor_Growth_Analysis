function plot_u2D(u2D,app)
    % Plot of tumor growth
    imagesc(u2D,'Interpolation','bilinear','parent',app);
end