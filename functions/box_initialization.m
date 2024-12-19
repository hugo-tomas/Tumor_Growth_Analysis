function box_initialization(dim,app_graph1,app_graph2,app_graph3)
    % Start first and second graph with zeros matrix
    imagesc(zeros(dim),'Interpolation','bilinear','parent',app_graph1);
    colorbar(app_graph1,'northoutside')
    imagesc(zeros(dim),'Interpolation','bilinear','parent',app_graph2);
    colorbar(app_graph2,'northoutside')
    imagesc(ones(dim),'Interpolation','bilinear','parent',app_graph3);
    colorbar(app_graph3,'northoutside')
    clim(app_graph3,[0,1]);
end

