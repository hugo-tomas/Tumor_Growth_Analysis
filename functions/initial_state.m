function [dim,dxy,dp,u,distance]=initial_state(dim_input,xymax,pmax,shape_i,radius_i,cells_i,options_flowers)
    % Correction of initial dimension
    [dim,C]=correction_dimensions(dim_input);

    % 3D matrix: u(x,y,p)
    u=zeros(dim,dim,dim);
    
    % Spacial and proliferation step
    % xymax --- dim squares
    %  dxy  ---  1 square
    dxy=xymax/dim;
    dp=pmax/dim;
    
    % Introduce the initial state of u(x,y,p)
    % 1 square --- dp=pmax/dim (/days)
    %    z     ---  p_i (/days)
    [r,c,distance,tetas]=centroide_distance(dim,C,xymax,shape_i,radius_i,options_flowers);
    volume=(4/3)*pi()*radius_i^3;
    density=cells_i/volume;

    prol_mean=1.7e-2+3e-3*randn();
    
    % Define 2D Gaussian function
    if strcmp(shape_i,"Flowers")
        f = @(x,y,A,C,sigma) A * exp(-((x-C).^2+(y-C).^2)./(2*(sigma+options_flowers(1)*cos(options_flowers(2)*tetas)).^2));
    elseif strcmp(shape_i,"Spheric")
        f = @(x,y,A,C,sigma) A * exp(-(((x-C).^2)/(2*(sigma^2)) + ((y-C).^2)/(2*(sigma^2))));
    elseif strcmp(shape_i,"Random")
        f= @(x,y,A,C,sigma) A*rand(size(x))+((rand(size(x))*2*sigma)-sigma);
    end
    
    % Evaluate function over a meshgrid
    [x,y] = meshgrid(1:dim);
    densities=f(x,y,density,C,radius_i/2);

    % Set the minimum and maximum cell density values
    min_density = 0;
    max_density = 1e5;

    densities = (densities - min(densities(:))) ./ (max(densities(:)) - min(densities(:)));
    densities = densities .* (max_density - min_density) + min_density;

    densities(distance > radius_i) = 0;
    densities = densities.*(density / (sum(sum(densities))/sum(sum(densities~=0))));

    for i=1:length(r)
        u(r(i),c(i),round(prol_mean/dp))=densities(r(i),c(i));
    end
end