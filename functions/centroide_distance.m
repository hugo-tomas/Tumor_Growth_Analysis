function [r,c,distance,tetas]=centroide_distance(dim,C,xymax,shape,radius,options_flowers)
    % Distance between the centroide and the (i,j)
    distance=zeros(dim); 
    tetas=zeros(dim);
    for i=1:dim
        for j=1:dim
            distance(i,j)=sqrt(((i-C)*(xymax/dim))^2+((j-C)*(xymax/dim))^2);
            tetas(i,j)=atan((i-C)/(j-C));
        end
    end
    tetas(C,C)=0;
    % 1 square --- dxy=xymax/dim (cm)
    %    z     ---  radius (cm)
    if strcmp(shape,"Spheric")
        % Find indexs with lower distance to the centroid than radius
        [r,c]=find(distance<=radius);
    elseif strcmp(shape,"Random")
        % Find indexs with lower distance to the centroid than
        % radius*abs(sin(rand(size(distance))))
        [r,c]=find(distance<=radius*abs(cos(tetas*rand(size(tetas)))));
    elseif strcmp(shape,"Flowers")
%         %[r,c]=find(distance<=radius*exp(-(distance-(radius+1*cos(6*tetas))).^2));
        [r,c]=find(distance<=radius+options_flowers(1)*cos(options_flowers(2)*tetas));
    else
        fprintf("Initial Tumor Shape not defined!");
        warndlg("Initial Tumor Shape not defined!","ERROR");
    end

end