function [new_dim,C]=correction_dimensions(dim) 
    % In order to have a well define centroide, dim must be odd!
    if mod(dim,2) == 0
         C=dim/2+1;
         new_dim=dim+1;
    else
        C=dim/2+0.5;
        new_dim=dim;
    end
end