function [a,b]=derivative_borders(dim,i)
    % Using null derivate borders
    % x=-1 -> x=0
    % x=lim_max+1 -> x=lim_max
    a=i+1; b=i-1;
    if a==dim+1
        a=dim;
    elseif b==0
        b=1;
    end
end