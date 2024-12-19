function [u_new,u2D,prol_mean,R_met,p_average,hoc,nhoc]=newSimulation(dim,u,Dc,Dp,miu,K,dxy,dp,dt,distance,u2d_options,gamma_options,shape_options,apply_delta,delta)
    % New state of u
    u_new=u;
    % u(x,y)
    u2D=zeros(dim);
    gamma=ones(dim);

    % p(x,y)
    prol_mean=zeros(dim);
    M=zeros(dim);
    
    M_num=0; M_den=0;
    avg_prol_num=0;
    avg_prol_den=0;
    
    for i=1:size(u,1) % x
        [a,b]=derivative_borders(dim,i);
        for j=1:size(u,2) % y
            [c,d]=derivative_borders(dim,j);

            u2D(i,j)=sum(u_new(i,j,:));
            mean_u2d=mean(u2d_options);
            if (~isempty(u2d_options) && ~isempty(gamma_options))
                if strcmp(shape_options,"Retangular")
                    if u2D(i,j)<=u2d_options(1)
                        gamma(i,j)=gamma_options(1);
                    elseif u2D(i,j)>u2d_options(1) && u2D(i,j)<=u2d_options(2)
                        gamma(i,j)=gamma_options(2);
                    else
                        gamma(i,j)=gamma_options(3);
                    end
                elseif strcmp(shape_options,"Hyperbolic Tangent")
                   gamma(i,j)=(u2D(i,j)<=mean_u2d)*(((gamma_options(2)-gamma_options(1))/2)*(tanh(u2D(i,j)-u2d_options(1))+1)+gamma_options(1))+(u2D(i,j)>mean_u2d)*(((gamma_options(3)-gamma_options(2))/2)*(tanh(u2D(i,j)-u2d_options(2))+1)+gamma_options(2));
                end
            end

            sum_prol=0; % Used in Proliferation Mean
            sum_M=0; % Used in Metabolic Radius
            for k=1:size(u,3) % p
                [e,f]=derivative_borders(dim,k);
                % u_new construction
                u_new(i,j,k)=u(i,j,k)+dt*Dc*(((u(a,j,k)+u(b,j,k)+u(i,c,k)+u(i,d,k)-4*u(i,j,k))/(dxy^2)));
                u_new(i,j,k)=u_new(i,j,k)+dt*Dp*((u(i,j,e)+u(i,j,f)-2*u(i,j,k))/(dp^2));
                if ((~isempty(u2d_options) && ~isempty(gamma_options))|| apply_delta==1)
                    u_new(i,j,k)=u_new(i,j,k)+dt*(delta(i,j)*gamma(i,j)*k*dp-miu)*(1-((1/K)*u2D(i,j)))*u(i,j,k);
                else
                     u_new(i,j,k)=u_new(i,j,k)+dt*(k*dp-miu)*(1-((1/K)*u2D(i,j)))*u(i,j,k);
                end
    
                % Sum of proliferation state*density of cells
                sum_prol=sum_prol+(u_new(i,j,k)*k*dp);
                % Parameters to calculate Proliferation Average (Base Article formula)
                avg_prol_num=avg_prol_num+(u_new(i,j,k)*k*dp)*((distance(i,j))^2);
                % Sum to calculate M(x,y)
                if (~isempty(u2d_options) && ~isempty(gamma_options))
                    sum_M=sum_M+(gamma(i,j)*k*dp-miu)*(1-((1/K)*u2D(i,j)))*u(i,j,k);
                else
                     sum_M=sum_M+(k*dp-miu)*(1-((1/K)*u2D(i,j)))*u(i,j,k);
                end
            end

            u2D(i,j)=sum(u_new(i,j,:)); % sum all u(i,j) to several proliferate taxes
            avg_prol_den=avg_prol_den+u2D(i,j)*((distance(i,j))^2); % Proliferation Average (Base Article formula)
    
            % Cell density treshold to consider prol_mean
            prol_mean(i,j)=sum_prol/u2D(i,j);
    
            M(i,j)=sum_M;
            % Parameters to obtein Metabolic Radius
            M_num=M_num+(sum_M*((distance(i,j))^3));
            M_den=M_den+(sum_M*((distance(i,j))^2));
        end
    end
    
    
    % Metabolic Radius
    R_met=M_num/M_den;
    
    % Proliferation Average
    p_average=avg_prol_num/avg_prol_den;
    
    % HOC
    % Index of maximum of proliferation
    if sum(sum(M))==0
        hoc=0;
    else
        [r_hotspot,c_hotspot]=find(M==max(max(M)));
        % Mininum distance with higher proliferation capacity
        distances=zeros(size(r_hotspot));
        for idx=1:length(r_hotspot)
            distances(idx)=distance(r_hotspot(idx),c_hotspot(idx));
        end
        hoc=mean(distances);
    end
    %NHOC
    nhoc=hoc/R_met;
end