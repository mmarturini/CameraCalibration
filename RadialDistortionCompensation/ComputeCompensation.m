function imageData =  ComputeCompensation(imageData, k, K, iimage)
% function that compute the compensated u and v and update the XYpixel
% to apply the iterative distortion

    u0 = K(1,3);
    v0 = K(2,3);
    alpha_u = K(1,1);
    alpha_v = K(2,2);

    for ii=1:length(iimage)
    
        XYmm = imageData(ii).XYmm;
        XYpixel = imageData(ii).XYpixel;
        P = imageData(ii).P;
        
        for jj=1:length(XYpixel)
            u_hat = XYpixel(jj,1); % u
            v_hat = XYpixel(jj,2); % v
            Xmm = XYmm(jj,1);      % x
            Ymm = XYmm(jj,2);      % y
            homog_coords = [Xmm; Ymm; 0; 1];
            proj = P * homog_coords;
            u = proj(1)/proj(3);
            v = proj(1)/proj(3);
        
            % normalized coordinates
            x_hat = (u_hat - u0)/alpha_u;
            y_hat = (v_hat - v0)/alpha_v;
            x0 = (u-u0)/alpha_u;
            y0 = (v-v0)/alpha_v; % initial guess
        
            options = optimoptions('fsolve','Display','none');
        
            fun = @(x)norm_nonlin_system(x, x_hat, y_hat, k);
            X0 = [x0 y0];
            x = fsolve(fun, X0, options);
        
            % unnormalize
            compensated_u = x(1)*alpha_u + u0;
            compensated_v = x(2)*alpha_v + v0;
        
            imageData(ii).XYpixel(jj,:) = [compensated_u,compensated_v]; % I store in XY pixel the compensated values of u and v
        end
    end
end
    

function F = norm_nonlin_system(x, x_hat, y_hat, k)
    F(1) = x_hat - x(1)*(1 + k(1)*(x(1)^2 + x(2)^2) + k(2)*(x(1)^4 + 2*x(1)^2 * x(2)^2 + x(2)^4) );
    F(2) = y_hat - x(2)*(1 + k(1)*(x(1)^2 + x(2)^2) + k(2)*(x(1)^4 + 2*x(1)^2 * x(2)^2 + x(2)^4) );
end
