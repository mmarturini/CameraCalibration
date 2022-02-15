%% Estimating radial distortion parameters k1 and k2
% Assuming the points detected by the detectCheckerboardPoints() are the
% measured one on the image without taking into consideration the radial
% distortion compensation, hence are the distorted ones, cause they are the
% one taken on the distorted image

% K is just one, assuming the skew angle is irrelevant 
function k = kEstimation(imageData, K, iimage)
    u0 = K(1,3);
    v0 = K(2,3);
    alpha_u = K(1,1);
    alpha_v = K(2,2);
    
    A = [];
    b = [];
    for ii=1:length(iimage)
        XYpixel=imageData(ii).XYpixel;
        XYmm=imageData(ii).XYmm;
        P = imageData(ii).P;
    
        for jj=1:length(XYpixel)
            u_hat = XYpixel(jj,1);  % u
            v_hat = XYpixel(jj,2);  % v
            Xmm = XYmm(jj,1);       % x
            Ymm = XYmm(jj,2);       % y
    
            homog_coords = [Xmm; Ymm; 0; 1];
            proj = P * homog_coords;
            u = proj(1)/proj(3);
            v = proj(1)/proj(3);
    
            rd2 = ((u-u0)/alpha_u)^2 + ((v-v0)/alpha_v)^2;
    
            A = [A; ...
                 (u-u0)*rd2 (u-u0)*rd2^2; ...
                 (v-v0)*rd2 (v-v0)*rd2^2;];
            b = [b; ...
                 u_hat-u; ...
                 v_hat-v];
        end
    end
    
    k = inv((A')*A)*(A')*b;
    
end
