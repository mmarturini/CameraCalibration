function [imageData, K] = ZhangCalibration(imageData, iimage)
% Function to estimate P using Zhang calibration procedure

    for ii = 1:length(iimage)
        XYpixel=imageData(ii).XYpixel;
        XYmm=imageData(ii).XYmm;
        A=[];
        b=[];
    
        for jj=1:length(XYpixel)
    
            Xpixel = XYpixel(jj,1); % u
            Ypixel = XYpixel(jj,2); % v
            Xmm = XYmm(jj,1);       % x
            Ymm = XYmm(jj,2);       % y
    
            % building homogeneous coordinates
            m = [Xmm; Ymm; 1];
            O = [0;0;0];
            % A is build upon previous values of A, starting from an empty matrix
            A = [A;...
                 m' O' -Xpixel*m'; ...
                 O' m' -Ypixel*m'];
            b = [b;...
                 0;
                 0];
        end
        [~,~,V] = svd(A);
        h=V(:,end); % last svd vector associated to V
    
        imageData(ii).H = reshape(h, [3 3])';    
    end
    
    V = [];
    
    for ii=1:length(iimage)
        H = imageData(ii).H; 
    
        v12 = [H(1,1)*H(1,2); ...
               H(1,1)*H(2,2) + H(2,1)*H(1,2); ... 
               H(2,1)*H(2,2); ...
               H(3,1)*H(1,2) + H(1,1)*H(3,2); ...
               H(3,1)*H(2,2) + H(2,1)*H(3,2); ...
               H(3,1)*H(3,2)];
    
        v11 = [H(1,1)*H(1,1); ...
               H(1,1)*H(2,1) + H(2,1)*H(1,1); ... 
               H(2,1)*H(2,1); ...
               H(3,1)*H(1,1) + H(1,1)*H(3,1); ... 
               H(3,1)*H(2,1) + H(2,1)*H(3,1); ...
               H(3,1)*H(3,1)];
    
        v22 = [H(1,2)*H(1,2); ...
               H(1,2)*H(2,2) + H(2,2)*H(1,2); ... 
               H(2,2)*H(2,2); ...
               H(3,2)*H(1,2) + H(1,2)*H(3,2); ... 
               H(3,2)*H(2,2) + H(2,2)*H(3,2); ...
               H(3,2)*H(3,2)];
    
        V = [V; ...
             v12'; ...
             (v11-v22)'];
    end
    [~,~,S] = svd(V);
    b = S(:,end);
    
    B = [b(1) b(2) b(4); ...
         b(2) b(3) b(5); ...
         b(4) b(5) b(6)];
    
    
    d = eig(B);
    isposdef = all(d > 0);
    if isposdef == false 
        B = -B;
    end
    
    R = chol(B);  % by default matlab returns the upper triangular
    L = R';       % L is lower triangular
    K = inv(L');  % given Cholesky factorization, the calibration matrix is just the inverse of the transpose
    K = K/K(3,3); 
    
    for ii=1:length(iimage)
        H = imageData(ii).H;  %h1, h2 and h3 are the columns
        h1  = H(:,1);
        h2 = H(:,2);
        h3 = H(:,3);
        lambda = 1/norm(K\h1); % inv(K)*h1 same as K\h1 but faster and more accurate, x = A\B solves the system of linear equations A*x = B
        r1 = lambda*(K\h1);
        r2 = lambda*(K\h2);
        r3 = cross(r1,r2);
        t = lambda*(K\h3);
        P = K*[r1 r2 r3 t];
        imageData(ii).R = [r1 r2 r3];   % extrinsics R and t
        imageData(ii).t = t;            
        imageData(ii).P = P;            % perspective projection matrix P
    end

end
