%% Loading images and getting correspondences

clear imageData

iimage=[4 7 10 16]; 

for ii=1:length(iimage)
    imageFileName = fullfile('images',['image' num2str(iimage(ii)) '.tif']);
    imageData(ii).I = imread(imageFileName);
    imageData(ii).XYpixel = detectCheckerboardPoints(imageData(ii).I);
end

squaresize = 30; %mm

for ii=1:length(iimage)

    XYpixel=imageData(ii).XYpixel;
    clear Xmm Ymm;
    for jj=1:length(XYpixel)
        [row,col] = ind2sub([12,13],jj); % linear index to row,col
        
        Xmm = (col-1)*squaresize;
        Ymm = (row-1)*squaresize;

        imageData(ii).XYmm(jj,:)=[Xmm,Ymm];
    end
end


%% Radial distortion Calibration 
% Trying to apply the calibration procedure compensating radial distortion.
% The code give an error at the second iteration of the while loop because
% in the Zhang calibration procedure the matrix B is not positive definite.

eps = 0.001;
diff_P = 10000;
diff_k = 10000;
old_k = [100; 100];
old_P = [1e6 1e6 1e6 1e6;
         1e6 1e6 1e6 1e6;
         1e6 1e6 1e6 1e6];

while (diff_P > eps && diff_k > eps)
    
    [imageData, K] = ZhangCalibration(imageData, iimage);
    k = kEstimation(imageData, K, iimage);
    imageData = ComputeCompensation(imageData, k, K, iimage);

    diff_k = norm(k-old_k);
    diff_P = norm(imageData(1).P - old_P);
    old_k = k;
    old_P = imageData(1).P;
end






