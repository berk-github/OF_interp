function [mvcol, mvrow] = myHornSchunck_analytical(im1, im2, alpha, iter, mvcolInitial, mvrowInitial)

% Horn-Schunck optical flow method 
% Horn, B.K.P., and Schunck, B.G., Determining Optical Flow, AI(17), No.
% 1-3, August 1981, pp. 185-203 http://dspace.mit.edu/handle/1721.1/6337
%
% Usage:
% [mvcol, mvrow] = myHornSchunck(im1, im2, alpha, ite, uInitial, vInitial, displayFlow)
%
% Inputs:
% -im1,im2 : two consecutive images.
% -alpha : coefficient for the smoothness term. for smoother flow, use
%          bigger alpha.
% -iter : number of iterations.
% -mvcolInitial, mvrowInitial : initial values for the flow. by default
%                               they're zero.
%
% Outputs:
% mvcol : flow vectors across columns (horizontal x-direction). a matrix size of the input images
% mvrow : flow vectors across rows (vertical y-direction). a matrix size of the input images
% 
% Flow field is from image 1 to image 2. 
%
% You can display the flow field: 
% figure, imshow(im1,[])
% hold on
% quiver(mvcol, mvrow, 3, 'linewidth', 2);
% set(gca,'YDir','reverse');
% hold off
% Note: You may want to downsample the flow field for better visualization.

%% Default parameters
if nargin<1 || nargin<2
    disp('You need at least 2 images for input!')
end
if nargin<3
    alpha = 1;  % Default smoothing factor
end
if nargin<4
    iter = 100; % Default number of iterations
end
if nargin<5 || nargin<6
    mvcolInitial = zeros(size(im1(:,:,1)));
    mvrowInitial = zeros(size(im2(:,:,1)));
elseif size(mvcolInitial,1) ==0 || size(mvrowInitial,1)==0
    mvcolInitial = zeros(size(im1(:,:,1)));
    mvrowInitial = zeros(size(im2(:,:,1)));
end
 
%% Convert images to grayscale
if size(size(im1),2) == 3
    im1 = rgb2gray(im1);
end
if size(size(im2),2)==3
    im2 = rgb2gray(im2);
end
im1 = double(im1);
im2 = double(im2);

%% Smoothing images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 8;  % standard deviation for the gaussian filter
kernelSize = 2*(sigma*3);

xpts = -(kernelSize/2):(1+1/kernelSize):(kernelSize/2);
gaussFilter = (1/(sqrt(2*pi)*sigma)) * exp (-(xpts.^2)/(2*sigma^2));

% use 1D kernel twice for 2D smoothing
im1 = conv2(conv2(im1,gaussFilter,'same'),gaussFilter','same');
im2 = conv2(conv2(im2,gaussFilter,'same'),gaussFilter','same');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;

%%
% Set initial value for the flow vectors
mvcol = mvcolInitial;
mvrow = mvrowInitial;

% Estimate spatiotemporal derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Horn-Schunck original method
fx = conv2(im1, 0.25*[-1 1; -1 1],'same') + conv2(im2, 0.25*[-1 1; -1 1],'same');
fy = conv2(im1, 0.25*[-1 -1; 1 1],'same') + conv2(im2, 0.25*[-1 -1; 1 1],'same');
ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');

% using Simonchelli filters for the derivatives
% p = [0.037659  0.249153  0.426375  0.249153  0.037659];
% d1 = [0.109604  0.276691  0.000000 -0.276691 -0.109604];
% fx = 0.5*(conv2(p, d1, im1, 'same') + conv2(p, d1, im2, 'same'));
% fy = 0.5*(conv2(d1, p, im1, 'same') + conv2(d1, p, im2, 'same'));
% ft = conv2(im2, (1/9)*ones(3),'same') - conv2(im1, (1/9)*ones(3),'same');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Averaging kernel
ave_kernel=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];

% Iterations
for i=1:iter
    % Compute local averages of the flow vectors
    mvcolAvg=conv2(mvcol, ave_kernel ,'same');
    mvrowAvg=conv2(mvrow, ave_kernel ,'same');
    % Compute flow vectors constrained by its local average and the optical flow constraints
    mvcol = mvcolAvg - (fx .* ((fx .* mvcolAvg) + (fy .* mvrowAvg) + ft)) ./ (alpha^2 + fx.^2 + fy.^2); 
    mvrow = mvrowAvg - (fy .* ((fx .* mvcolAvg) + (fy .* mvrowAvg) + ft)) ./ (alpha^2 + fx.^2 + fy.^2);
end

mvcol(isnan(mvcol))=0;
mvrow(isnan(mvrow))=0;

% toc