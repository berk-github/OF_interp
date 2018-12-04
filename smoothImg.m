function smoothIM = smoothImg(im1,sm)

% sigma = 1;  % standard deviation for the gaussian filter
sigma = sm;
kernelSize = 2*(sigma*3);

xpts = -(kernelSize/2):(1+1/kernelSize):(kernelSize/2);
gaussFilter = (1/(sqrt(2*pi)*sigma)) * exp (-(xpts.^2)/(2*sigma^2));

% use 1D kernel twice for 2D smoothing
smoothIM = conv2(conv2(im1,gaussFilter,'same'),gaussFilter','same');
