function newIntMat = buildIntMV(im1, im2, mvrow1, mvcol1, mvrow2, mvcol2, numadd, intMethod)

newIntMat = squeeze(zeros([size(im1),numadd]));

for z = 1 : numadd
    
    fract = ( z / (numadd+1) );
    backfract = 1 - fract;
    
    if strcmp(intMethod,'linear')
        im_interp = fract * im1 + backfract * im2;
    elseif strcmp(intMethod,'acgi')
        [ row , col ] = size( im1 );
        [ c , r ] = meshgrid( 1 : col , 1 : row );
        
        % im1est = fast_interp_general( im2 , mvcol1 + c , mvrow1 + r ) ;
        im1est = interp2( im2 , mvcol1 + c , mvrow1 + r );
        im1samp = backfract * im1 + fract * im1est;
        
        % im2est = fast_interp_general( im1 , mvcol2 + c , mvrow2 + r ) ;
        im2est = interp2( im1 , mvcol2 + c , mvrow2 + r );
        im2samp = fract * im2 + backfract * im2est;
        
        im2_inter = griddata(mvcol2*backfract+c,mvrow2*backfract+r,im2samp,c,r);
        im1_inter = griddata(mvcol1*fract+c,mvrow1*fract+r,im1samp,c,r);
        im_interp = im2_inter*fract + im1_inter*backfract;   
    elseif strcmp(intMethod,'foo')
        [ row , col ] = size( im1 );
        [ c , r ] = meshgrid( 1 : col , 1 : row );
        
        % im1est = fast_interp_general( im2 , mvcol1 + c , mvrow1 + r ) ;
        % im1est = interp2( im2 , mvcol1 + c , mvrow1 + r );
        % im1samp = backfract * im1 + fract * im1est;
        im1samp = im1;
        
        % im2est = fast_interp_general( im1 , mvcol2 + c , mvrow2 + r ) ;
        % im2est = interp2( im1 , mvcol2 + c , mvrow2 + r );
        % im2samp = fract * im2 + backfract * im2est;
        im2samp = im2;
        
        im2_inter = griddata(mvcol2+c,mvrow2+r,im2samp,c,r);
        im1_inter = griddata(mvcol1+c,mvrow1+r,im1samp,c,r);
        im_interp = im2_inter*fract + im1_inter*backfract;
    elseif strcmp(intMethod,'single')
        [ row , col ] = size( im1 );
        [ c , r ] = meshgrid( 1 : col , 1 : row );
        backfract=1;
        fract=1;
%         im1est = interp2( im1 , mvcol1 + c , mvrow1 + r );
%         im2est = interp2( im2 , mvcol1 - c , mvrow1 - r );
        im1_inter = griddata(-mvcol1*fract+c,-mvrow1*fract+r,im1,c,r);
        im2_inter = griddata(mvcol1*fract+c,mvrow1*fract+r,im2,c,r);
        im_interp = .5*(im2_inter*fract + im1_inter*backfract);
%         im_interp = .25*(im2_inter + im1_inter + im1est + im2est);
    elseif strcmp(intMethod,'2foo')
        [ row , col ] = size( im1 );
        [ c , r ] = meshgrid( 1 : col , 1 : row );
        backfract=1;
        fract=1;
%         im1est = interp2( im2 , mvcol1 + c , mvrow1 + r );
%         im1samp = .5*(backfract * im1 + fract * im1est);
%         im2est = interp2( im1 , mvcol2 + c , mvrow2 + r );
%         im2samp = .5*(fract * im2 + backfract * im2est);
        im1samp = im1;
        im2samp = im2;
%         im2_inter = griddata(mvcol2*backfract+c,mvrow2*backfract+r,im2samp,c,r);
%         im1_inter = griddata(mvcol1*fract+c,mvrow1*fract+r,im1samp,c,r);
        im2_inter = griddata(-mvcol2*backfract+c,-mvrow2*backfract+r,im2samp,c,r);
        im1_inter = griddata(-mvcol1*fract+c,-mvrow1*fract+r,im1samp,c,r);
        im_interp = .5*(im2_inter*fract + im1_inter*backfract);
    elseif strcmp(intMethod,'3foo')
        [ row , col ] = size( im1 );
        [ c , r ] = meshgrid( 1 : col , 1 : row );
        backfract=1;
        fract=1;
        im1est = interp2( im2 , 2*mvcol1 + c , 2*mvrow1 + r );
        im1samp = .5*(backfract * im1 + fract * im1est);
        im2est = interp2( im1 , 2*mvcol2 + c , 2*mvrow2 + r );
        im2samp = .5*(fract * im2 + backfract * im2est);
%         im1samp = im1;
%         im2samp = im2;
        im2_inter = griddata(2*mvcol2*backfract+c,2*mvrow2*backfract+r,im2samp,c,r);
        im1_inter = griddata(2*mvcol1*fract+c,2*mvrow1*fract+r,im1samp,c,r);
%         im2_inter = griddata(-mvcol2*backfract+c,-mvrow2*backfract+r,im2samp,c,r);
%         im1_inter = griddata(-mvcol1*fract+c,-mvrow1*fract+r,im1samp,c,r);
        im_interp = .5*(im2_inter*fract + im1_inter*backfract);
    elseif strcmp(intMethod,'4foo')
        [ row , col ] = size( im1 );
        [ c , r ] = meshgrid( 1 : col , 1 : row );
        
        estIm1_1 = interp2( c, r, im1 , -mvcol1 + c , -mvrow1 + r );
        estIm2_1 = interp2( c, r, im2 , mvcol1 + c , mvrow1 + r );
        estIm_1 = .5*(estIm1_1 + estIm2_1);
        
        estIm2_2 = interp2( c, r, im2 , -mvcol2 + c , -mvrow2 + r );
        estIm1_2 = interp2( c,r, im1 , mvcol2 + c , mvrow2 + r );
        estIm_2 = .5*(estIm1_2 + estIm2_2);
        
        im_interp = .5*(estIm_1 + estIm_2);
    else
        [ row , col ] = size( im1 );
        [ c , r ] = meshgrid( 1 : col , 1 : row );
        backfract=1;
        fract=1;
        im1est = interp2( im2 , mvcol1 + c , mvrow1 + r );
        im1samp = .5*(backfract * im1 + fract * im1est);
        im2est = interp2( im1 , mvcol2 + c , mvrow2 + r );
        im2samp = .5*(fract * im2 + backfract * im2est);
        im2_inter = griddata(-mvcol2*backfract+c,-mvrow2*backfract+r,im2samp,c,r);
        im1_inter = griddata(-mvcol1*fract+c,-mvrow1*fract+r,im1samp,c,r);
        im_interp = .5*(im2_inter*fract + im1_inter*backfract);
    end
    im_interp(isnan(im_interp))=0;
    newIntMat(:,:,z) = im_interp;
end
