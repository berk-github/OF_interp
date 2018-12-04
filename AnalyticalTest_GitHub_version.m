clear
%%%% Generate data - Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inRange=[-10 10];  
inSpac=.1; %.25;
zRange=[-20 20];
zSpac=8;
zStep=round(zSpac/inSpac);

%% New set
[x,y,z]=meshgrid( inRange(1):inSpac:inRange(2), inRange(1):inSpac:inRange(2), ...
    zRange(1):inSpac:zRange(2) ) ;  % define the interval and sampling rate

xi = .3 * (y.^2) + .15*x.^2;
yj = .3 * (1-x.^2) .* (y-1) - .3*y.*x; 
zk = .3 * -z .* (1-x.^2); 

xi_orig = xi;
yj_orig = yj;
zk_orig = zk;
temp1 = rand(size(xi));
temp1(temp1 < .7) = 0;
temp1 = temp1 * 20;
temp2 = rand(size(yj));
temp2(temp2 < .7) = 0;
temp2 = temp2 * 20;
temp3 = rand(size(zk));
temp3(temp3 < .7) = 0;
temp3 = temp3 * 10;
xi = xi + temp1;
yj = yj + temp2;
zk = zk + temp3;

syms x y z
divergence([x*sin(z),y*x,cos(z)-z*x],[x y z])

%%%% Choose a dataset - End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate data - End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%

imReg=20:(size(xi,2)-19);
m=length(imReg);
n=length(imReg);
% zSlice1=90;
% zSlice2=zSlice1+zStep*2;
% zSlice_original=zSlice1+zStep;
imRange=90;%40:160;
mseLst=zeros(length(imRange), 8);
divLst=zeros(length(imRange), 5);
mseLst2=zeros(length(imRange), 8);
divLst2=zeros(length(imRange), 5);
msedivLstLin=zeros(length(imRange), 5);
cnt=1;
for imNum=imRange
    
    zSlice1=imNum;
    zSlice2=zSlice1+zStep*2;
    zSlice_original=zSlice1+zStep;
    
    Vx0 = xi(:,:,zSlice1-1);
    Vx1 = xi(:,:,zSlice1);
    Vx2 = xi(:,:,zSlice2);
    Vx3 = xi(:,:,zSlice2+1);
    
    Vy0 = yj(:,:,zSlice1-1);
    Vy1 = yj(:,:,zSlice1);
    Vy2 = yj(:,:,zSlice2);
    Vy3 = yj(:,:,zSlice2+1);
    
    Vz0 = zk(:, :,zSlice1-1);
    Vz1 = zk(:, :,zSlice1);
    Vz2 = zk(:, :,zSlice2);
    Vz3 = zk(:, :,zSlice2+1);
    
    In0 = sqrt(Vx0.^2 + Vy0.^2 + Vz0.^2);
    In1 = sqrt(Vx1.^2 + Vy1.^2 + Vz1.^2);
    In2 = sqrt(Vx2.^2 + Vy2.^2 + Vz2.^2);
    In3 = sqrt(Vx3.^2 + Vy3.^2 + Vz3.^2);
    
    Vx_original = xi(:,:,zSlice_original);
    Vy_original = yj(:,:,zSlice_original);
    Vz_original = zk(:,:,zSlice_original);
    Vz_original_up = zk(:,:,zSlice_original-1);
    Vz_original_down = zk(:,:,zSlice_original+1);
    In_original = sqrt(Vx_original.^2 + Vy_original.^2 + Vz_original.^2);  
    
    spacVec = [inSpac, inSpac, zSpac];
    
    %%
    numadd=1;
    % interpolate and reconstruct the new image using linear interpolation.
    vLin1=0;
    uLin1=0;
    vLin2=0;
    uLin2=0;
    
    In_lin = buildIntMV(In1, In2, vLin1, uLin1, vLin2, uLin2, numadd, 'linear');
    Vx_lin = buildIntMV(Vx1, Vx2, vLin1, uLin1, vLin2, uLin2, numadd, 'linear');
    Vy_lin = buildIntMV(Vy1, Vy2, vLin1, uLin1, vLin2, uLin2, numadd, 'linear');
    Vz_lin = buildIntMV(Vz1, Vz2, vLin1, uLin1, vLin2, uLin2, numadd, 'linear');
    
    mseIn_lin = sum(sum((In_original(imReg,imReg) - In_lin(imReg,imReg)).^2))/(m*n);
    mseVx_lin = sum(sum((Vx_original(imReg,imReg) - Vx_lin(imReg,imReg)).^2))/(m*n);
    mseVy_lin = sum(sum((Vy_original(imReg,imReg) - Vy_lin(imReg,imReg)).^2))/(m*n);
    mseVz_lin = sum(sum((Vz_original(imReg,imReg) - Vz_lin(imReg,imReg)).^2))/(m*n);
    
    %%
    alpha1=1;
    ite1=200;
    uInitial = zeros(size(In1));
    vInitial = zeros(size(In2));
    
 
    [uHS1,vHS1]=myHornSchunck_analytical(In1,In2,alpha1,ite1,uInitial,vInitial);
    [uHS2,vHS2]=myHornSchunck_analytical(In2,In1,alpha1,ite1,uInitial,vInitial);
    
    % interpolate and reconstruct the new image using the motion vectors
    % u : column vectors
    % v : row vectors
    In_HS = buildIntMV(In1, In2, vHS1, uHS1, vHS2, uHS2, numadd, 'foo');
    Vx_HS = buildIntMV(Vx1, Vx2, vHS1, uHS1, vHS2, uHS2, numadd, 'foo');
    Vy_HS = buildIntMV(Vy1, Vy2, vHS1, uHS1, vHS2, uHS2, numadd, 'foo');
    Vz_HS = buildIntMV(Vz1, Vz2, vHS1, uHS1, vHS2, uHS2, numadd, 'foo');
    
    mseIn_HS = sum(sum((In_original(imReg,imReg) - In_HS(imReg,imReg)).^2))/(m*n);
    mseVx_HS = sum(sum((Vx_original(imReg,imReg) - Vx_HS(imReg,imReg)).^2))/(m*n);
    mseVy_HS = sum(sum((Vy_original(imReg,imReg) - Vy_HS(imReg,imReg)).^2))/(m*n);
    mseVz_HS = sum(sum((Vz_original(imReg,imReg) - Vz_HS(imReg,imReg)).^2))/(m*n);
    
    %%
    gammaS=1;
    alphaS=1;
    mvCoeffS=1;
    iterS=200;
    
    nRuns = length(gammaS)*length(alphaS)*length(mvCoeffS)*length(iterS);
    disp(['Total ' int2str(nRuns) ' runs.'])
    
    mseIn_DF = zeros(length(gammaS),length(alphaS),length(mvCoeffS),length(iterS));
    mseVx_DF = zeros(length(gammaS),length(alphaS),length(mvCoeffS),length(iterS));
    mseVy_DF = zeros(length(gammaS),length(alphaS),length(mvCoeffS),length(iterS));
    mseVz_DF = zeros(length(gammaS),length(alphaS),length(mvCoeffS),length(iterS));
    errDivDF = zeros(length(gammaS),length(alphaS),length(mvCoeffS),length(iterS));
    %sliceDiv_DF = zeros(length(gammaS),length(alphaS),length(mvCoeffS),length(iterS));
    sliceDiv2_DF = zeros(length(gammaS),length(alphaS),length(mvCoeffS),length(iterS));
    
    I_set(:,:,1)=In0;
    I_set(:,:,2)=In1;
    I_set(:,:,3)=In2;
    I_set(:,:,4)=In3;
    
    Vx_set(:,:,1)=Vx0;
    Vx_set(:,:,2)=Vx1;
    Vx_set(:,:,3)=Vx2;
    Vx_set(:,:,4)=Vx3;
    
    Vy_set(:,:,1)=Vy0;
    Vy_set(:,:,2)=Vy1;
    Vy_set(:,:,3)=Vy2;
    Vy_set(:,:,4)=Vy3;
    
    Vz_set(:,:,1)=Vz0;
    Vz_set(:,:,2)=Vz1;
    Vz_set(:,:,3)=Vz2;
    Vz_set(:,:,4)=Vz3;
    
    minMSE=[100 0 0 0 0 0 0 0];
    minDIV=[10000 0 0 0 0];
    minMSE2=[100 0 0 0 0 0 0 0];
    minDIV2=[10000 0 0 0 0];
    nR = 1;
    for ii = 1 : length(gammaS)
        for jj = 1 : length(alphaS)
            for kk = 1 : length(mvCoeffS)
                for ll = 1 : length(iterS)
                    
                    tic
                    gamma = gammaS(ii);
                    alpha2 = alphaS(jj);
                    mvCoeff = mvCoeffS(kk);
                    iter = iterS(ll);
                    disp(['Parameters-> gamma: ' num2str(gamma) ' alpha: ' num2str(alpha2) ' mvCoefficient: ' num2str(mvCoeff) ' iter: ' num2str(iter)])
                    [mvrowDF1,mvcolDF1,errDivDF(ii,jj,kk,ll)] = getMVdivfree3_actualdVz_tester_analytical(I_set, Vx_set, Vy_set, Vz_set, iter, gamma, alpha2,1,[1 1 16]);
                    [mvrowDF2,mvcolDF2,errDivDF(ii,jj,kk,ll)] = getMVdivfree3_actualdVz_tester_analytical(I_set(:,:,end:-1:1), Vx_set(:,:,end:-1:1), Vy_set(:,:,end:-1:1), Vz_set(:,:,end:-1:1), iter, gamma, alpha2,1,[1 1 16]);
                    
                    % interpolate and reconstruct the new image using the motion vectors
                    if (min(mvrowDF1(:))+max(mvrowDF1(:)))==0 || (min(mvcolDF1(:))+max(mvcolDF1(:)))==0
                        if (min(mvrowDF1(:))+max(mvrowDF1(:))+min(mvcolDF1(:))+max(mvcolDF1(:)))==0
                            disp('Flow vectors (both) are 0! Try different parameters.')
                        elseif (min(mvrowDF1(:))+max(mvrowDF1(:)))==0
                            disp('Flow vectors (row) are 0! Try different parameters.')
                        elseif (min(mvcolDF1(:))+max(mvcolDF1(:)))==0
                            disp('Flow vectors (column) are 0! Try different parameters.')
                        end
                        mseIn_DF(ii,jj,kk,ll)=500;
                        mseVx_DF(ii,jj,kk,ll)=500;
                        mseVy_DF(ii,jj,kk,ll)=500;
                        mseVz_DF(ii,jj,kk,ll)=500;
                        %sliceDiv_DF(ii,jj,kk,ll,:)=10000;
                        sliceDiv2_DF(ii,jj,kk,ll,:)=10000;
                        %                     continue
                    else
                        In_DF = buildIntMV(In1, In2, mvrowDF1, mvcolDF1, mvrowDF2, mvcolDF2, numadd, '4foo');
                        Vx_DF = buildIntMV(Vx1, Vx2, mvrowDF1, mvcolDF1, mvrowDF2, mvcolDF2, numadd, '4foo');
                        Vy_DF = buildIntMV(Vy1, Vy2, mvrowDF1, mvcolDF1, mvrowDF2, mvcolDF2, numadd, '4foo');
                        Vz_DF = buildIntMV(Vz1, Vz2, mvrowDF1, mvcolDF1, mvrowDF2, mvcolDF2, numadd, '4foo');
                        mseIn_DF(ii,jj,kk,ll) = sum(sum((In_original(imReg,imReg) - In_DF(imReg,imReg)).^2))/(m*n);
                        mseVx_DF(ii,jj,kk,ll) = sum(sum((Vx_original(imReg,imReg) - Vx_DF(imReg,imReg)).^2))/(m*n);
                        mseVy_DF(ii,jj,kk,ll) = sum(sum((Vy_original(imReg,imReg) - Vy_DF(imReg,imReg)).^2))/(m*n);
                        mseVz_DF(ii,jj,kk,ll) = sum(sum((Vz_original(imReg,imReg) - Vz_DF(imReg,imReg)).^2))/(m*n);
                    end                

                    disp(['Run ' int2str(nR) ' complete.'])
                    
                    nR = nR+1;
                    if mseIn_DF(ii,jj,kk,ll) < minMSE(1)
                        minMSE(1)=mseIn_DF(ii,jj,kk,ll);
                        minMSE(2)=mseVx_DF(ii,jj,kk,ll);
                        minMSE(3)=mseVy_DF(ii,jj,kk,ll);
                        minMSE(4)=mseVz_DF(ii,jj,kk,ll);
                        minMSE(5)=gamma;
                        minMSE(6)=alpha2;
                        minMSE(7)=mvCoeff;
                        minMSE(8)=iter;
                    end
                    if abs(sliceDiv2_DF(ii,jj,kk,ll)) < abs(minDIV(1))
                        minDIV(1)=sliceDiv2_DF(ii,jj,kk,ll);
                        minDIV(2)=gamma;
                        minDIV(3)=alpha2;
                        minDIV(4)=mvCoeff;
                        minDIV(5)=iter;
                    end
                    
                end
            end
        end
    end
    
    mseLst(cnt,:)=minMSE(:);
    divLst(cnt,:)=minDIV(:);
    mseLst2(cnt,:)=minMSE2(:);
    divLst2(cnt,:)=minDIV2(:);
    
    cnt=cnt+1;
    disp(num2str(imNum))
end


%%
Vx_original_NoNoise = xi_orig(:,:,zSlice_original);
Vy_original_NoNoise = yj_orig(:,:,zSlice_original);

figure, myquiver(Vx_original_NoNoise,Vy_original_NoNoise,5,'simple')
title('Original')
figure, myquiver(Vx_original,Vy_original,5,'simple')
title('Original with Gaussian Noise')
figure, myquiver(Vx_lin,Vy_lin,5,'simple')
title('Linear Interpolation')
figure, myquiver(Vx_HS,Vy_HS,5,'simple')
title('Horn-Schunck Based Interpolation')
figure, myquiver(Vx_DF,Vy_DF,5,'simple')
title('Proposed Method')

