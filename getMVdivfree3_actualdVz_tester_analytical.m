function [mvrow,mvcol,errDiv] = getMVdivfree3_actualdVz_tester_analytical(I, Vx, Vy, Vz, ite, gamma, alpha, mvCoeff, spacVec)

% I_n0 = I(:,:,1);
I_n1 = I(:,:,2);
I_n2 = I(:,:,3);
% I_n3 = I(:,:,4);

% Vx0 = Vx(:,:,1);
Vx1 = Vx(:,:,2);
Vx2 = Vx(:,:,3);
% Vx3 = Vx(:,:,4);

% Vy0 = Vy(:,:,1);
Vy1 = Vy(:,:,2);
Vy2 = Vy(:,:,3);
% Vy3 = Vy(:,:,4);

% Vz0 = Vz(:,:,1);
% Vz1 = Vz(:,:,2);
% Vz1_down = Vz(:,:,3);
% Vz2_up = Vz(:,:,4);
% Vz2 = Vz(:,:,5);
% Vz3 = Vz(:,:,6);

% Vz1 = Vz(:,:,1);
% % Vz1_down = Vz(:,:,2);
% % Vz2_up = Vz(:,:,3);
% Vz2 = Vz(:,:,4);
Vz1 = Vz(:,:,2);
Vz2 = Vz(:,:,3);

% I_n1 = sqrt(Vx1.^2 + Vy1.^2 + Vz1.^2);
% ind1 = find( I_n1 < 0 );
% I_n1( ind1 )=0;

%%
sm = 8;
I_n1 = smoothImg(I_n1,sm);
I_n2 = smoothImg(I_n2,sm);
Vx1 = smoothImg(Vx1,sm);
Vy1 = smoothImg(Vy1,sm);
Vz1 = smoothImg(Vz1,sm);
Vx2 = smoothImg(Vx2,sm);
Vy2 = smoothImg(Vy2,sm);
Vz2 = smoothImg(Vz2,sm);

%%
% Iz = conv2(I_n2, (1/9)*ones(3),'same') - conv2(I_n1, (1/9)*ones(3),'same');
% Iz = conv2(I_n3, (1/9)*ones(3),'same') - conv2(I_n0, (1/9)*ones(3),'same');
% Iz = conv2(I_n2, 0.25*ones(2), 'same') + conv2(I_n1, -0.25*ones(2), 'same');

% Ix(:,1)=0;
% Ix(:,end)=0;
% Iy(1,:)=0;
% Iy(end,:)=0;
% Iz(:,1)=0;
% Iz(:,end)=0;
% Iz(1,:)=0;
% Iz(end,:)=0;
simoncelli=0;

if simoncelli==1
    
    % Initial first term
    % Simoncelli matched balanced filters
    p  = [0.030320  0.249724  0.439911  0.249724  0.030320];
    d1 = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
    d2 = [0.232905  0.002668 -0.471147  0.002668  0.232905];
    Ix = conv2(p, d1, I_n1, 'same') + conv2(p, d1, I_n2, 'same');
    Iy = conv2(d1, p, I_n1, 'same') + conv2(d1, p, I_n2, 'same');
    Iz = .5*(I_n2 - I_n1)/spacVec(3);
    
    % Defining the first term
    dVx_dx_n1 = conv2(p, d1, Vx1, 'same')/spacVec(1);
    dVy_dy_n1 = conv2(d1, p, Vy1, 'same')/spacVec(2);
    dVx_dx_n2 = conv2(p, d1, Vx2, 'same')/spacVec(1);
    dVy_dy_n2 = conv2(d1, p, Vy2, 'same')/spacVec(2);
    dVz_dz = 1*(Vz2 - Vz1)/spacVec(3);
    
    Term1 = dVx_dx_n1 + dVy_dy_n1 + dVz_dz + dVx_dx_n2 + dVy_dy_n2 + dVz_dz;
    
    % Defining the second term
    d2Vx_dx2_n1 = conv2(p, d1, dVx_dx_n1, 'same')/spacVec(1);
    
    dVy_dx_n1 = conv2(p, d1, Vy1, 'same')/spacVec(1);
    d2Vy_dydx_n1 = conv2(d1, p, dVy_dx_n1, 'same')/spacVec(2);
    
    d2Vx_dx2_n2 = conv2(p, d1, dVx_dx_n2, 'same')/spacVec(1);
    
    dVy_dx_n2 = conv2(p, d1, Vy2, 'same')/spacVec(1);
    d2Vy_dydx_n2 = conv2(d1, p, dVy_dx_n2, 'same')/spacVec(2);
    
    Term2 = d2Vx_dx2_n1 + d2Vy_dydx_n1 - d2Vx_dx2_n2 - d2Vy_dydx_n2;
    
    % Defining Term 3
    dVx_dy_n1 = conv2(d1, p, Vx1, 'same')/spacVec(2);
    d2Vx_dxdy_n1 = conv2(p, d1, dVx_dy_n1, 'same')/spacVec(1);
    
    d2Vy_dy2_n1 = conv2(d1, p, dVy_dy_n1, 'same')/spacVec(2);
    
    dVx_dy_n2 = conv2(d1, p, Vx2, 'same')/spacVec(2);
    d2Vx_dxdy_n2 = conv2(p, d1, dVx_dy_n2, 'same')/spacVec(1);
    
    d2Vy_dy2_n2 = conv2(d1, p, dVy_dy_n2, 'same')/spacVec(2);
    
    Term3 = d2Vx_dxdy_n1 + d2Vy_dy2_n1 - d2Vx_dxdy_n2 - d2Vy_dy2_n2;

else 
    
%     kernel_x = .5*[1 0 -1];
%     kernel_y = .5*[-1 0 1]';
    kernel_x = .25*[1 0 -1; 1 0 -1];
%     kernel_x = .25*[-1 0 1; -1 0 1];
    kernel_y = .25*[-1 0 1; -1 0 1]';
%     kernel_y = kernel_x';
%     kernel_x = (1/6)*[1 0 -1; 1 0 -1; 1 0 -1];
%     kernel_y = (1/6)*[-1 0 1; -1 0 1; -1 0 1]';
%     kernel_x = [1/8 0 -1/8; 1/4 0 -1/4; 1/8 0 -1/8];
%     kernel_y = [-1/8 0 1/8; -1/4 0 1/4; -1/8 0 1/8]';
    
    kernel_x = (0.25)*[-1 1; -1 1];
%     kernel_x = .25*[1 -1; 1 -1];
%     kernel_y = kernel_x';
    kernel_y = (0.25)*[-1 -1; 1 1];
%     kernel_y = .25*[1 1; -1 -1];

    % Initial first term
    C=1;
    Ix = C*(conv2(I_n1, kernel_x, 'same')/spacVec(1) + conv2(I_n2, kernel_x, 'same')/spacVec(1));
    Iy = C*(conv2(I_n1, kernel_y, 'same')/spacVec(2) + conv2(I_n2, kernel_y, 'same')/spacVec(2));
    Iz = C*(.5*(I_n2 - I_n1)/spacVec(3));
%     Iz = conv2(I_n1, 0.25*ones(2),'same') + conv2(I_n2, -0.25*ones(2),'same');
    
    % Defining the first term
    dVx_dx_n1 = conv2(Vx1, kernel_x, 'same')/spacVec(1);
    dVy_dy_n1 = conv2(Vy1, kernel_y, 'same')/spacVec(2);
    dVx_dx_n2 = conv2(Vx2, kernel_x, 'same')/spacVec(1);
    dVy_dy_n2 = conv2(Vy2, kernel_y, 'same')/spacVec(2);
    dVz_dz = 1.*(Vz2 - Vz1)/spacVec(3);
%     dVz_dz = conv2(Vz1, 0.25*ones(2),'same') + conv2(Vz2, -0.25*ones(2),'same');
    
    Term1 = dVx_dx_n1 + dVy_dy_n1 + dVz_dz + dVx_dx_n2 + dVy_dy_n2 + dVz_dz;
    
    % Defining the second term
    d2Vx_dx2_n1 = conv2(dVx_dx_n1, kernel_x, 'same')/spacVec(1);
    
%     dVy_dx_n1 = conv2(Vy1, kernel_x, 'same')/spacVec(1);
%     d2Vy_dydx_n1 = conv2(dVy_dx_n1, kernel_y, 'same')/spacVec(2);
    d2Vy_dydx_n1 = conv2(dVy_dy_n1, kernel_x, 'same');

    d2Vx_dx2_n2 = conv2(dVx_dx_n2, kernel_x, 'same')/spacVec(1);
    
%     dVy_dx_n2 = conv2(Vy2, kernel_x, 'same')/spacVec(1);
%     d2Vy_dydx_n2 = conv2(dVy_dx_n2, kernel_y, 'same')/spacVec(2);
    d2Vy_dydx_n2 = conv2(dVy_dy_n2, kernel_x, 'same');
    
    Term2 = d2Vx_dx2_n1 + d2Vy_dydx_n1 - d2Vx_dx2_n2 - d2Vy_dydx_n2;
    
    % Defining Term 3
%     dVx_dy_n1 = conv2(Vx1, kernel_y, 'same')/spacVec(2);
%     d2Vx_dxdy_n1 = conv2(dVx_dy_n1, kernel_x, 'same')/spacVec(1);
    d2Vx_dxdy_n1 = conv2(dVx_dx_n1, kernel_y, 'same');

    d2Vy_dy2_n1 = conv2(dVy_dy_n1, kernel_y, 'same')/spacVec(2);
    
%     dVx_dy_n2 = conv2(Vx2, kernel_y, 'same')/spacVec(2);
%     d2Vx_dxdy_n2 = conv2(dVx_dy_n2, kernel_x, 'same')/spacVec(1);
    d2Vx_dxdy_n2 = conv2(dVx_dx_n2, kernel_y, 'same');
    
    d2Vy_dy2_n2 = conv2(dVy_dy_n2, kernel_y, 'same')/spacVec(2);
    
    Term3 = d2Vx_dxdy_n1 + d2Vy_dy2_n1 - d2Vx_dxdy_n2 - d2Vy_dy2_n2;
    
end

% sum(sum(Term1))
% sum(sum(Term2))
% sum(sum(Term3))
% sum(sum(Iz))
% sum(sum(Ix))
% sum(sum(Iy))
%%
denomDet = ((gamma.^2) .* (( ( Ix .* Term3 ) - ( Iy .* Term2 ) ).^2)) + (alpha .* ( Ix.^2 + Iy.^2 + (gamma.^2) .*( Term2.^2 + Term3.^2 ) + alpha ));

A = alpha .* ( ( Ix .* Iz ) + (gamma.^2) .* ( Term1 .* Term2 ) ) + (gamma.^2) .* ( ( Ix .* Iz .* Term3.^2 ) + ...
    ( Iy.^2 .* Term1 .* Term2 ) - ( Iy .* Iz .* Term2 .* Term3 ) - ( Ix .* Iy .* Term1 .* Term3 ) );

B = alpha .* ( ( Iy .* Iz ) + (gamma.^2) .* ( Term1 .* Term3 ) ) + (gamma.^2) .* ( ( Iy .* Iz .* Term2.^2 ) + ...
    ( Ix.^2 .* Term1 .* Term3 ) - ( Ix .* Iz .* Term2 .* Term3 ) - ( Ix .* Iy .* Term1 .* Term2 ) );

coefu1 = (gamma.^2) * ( ( Ix .* Term3 ) - ( Iy .* Term2 ) ).^2 + alpha * ( Ix.^2 + (gamma.^2)*(Term2.^2) );
coefu2 = alpha * ( ( Ix .* Iy ) + (gamma.^2)*( Term2 .* Term3 ) );
coefv1 = coefu2;    % coefv1 = alpha * ( ( Ix .* Iy ) + (gamma.^2)*( Term2 .* Term3 ) );
coefv2 = (gamma.^2) * ( ( Ix .* Term3 ) - ( Iy .* Term2 ) ).^2 + alpha * ( Iy.^2 + (gamma.^2)*(Term3.^2) );

% Averaging kernel
kernel_ave=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
% kernel_ave=[1/8 1/8 1/8;1/8 0 1/8;1/8 1/8 1/8];

uInitial = zeros(size(I_n1));
vInitial = zeros(size(I_n1));
% uInitial = Vx2-Vx1;
% vInitial = Vy2-Vy1;

u = uInitial;
v = vInitial;

% Iterations
for i=1:ite

    uAvg = conv2(u,kernel_ave,'same');
    vAvg = conv2(v,kernel_ave,'same');
	
    u = uAvg - (( ( coefu1 .* uAvg ) + ( coefv1 .* vAvg ) + A )  ./ denomDet);
    v = vAvg - (( ( coefu2 .* uAvg ) + ( coefv2 .* vAvg ) + B )  ./ denomDet);

end

u(isnan(u))=0;
v(isnan(v))=0;
% u(isinf(u))=0;
% v(isinf(v))=0;

%%
% mvCoeff=.1;   %10;
if mvCoeff == 0
    mvCoeffu=max(u(:));
    mvCoeffv=max(v(:));
else
    mvCoeffu=mvCoeff;
    mvCoeffv=mvCoeff;
end
mvcol=u/mvCoeffu;  %best result
mvrow=v/mvCoeffv;

mvcol(isnan(mvcol))=0;
mvrow(isnan(mvrow))=0;
mvcol(isinf(mvcol))=0;
mvrow(isinf(mvrow))=0;

eDiv = Term1 + u.*Term2 + v.*Term3;
errDiv = sum(sum(eDiv.^2));

%%
recur=0;
if recur==1
    uInitial = mvcol;
    vInitial = mvrow;
    u = uInitial;
    v = vInitial;
    
    % Iterations
    for i=1:ite
        
        uAvg = conv2(u,kernel_ave,'same');
        vAvg = conv2(v,kernel_ave,'same');
        
        u = uAvg - ( ( coefu1 .* uAvg ) + ( coefv1 .* vAvg ) + A )  ./ denomDet;
        v = vAvg - ( ( coefu2 .* uAvg ) + ( coefv2 .* vAvg ) + B )  ./ denomDet;
        
    end
    
    u(isnan(u))=0;
    v(isnan(v))=0;
    u(isinf(u))=0;
    v(isinf(v))=0;
    
    % mvCoeff=.1;   %10;
    if mvCoeff == 0
        mvCoeffu=max(u(:));
        mvCoeffv=max(v(:));
    else
        mvCoeffu=mvCoeff;
        mvCoeffv=mvCoeff;
    end
    mvcol=u/mvCoeffu;  %best result
    mvrow=v/mvCoeffv;
    
    mvcol(isnan(mvcol))=0;
    mvrow(isnan(mvrow))=0;
    mvcol(isinf(mvcol))=0;
    mvrow(isinf(mvrow))=0;
    
    eDiv = Term1 + u.*Term2 + v.*Term3;
    errDiv = sum(sum(eDiv.^2));
end

