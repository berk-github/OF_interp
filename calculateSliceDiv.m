function [calcDiv, derivativeList] = calculateSliceDiv(Vx,Vy,Vz_up,Vz1,Vz_down,spacVec)


kernel_x = .5*[1 0 -1];
kernel_y = .5*[-1 0 1]';
kernel_y = .5*[1 0 -1]';

% kernel_x = (0.25)*[-1 1; -1 1];
% kernel_y = (0.25)*[-1 -1; 1 1];

%%%%%%%%%%%%%%%%%%%%%%
% Vx1=[zeros(size(Vx,1),2) Vx];
% Vx2=[Vx zeros(size(Vx,1),2)];
% temp = .5*(Vx2-Vx1);
% % dVx_dx = temp(2:end-1,3:end-2);
% dVx_dx = temp(:,2:end-1);
% 
% Vy1=[zeros(2,size(Vy,2)); Vy];
% Vy2=[Vy; zeros(2,size(Vy,2))];
% temp = .5*(Vy2-Vy1);
% % dVy_dy = temp(3:end-2,2:end-1);
% dVy_dy = temp(2:end-1,:);
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Vx1=[Vx(:,2:-1:1) Vx];
% Vx2=[Vx Vx(:,end:-1:end-1)];
% temp = .5*(Vx2-Vx1);
% dVx_dx = temp(:,2:end-1);
% dVx_dx = temp(2:end-1,3:end-2);

% Vy1=[Vy(2:-1:1,:); Vy];
% Vy2=[Vy; Vy(end:-1:end-1,:)];
% temp = .5*(Vy2-Vy1);
% dVy_dy = temp(2:end-1,:);
% dVy_dy = temp(3:end-2,2:end-1);
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% dVx_dx = conv2(Vx, [-1 0 1], 'same');  % *10e15
% dVx_dx = dVx_dx(2:end-1,2:end-1);
% dVy_dy = conv2(Vy, [-1 0 1]', 'same');  % *10e15
% dVy_dy = dVy_dy(2:end-1,2:end-1);
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% dVx_dx = .5*conv2(Vx, [-1 0 1], 'same');  % *10e15
% dVy_dy = .5*conv2(Vy, [-1 0 1]', 'same');  % *10e15
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
% temp = .5*(Vz2-Vz0);
% dVz_dz1 = temp(2:end-1,2:end-1);
% 
% temp = .5*(Vz0-Vz2);
% dVz_dz2 = temp(2:end-1,2:end-1);
%%%%%%%%%%%%%%%%%%%%%%

% dVz_dz1 = .05*(Vz2-Vz0);
% dVz_dz2 = .05*(Vz0-Vz2);
% dVz_dz3 = Vz1-Vz0;
% dVz_dz4 = Vz2-Vz1;
% dVz_dz5 = Vz1-Vz2;
% dVz_dz6 = Vz0-Vz1;

% div1 = dVx_dx + dVy_dy + dVz_dz1;
% div1 = sum(div1(:));
% div2 = dVx_dx + dVy_dy + dVz_dz2;
% div2 = sum(div2(:));
% div3 = dVx_dx + dVy_dy + dVz_dz3;
% div3 = sum(div3(:));
% div4 = dVx_dx + dVy_dy + dVz_dz4;
% div4 = sum(div4(:));
% div5 = dVx_dx + dVy_dy + dVz_dz5;
% div5 = sum(div5(:));
% div6 = dVx_dx + dVy_dy + dVz_dz6;
% div6 = sum(div6(:));

Vx=[Vx(:,2) Vx Vx(:,end-1)];
Vy=[Vy(2,:); Vy; Vy(end-1,:)];

dVx_dx = conv2(Vx, kernel_x, 'same')/spacVec(1);
dVy_dy = conv2(Vy, kernel_y, 'same')/spacVec(2);
dVz_dz = .5*(Vz_down-Vz_up)/spacVec(3);

% %%%%%%%%%%%%%%%%%%%%%%
% dVx_dx = dVx_dx(2:end-1,2:end-1);
% dVy_dy = dVy_dy(2:end-1,2:end-1);
% dVz_dz = dVz_dz(2:end-1,2:end-1);
% %%%%%%%%%%%%%%%%%%%%%%

dVx_dx = dVx_dx(:,2:end-1);
dVy_dy = dVy_dy(2:end-1,:);
% dVz_dz=[zeros(size(dVz_dz,1),2) dVz_dz zeros(size(dVz_dz,1),2)];
% dVz_dz=[zeros(size(2,dVz_dz),2); dVz_dz; zeros(2,size(dVz_dz,2))];

inOffset=2;
dVx_dx1 = dVx_dx(inOffset:end-inOffset+1,inOffset:end-inOffset+1);
dVy_dy1 = dVy_dy(inOffset:end-inOffset+1,inOffset:end-inOffset+1);
dVz_dz1 = dVz_dz(inOffset:end-inOffset+1,inOffset:end-inOffset+1);

div1 = dVx_dx1 + dVy_dy1 + dVz_dz1;
% div1 = sum(div1(:));
div1 = sum(abs(div1(:)));

% calcDiv = [div1 div2];
calcDiv = div1;
derivativeList = [sum(abs(dVx_dx(:))) sum(dVy_dy(:)) sum(abs(dVz_dz(:)))];




