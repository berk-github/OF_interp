function myquiver(xi,yj,rSize,Vopt,varargin)

% rSize : region size.
% scale : vector length scale.

if nargin>6
    error('Too many input arguments!')
elseif nargin==6
    x_cor=varargin{1};
    y_cor=varargin{2};
    [row,col,~,~] = size(xi); 
    xMin=min(x_cor(:));
    xMax=max(x_cor(:));
    yMin=min(y_cor(:));
    yMax=max(y_cor(:));
    maxspac=2*(abs(xMax-xMin)/(col-1));
%     maxspac=abs(yMax-yMin)/(row-1);
else
    disp('Not enough grid info. Regular meshgrid will be used.')
    [row,col,~,~] = size(xi);     
    [x_cor,y_cor] = meshgrid( 1:col, 1:row ); 
    % this will be the maximum arrow length in one direction.
    % maxspac = 4 * abs(int_start-int_end)/(row-1);
    maxspac=rSize*2;
end

if length(size(xi))==3
    disp('First argument has 3 dimensions. selecting first slice for display.')
    xi=xi(:,:,1);
end
if length(size(yj))==3
    disp('Second argument has 3 dimensions. selecting first slice for display.')
    yj=yj(:,:,1);
end  

lineWidth=1;    %2
arrowSize=1.5;    %1.5

% u1 = xi;
% v1 = yj;

c_x = maxspac/max(abs(xi(:)));  % scale factor for the maximum allowable arrow length
c_y = maxspac/max(abs(yj(:)));

u1=c_x * xi; 
v1=c_y * yj;

rRange=rSize:rSize:row;
cRange=rSize:rSize:col;
    
if strcmp(Vopt,'simple')
%     % imshow(displayImg,[], 'Border', 'tight');
%     % hold on;
%     % rSize=3; %region size
%     % scale=5;
%     % Enhance the quiver plot visually by showing one vector per region
%     for i=1:size(x_cor,1)
%         for j=1:size(y_cor,2)
%             if floor(i/rSize)~=i/rSize || floor(j/rSize)~=j/rSize
%                 u1(i,j)=0;
%                 v1(i,j)=0;
%             end
%         end
%     end
%     % figure
%     quiver(x_cor, y_cor, u1, v1, 'color', 'b', 'linewidth', lineWidth);
%     set(gca,'YDir','reverse');
    cnti=1;
    for i=1:size(x_cor,1)
        cntj=1;
        for j=1:size(y_cor,2)
            
            if floor(i/rSize)~=i/rSize || floor(j/rSize)~=j/rSize
                u1(i,j)=0;
                v1(i,j)=0;
            else
                u2(cnti,cntj)=u1(i,j);
                v2(cnti,cntj)=v1(i,j);
                cntj=cntj+1;
            end
        end
        if floor(i/rSize)==i/rSize
            cnti=cnti+1;
        end
            
    end
    hq = quiver(x_cor(rRange,cRange), y_cor(rRange,cRange), u2, v2, 'Autoscale','off', 'linewidth',lineWidth);
%     hq = quiver(x_cor(rRange,cRange), y_cor(rRange,cRange), u1(rRange,cRange), v1(rRange,cRange),'AutoScale','off', 'color', 'b', 'linewidth', lineWidth);
    hq.MaxHeadSize = .8;
    hq.AutoScale = 'off';
%     adjust_quiver_arrowhead_size(hq,arrowSize);
    axis([0 inf 0 inf])
    axis off
else
    rRange=rSize:rSize:row;
    cRange=rSize:rSize:col;
    hq = quiver(x_cor(rRange,cRange), y_cor(rRange,cRange), u1(rRange,cRange), v1(rRange,cRange),'AutoScale','off', 'color', 'b', 'linewidth', lineWidth);
    hq.MaxHeadSize = .8;
%     adjust_quiver_arrowhead_size(hq,arrowSize);
    
    %get the line position (first handle)
%     hkid = get(hq,'children');
%     X = get(hkid(1),'XData');
%     Y = get(hkid(1),'YData');
    X = get(hq,'XData');
    Y = get(hq,'YData');
    
%     delete(hq);
    % %right version (with annotation)
    % hax_2 = subplot(1,2,2);
    cmap = jet(361); %colormap, 116 because angles goes up to 115 degrees
    
    for ii = 1:3:length(X)-1
        
        headWidth = 200 * sqrt((X(ii+1)-X(ii)).^2 + (Y(ii+1)-Y(ii)).^2); % set the headWidth, function of length of arrow
        angled = floor(atan2(Y(ii+1)-Y(ii),X(ii+1)-X(ii))*180/pi); %get the angle
        if angled<0
            angled=angled+360;
        end
        angled=angled+1;
        ah = annotation('arrow',...
            'Color', cmap(angled,:),'LineWidth',lineWidth,...
            'headStyle','vback2','HeadLength',headWidth,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[X(ii) Y(ii) X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
    end
    axis off;
    title('Quiver - annotations ','FontSize',16);
end