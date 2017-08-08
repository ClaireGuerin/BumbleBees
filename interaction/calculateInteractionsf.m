function fInteractions = calculateInteractions(tag,shape,expansion,im,vis)
% frm = 3560;
% 
% fImage = read(vid,frm);
% im = fImage;
% 
% xcenter = coordFixed(frm,:,1);
% ycenter = coordFixed(frm,:,2);
% xfront = coordFixed(frm,:,3);
% yfront = coordFixed(frm,:,4);
% 
% fHeight = size(fImage,1);
% fWidth = size(fImage,2);
% 
% expansion = 5*1.2;
% shape = [len;width];
% vis = 1;

xcenter = tag(1,:);
ycenter = tag(2,:);
xfront = tag(3,:);
yfront = tag(4,:);

fHeight = size(im,1);
fWidth = size(im,2);


% ELLIPSES PARAMETERS

a = expansion * shape(1,:)/2;
b = expansion * shape(2,:)/2;
angle = atan2(yfront - ycenter,xfront - xcenter);
%         x1 = xcenter + (1/2*a.*cos(angle));
%         x2 = xcenter + (3/2*a.*cos(angle + pi));
%         y1 = ycenter + (1/2*a.*sin(angle));
%         y2 = ycenter + (3/2*a.*sin(angle + pi));
x0 = xcenter - 1/2*a.*cos(angle);
y0 = ycenter - 1/2*a.*sin(angle);

t = linspace(0,2*pi,10);
X = a.*cos(t.');
Y = b.*sin(t.');
%         w = atan2(y1-y2,x1-x2);
%         x = (x1+x2)/2 + X.*cos(w) - Y.*sin(w);
%         y = (y1+y2)/2 + X.*sin(w) + Y.*cos(w);
%         x = x0 + X.*cos(w) - Y.*sin(w);
%         y = y0 + X.*sin(w) + Y.*cos(w);
x = x0 + X.*cos(angle) - Y.*sin(angle);
y = y0 + X.*sin(angle) + Y.*cos(angle);

bHere = find(any(x));
nPresentB = size(bHere,2);

nPop = size(xcenter,2);
fInteractions = zeros(nPop,nPop);

if vis == 1
    imshow(im)
    hold on
end

if nPresentB > 1
    BW = nan(fHeight,fWidth,nPresentB);
    
    for bee = 1:nPresentB
        BW(:,:,bee) = poly2mask(x(:,bHere(bee)), y(:,bHere(bee)), fHeight, fWidth);
    end
    
    combos = nchoosek(1:nPresentB,2);
    
    if vis == 1
        plot(x(:,bHere),y(:,bHere), 'y', 'LineWidth', 2)
    end
    
    for pair = 1:size(combos,1)
        bee1 = combos(pair,1);
        bee2 = combos(pair,2);
        interbee = any(any(BW(:,:,bee1)+BW(:,:,bee2)>1));
        
        intB1 = bHere(bee1);
        intB2 = bHere(bee2);
        fInteractions(intB1,intB2) = interbee;
        
        if vis == 1 && interbee
            plot(x(:,[bHere(bee1),bHere(bee2)]),y(:,[bHere(bee1),bHere(bee2)]),'r', 'LineWidth', 2)
        end
        
        %interactions(bee1,bee2,frm) = interbee;
    end
end
%%

%[x, y] = ind2sub([frameWidth,frameHeight],px);
%presence = x^2./a.^2 + y^2./b.^2 <= 1 + interaction_range;
%temp = full(pxVal,px);
%pxVal(px,:) = presence;
%pxVal = sparse(pxVal);
%clearvars im xcenter ycenter xfront yfront a b slope angle x1 x2 y1 y2 t x y w x y bHere nPresentB BW combos interbee