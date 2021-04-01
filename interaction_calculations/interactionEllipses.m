function fInteractions = interactionEllipses(tag,shape,expansion,imSize)

xcenter = tag(1,:);
ycenter = tag(2,:);
xfront = tag(3,:);
yfront = tag(4,:);

fWidth = imSize(1);
fHeight = imSize(2);

% ELLIPSES PARAMETERS

a = expansion * shape(1,:)/2;
b = expansion * shape(2,:)/2;
angle = atan2(yfront - ycenter,xfront - xcenter);
x0 = xcenter - 1/2*a.*cos(angle);
y0 = ycenter - 1/2*a.*sin(angle);

t = linspace(0,2*pi,10);
X = a.*cos(t.');
Y = b.*sin(t.');
x = x0 + X.*cos(angle) - Y.*sin(angle);
y = y0 + X.*sin(angle) + Y.*cos(angle);

% INTERACTIONS BEES PRESENT ON THE FRAME

bHere = find(any(x));
nPresentB = size(bHere,2);

nPop = size(xcenter,2);
fInteractions = zeros(nPop,nPop);

if nPresentB > 1
    BW = nan(fHeight,fWidth,nPresentB);
    
    for bee = 1:nPresentB
        BW(:,:,bee) = poly2mask(x(:,bHere(bee)), y(:,bHere(bee)), fHeight, fWidth);
    end
    
    combos = nchoosek(1:nPresentB,2);
    
    for pair = 1:size(combos,1)
        bee1 = combos(pair,1);
        bee2 = combos(pair,2);
        interbee = any(any(BW(:,:,bee1)+BW(:,:,bee2)>1));
        
        intB1 = bHere(bee1);
        intB2 = bHere(bee2);
        fInteractions(intB1,intB2) = interbee;
        
    end
end