clear
xdim = 13;
ydim = 13;
    xcoord = repmat(1:xdim,ydim,1);
    xcoord = xcoord(:);                % vector of x-coordinates for each patch
    ycoord = repmat([1:ydim]',1,xdim);
    ycoord = ycoord(:);                % vector of y-coordinates for each patch

dists=zeros(length(xcoord),length(xcoord));
    for ii = 1:length(xcoord)
        for i = 1:length(xcoord)
            dists(ii,i)=((xcoord(ii)-xcoord(i))^2+(ycoord(ii)-ycoord(i))^2)^.5;
        end
    end
    dists=floor(dists);

inds = 1:length(xcoord);
coords = [xcoord ycoord inds'];
center = 85;    

ones = coords(dists(:,center)==1,:);
twos = coords(dists(:,center)==2,:);
threes = coords(dists(:,center)==3,:);
fours = coords(dists(:,center)==4,:);
fives = coords(dists(:,center)==5,:);

plot(ones(:,1),ones(:,2),'.')
axis([1 xdim 1 ydim])
grid on

plot(twos(:,1),twos(:,2),'.')
axis([1 xdim 1 ydim])
grid on

plot(threes(:,1),threes(:,2),'.')
axis([1 xdim 1 ydim])
grid on

plot(fours(:,1),fours(:,2),'.')
axis([1 xdim 1 ydim])
grid on

plot(fives(:,1),fives(:,2),'.')
axis([1 xdim 1 ydim])
grid on