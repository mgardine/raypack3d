function ray_plotcontour(X,Y,Z,w,direction,line)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_plotcontour
%
% This function is internal to the ray_plotmodel function. It generates a
% contour plot of a volumetric slice given by the direction and line
% inputs.
%
% Author: 
% Matt Gardine
% May 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_x = min(min(min(X)));
max_x = max(max(max(X)));
min_y = min(min(min(Y)));
max_y = max(max(max(Y)));
min_z = min(min(min(Z)));
max_z = max(max(max(Z)));

max_w=ceil(max(max(max(w))));
min_w=floor(min(min(min(w))));

x = unique(X);
y = unique(Y);
z = unique(Z);

num_x=size(X,2);
num_y=size(X,1);
num_z=size(X,3);

if strcmp(direction,'xz')
    hslice = surf(x,z,zeros(num_z,num_x)+line);
    rotate(hslice,[-1,0,0],0);
    xd = get(hslice,'XData');
    zd = get(hslice,'YData');
    yd = get(hslice,'ZData');
    delete(hslice);
    hslice = slice(X,Y,Z,w,xd,yd,zd);
    xd = get(hslice,'XData');
    yd = get(hslice,'YData');
    zd = get(hslice,'CData');
    delete(hslice);
    contour(xd,yd,zd,linspace(min_w,max_w,(max_w-min_w)*5+1));
    view(0,180)
elseif strcmp(direction,'yz')
    hslice = surf(y,z,zeros(num_z,num_y)+line);
    rotate(hslice,[0,0,1],0);
    yd = get(hslice,'XData');
    zd = get(hslice,'YData');
    xd = get(hslice,'ZData');
    delete(hslice);
    hslice = slice(X,Y,Z,w,xd,yd,zd);
    xd = get(hslice,'XData');
    yd = get(hslice,'YData');
    zd = get(hslice,'CData');
    delete(hslice);
    contour(xd,yd,zd,linspace(min_w,max_w,(max_w-min_w)*5+1));
    view(90,180)
elseif strcmp(direction,'xy')
    hslice = slice(X,Y,Z,w,[],[],line);
    xd = get(hslice,'XData');
    yd = get(hslice,'YData');
    zd = get(hslice,'CData');
    delete(hslice);
    contour(xd,yd,zd,linspace(min_w,max_w,(max_w-min_w)*5+1));
    view(90,-90)
end

colormap(flipud(jet));
xlim([min_x max_x]);
ylim([min_y max_y]);
zlim([min_z max_z]);
xlabel('Northing (km)')
ylabel('Easting (km)')
zlabel('Depth (km)')

colorbar