function [x,y,z]=ray_latlon2xyz_flat(lat,lon,alt,ref_lat,ref_lon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_latlon2xyz_flat
%
% This function takes latitude and longitude coordinates and converts them
% into a NED (North [positve x], East [positive y], Down [positive z])
% cartesian coordinate system, assuming a flat earth. In this tranformation, 
% x, y, and z are in kilometers. For a spherical earth
% projection, use ray_lonlon2xyz.
%
% Usage: [x,y,z]=ray_latlon2xyz_flat(lat,lon,alt,ref_lat,ref_lon)
%
% Inputs:
%   lat:        Input latitude coordinate
%   lon:        Input longitude coordinate
%   alt:        Input altitude (in kilometers)
%   ref_lat:    Reference latitude where x will be equal to 0
%   ref_lon:    Reference longitude where y will be equal to 0
%
% Outputs:
%   x:          X coordinate of input
%   y:          Y coordinate of input
%   z:          Z coordinate of input
%
% Author:
% Matt Gardine
% February 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 6367;
%convert to radians
lat = lat*pi/180;
lon = lon*pi/180;
ref_lat = ref_lat*pi/180;
ref_lon = ref_lon*pi/180;

enu(1) = (lon-ref_lon)*cos(ref_lat)*R;
enu(2) = (lat-ref_lat)*R;
enu(3) = alt;

rotation = [0 1 0; 1 0 0; 0 0 -1];
ned = enu*rotation;

x = ned(1);
y = ned(2);
z = ned(3);
