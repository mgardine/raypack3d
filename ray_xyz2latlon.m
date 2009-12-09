function [lat,lon,alt] = ray_xyz2latlon(x,y,z,ref_lat,ref_lon,ref_alt,projection)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function ray_xyz2latlon
%
% This function takes NED (North [positve x], East [positive y], Down 
% [positive z]) cartesian coordinate system and converts them to latitude & 
% longitude coordinates, for either a flat earth or a spherical earth. 
% In this tranformation, x, y, and z are in kilometers.
%
% Usage: [lat,lon,elev]=ray_latlon2xyz(x,y,z,ref_lat,ref_lon,ref_alt)
%
% Inputs:
%   x:          X coordinate of data point
%   y:          Y coordinate of data point
%   z:          Z coordinate of data point
%   ref_lat:    Reference latitude where x will be equal to 0
%   ref_lon:    Reference longitude where y will be equal to 0
%   ref_alt:    Reference altitude where z will be equal to 0
%   projection: Desired projection type.  Valid options are 'flat' for a
%                   flat-earth projection, or 'spherical' for a
%                   spherical-earth projection.
%
% Outputs:
%   lat:        Output latitude coordinate
%   lon:        Output longitude coordinate
%   elev:       Output elevation (in kilometers)
%
% Author:
% Matt Gardine
% February 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(projection,'spherical')
  % Convert north, east, down coordinates (labeled x, y, z) to 
  % east, north, up (ENU) coordinates
  %ned = [x; y; z+1.169149329402717e+04];
  ned = [x; y; z];
  rotation = [0 1 0; 1 0 0; 0 0 -1];
  
  enu = rotation*ned;
  
  x = enu(1);
  y = enu(2);
  z = enu(3);
  
  % Convert east, north, up coordinates (labeled x, y, z) to ECEF
  % coordinates. The reference point (phi, lambda, h) must be given. All distances are in metres
  
  ref_lat = ref_lat/180*pi;
  ref_lon = ref_lon/180*pi;
  a = 6378137.0; % earth semimajor axis in meters 
  f = 1/298.257223563; % reciprocal flattening 
  e2 = 2*f -f^2; % eccentricity squared 
  b = a*(1-f);% semi-minor axis
  ep2 = f*(2-f)/((1-f)^2); % second eccentricity squared
  
  % Convert reference lat, lon, altitude to ECEF X,Y,Z
  ref_chi = sqrt(1-e2*(sin(ref_lat)).^2); 
  Xr = (a./ref_chi +ref_alt).*cos(ref_lat).*cos(ref_lon); 
  Yr = (a./ref_chi +ref_alt).*cos(ref_lat).*sin(ref_lon); 
  Zr = (a*(1-e2)./ref_chi + ref_alt).*sin(ref_lat);
  
  phiP = atan2(Zr,sqrt(Xr^2+Yr^2)); % Geocentric latitude
 
  X = -sin(ref_lon)*x - cos(ref_lon)*sin(phiP)*y + cos(ref_lon)*cos(phiP)*z + Xr;
  Y =  cos(ref_lon)*x - sin(ref_lon)*sin(phiP)*y + cos(phiP)*sin(ref_lon)*z + Yr;
  Z = cos(phiP)*y + sin(phiP)*z + Zr;

  % Convert ECEF to Lat,Lon,altitude
  r2 = X.^2+Y.^2;
  r = sqrt(r2);
  E2 = a^2 - b^2;
  F = 54*b^2*Z.^2;
  G = r2 + (1-e2)*Z.^2 - e2*E2;
  c = (e2*e2*F.*r2)./(G.*G.*G);
  s = ( 1 + c + sqrt(c.*c + 2*c) ).^(1/3);
  P = F./(3*(s+1./s+1).^2.*G.*G);
  Q = sqrt(1+2*e2*e2*P);
  ro = -(e2*P.*r)./(1+Q) + sqrt((a*a/2)*(1+1./Q) - ((1-e2)*P.*Z.^2)./(Q.*(1+Q)) - P.*r2/2);
  tmp = (r - e2*ro).^2;
  U = sqrt( tmp + Z.^2 );
  V = sqrt( tmp + (1-e2)*Z.^2 );
  zo = (b^2*Z)./(a*V);
 
  alt = U.*( 1 - b^2./(a*V));
  lat = atan( (Z + ep2*zo)./r )*180/pi;
  lon = atan2(Y,X)*180/pi;
  
elseif strcmp(projection,'flat')
    R = 6367;
    ned = [x; y; z];
    rotation = [0 1 0; 1 0 0; 0 0 -1];
  
    enu = rotation*ned;
  
    x = enu(1);
    y = enu(2);
    z = enu(3);

    lat = ref_lat + y*180/(pi*R);
    lon = ref_lon + x*180/(pi*R*cos(lat*pi/180));
    elev = z;
end
