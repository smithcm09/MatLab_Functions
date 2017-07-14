function [Dist, horizRadius, totalRadius] = distCalc2 ( sx, sy, sz, rx, ry, rz )
% distCalc.m
% Given a two coordinates (lat, lon, depth), calculate the distance
% between.  Not quite right because a triangle geometry is assumed to add
% the depth, but should be good enough for distances less than 100 km.
% DOES NOT REQUIRE MAPPING TOOLBOX.
% ------------------------------------------------------------------------
% 
%   [Dist, horizRadius, totalRadius] = distCalc ( sx, sy, sz, rx, ry, rz );
%   INPUT:
%   	sx: source longitude in decimal degrees (negative is W longitude)
%  	sy: source latitude in decimal degrees
%  	sz: source depth in km (from velocity model datum, positive down)
%  	rx: reciever longitude in decimal degrees
%  	ry: reciever latitude in decimal degrees
%  	rz: reciever depth in km (from velocity model datum, positive down)
%   OUTPUT:
%  	Dist: Structure of distances from the source
%  		Dist.x: Distance east from source (km)
%  		Dist.y: Distance north from source (km)
%  		Dist.z: Depth from source (km) (same reference datum as sz)
%  	horizRadius: Horizontal radius between source and reciever
%  	totalRadius: Total radius between source and reciever
% -------------------------------------------------------------------------

test = 0;
if test == 1
        sx = -121.695620;
        sy = 45.372779;
        sz = 7;
        rx = -121.791520;
        %rx = -121.695620;
        ry =  45.289670;
        %ry =  45.372779;
        rz = 0;
end
Dist = struct('x', [], 'y', [], 'z', []);
totalRadius = [];
horizRadius = [];
% Use the haversine formula (stable at small distances)
R = 6371; % km
radConv = pi/180;
dLat = (ry-sy)*radConv;
dLon = (rx-sx)*radConv; 
% Calculate total horizontal distance
a = sin(dLat/2) * sin(dLat/2) +...
        cos(sy*radConv) * cos(ry*radConv) *... 
        sin(dLon/2) * sin(dLon/2); 
c = 2 * atan2(sqrt(a), sqrt(1-a)); 
horizRadius = R * c;
% Calculate total west-east distance
dLat = 0;
a = sin(dLat/2) * sin(dLat/2) +...
        cos(sy*radConv) * cos(ry*radConv) *... 
        sin(dLon/2) * sin(dLon/2); 
c = 2 * atan2(sqrt(a), sqrt(1-a)); 
Dist.x = R * c;
% Calculate total north-south distance
dLon = 0;
dLat = (ry-sy)*radConv;
a = sin(dLat/2) * sin(dLat/2) +...
        cos(sy*radConv) * cos(ry*radConv) *... 
        sin(dLon/2) * sin(dLon/2); 
c = 2 * atan2(sqrt(a), sqrt(1-a)); 
Dist.y = R * c;
%-------------------------------------------------------

% Now assume a triangle and add the depth coordinate
Dist.z = abs(rz+sz);

totalRadius = sqrt(Dist.z^2+horizRadius^2);