function CWT_ICWT_MHW(id, imseries_index, protein, timeSlice, aScale, aMax, da, mesh)
%
% CWT_ICWT_MHW.m
% 2D, MHW

% timeSlice is the frame to deconstruct
% mesh: number of mesh points per pixel

load_data_info

datafile = sprintf('%s/%s_f%03d.mat', original_data_location, protein, timeSlice);
load(datafile)

imageSlice = double(im);

ascale = aScale;
a_max = aMax;

pad = 2 * 3 * a_max; % needs to be >= 2 * 3 * a_max

X = size(imageSlice,2);
Y = size(imageSlice,1);

Lx = mesh * X;
xLength = Lx;
bxLength = xLength; % xLength = bLength
LLx = Lx + 2 * pad;

Ly = mesh * Y; 
yLength = Ly;
byLength = yLength; % xLength = bLength
LLy = Ly + 2 * pad;

c_psi = 2*pi^2;

cwt = [];

cwt_pad = zeros(LLy,LLx);
icwt_pad = zeros(LLy,LLx);

if mesh > 1,
    imZoom = zoomImage(imageSlice, mesh); % creates imageZoom
else
    imZoom = double(imageSlice); % convert to double
end

clear imageSlice
imZoom = padarray(imZoom, [pad pad]);

aWin = 3 * ceil(ascale); % 3 x 'a' window size

[lengthX, lengthY] = meshgrid(1:2*aWin,1:2*aWin);

winD = 1/ascale ...
    * (2 - ((lengthX-aWin)/ascale).^2 - ((lengthY-aWin)/ascale).^2) ... 
    .* exp(-((lengthX-aWin)/ascale).^2/2 - ((lengthY-aWin)/ascale).^2/2);

for bY = 1:LLy-2*aWin,
    for bX = 1:LLx-2*aWin,

        cwt_part = imZoom(bY:bY+2*aWin-1,bX:bX+2*aWin-1) .* winD ;
        cwt_pad(bY+aWin,bX+aWin) = sum(cwt_part(:));

    end
end

for bY = 1:LLy-2*aWin, 
    for bX = 1:LLx-2*aWin,

        icwt_part = zeros(LLy,LLx);
        icwt_part(1+bY:2*aWin+bY,1+bX:2*aWin+bX) = 1/c_psi*cwt_pad(bY+aWin,bX+aWin).*winD/(ascale^3)*da;

        icwt_pad = icwt_pad + icwt_part;
    end
end

cwt = cwt_pad(pad+1:end-pad,pad+1:end-pad);
icwt = icwt_pad(pad+1:end-pad,pad+1:end-pad);

icwt_file = sprintf('%s/%s/cwt_icwt_%s_t%03d_a%03d.mat', destination, id, id, timeSlice, aScale);
save( icwt_file, 'cwt', 'icwt', '-v7.3')
