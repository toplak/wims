function CWT_BGW(id, imseries_index, protein, timeSlice, angle_index, total_angles, th_cluster, a_value, xi, mesh)
% CWT_BGW.m
% 2D, BGW, Guillimin
% 
% th_cluster - number of angles to group together for the run
% mesh - number of mesh points per pixel

datafile = sprintf('%s_f%03d.mat', protein, timeSlice);
load(sprintf('%s/%s',original_data_location,datafile)); % matrix called 'im' or 'im2'

imageSlice = double(im);

a_str = decimal2pt(a_value);
xi_str = decimal2pt(xi);

% directory to save results
ddir = sprintf('%s_%03d_%02d_%s', id, timeSlice, a_value, xi_str);

pad = 2.0 * 3 * a_value; % needs to be >= 2 * 3 * a_max

Lx = mesh * size(imageSlice, 2);
xLength = Lx;
bxLength = xLength; % xLength = bLength
LLx = Lx + 2 * pad;
Ly = mesh * size(imageSlice, 1);
yLength = Ly;
byLength = yLength; % xLength = bLength
LLy = Ly + 2 * pad;

imageSlice = padarray(imageSlice, [pad pad], mean(imageSlice(:)));

for the_angle = angle_index:angle_index + (th_cluster - 1),

    cwtPadded = zeros(LLy,LLx);

    theta = (the_angle-1)*pi/total_angles;

    ascale = a_value; 
    aWin = 2.0 * 3 * ceil(ascale) + 1; % 3 x 'a' window size

    [lengthX, lengthY] = meshgrid(1:2*aWin+1,1:2*aWin+1);

    winD = 1/ascale * (2 - (((lengthX-aWin)*cos(theta) - (lengthY-aWin)*sin(theta))/ascale).^2 - (((lengthX-aWin)*sin(theta) + (lengthY-aWin)*cos(theta))).^2/(xi^2)).*exp(-(((lengthX-aWin)*cos(theta) - (lengthY-aWin)*sin(theta))/ascale).^2/2 - (((lengthX-aWin)*sin(theta) + (lengthY-aWin)*cos(theta))/(xi)).^2/2);

    clear lengthX lengthY

    for bY = 1:LLy-2*aWin,
        for bX = 1:LLx-2*aWin,

            cwtPadded(bY+aWin,bX+aWin) = sum(sum( imageSlice(bY:bY+2*aWin,bX:bX+2*aWin) .* winD ));
        end
    end

    clear winD

    cwt = cwtPadded(pad+1:end-pad,pad+1:end-pad);

    cwt_file = sprintf('%s/%s/%s/cwt_%s_t%03d_th%03d.mat', destination, id, ddir, id, timeSlice, the_angle);

    save(cwt_file, 'cwt')

    % save max and min wavelet coefficients
    max_cwt = max(cwt(:));
    min_cwt = min(cwt(:));
    cwt_max_file = sprintf('%s/%s/max/max_%s_f%03d_a%s_xi%s_th%03d.mat', destination, id, id, timeSlice, a_str, xi_str, the_angle);
    save(cwt_max_file, 'max_cwt', 'min_cwt')

    clear cwt cwtPadded
end

