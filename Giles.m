%% Giles: Automated fracture tracing script
% Fracture tracing script based on contrast difference, currently only 
% tested on SEM BSE images. To use: run script. Change parameters according
% to results. Details on the parameters are found in the pdf manual. The
% script requires the image processing toolbox. 
%
% This is an early version, so look out for future updates! 
%
% Giles stands for:
% = fracture tracinG by medIan filter, skeLetonisation, and targEted cloSure =
% and also honours Giles O., tracer of fractures. 
% 
% ========================================================================
% !! When using this script, please cite the following two papers !!
%
% - Aben, F.M., Brantut, N., and Mitchell, T. [submitted], Off-fault damage
% characterisation during and after experimental quasi-static and dynamic
% rupture in crustal rock from laboratory P-wave tomography and
% microstructures.
%
% - Griffiths, L., Heap, M. J., Baud, P., Schmittbuhl, J. [2017],
% Quantification of microcrack characteristics and implications for
% stiffness and strength of granite. International Journal of Rock
% Mechanics and Mining Sciences, 100, 138-150
% ========================================================================
% 
%
% 
%
% --------- F.M. Aben, 03-2020 --------

clear all;
close all;

%% --------------------------- SET PARAMETERS ----------------------------
% MEDIAN FILTER
aperture       = [25 25];

% BINARIZATION AND FILTER PARAMETERS
bw_threshold   = 14;                                                        % binarization threshold (1-100)
wiener_area    = [10 10];                                                   % neighbourhoodsize (x,y) for wiener filter
rescale_range  = [0 1];                                                    % min-max range to rescale grayscale values

% MORPHOLOGICAL PARAMETERS
frac_length     = 8;                                                       % minimum fracture length (given by the major axes length, regionprops function), in pixels
object_limit    = round((2/3).*frac_length);                                % largest object size (nr of pixels) to be removed
gap_length1     = round((1/3).*frac_length);                                % gap length closed by morphological dilation-erosion

% TARGETED CLOSING
ecc_threshold   = 0.85;                                                     % eccentricity threshold (0 is circular, 1 is line)
angle_var       = 15;                                                       % variation in alignment between major axis of component and end-point of a second component (-angle_var to + angle_var, in degrees)
gap_length2     = 20;                                                       % gap length between line segments closed by targeted connections

% PRUNING
pruning_length  = frac_length;                                              % length of branches to be pruned
frac_length2 = 4.*frac_length;                                              % last size-based filter

%% ============================ READ DATA ================================

[FILENAME, PATHNAME] = uigetfile({'*.jpg'},'Select image');
Im_original = imread([PATHNAME FILENAME]);


%% ========================== BINARIZATION ===============================
% FILTER
tic;
Im_filt = medfilt2(Im_original, aperture,'symmetric');                      % median filter
Im_filt = Im_filt - Im_original;                                            % difference image
Im_filt = wiener2(Im_filt, wiener_area);                                    % remove noise
figure; imshow(Im_filt);
toc;

% BINARIZE
tic;
Im_bin                          = mat2gray(Im_filt,rescale_range);          % rescale grayscale values
Im_bin(Im_filt < bw_threshold)  = 0;                                        % set values to zero, based on non-rescaled image
Im_bin(Im_filt >= bw_threshold) = 1;                                        % set values to one, based on non-rescaled image
Im_bin                          = im2bw(Im_bin);                            % binarise image (convert to logical array)
toc;
figure; imshow(Im_bin);

%% =================== CONNECTION OF LINEAR OBJECTS ======================
% ------------------------ MORPHOLOGICAL CLOSING -------------------------

Im_bin      = bwareaopen(Im_bin,object_limit);                              % remove smallest objects
tic;
for alpha   = 1:181                                                         % angle of the line object
    se      = strel('line' ,gap_length1 ,alpha);                             % set orientation of line object
    Im_bin  = imclose(Im_bin, se);                                          % close gaps between line segments by dilation and erosion
end

Im_bin      = bwmorph(Im_bin,'thin', Inf);                                  % skeletonize binary image
toc;


% ------------------------- TARGETED CLOSING -----------------------------

% GET CONNECTED COMPONENT PROPERTIES
STATS           = regionprops(Im_bin,'Eccentricity',...
                  'MajorAxisLength','Orientation');                         % get properties of connected components
STATS           = struct2table(STATS);
CC              = bwconncomp(Im_bin);                                       % get indices for all connected components
ENDS            = bwmorph(Im_bin, 'endpoints');                             % get end-points for ALL connected components
[y_ends,x_ends] = find(ENDS == 1);                                          % x,y-coordinates of end-points

I               = find(STATS.MajorAxisLength >= frac_length ...
                  & STATS.Eccentricity > ecc_threshold);                    % filter connected components

Im_bin_conn     = Im_bin;                                                   % set image that will be connected

disp(['Found ' num2str(length(I)) ' suited connected objects']);

% LOOP OVER FILTERED CONNECTED COMPONENTS
for j = 1:length(I);
    
    disp(['analysing object ' num2str(j) ' of ' num2str(length(I))]);
    
    Im_work                        = Im_bin;
    Im_work(CC.PixelIdxList{I(j)}) = 0;                                     % remove component under consideration from work image
    ENDS_obj                       = bwmorph(Im_bin - Im_work, 'endpoints');% get end-notes of component under consideration
    [y_obj,x_obj]                  = find(ENDS_obj == 1);                   % get x-y coordinates of end-nodes
    
    % END-POINTS OF CONSIDERED COMPONENT
    for i = 1:length(y_obj)
        
        J         = find(y_ends == y_obj(i) & x_ends == x_obj(i));          % set end-points of considered component to zero
        y_ends(J) = NaN;
        x_ends(J) = NaN;
        
    end
    
    % FIND AND CONNECT COMPONENTS FOR EACH END-POINT
    for i = 1:length(y_obj)
        
        dist    = sqrt(abs(x_obj(i) - x_ends).^2 + ...
                  abs(y_obj(i) - y_ends).^2);                               % obtain distances between considered end-point and all end-points
        angle   = atand(-(y_obj(i) - y_ends)./(x_obj(i) - x_ends)) - ...    % obtain angles between above
                  STATS.Orientation(I(j));
        K       = find(dist <= gap_length2 & abs(angle) <= angle_var);       % apply filter based on angle and gap distance
        
        if isempty(K)
            
            continue
            
        elseif length(K) >= 2
            
            [~,Y]   = sort(dist(K));                                        % if multiple K, sort on shortest distance and only use that one
            K       = K(Y(1));
            
        end
        
        dx = x_obj(i) - x_ends(K);                                          % x-distance between considered and target end-point
        dy = y_obj(i) - y_ends(K);                                          % y-distance between considered and target end-point
        
        % GET X,Y-COORDINATES OF POINTS TO BE FILLED
        if abs(dx) >= abs(dy)                                               % for x-distance larger than y-distance
            x_fill     = 0:1:(abs(dx));                                     % x-indices
            if dx < 0                                                       % for target end-point 'in front of' considered end-point
                y_fill = x_fill .* (dy./dx);
                x_fill = x_fill + x_obj(i);
                y_fill = y_fill + y_obj(i);
            else                                                            % for target end-point 'behind' considered end-point
                y_fill = x_fill .* (dy./dx);
                x_fill = x_fill + x_ends(K);
                y_fill = y_fill + y_ends(K);
            end
            y_fill     = round(y_fill);                                     % make integer y-indices
        else                                                                % for y-distance larger than x-distance
            y_fill     = 0:1:(abs(dy));                                     % y-indices
            if dy < 0                                                       % for target end-point 'in front of' considered end-point
                x_fill = -(y_fill .* (dx./dy));
                y_fill = y_fill + y_obj(i);
                y_fill = flip(y_fill);
                x_fill = x_fill + x_ends(K);
            else                                                            % for target end-point 'behind' considered end-point
                x_fill = -(y_fill .* (dx./dy));
                y_fill = y_fill + y_ends(K);
                y_fill = flip(y_fill);
                x_fill = x_fill + x_obj(i);
            end
            x_fill     = round(x_fill);                                     % make integer y-indices
        end
        
        fill              = sub2ind(size(Im_bin),y_fill,x_fill);            % fill connection on image
        
        clearvars y_fill x_fill
        Im_bin_conn(fill) = 1;
    end
end

%% ==================== FILTER CIRCULAR COMPONENTS =======================
Im_bin_conn      = bwmorph(Im_bin_conn,'thin', Inf);                        % skeletonize binary image

STATS   = regionprops(Im_bin_conn,'MajorAxisLength','Eccentricity');        % get new stats
STATS   = struct2table(STATS);
CC      = bwconncomp(Im_bin_conn);                                          % get indices new connected components
I       = find(STATS.MajorAxisLength <= frac_length.*2 ...
          & STATS.Eccentricity < ecc_threshold);                            % find objects smaller than frac_length .* 2

for i = 1:length(I)
    Im_bin_conn(CC.PixelIdxList{I(i)}) = 0;                                 % remove smaller objects
end
% figure; imshow(Im_bin_conn);

%% ============================= PRUNING =================================
close all;

disp('pruning...');tic;
% EROSION
BRANCH = bwmorph(Im_bin_conn,'branchpoints');

for i = 1:pruning_length
    ENDS            = bwmorph(Im_bin_conn, 'endpoints');                    % get end-points for connected components
    OVERLAP = ENDS & BRANCH;
    ENDS = ENDS - OVERLAP;
    END_OLD{i} = ENDS;
    Im_bin_conn = Im_bin_conn - ENDS;
    Im_bin_conn = logical(Im_bin_conn);
end

% DILATION
nhood = ones(3,3);
for i = 1:pruning_length
    ENDS            = bwmorph(Im_bin_conn, 'endpoints');                    % get end-points for connected components
    OVERLAP = ENDS & BRANCH;
    ENDS = ENDS - OVERLAP;
    ENDS = imdilate(ENDS,nhood);
    ADD = ENDS & END_OLD{(pruning_length+1 - i)};
    Im_bin_conn = Im_bin_conn + ADD;
    Im_bin_conn = logical(Im_bin_conn);
end
toc;

x = 4;
STATS   = regionprops(Im_bin_conn,'MajorAxisLength');                       % get new stats
STATS   = struct2table(STATS);
CC      = bwconncomp(Im_bin_conn);                                          % get indices new connected components
I       = find(STATS.MajorAxisLength <= frac_length2);                    % find objects smaller than frac_length .* x
for i = 1:length(I)
    Im_bin_conn(CC.PixelIdxList{I(i)}) = 0;                                 % remove smaller objects
end

%% ============================== FIGURE =================================

[B,L]   = bwboundaries(Im_bin_conn);

figure; 
hold on;
imshow(Im_original);
hold on;
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end

%% ========================== SAVE DATA ==================================
save([FILENAME '_traces.mat'], 'Im_bin_conn');

param.bw_threshold = bw_threshold;
param.wiener_area = wiener_area;
param.rescale_range = rescale_range;
param.aperture = aperture;
param.frac_length = frac_length;
param.frac_length2 = frac_length2;
param.ecc_threshold = ecc_threshold;
param.angle_var = angle_var;
param.gap_length2 = gap_length2;
param.gap_length1 = gap_length1;
param.object_limit = object_limit;

save([FILENAME '_params.mat'], 'param');
