% Perform depth-based blurring of input (source) image based on its depth
% map (input), to create an output image with only image regions lying at
% depths of user's interest (input) in focus, and everything else blurred
function im_db = depth_blur(im_src, depth_map, depth_focus, blur_rate)
% INPUT ARGUMENTS:
% - im_src      : Input (source) image on which depth-blurring will be done
% - depth_map   : Depth map of the input image
% - depth_focus : Depth range(s) of the image regions to be kept in focus
% - blur_rate   : Extent to which regions of non-interest should be blurred
%
%       CAUTION : If you don't supply any value for the last input argument,
%                 I will assume it as 10, and perform blurring accordingly!
%
% OUTPUT:
% -       im_db : Output image; To view it, use imshow(im_db). Only image
%                 regions lying in user's depths of interest are in focus
%
%      AUTHOR   : Subhayan Mukherjee (subhayan001@gmail.com)
%      DATE     : Apr.2014
%
% Basic error checking
if ~exist('im_src', 'var') || ~exist('depth_map', 'var') || ~exist('depth_focus', 'var'), disp('First three parameters must be specified!'); end
%
% Default value for the last input argument
if ~exist('blur_rate', 'var'), blur_rate = 10; end

I = imread('3.jpg');
D = imread('2.jpg');
im_db = imfilter(I, fspecial('gaussian', blur_rate, blur_rate), 'replicate');

for i = 1: 1: size(I, 1)
    for j = 1: 1: size(I, 2)
        if numel(find(depth_focus == D(i, j), 1, 'first')) > 0
            im_db(i, j, :) = I(i, j, :);
        end
    end
end

im_db = uint8(im_db);