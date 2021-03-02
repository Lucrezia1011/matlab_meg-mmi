function [tfced] = matlab_tfce_transform_gridtime(img,H,E,dh)
%MATLAB_TFCE_TRANSFORM performs threshold free cluster enhancement
%   [tfced] = matlab_tfce_transform(img,H,E,C,ndh) performs threshold
%   free cluster enhancement on 'img' as per Smith & Nichols (2009).
%   -- img the 3D image to be transformed
%   -- H height exponent,   H = 2
%   -- E extent exponent,   E = 0.5
%   -- C connectivity,      C = 80
%   -- dh size of steps for cluster formation   dh = 0.1
% https://github.com/markallenthornton/MatlabTFCE/blob/master/matlab_tfce_transform.m

% set cluster thresholds
threshs = 0:dh:max(img(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(img(:));

% find connected components
vals = zeros(nvox,1);

% bwconncomp_nd is in private Matlab directory not allowed in path
% /usr/local/apps/Matlab/R2020a/toolbox/images/images/private 
% copied to ~/matlab/

pixelIdxList = arrayfun(@(x) builtin('_pixelIdxListsn', bsxfun(@ge,img,x), 3^4-1), threshs,'UniformOutput', false);

for h = 1:ndh
    clustsize = zeros(nvox,1);
%     ccc = cc(h);
%     voxpercc = cellfun(@numel,ccc.PixelIdxList);
%     
%     for c = 1:ccc.NumObjects
%         clustsize(ccc.PixelIdxList{c}) = voxpercc(c);
%     end
%     
    ccc = pixelIdxList{h};
    voxpercc = cellfun(@numel,ccc);
    for c = 1:length(pixelIdxList{h})
        clustsize(ccc{c}) = voxpercc(c);
    end
    % calculate transform
    curvals = (clustsize.^E).*(threshs(h)^H);
    vals = vals + curvals;
end
tfced = NaN(size(img));
tfced(:) = vals.*dh;

end