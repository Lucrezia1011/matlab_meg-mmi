function [tfced] = matlab_tfce_transform_MEG(T,H,E,dh,neighbours)
%MATLAB_TFCE_TRANSFORM_MEG performs threshold free cluster enhancement 
%   [tfced] = matlab_tfce_transform_MEG(img,H,E,C,dh,neighbours) performs threshold
%   free cluster enhancement on a free spatial arrangement as per Smith & Nichols (2009).
%   -- T the 3D image to be transformed
%   -- H height exponent,   H = 2
%   -- E extent exponent,   E = 0.5
%   -- dh size of steps for cluster formation   dh = 0.1
%   -- neighbours, structure with information of channel neighbours
%       cfg_neighb          = [];
%       cfg_neighb.method   = 'template';
%       cfg_neighb.template = 'CTF275_neighb.mat';
%       cfg_neighb.channel  = channels;
%       neighbours          = ft_prepare_neighbours(cfg_neighb, hdr);
%       channels = {neighbours.label};
%       for n = 1:length(neighbours)
%           [~,~,iB] = intersect(neighbours(n).neighblabel, channels );
%           neighbours(n).neighbnum =iB;
%       end
% Based on grid TFCE algorithm:
% https://github.com/markallenthornton/MatlabTFCE/blob/master/matlab_tfce_transform.m

% set cluster thresholds
threshs = 0:dh:max(T(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(T(:));

% find connected components
vals = zeros(nvox,1);

for h = 1:ndh
    
    Ttresh = T>=threshs(h);
    
    idx = find(Ttresh);
    neighboursT = neighbours(Ttresh);
    A = zeros(size(T,2)); % adjancency matrix
    for n = 1:length(neighboursT)
        B = zeros(1,size(T,2));
        B(neighboursT(n).neighbnum) = 1;
        A(idx(n),Ttresh & B)  =1;
    end
    G = graph(A);
    bins = conncomp(G);
    binnodes = accumarray(bins', 1:numel(bins), [], @(v) {sort(v')});
    
    for binidx = numel(binnodes):-1:1
        if ismember(binnodes{binidx}, find(Ttresh)) == 0
            binnodes(binidx)  =[];
        end
    end
    
    clustsize = zeros(nvox,1);
   
    voxpercc = cellfun(@numel,binnodes);
    for c = 1:length(binnodes)
        clustsize(binnodes{c}) = voxpercc(c);
    end
    % calculate transform
    curvals = (clustsize.^E).*(threshs(h)^H);
    vals = vals + curvals;
    
end
tfced = NaN(size(T));
tfced(:) = vals.*dh;

end