function [tfced] = matlab_tfce_transform_MEGtime(T,H,E,dh,neighbours)
%MATLAB_TFCE_TRANSFORM performs threshold free cluster enhancement
%   [tfced] = matlab_tfce_transform(img,H,E,C,ndh) performs threshold
%   free cluster enhancement on 'img' as per Smith & Nichols (2009).
%   -- img the 3D image to be transformed
%   -- H height exponent,   H = 2
%   -- E extent exponent,   E = 0.5
%   -- A adjecency matrix of neighbouring channels
%   -- dh size of steps for cluster formation   dh = 0.1
% https://github.com/markallenthornton/MatlabTFCE/blob/master/matlab_tfce_transform.m

% set cluster thresholds
threshs = 0:dh:max(T(:));
threshs = threshs(2:end);
ndh = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(T(:));

% find connected components
vals = zeros(nvox,1);
if size(T,1)~=length(neighbours)
    T = T';
end
% cfg_neighb        = [];
% cfg_neighb.method = 'template';%'distance';
% cfg_neighb.template = 'CTF275_neighb.mat';
% cfg_neighb.channel = channels;
% neighbours        = ft_prepare_neighbours(cfg_neighb, hdr);

% channels = {neighbours.label};
% for n = 1:length(neighbours)
%     [~,~,iB] = intersect(neighbours(n).neighblabel, channels );
%     neighbours(n).neighbnum =iB;
% end

N = size(T,1);
for h = 1:ndh
    
      
    edgeA = cell(size(T,2),1);
    
    parfor tt = 1:size(T,2)-1
%      
        % Add spatial edges for next time points
        
        % Only need to check t+1 (t-1 is equivalently calculated in previous loop iteration)
%         Ttresh_m = T(:,tt-1)>=threshs(h);
        Ttresh = T(:,tt)>=threshs(h);
        Ttresh_p = T(:,tt+1)>=threshs(h);
        
        idx = find(Ttresh);
      
        neighboursT = neighbours(Ttresh);    
        edge_Add = cell(length(neighboursT),1);
        A = zeros(N);
        for n = 1:length(neighboursT)
            B = zeros(N,1);
            B(neighboursT(n).neighbnum) = 1;
            B(idx(n)) = 1;
            edge_p = find(Ttresh_p & B);
%             edge_m = find(Ttresh_m & B);
%             edge_Add{n,1} = [ones(nnz(edge_p)+nnz(edge_m),1)*idx(n) + N*(tt-1), ...
%                 [edge_p + N*(tt); edge_m + N*(tt-2)] ];
            edge_Add{n,1} = [ones(nnz(edge_p),1)*idx(n) + N*(tt-1), edge_p + N*(tt)];
            A(idx(n),Ttresh & B)  =1;                
        end
        Gt = graph(A);
        
        edgeA{tt} = [Gt.Edges.EndNodes + N*(tt-1); cell2mat(edge_Add)];      
        
    end
    
    tt = size(T,2);
    Ttresh = T(:,tt)>=threshs(h);
    
    idx = find(Ttresh);
    neighboursT = neighbours(Ttresh);
    A = zeros(N);
    for n = 1:length(neighboursT)
        B = zeros(N,1);
        B(neighboursT(n).neighbnum) = 1;
        B(idx(n)) = 1;
        A(idx(n),Ttresh & B)  =1;
    end
    Gt = graph(A);
    
    edgeA{tt} =Gt.Edges.EndNodes + N*(tt-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    edgeA = cell2mat(edgeA);
    
        G = graph(edgeA(:,1),edgeA(:,2),1);
        nodes = unique(G.Edges.EndNodes);
        bins = conncomp(G);
        binnodes = accumarray(bins', 1:numel(bins), [], @(v) {sort(v')});
        binnodes1 = [];
        n = 0; % only include with edges
        for binidx = 1:numel(binnodes)    
            if ismember(binnodes{binidx}, nodes) 
                n = n+1;
                binnodes1{n,1}  = binnodes{binidx};
            end
        end
    
        clustsize = zeros(nvox,1);

        voxpercc = cellfun(@numel,binnodes1);
        for c = 1:length(binnodes1)
            clustsize(binnodes1{c}) = voxpercc(c);
        end
        % calculate transform
        curvals = (clustsize.^E).*(threshs(h)^H);
        vals = vals + curvals;
    
    clc; fprintf('Calculating %.0f/%.0f\n',h,ndh)
end
tfced = NaN(size(T));
tfced(:) = vals.*dh;

end