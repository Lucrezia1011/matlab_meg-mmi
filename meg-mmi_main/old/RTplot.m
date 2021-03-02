function RTplot(RT,bins)
[n,edges] = histcounts(RT,bins);
hold all

m = max(n);

for ii = 1:bins
    fill([edges(ii) edges(ii+1) edges(ii+1) edges(ii)] ,[-1 -1 1 1], [1 0 0],'edgecolor','none','facealpha',n(ii)*0.3/m)
end