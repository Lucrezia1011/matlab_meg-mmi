function z = fisher(r)
% Fisher transformation of correlation coefficient 
% USE: z = fisher(r)

z = 0.5 * log((1+r)./(1-r));
end