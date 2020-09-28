function z = fisher(r)
% Fisher transformation of correlation coefficient 
% USE: z = fisher(r)

if ischar(r)
    r=str2double(r);
end
z = 0.5 * log((1+r)./(1-r));

