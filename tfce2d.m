function TFCE = tfce2d(Tst)

E = 0.5;
H = 2;
S=regionprops(Tst>0,'PixelIdxList','PixelList');
TFCE = tfextent(Tst,S,E,H);

S=regionprops(Tst<0,'PixelIdxList','PixelList');
TFCEn = tfextent(-Tst,S,E,H);

TFCE(TFCEn>0) = -TFCEn(TFCEn>0);

function TFCE = tfextent(Tst,S,E,H)
[Nx,Ny] = size(Tst); % time dimension
TFCE = zeros(size(Tst));
% convert into area
dx = 1;

if Nx == 1 || Ny == 1 % 1D
    
    for ii = 1:length(S)
        a = S(ii).PixelIdxList;
        
        h = Tst(a);
        
        dh = 0.1;
        for p = 1:length(h)
            A = 0;
            hp = dh;
            while hp < h(p)
                e = hp <= h;
                de = diff(e);
                de(end+1,:) = 0;
                e = e & ~de;
                for t = 1:length(e)
                    if e(t) == 0 && t<p
                        e(1:t) = 0;
                    elseif e(t) == 0 && t>p
                        e(t:end) = 0;
                    end
                end
                % e =  time extent
                
                eE = (sum(e)*dx)^E;
                A = A + eE*(hp^H)*dh;
                hp = hp + dh;
            end
            TFCE(a(p)) = A;
        end
    end
    
    
else % 2D
    
    for ii = 1:length(S)
        a = S(ii).PixelList;
        for ff = 1:Ny
            indt = a(a(:,1) == ff,2);
            h = Tst(indt,ff);
            
            dh = 0.1;
            for p = 1:length(h)
                A = 0;
                hp = dh;
                while hp < h(p)
                    e = hp <= h;
                    de = diff(e);
                    de(end+1,:) = 0;
                    e = e & ~de;
                    for t = 1:length(e)
                        if e(t) == 0 && t<p
                            e(1:t) = 0;
                        elseif e(t) == 0 && t>p
                            e(t:end) = 0;
                        end
                    end
                    % e =  time extent
                    
                    eE = (sum(e)*dx)^E;
                    A = A + eE*(hp^H)*dh;
                    hp = hp + dh;
                end
                TFCE(indt(p),ff) = A;
            end
        end
        
        tfce1 = TFCE;
        for tt = 1:Nx
            indf = a(a(:,2) == tt,1);
            h = Tst(tt,indf);
            
            for p = 1:length(h)
                df = 1;
                ind = (indf - indf(p)) == df;
                while h(p) > h(ind)
                    TFCE(tt,indf(p)) = TFCE(tt,indf(p)) + tfce1(tt,indf(ind))^E;
                    df = df+1;
                    ind = (indf - indf(p)) == df;
                end
                
                df = -1;
                ind = (indf - indf(p)) == df;
                while h(p) > h(ind)
                    TFCE(tt,indf(p)) = TFCE(tt,indf(p)) + tfce1(tt,indf(ind))^E;
                    df = df-1;
                    ind = (indf - indf(p)) == df;
                end
            end
            
        end
    end
    
end