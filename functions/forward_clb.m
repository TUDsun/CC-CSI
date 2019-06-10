function V = forward_clb(v1, N, src_n, ptype)
if nargin<4, ptype = PT.FULL; end
Np = 1:3;
V = reshape(full(v1), Np(ptype), N(1), N(2), N(3), src_n);
switch ptype
    case PT.FULL
        V(1, :, :, :, :)    = (V(1, :, :, :, :) + circshift(V(1, :, :, :, :), 1, 2)) / 2;
        V(2, :, :, :, :)    = (V(2, :, :, :, :) + circshift(V(2, :, :, :, :), 1, 3)) / 2;
        V(3, :, :, :, :)    = (V(3, :, :, :, :) + circshift(V(3, :, :, :, :), 1, 4)) / 2;
    case PT.TE
        %         V(1, :, :, :, :)    = (V(1, :, :, :, :) + circshift(V(1, :, :, :, :), 1, 2)) / 2;
        %         V(2, :, :, :, :)    = (V(2, :, :, :, :) + circshift(V(2, :, :, :, :), 1, 3)) / 2;
        V1          = V(1, :, :, :, :);
        V2          = V(2, :, :, :, :);
        
        V12         = circshift(V1, -1, 2 + 1);
        V13         = circshift(V1, -1, 3 + 1);
        
        V123        = circshift(V12, -1, 3 + 1);
        
        V21         = circshift(V2, -1, 1 + 1);
        V23         = circshift(V2, -1, 3 + 1);
        
        V213        = circshift(V21, -1, 3 + 1);
        
        V(1, :, :, :, :)    = (V1 + V12 + V13 + V123) / 4;
        V(2, :, :, :, :)    = (V2 + V21 + V23 + V213) / 4;
    case PT.TM
end
V   = reshape(V, Np(ptype), prod(N), src_n);
end