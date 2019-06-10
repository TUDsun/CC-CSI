function X = back_clb(Y, N, datype)

Y = reshape(full(Y), N(1), N(2), N(3));
switch datype
    case PT.FULL
        Y1 = (Y + circshift(Y, -1, 1)) / 2;
        Y2 = (Y + circshift(Y, -1, 2)) / 2;
        Y3 = (Y + circshift(Y, -1, 3)) / 2;
        Z  = [Y1(:).'; Y2(:).'; Y3(:).'];
    case PT.TE
        %         Z(1, :, :, :) = (Y + circshift(Y, -1, 1)) / 2;
        %         Z(2, :, :, :) = (Y + circshift(Y, -1, 2)) / 2;
        Y1          = circshift(Y, +1, 1);
        Y2          = circshift(Y, +1, 2);
        Y3          = circshift(Y, +1, 3);
        Y13         = circshift(Y1, +1, 3);
        Y23         = circshift(Y2, +1, 3);
        Z(1, :, :, :)  = (Y + Y2 + Y3 + Y23) / 4;
        Z(2, :, :, :)  = (Y + Y1 + Y3 + Y13) / 4;
    case PT.TM
        Z = Y;
end
num_MD  = numel(Z);
X       = sparse(1:num_MD, 1:num_MD, Z(:), num_MD, num_MD);
end