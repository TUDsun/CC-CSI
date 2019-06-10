function ro = InvInd(Mo, ptype)
switch ptype
    case PT.TE
        Mo   = repmat(Mo.', 2, 1);
        Mo   = Mo(:);
        ro   = find(Mo == 1);
    case PT.TM
        ro   = find(Mo == 1);
    otherwise
        Mo   = repmat(Mo.', 3, 1);
        Mo   = Mo(:);
        ro   = find(Mo == 1);
end
end