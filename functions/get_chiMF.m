function x = get_chiMF(v1, v2, N, src_n, omega, datype)

switch datype
    case PT.FULL
        V1      = forward_clb(v1, N, src_n);
        V2      = forward_clb(v2, N, src_n);
        x       = sum(sum(V1 .* conj(V2), 3), 1) ./ (sum(sum(V2 .* conj(V2), 3), 1) + eps);
        x       = x(:);
    case PT.TM
        V1R     = cellfun(@(x, y) real(sum(x .* conj(y), 2)), v1, v2, ...
            'UniformOutput', false);
        
        V2R 	= cellfun(@(x) sum(x .* conj(x), 2), v2, ...
            'UniformOutput', false);
        
        xR      = sum(cell2mat(V1R), 2) ./ sum(cell2mat(V2R), 2);
        
        V1I     = cellfun(@(x, y, z) imag(sum(x .* conj(y), 2)) * omega{1} / z, v1, v2, omega, ...
            'UniformOutput', false);
        
        V2I 	= cellfun(@(x, y) sum(x .* conj(x), 2) * (omega{1} / y) ^ 2, v2, omega, ...
            'UniformOutput', false);
        
        xI      = sum(cell2mat(V1I), 2) ./ sum(cell2mat(V2I), 2);
        x       = xR + 1j * xI;
    case PT.TE
        V1      = cellfun(@(x) forward_clb(x, N, src_n, datype), v1, ...
            'UniformOutput', false);
        
        V2      = cellfun(@(x) forward_clb(x, N, src_n, datype), v2, ...
            'UniformOutput', false);
        
        V1R     = cellfun(@(x, y) real(sum(sum(x.*conj(y), 3), 1)).', V1, V2, ...
            'UniformOutput', false);
        
        V2R 	= cellfun(@(x) sum(sum(x.*conj(x), 3), 1).', V2, ...
            'UniformOutput', false);
        
        xR      = sum(cell2mat(V1R), 2) ./ sum(cell2mat(V2R), 2);
        
        V1I     = cellfun(@(x, y, z) imag(sum(sum(x .* conj(y), 3), 1)).' * omega{1} / z, V1, V2, omega, ...
            'UniformOutput', false);
        
        V2I 	= cellfun(@(x, y) sum(sum(x .* conj(x), 3), 1).'*(omega{1} / y) ^ 2, V2, omega, ...
            'UniformOutput', false);
        
        xI      = sum(cell2mat(V1I), 2) ./ sum(cell2mat(V2I), 2);
        x       = xR + 1j * xI;
    otherwise
        V1TE    = cellfun(@(x) forward_clb(x, N, src_n, PT.TE), v1(2, :), ...
            'UniformOutput', false);
        
        V2TE    = cellfun(@(x) forward_clb(x, N, src_n, PT.TE), v2(2, :), ...
            'UniformOutput', false);
        
        V1RTE   = cellfun(@(x, y) real(sum(sum(x .* conj(y), 3), 1)).', V1TE, V2TE, ...
            'UniformOutput', false);
        
        V2RTE	= cellfun(@(x) sum(sum(x .* conj(x), 3), 1).', V2TE, ...
            'UniformOutput', false);
        
        V1ITE   = cellfun(@(x, y, z) imag(sum(sum(x .* conj(y), 3), 1)).' * omega{1} / z, V1TE, V2TE, omega, ...
            'UniformOutput', false);
        
        V2ITE	= cellfun(@(x, y) sum(sum(x .* conj(x), 3), 1).'*(omega{1} / y)^2, V2TE, omega, ...
            'UniformOutput', false);
        
        
        V1RTM   = cellfun(@(x, y) real(sum(x .* conj(y), 2)), v1(1, :), v2(1, :), ...
            'UniformOutput', false);
        
        V2RTM 	= cellfun(@(x) sum(x .* conj(x), 2), v2(1, :), ...
            'UniformOutput', false);
        
        V1ITM   = cellfun(@(x, y, z) imag(sum(x .* conj(y), 2))*omega{1} / z, v1(1, :), v2(1, :), omega, ...
            'UniformOutput', false);
        
        V2ITM	= cellfun(@(x, y) sum(x .* conj(x), 2) * (omega{1} / y) ^ 2, v2(1, :), omega, ...
            'UniformOutput', false);
        
        xR      = (sum(cell2mat(V1RTM), 2) + sum(cell2mat(V1RTE), 2)) ./ (sum(cell2mat(V2RTM), 2) + sum(cell2mat(V2RTE), 2));
        xI      = (sum(cell2mat(V1ITM), 2) + sum(cell2mat(V1ITE), 2)) ./ (sum(cell2mat(V2ITM), 2) + sum(cell2mat(V2ITE), 2));
        x       = xR + 1j * xI;
end
end