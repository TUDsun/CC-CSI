function [J_cell, M_cell, Ms] = myassign_source(grid3d, srcj_array, srcm_array)


chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid.');
chkarg(istypesizeof(srcj_array, 'Source', [1 0]), ...
    '"srcj_array" should be row vector of instances of Source.');
srcn = max(length(srcj_array),length(srcm_array));
ind_Ms = zeros(1,srcn);
if ~isempty(srcj_array)
    J_cell = cell(length(srcj_array), Axis.count);
    ii = 1;
    for src = srcj_array
        for w = Axis.elems
            %             J_cell{ii, w} = sparse3d(grid3d.N);
            J_cell{ii, w} = ndSparse.build(grid3d.N);
        end
        ii = ii+1;
    end
    ii = 1;
    
    for src = srcj_array
        for w = Axis.elems
            [ind, JMw_patch] = src.generate(w, grid3d);
            if ~isempty(JMw_patch)
                J_cell{ii, w}(ind{:}) = J_cell{ii, w}(ind{:}) + JMw_patch;  % superpose sources
                ind_Ms(ii) = ((ind{3} - 1) * grid3d.N(2) + ind{2} - 1) * grid3d.N(1) + ind{1};
            end
        end
        ii = ii + 1;
    end
else
    J_cell = cell(length(srcm_array), Axis.count);
    ii = 1;
    for src = srcm_array
        for w = Axis.elems
            J_cell{ii, w} = ndSparse.build(grid3d.N);
        end
        ii = ii+1;
    end
end
if ~isempty(srcm_array)
    M_cell = cell(length(srcm_array), Axis.count);
    ii = 1;
    for src = srcm_array
        for w = Axis.elems
            M_cell{ii, w} = ndSparse.build(grid3d.N);
        end
        ii = ii+1;
    end
    ii = 1;
    for src = srcm_array
        for w = Axis.elems
            [ind, JMw_patch] = src.generate(w, grid3d);
            if ~isempty(JMw_patch)
                M_cell{ii, w}(ind{:}) = M_cell{ii, w}(ind{:}) + JMw_patch;  % superpose sources
                ind_Ms(ii) = ((ind{3} - 1) * grid3d.N(2) + ind{2} - 1) * grid3d.N(1) + ind{1};
            end
        end
        ii = ii+1;
    end
else
    M_cell = cell(length(srcj_array), Axis.count);
    ii = 1;
    for src = srcj_array
        for w = Axis.elems
            %             M_cell{ii, w} = sparse3d(grid3d.N);
            M_cell{ii, w} = ndSparse.build(grid3d.N);
        end
        ii = ii+1;
    end
end
ind_Ms = Axis.count*(ind_Ms - 1);
ind_Ms = [ind_Ms + 1; ind_Ms + 2; ind_Ms + 3];
ind_Ms = ind_Ms(:).';
Ms = sparse(1:Axis.count*srcn, ind_Ms, ones(1,Axis.count*srcn), Axis.count*srcn, Axis.count*prod(grid3d.N));
end


