%% ST 
% Enumeration class for the Cartesian axes x, y, z.

%%% Description
% The only two instances |ST.inc|, and |ST.sct| of |ST| are used
% extensively in MaxwellFDFD to represent and specify the simulations of the incident 
% field and the scattering field.
%
% |ST| is a subclass of |Enumerated|, so it inherits all the methods of
% |Enumerated|.  See <Enumerated.html |Enumerated|> for more details.

%%% Instances
% * |ST.inc|: instance representing the simulation of the incident field 
% * |ST.sct|: instance representing the simulation of the scattering field 

%%% See Also
% <Enumerated.html |Enumerated|>

classdef ST < Enumerated
	enumeration
		inc('inc')
		sct('sct')
	end

	methods (Static)
		function elems = elems(ind)
			elems = [ST.inc, ST.sct];
			if nargin > 0  % ind
				elems = elems(ind);
			end
        end
		
		function count = count()
			count = length(ST.elems);
		end
    end
end
