%% Axis
% Enumeration class for the Cartesian axes x, y, z.

%%% Description
% The only three instances |Axis.x|, |Axis.y|, and |Axis.z| of |Axis| are used
% extensively in MaxwellFDFD to represent and specify the x-, y-, and z-axis of
% the Cartesian coordinate system.
%
% |Axis| is a subclass of |Enumerated|, so it inherits all the methods of
% |Enumerated|.  See <Enumerated.html |Enumerated|> for more details.

%%% Instances
% * |Axis.x|: instance representing the x-axis
% * |Axis.y|: instance representing the y-axis
% * |Axis.z|: instance representing the z-axis

%%% Methods
% Below, |axis| is an instance of |Axis|.
%
% * |[p, q, r] = cycle(axis)|: |[p, q, r]| is a cyclic permutation of |[Axis.x,
% Axis.y, Axis.z]| satisfying |r == axis|.

%%% Example
%   % Test conversion to integers and strings.
%   fprintf('# of instances of Axis: %d\n', Axis.count);
%   for w = Axis.elems
%       fprintf('The integer value of Axis.%s is %d\n', char(w), int(w));
%   end
%
%   % Test cyclic permutation.
%   [p q r] = cycle(Axis.y);
%   fprintf('The cyclic permutation of [x, y, z] beginning with y is [%s, %s, %s]\n', char(r), char(p), char(q));

%%% See Also
% <Enumerated.html |Enumerated|>

classdef OT < Enumerated
    enumeration
        So('So')
        L('L')
        LB('LB')
        U('U')
        Sd('Sd')
        LO('LO')
    end
    
    methods (Static)
        function elems = elems(ind)
            elems = [OT.So, OT.L, OT.LB, OT.U, OT.Sd, OT.LO];
            if nargin > 0  % ind
                elems = elems(ind);
            end
        end
        
        function count = count()
            count = length(OT.elems);
        end
        function name = name(x)
            switch x
                case OT.So
                    name = 'So';
                case OT.L
                    name = 'L';
                case OT.LB
                    name = 'LB';
                case OT.U
                    name = 'U';
                case OT.Sd
                    name = 'Sd';
                case OT.LO
                    name = 'LO';
            end
        end
    end
end
