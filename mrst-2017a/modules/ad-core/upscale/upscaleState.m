function state = upscaleState(coarsemodel, model, state)
%Create a upscaled state by simple processing of values
%
% SYNOPSIS:
%   state_coarse = upscaleState(coarsemodel, model, state_fine)
%
% DESCRIPTION:
%   Convert a state for a fine model into a realization of the same state
%   for a coarse model.
%
% REQUIRED PARAMETERS:
%   coarsemodel - A coarse model derived from the fine model.
%
%   model       - The fine model. Subclass of ReservoirModel.
%
%   state       - State to be converted. Should correspond to model.
%
%
% RETURNS:
%   state       - Coarse state suitable for the coarsemodel.
%

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    p = coarsemodel.G.partition;
    CG = coarsemodel.G;
    
    pvc = coarsemodel.operators.pv;
    pvf = model.operators.pv;
    
    counts = accumarray(p, 1);
    state_f = state;
    
    nph = size(state.s, 2);
    pvs = bsxfun(@times, state.s, pvf);
    
    % Calculate saturations based on new and old pore volume
    s = zeros(CG.cells.num, nph);
    for i = 1:nph
        s(:, i) = accumarray(p, pvs(:, i))./pvc;
    end
    state.s = s;
    if isfield(state, 'components')
        z = [state.components{:}];
        pvz = bsxfun(@times, z, pvf);
        for i = 1:numel(state.components)
            state.components{i} = accumarray(p, pvz(:, i))./pvc;
        end
    end
    
    if isfield(state, 'rs')
        sO = model.getProp(state_f, 'sO');
        sO_c = coarsemodel.getProp(state, 'sO');
        state.rs = accumarray(p, sO.*state.rs.*pvf)./(sO_c.*pvc);
        state.rs(~isfinite(state.rs)) = 0;
    end
    if isfield(state, 'rv')
        sG = model.getProp(state_f, 'sG');
        sG_c = coarsemodel.getProp(state, 'sG');
        state.rv = accumarray(p, sG.*state.rv.*pvf)./(sG_c.*pvc);
        state.rv(~isfinite(state.rv)) = 0;
    end
    % Average the pressure (not entirely correct for compressible systems,
    % but we won't start evaluating properties in here).
    state.pressure = accumarray(p, state.pressure)./counts;
    if isfield(state, 'T')
        state.T = accumarray(p, state.T)./counts;
    end
    
    state.flux = zeros(CG.faces.num, nph);
end
