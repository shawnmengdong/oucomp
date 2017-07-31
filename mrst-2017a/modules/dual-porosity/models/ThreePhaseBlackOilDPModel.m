classdef ThreePhaseBlackOilDPModel < DualPorosityReservoirModel
% Three phase with optional dissolved gas and vaporized oil
properties
end

methods
    function model = ThreePhaseBlackOilDPModel(G, rock, fluid, rock_matrix, fluid_matrix, dp_info, varargin)
        model = model@DualPorosityReservoirModel(G, rock, fluid, rock_matrix, fluid_matrix, dp_info, varargin{:});

        % Typical black oil is disgas / dead oil, but all combinations
        % are supported
        model.vapoil = false;
        model.disgas = false;

        % Max increments
        model.drsMaxAbs = inf;
        model.drsMaxRel = inf;

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        % All phases are present
        model.oil = true;
        model.gas = true;
        model.water = true;
        model.saturationVarNames = {'sw', 'so', 'sg'};
		model.matrixVarNames = {'pwm','pom','pgm','swm','som','sgm'};

        model = merge_options(model, varargin{:});

        d = model.inputdata;
        if ~isempty(d)
            % Assume ECL-style input deck, as this is the only
            % supported format at the moment.
            if isfield(d, 'RUNSPEC')
                if isfield(d.RUNSPEC, 'VAPOIL')
                    model.vapoil = d.RUNSPEC.VAPOIL;
                end
                if isfield(d.RUNSPEC, 'DISGAS')
                    model.disgas = d.RUNSPEC.DISGAS;
                end
            else
                error('Unknown dataset format!')
            end
        end
    end
    
    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)
        switch(lower(name))
            case {'rsm', 'rvm'}
                % RS and RV for gas dissolving into the oil phase and oil
                % components vaporizing into the gas phase respectively.
                fn = lower(name);
                index = 1;
            otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@DualPorosityReservoirModel(model, name);
        end
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsBlackOilDP(state0, state, model, dt, ...
                        drivingForces, varargin{:});

    end

    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Check parent class
        state = validateState@DualPorosityReservoirModel(model, state);
        nc = model.G.cells.num;
        if model.disgas
            % RS must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'rsm', nc, 1);
        else
            % RS does not really matter. Assign single value.
            fn = model.getVariableField('rs');
            if ~isfield(state, fn)
                dispif(model.verbose, ...
                    ['Missing field "', fn, '" added since disgas is not enabled.\n']);
                state.(fn) = 0;
            end
            clear fn
            
            % Doing the same for the matrix
            fn = model.getVariableField('rsm');
            if ~isfield(state, fn)
                dispif(model.verbose, ...
                    ['Missing field "', fn, '" added since disgas is not enabled.\n']);
                state.(fn) = 0;
            end
            clear fn
        end
        if model.vapoil
            % RV must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'rv', nc, 1);
            
            % RV must be supplied for all cells. This may cause an error.
            model.checkProperty(state, 'rvm', nc, 1);
        else
            % RS does not really matter. Assign single value.
            fn = model.getVariableField('rv');
            if ~isfield(state, fn)
                dispif(model.verbose, ...
                    ['Missing field "', fn, '" added since vapoil is not enabled.\n']);
                state.(fn) = 0;
            end
            clear fn
            
            % Doing the same for the matrix
            fn = model.getVariableField('rvm');
            if ~isfield(state, fn)
                dispif(model.verbose, ...
                    ['Missing field "', fn, '" added since disgas is not enabled.\n']);
                state.(fn) = 0;
            end
            clear fn
        end
    end
    
    % --------------------------------------------------------------------%
    function scaling = getScalingFactorsCPR(model, problem, names)
        % Get approximate, impes-like pressure scaling factors
        nNames = numel(names);
        
        scaling = cell(nNames, 1);
        handled = false(nNames, 1);
        
        % Take averaged pressure for scaling factors
        state = problem.state;
        fluid = model.fluid;
        p = mean(state.pressure);
        
        for iter = 1:nNames
            name = lower(names{iter});
            switch name
                case 'oil'
                    if model.disgas
                       rs = fluid.rsSat(p);
                       bO = fluid.bO(p, rs, true);
                    else
                       bO = fluid.bO(p);
                    end
                    s = 1./bO;
                case 'water'
                    bW = fluid.bW(p);
                    s = 1./bW;
                case 'gas'
                    if model.vapoil
                        rv = fluid.rvSat(p);
                        bG = fluid.bG(p, rv, true);
                    elseif model.gas
                        bG = fluid.bG(p);
                    end
                    s = 1./bG;
                otherwise
                    continue
            end
            sub = strcmpi(problem.equationNames, name);
            
            scaling{iter} = s;
            handled(sub) = true;
        end
        if ~all(handled)
            % Get rest of scaling factors from parent class
            other = getScalingFactorsCPR@DualPorosityReservoirModel(model, problem, names(~handled));
            [scaling{~handled}] = other{:};
        end
    end
       
end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
