function [up, theta, r] = multiphaseUpwindIndices(G, vT, T, K, upstr)
    nPh = numel(G);
    nF  = numel(T);
    G = flattenAndMakeDouble(G);
    K = flattenAndMakeDouble(K);

    [G, sortInd] = sort(G, 2);
    theta = repmat(vT, 1, nPh);
    for l = 1:nPh
        for j = 1:nPh
            if j == l
                continue
            end
            flag = repmat(j > l, nF, 1);
            kj = nan(nF, 1);
            for k = 1:nPh
                % Upstream weight mobility
                kloc = upstr(flag, K(:, k));
                % Lookup sort index for this phase (G is sorted by weight,
                % K is not)
                ind = sortInd(:, j) == k;
  
                kj(ind) = kloc(ind);
            end
            theta(:, l) = theta(:, l) + T.*(G(:, l) - G(:, j)).*kj;
        end
    end
    subs = theta < 0;
    [ix, r] = max(~subs, [], 2);
    r = r - 1;
    r(all(subs, 2)) = nPh;
    
    up = false(nF, nPh);
    for l = 1:nPh
        up(:, l) = (sortInd(:, l) > r);
    end
end

function d = flattenAndMakeDouble(d)
    d = cellfun(@double, d, 'UniformOutput', false);
    d = [d{:}];
end

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
