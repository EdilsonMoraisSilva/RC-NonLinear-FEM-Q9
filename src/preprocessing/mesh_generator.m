% =========================================================================
% mesh_generator
%
% Universidade Federal do Pará (UFPA)
% Instituto de Tecnologia (ITEC)
% Faculdade de Engenharia Civil
% Author: Edilson Morais Lima e Silva
%
% Generates a structured 2D mesh of 9-node quadratic quadrilateral (Q9)
% elements by sweeping through one or more rectangular segments stacked
% along the y direction.
%
% Each segment defines an independent block of the structure. The
% generator places (2*Nx + 1) x (2*Ny + 1) nodes per segment (the extra
% factor of 2 accommodates the mid-side and centre nodes of the Q9
% element) and accumulates them sequentially in Dados.coords.
%
% Consecutive segments share their common edge: when the segment
% dimensions change between t-1 and t, the offset flags 'm' and 'n'
% prevent duplicated nodes along that shared boundary.
%
% INPUT:
%   Dados   - global problem data struct (will receive Dados.coords)
%   Segment - K x 4 matrix, one row per segment, with columns:
%               col 1: segment length along x  (m)
%               col 2: segment length along y  (m)
%               col 3: number of elements in x (Nx)
%               col 4: number of elements in y (Ny)
%
% OUTPUT:
%   Dados.coords - N x 2 array of nodal coordinates [x, y]
% =========================================================================

function [Dados] = mesh_generator(Dados, Segment)

% Number of element divisions per segment
Nx = Segment(:, 3);
Ny = Segment(:, 4);

% --- Initialise auxiliary counters ---------------------------------------
k  = 1;     % running node index in Dados.coords
Tx = 0;     % cumulative offset in x (kept for future multi-column meshes)
Ty = 0;     % cumulative offset in y between segments
c  = 0;     % spare counter (reserved)

% Skip flags: when a segment shares an edge with the previous one,
% the first row/column is skipped to avoid duplicated nodes.
m = 0;      % column-skip flag (x direction)
n = 0;      % row-skip flag    (y direction)

% =========================================================================
% Loop over segments (stacked in the y direction)
% =========================================================================
for t = 1:size(Segment, 1)

    % Sub-element grid spacing inside this segment.
    % Note the factor 2: Q9 elements have nodes at corners, mid-sides
    % and centre, so the nodal grid is twice as fine as the element grid.
    dx = Segment(t, 1) / Nx(t);
    dy = Segment(t, 2) / Ny(t);

    % --- Decide whether to skip the shared boundary row/column -----------
    % Only relevant from the second segment onward.
    if t ~= 1
        if Segment(t, 1) - Segment(t-1, 1) ~= 0
            m = 1;   % x length changed -> skip first column of nodes
        else
            m = 0;
        end

        if Segment(t, 2) - Segment(t-1, 2) ~= 0
            n = 1;   % y length changed -> skip first row of nodes
        else
            n = 0;
        end
    end

    % --- Generate nodes inside this segment ------------------------------
    for i = n:2*Ny(t)
        for j = m:2*Nx(t)
            Dados.coords(k, 1) = (dx/2) * j;
            Dados.coords(k, 2) = (dy/2) * i + Ty;
            k = k + 1;
        end
    end

    % --- Legacy block (unconditional generation) -------------------------
    % Kept commented for reference. The current implementation above
    % handles the first segment naturally via the (m,n) = (0,0) defaults.
    %
    % else
    %     for i = 0:2*Ny(t)
    %         for j = 0:2*Nx(t)
    %             Dados.coords(k,1) = (dx/2)*j;
    %             Dados.coords(k,2) = (dy/2)*i + Ty;
    %             k = k + 1;
    %         end
    %     end
    % end

    % Advance the y offset for the next segment
    Ty = Ty + Segment(t, 2);
end

end