% =========================================================================
% obtemgdlelem
%
% Returns the global DOF indices of an element, in the same node-then-DOF
% ordering used by the element stiffness and force routines.
% =========================================================================

function [GDLelem] = obtemgdlelem(Dados, Nosdoelem, GDL, Numnoselem)

cont = 0;
for no = 1:Numnoselem
    for i = Nosdoelem(no)
        for j = 1:Dados.NumGDLNo
            cont = cont + 1;
            GDLelem(cont, 1) = GDL(i, j);
        end
    end
end

return
