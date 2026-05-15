
%Universidade Federal do Pará
%Intituto de Tecnologia
%Faculdade de Engenharia
%Author 1: Remo Magalhães de Souza
%Author 2: Edilson Morais

%Desenha Malha indeformed considerando  o element de Membrana
function ElementMembranaDeformada(Dados,Resultado,EscalaDesloc)

%Elemento Membrana
% calcula as coordenadas deformeds
figure(2)
axis equal  
for i = 1: Dados.NumNos
   coordsDef(i,1) = Dados.coords(i,1) + Resultado.DeslocNos(i,1)*EscalaDesloc;
   coordsDef(i,2) = Dados.coords(i,2) + Resultado.DeslocNos(i,2)*EscalaDesloc;
end

% determina dimensoes minimas e maximas da structure
Dados.coordsTot = [Dados.coords;
                   coordsDef];

xmin = min(Dados.coordsTot(:,1));
ymin = min(Dados.coordsTot(:,2));
xmax = max(Dados.coordsTot(:,1));
ymax = max(Dados.coordsTot(:,2));

dx = xmax - xmin;
dy = ymax - ymin;

lado = max([dx,dy]);
delta = lado/10;

axis ([xmin-delta,xmax+delta,ymin-delta,ymax+delta]);

delta = lado/40;


%Create the patch by specifying the Faces, Vertices, and FaceVertexCData properties as well as the FaceColor property.
for i=1:Dados.NumElem
    faces(i,:) = [Dados.conect(i,1) Dados.conect(i,5) Dados.conect(i,2) Dados.conect(i,6) Dados.conect(i,3) Dados.conect(i,7) Dados.conect(i,4) Dados.conect(i,8)];
end

patch('faces',faces,'vertices',Dados.coords,'FaceColor','none', 'Linestyle', ':', 'EdgeColor', [0.1 0.1 0.1])
patch('faces',faces,'vertices',coordsDef,'FaceColor','none', 'EdgeColor', [1.0 0.0 0.0])



%Generates boundary conditions
for no = 1:Dados.NumNos
                   if (Dados.restr(no,1) == 1)
                         d = delta;
                         %cx = coordsDef(no,1);
                         %cy = coordsDef(no,2);
                         cx = Dados.coords(no,1);
                         cy = Dados.coords(no,2);

                         %Creates the triangle representing the DOF restraint
                         %relative to x-translation
                         ux =  [cx, cx-d,   cx-d, cx];   %Vertices do triangle x
                         uy =  [cy, cy+d/2, cy-d/2, cy]; %Vertices do triangle y

                         H=line(ux,uy);
                         set(H,'color',[0,0.7,0]);
                         set(H,'LineWidth',[2]);
                      end

                      if (Dados.restr(no,2) == 1)
                         cx = coordsDef(no,1);
                         cy = coordsDef(no,2);

                         ux =  [cx, cx-d/2, cx+d/2, cx];
                         uy =  [cy, cy-d,   cy-d, cy];

                         H=line(ux,uy);
                         set(H,'color',[0,0.7,0]);
                         set(H,'LineWidth',[2]);
                      end

% 
%                   if (normF(no) ~= 0)
%                      [ux, uy, uz] = criaSeta(Dados.CargaNo(no,:), Dados.coords(no,:), escalaVetor, 0);
% 
%                      H=line(ux,uy,uz);
%                      set(H,'color',[0,0.7,0]);
%                      set(H,'LineWidth',[2]);
                    %end
end

title ('Configuracao Deformada');
xlabel ('X');
ylabel ('Y');

set(gcf,'DefaultLineColor','red') 
return