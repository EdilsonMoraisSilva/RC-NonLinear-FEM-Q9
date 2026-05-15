% Universidade Federal do Pará (UFPA)
% Technology Center (CT)
% Civil Construction Department (DCC)
% Civil Engineering Program
% Author: Remo Magalhaes de Souza       (remo@ufpa.br)
% Modified by: Edilson Morais


function DesenhaMalhaIndefQ9(Dados,escala,ArmaduraDirx,ArmaduraDiry)


figure (1);

set(gcf,'DefaultLineColor','red')  % define a cor "default" das linhas
axis equal                         % mesma escala em x e y


% Determine the minimum and maximum dimensions of the structure
xmin = min(Dados.coords(:,1));
ymin = min(Dados.coords(:,2));
xmax = max(Dados.coords(:,1));
ymax = max(Dados.coords(:,2));

dx = xmax - xmin;
dy = ymax - ymin;

lado = max([dx,dy]);
delta = lado/10;

axis ([xmin-delta,xmax+delta,ymin-delta,ymax+delta]);

delta = lado/120;

for no = 1: length(Dados.CargaNo(:,1))
   normF(no) = norm(Dados.CargaNo(no,:));
end    

forcaMax = max(normF);

escalaVetor = escala*lado/forcaMax;


%Creates region representing the quadrilateral finite elements
Dados.coords;
for i=1:Dados.NumElem
    faces(i,:) = [Dados.conect(i,1) Dados.conect(i,5) Dados.conect(i,2) Dados.conect(i,6) Dados.conect(i,3) Dados.conect(i,7) Dados.conect(i,4) Dados.conect(i,8)];
end
patch('faces',faces,'vertices',Dados.coords,'FaceColor',[0.7 0.7 0.7], 'EdgeColor', [0.1 0.1 0.1])

NumElem = length(Dados.conect(:,1));
%
for el = 1:NumElem 
   
   % determina os nos do element
   noI = Dados.conect(el,1);
   noJ = Dados.conect(el,2);
   noK = Dados.conect(el,3);
   noL = Dados.conect(el,4);
   noM = Dados.conect(el,5);
   noN = Dados.conect(el,6);
   noO = Dados.conect(el,7);
   noP = Dados.conect(el,8);
   noQ = Dados.conect(el,9);
   
   % determina as coordenadas dos nos do element
   coordsElem(1,:) = Dados.coords(noI,:);  
   coordsElem(2,:) = Dados.coords(noJ,:);  
   coordsElem(3,:) = Dados.coords(noK,:);  
   coordsElem(4,:) = Dados.coords(noL,:);  
   coordsElem(5,:) = Dados.coords(noM,:);  
   coordsElem(6,:) = Dados.coords(noN,:);  
   coordsElem(7,:) = Dados.coords(noO,:);  
   coordsElem(8,:) = Dados.coords(noP,:);
   coordsElem(9,:) = Dados.coords(noQ,:);
   % Encontra o centroide do element de forma a centralizar o seu numero
   % (por exemplo element "1")
   coordCentro = (coordsElem(1,:) + coordsElem(2,:) + coordsElem(3,:)+coordsElem(4,:))/4;
  
   switch Dados.Desenho.LabelEl
       case 'sim'
           %Returns o numero do element
           elemNum = sprintf('%d',el);
           Ht=text(coordCentro(1),coordCentro(2),elemNum);
           set(Ht,'color',[0,0,1]);
       case 'nao'
   end

end    

%Generates boundary conditions and node labels of the structure
hold on;
for no=1: Dados.NumNos
  xNo = Dados.coords(no,1);
  yNo = Dados.coords(no,2);
  
  plot(xNo,yNo,'bo','Markersize',6);
  
    switch Dados.Desenho.LabelNo
        case 'sim'
            noNum = sprintf('%d',no);
            Ht=text(xNo+delta,yNo+delta,noNum);
            set(Ht,'color',[0.7,0,0]);
        case 'nao'
    end
  
  elta = lado/40;
  for cont = 1:length(Dados.conect(:,1))
        switch Dados.TipoElemento(cont,:)
            case 'ElemMembrana'
                    if (Dados.restr(no,1) == 1)
                         d = elta;
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
                         d = elta; 
                         cx = Dados.coords(no,1);
                         cy = Dados.coords(no,2);

                         ux =  [cx, cx-d/2, cx+d/2, cx];
                         uy =  [cy, cy-d,   cy-d, cy];

                         H=line(ux,uy);
                         set(H,'color',[0,0.7,0]);
                         set(H,'LineWidth',[2]);
                      end


                  if (normF(no) ~= 0)
                     [ux, uy] = criaSeta(Dados,Dados.CargaNo(no,:), Dados.coords(no,:), escalaVetor, 0);

                     H=line(ux,uy);
                     set(H,'color',[0,0.7,0]);
                     set(H,'LineWidth',[2]);
                  end

            case 'ElemPlacaMindlin'
                  if Dados.restr(no,:) == [1 1 1]
                         d = delta;
                         cx = Dados.coords(no,1);
                         cy = Dados.coords(no,2);

                         ux =  [cx,cx,cx-d,cx-d,cx,cx];
                         uy =  [cy,cy+d/2,cy+d/2,cy-d/2,cy-d/2,cy];

                         H=line(ux,uy);
                         set(H,'color',[0,0.7,0]);
                         set(H,'LineWidth',[2]);

                  elseif (Dados.restr(no,1) == 1) %Apenas o DOF a translation em z is restrained
                         d = delta;
                         cx = Dados.coords(no,1);
                         cy = Dados.coords(no,2);

                         ux =  [cx, cx-d,   cx-d, cx];
                         uy =  [cy, cy+d/2, cy-d/2, cy];

                         H=line(ux,uy);
                         set(H,'color',[0,0.7,0]);
                         set(H,'LineWidth',[2]);
                     

%                       if (Dados.restr(no,3) == 1)
%                          cx = Dados.coords(no,1);
%                          cy = Dados.coords(no,2);
% 
%                          ux =  [cx, cx-d/2, cx+d/2, cx];
%                          uy =  [cy, cy-d,   cy-d, cy];
% 
%                          H=line(ux,uy);
%                          set(H,'color',[0,0.7,0]);
%                          set(H,'LineWidth',[2]);
%                       end
                  end

      end
  end
end


switch Dados.TipoArm
    case 'Discretizada'
        if ArmaduraDirx == 'sim'
            %Desenha Armadura discretizada para x
            for i=1:length(Dados.Arm.NosArm.x(:,1))
            plot(Dados.Arm.CoordArm.x(i).Coord(:,1),Dados.Arm.CoordArm.x(i).Coord(:,2),'-k','LineWidth',2)
            end
        end
        
        if ArmaduraDiry == 'sim'
            %Desenha Armadura discretizada para y
            for i=1:length(Dados.Arm.NosArm.y(:,1))
            plot(Dados.Arm.CoordArm.y(i).Coord(:,1),Dados.Arm.CoordArm.y(i).Coord(:,2),'-k','LineWidth',2)
            end
        end
end

title ('Malha de Elementos Finitos');
xlabel ('X');
ylabel ('Y');

set(gcf,'DefaultLineColor','red') 

return

