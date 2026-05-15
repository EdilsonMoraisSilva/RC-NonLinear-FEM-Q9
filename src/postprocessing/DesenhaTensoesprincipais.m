% =========================================================================
% DesenhaTensoesprincipais
%
% Postprocessing routine that draws the principal stress field on top of
% the deformed mesh. Two arrows are drawn at each Gauss point, one per
% principal direction; arrow length is proportional to |sigma_i|, color
% indicates sign (tension / compression).
%
% Universidade Federal do Pará (UFPA)
% Authors: Remo Magalhães de Souza, Edilson Morais Lima e Silva
% =========================================================================

function DesenhaTensoesprincipais(Dados, Resultado, escala)
% DesenhaTensoesprincipais - Visualises principal stresses, isolines and cracks
%
% TensoesPrinc(1,:) = s1       TensoesPrinc(5,:) = crack flag s1
% TensoesPrinc(2,:) = s2       TensoesPrinc(6,:) = crack flag s2
% TensoesPrinc(3,:) = TauMax   TensoesPrinc(7,:) = crack opening w_s1 [m]
% TensoesPrinc(4,:) = theta    TensoesPrinc(8,:) = crack opening w_s2 [m]
%
% Outputs:
%   Figure 3  - Principal stresses (arrows: red=tension, blue=compression)
%   Figure 4  - Smooth gradient of s1 (no isolines)
%   Figure 5  - Smooth gradient of s2 (no isolines)
%   Figure 6  - Crack map

%% -----------------------------------------------------------------------
%  GENERAL SETTINGS
% -----------------------------------------------------------------------
lwSeta   = 1.5;
lwMalha  = 0.5;
nCores   = 512;
nGrid    = 300;

corTracao     = [0.85 0.10 0.10];   % red - tension
corCompressao = [0.10 0.25 0.85];   % blue    - compression
corMalha      = [0.20 0.20 0.20];

%% -----------------------------------------------------------------------
%  STRUCTURE BOUNDS
% -----------------------------------------------------------------------
xmin = min(Dados.coords(:,1));
ymin = min(Dados.coords(:,2));
xmax = max(Dados.coords(:,1));
ymax = max(Dados.coords(:,2));
lado = max([xmax-xmin, ymax-ymin]);
margem  = lado / 10;
axisLim = [xmin-margem, xmax+margem, ymin-margem, ymax+margem];

%% -----------------------------------------------------------------------
%  FACES DOS ELEMENTOS
% -----------------------------------------------------------------------
faces = zeros(Dados.NumElem, 8);
for i = 1:Dados.NumElem
    c = Dados.conect(i,:);
    faces(i,:) = [c(1) c(5) c(2) c(6) c(3) c(7) c(4) c(8)];
end

%% -----------------------------------------------------------------------
%  COLETA DE PONTOS DE GAUSS, TENSOES E FISSURAS
% -----------------------------------------------------------------------
todosS1     = [];
todosS2     = [];
todasCoords = [];
todosFlag1  = [];
todosFlag2  = [];
todosW1     = [];
todosW2     = [];
todosTheta  = [];

for i = 1:length(Dados.conect(:,1))
    Ptsprinc = Resultado.el(i).TensoesPrinc;
    Ptcoords = Resultado.el(i).CoordsPtsGauus;
    NumPts   = 4;

    for j = 1:NumPts
        todosS1(end+1)       = Ptsprinc(1,j);
        todosS2(end+1)       = Ptsprinc(2,j);
        todasCoords(end+1,:) = Ptcoords(j,:);
        todosTheta(end+1)    = Ptsprinc(4,j);

        if size(Ptsprinc,1) >= 8
            todosFlag1(end+1) = Ptsprinc(5,j);
            todosFlag2(end+1) = Ptsprinc(6,j);
            todosW1(end+1)    = Ptsprinc(7,j);
            todosW2(end+1)    = Ptsprinc(8,j);
        else
            todosFlag1(end+1) = 0;
            todosFlag2(end+1) = 0;
            todosW1(end+1)    = 0;
            todosW2(end+1)    = 0;
        end
    end
end

tensaomax   = max(max(abs([todosS1; todosS2])));
escalaVetor = escala * lado / tensaomax;

cmapDiv = buildDivergentColormap(nCores);   % azul-branco-vermelho (com sinal)

%% -----------------------------------------------------------------------
%  GRADE DE INTERPOLACAO PARA PCOLOR (Figs 4 e 5)
% -----------------------------------------------------------------------
xg = linspace(xmin, xmax, nGrid);
yg = linspace(ymin, ymax, nGrid);
[Xg, Yg] = meshgrid(xg, yg);

F_s1 = scatteredInterpolant(todasCoords(:,1), todasCoords(:,2), todosS1', 'natural','nearest');
F_s2 = scatteredInterpolant(todasCoords(:,1), todasCoords(:,2), todosS2', 'natural','nearest');

Zs1 = mascararForaDominio(Xg, Yg, F_s1(Xg,Yg), Dados.coords);
Zs2 = mascararForaDominio(Xg, Yg, F_s2(Xg,Yg), Dados.coords);

%% =======================================================================
%  FIGURE 3 - TENSOES PRINCIPAIS (SETAS CLASSICAS)
%  Red = tension | Blue = compression
% =======================================================================
figure(3); clf;
set(gcf,'Name','Tensoes Principais','Color',[0.97 0.97 0.97],'Position',[80 100 920 660]);
ax1 = axes('Parent',gcf);
hold(ax1,'on'); axis(ax1,'equal'); axis(ax1,axisLim);

% Malha de fundo branca
patch('Parent',ax1,'Faces',faces,'Vertices',Dados.coords,...
      'FaceColor',[1 1 1],'EdgeColor',corMalha,'LineWidth',lwMalha);

% Arrows: color by sign (tension/compression), size by magnitude
for i = 1:length(Dados.conect(:,1))
    Ptsprinc = Resultado.el(i).TensoesPrinc;
    Ptcoords = Resultado.el(i).CoordsPtsGauus;
    for j = 1:4
        desenhaSetas(Ptsprinc(:,j), Ptcoords(j,:), escalaVetor, ...
                     corTracao, corCompressao, lwSeta, ax1);
    end
end

xlabel(ax1,'Comprimento (m)','FontSize',11,'FontWeight','bold');
ylabel(ax1,'Altura (m)','FontSize',11,'FontWeight','bold');
title(ax1,'Tensoes Principais  |  {\color{red}Tracao}  /  {\color{blue}Compressao}',...
          'FontSize',13,'FontWeight','bold');
grid(ax1,'on'); ax1.GridAlpha = 0.12; box(ax1,'on'); hold(ax1,'off');

%% =======================================================================
%  FIGURE 4 - SMOOTH GRADIENT OF s1
% =======================================================================
figure(4); clf;
set(gcf,'Name','Tensao Principal Maxima s1','Color',[0.97 0.97 0.97],'Position',[130 80 920 660]);
ax4 = axes('Parent',gcf);
hold(ax4,'on'); axis(ax4,'equal'); axis(ax4,axisLim);

hp4 = pcolor(ax4, Xg, Yg, Zs1);
set(hp4,'EdgeColor','none'); shading(ax4,'interp');
colormap(ax4, cmapDiv);
axes(ax4); caxis([-tensaomax, tensaomax]);
cb4 = colorbar(ax4);
cb4.Label.String = 's1  (Pa)';
cb4.Label.FontSize = 10; cb4.Label.FontWeight = 'bold';

% Malha por cima
patch('Parent',ax4,'Faces',faces,'Vertices',Dados.coords,...
      'FaceColor','none','EdgeColor',corMalha,'LineWidth',lwMalha);

xlabel(ax4,'Comprimento (m)','FontSize',11,'FontWeight','bold');
ylabel(ax4,'Altura (m)','FontSize',11,'FontWeight','bold');
title(ax4,'Tensao Principal Maxima  s1','FontSize',13,'FontWeight','bold');
grid(ax4,'on'); ax4.GridAlpha = 0.08; box(ax4,'on'); hold(ax4,'off');

%% =======================================================================
%  FIGURE 5 - SMOOTH GRADIENT OF s2
% =======================================================================
figure(5); clf;
set(gcf,'Name','Tensao Principal Minima s2','Color',[0.97 0.97 0.97],'Position',[180 60 920 660]);
ax5 = axes('Parent',gcf);
hold(ax5,'on'); axis(ax5,'equal'); axis(ax5,axisLim);

hp5 = pcolor(ax5, Xg, Yg, Zs2);
set(hp5,'EdgeColor','none'); shading(ax5,'interp');
colormap(ax5, cmapDiv);
axes(ax5); caxis([-tensaomax, tensaomax]);
cb5 = colorbar(ax5);
cb5.Label.String = 's2  (Pa)';
cb5.Label.FontSize = 10; cb5.Label.FontWeight = 'bold';

patch('Parent',ax5,'Faces',faces,'Vertices',Dados.coords,...
      'FaceColor','none','EdgeColor',corMalha,'LineWidth',lwMalha);

xlabel(ax5,'Comprimento (m)','FontSize',11,'FontWeight','bold');
ylabel(ax5,'Altura (m)','FontSize',11,'FontWeight','bold');
title(ax5,'Tensao Principal Minima  s2','FontSize',13,'FontWeight','bold');
grid(ax5,'on'); ax5.GridAlpha = 0.08; box(ax5,'on'); hold(ax5,'off');

%% =======================================================================
%  FIGURE 6 - MAPA DE FISSURAS
% =======================================================================
temFissura = any(todosFlag1 > 0) || any(todosFlag2 > 0);

figure(6); clf;
set(gcf,'Name','Mapa de Fissuras','Color',[0.97 0.97 0.97],'Position',[230 40 920 660]);
ax6 = axes('Parent',gcf);
hold(ax6,'on'); axis(ax6,'equal'); axis(ax6,axisLim);

patch('Parent',ax6,'Faces',faces,'Vertices',Dados.coords,...
      'FaceColor',[1 1 1],'EdgeColor',[0.75 0.75 0.75],'LineWidth',0.4);

if ~temFissura
    text(mean([xmin xmax]), mean([ymin ymax]), ...
         'Nenhuma crack detectada neste incremento', ...
         'HorizontalAlignment','center','FontSize',12,'Color',[0.5 0.5 0.5],'Parent',ax6);
else
    comprTraco = lado / 25;
    wmax = max([todosW1, todosW2]);
    if wmax < 1e-12, wmax = 1; end

    cmapFiss = buildFissuraColormap(256);

    for k = 1:length(todasCoords)
        xk      = todasCoords(k,1);
        yk      = todasCoords(k,2);
        theta_k = todosTheta(k);

        % Fissura perpendicular a s1 (plano de fratura abre na direcao de tracao maxima)
        if todosFlag1(k) > 0
            w   = todosW1(k);
            cor = magnitudePara_Cor(w, wmax, cmapFiss);
            lw  = 0.8 + 2.5*(w/wmax);
            desenhaTracoFissura(xk, yk, theta_k + pi/2, comprTraco, cor, lw, ax6);
        end

        % Fissura perpendicular a s2
        if todosFlag2(k) > 0
            w   = todosW2(k);
            cor = magnitudePara_Cor(w, wmax, cmapFiss);
            lw  = 0.8 + 2.5*(w/wmax);
            desenhaTracoFissura(xk, yk, theta_k, comprTraco, cor, lw, ax6);
        end
    end

    colormap(ax6, cmapFiss);
    axes(ax6); caxis([0, wmax]);
    cb6 = colorbar(ax6);
    cb6.Label.String = 'Abertura de Fissura  w  (m)';
    cb6.Label.FontSize = 10; cb6.Label.FontWeight = 'bold';

    hL1 = line(ax6, NaN, NaN,'Color',[0.9 0.5 0.1],'LineWidth',1.2);
    hL2 = line(ax6, NaN, NaN,'Color',[0.7 0.0 0.0],'LineWidth',3.0);
    legend(ax6,[hL1 hL2],{'Fissura ativa (softening)','Fissura saturada'},...
           'Location','northeast','FontSize',9);

    nFiss = sum(todosFlag1>0) + sum(todosFlag2>0);
    nSat  = sum(todosFlag1==2) + sum(todosFlag2==2);
    fprintf('[Fissuras] Pts crackdos: %d  |  Pts saturados: %d\n', nFiss, nSat);
end

xlabel(ax6,'Comprimento (m)','FontSize',11,'FontWeight','bold');
ylabel(ax6,'Altura (m)','FontSize',11,'FontWeight','bold');
title(ax6,'Mapa de Fissuras','FontSize',13,'FontWeight','bold');
grid(ax6,'on'); ax6.GridAlpha = 0.08; box(ax6,'on'); hold(ax6,'off');

fprintf('\n[DesenhaTensoesprincipais] Tensao maxima = %.4e Pa\n', tensaomax);
fprintf('[DesenhaTensoesprincipais] Figuras 3, 4, 5 e 6 geradas.\n');

end % function principal


%% =========================================================================
%  SUBFUNCOES LOCAIS
% =========================================================================

function desenhaSetas(sprinc, coordCentro, escalaVetor, cTrac, cComp, lw, ax)
% Setas classicas: cor pelo sinal de cada tensao principal individualmente
    s1    = sprinc(1);
    s2    = sprinc(2);
    theta = sprinc(4);

    v1 = criaSeta(s1);
    v2 = criaSeta(s2);

    if s1 >= 0, cor1 = cTrac; else cor1 = cComp; end
    if s2 >= 0, cor2 = cTrac; else cor2 = cComp; end

    plotaSeta(v1, theta,        escalaVetor, coordCentro, cor1, lw, ax);
    plotaSeta(v1, theta+pi,     escalaVetor, coordCentro, cor1, lw, ax);
    plotaSeta(v2, theta+pi/2,   escalaVetor, coordCentro, cor2, lw, ax);
    plotaSeta(v2, theta+3*pi/2, escalaVetor, coordCentro, cor2, lw, ax);
end

function v = criaSeta(s)
    as = abs(s);
    vx = [0, as, as,      1.2*as, as,      as];
    vy = [0,  0, -0.1*as, 0,      0.1*as,  0];
    if s < 0, vx = -vx + 1.2*as; end
    v = [vx; vy];
end

function plotaSeta(v, theta, escalaVetor, coordCentro, cor, lw, ax)
    R  = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    u  = R * v;
    ux = u(1,:)*escalaVetor + coordCentro(1);
    uy = u(2,:)*escalaVetor + coordCentro(2);
    H  = line(ax, ux, uy);
    set(H,'Color',cor,'LineWidth',lw);
end

function desenhaTracoFissura(xc, yc, angulo, comp, cor, lw, ax)
    dx = 0.5*comp*cos(angulo);
    dy = 0.5*comp*sin(angulo);
    H  = line(ax,[xc-dx, xc+dx],[yc-dy, yc+dy]);
    set(H,'Color',cor,'LineWidth',lw);
end

function cor = magnitudePara_Cor(mag, vmax, cmap)
    n   = size(cmap,1);
    idx = round(1 + (n-1)*min(mag/vmax, 1.0));
    idx = max(1, min(n, idx));
    cor = cmap(idx,:);
end

function Zout = mascararForaDominio(Xg, Yg, Zin, coords)
    k      = convhull(coords(:,1), coords(:,2));
    dentro = inpolygon(Xg(:), Yg(:), coords(k,1), coords(k,2));
    Zout   = Zin;
    Zout(~reshape(dentro, size(Xg))) = NaN;
end

function cmap = buildDivergentColormap(n)
% Azul escuro -> Branco -> Vermelho escuro
    half = floor(n/2);
    r1 = linspace(0.09,1.00,half)';  g1 = linspace(0.30,1.00,half)';  b1 = linspace(0.85,1.00,half)';
    r2 = linspace(1.00,0.80,n-half)'; g2 = linspace(1.00,0.08,n-half)'; b2 = linspace(1.00,0.08,n-half)';
    cmap = [[r1;r2],[g1;g2],[b1;b2]];
end

function cmap = buildFissuraColormap(n)
% Amarelo claro -> Laranja -> Vermelho escuro
    r = linspace(1.00, 0.65, n)';
    g = linspace(0.95, 0.00, n)';
    b = linspace(0.20, 0.00, n)';
    cmap = [r, g, b];
end