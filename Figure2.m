clear; clc; close all;

 

% ---------- Styling ----------

set(groot,'defaultTextInterpreter','latex');

set(groot,'defaultLegendInterpreter','latex');

set(groot,'defaultAxesTickLabelInterpreter','latex');

set(groot,'defaultAxesFontSize',25);

set(groot,'defaultTextFontSize',25);

set(groot,'defaultLineLineWidth',2.5);

set(groot,'defaultAxesLineWidth',1.2);

set(groot,'defaultAxesTickDir', 'out');

set(groot,'defaultAxesBox', 'off')
 

% ---------- Axes grid ----------

lamR = linspace(0, 0.6, 400);   % x-axis: \hat{\lambda}_R

f    = linspace(0, 1.0, 400);   % y-axis: \hat f

[LamR,F] = meshgrid(lamR,f);

 

% ---------- Choose delta overlays ----------

deltas = [0.40, 0.50, 0.60];  % add/remove values here

colors = lines(numel(deltas));

 

% ---------- Base figure ----------

fig = figure('Color','w','Position',[100 100 900 650]); hold on; box on;

 

% Soft background shading for a representative delta (first one)

deltaShade = deltas(1);

fcritShade = (deltaShade - LamR) ./ (2*deltaShade);

regionGrowth = F > fcritShade;

regionExt    = F < fcritShade;

 

% Light shading — keep subtle so lines pop

hExt = contourf(LamR, F, double(regionExt), [1 1], ...
    'FaceAlpha', 0.10, 'LineStyle','none','HandleVisibility','off');
hGro = contourf(LamR, F, double(regionGrowth), [1 1], ...
    'FaceAlpha', 0.10, 'LineStyle','none','HandleVisibility','off');

 

% ---------- Critical lines for each delta ----------

legLines = gobjects(0);

legLabs  = {};

 

for i = 1:numel(deltas)

    d = deltas(i);

    fcrit = (d - lamR) ./ (2*d);

    % clip to [0,1] for plotting prettiness

    fcrit = max(0, min(1, fcrit));

    h = plot(lamR, fcrit, '-', 'Color', colors(i,:), 'DisplayName', sprintf('$\\delta=%.2g$', d));

    legLines(end+1) = h; %#ok<*SAGROW>

    legLabs{end+1}  = get(h,'DisplayName');

end


% ---------- Labels and legend ----------

xlabel('$\hat{\lambda}_R$');

ylabel('$\hat{f}$');

title('Bifurcation diagram in the $(\hat{\lambda}_R,\hat f)$ plane');

 

% Create legend: delta lines + markers

legend('Location','northeastoutside'); grid on; ylim([0 1]); xlim([0 0.6]);

 

% ---------- Annotations ----------

text(0.05, 0.08, 'Extinction region', 'Color',[0.2 0.2 0.2]);

text(0.3, 0.80, 'Growth region',     'Color',[0.2 0.2 0.2]);

 

% ---------- Export ----------

exportgraphics(fig, 'bifurcation_lambdaR_f.svg', 'Resolution', 300);

exportgraphics(fig, 'bifurcation_lambdaR_f.pdf', 'Resolution', 300);