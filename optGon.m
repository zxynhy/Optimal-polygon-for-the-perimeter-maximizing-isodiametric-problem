function optGon(m)
% Input:
%     m — an integer m ≥ 3 with an odd factor
% Output:
%     all tight frames (m-vector sets)
%     that achieve the optimal phase retrieval stability constant, as well as the 
%     corresponding optimal convex m-polygons

clc

% ——————— 1. Computation of all results ————————
tic
% angleSets: the counterclockwise angles of all optimal polygons, 
%            s×m matrix. Each row corresponds to the set of 
%            the counterclockwise angles of one optimal polygon, 
%            from e_1 to e_i, i=1,...,m.
%
% vectorSets: all optimal tight frames, 
%             2×m×s matrix. Each slice (a 2×m matrix) corresponds to 
%             one optimal tight frame.
[angleSets,vectorSets] = allOptPolygons(m);
toc

% ——————— 2. Visualization (optimal convex polygons) ————————
s = size(angleSets,1);
cols = lines(s);
ncol = ceil(sqrt(s));
nrow = ceil(s/ncol);

figure('Name','All optimal convex m-gons','Color','w');
for k = 1:s
    th = angleSets(k,:);
    E = [cos(th(:)), sin(th(:))];   % m×2
    V = [0,0];
    for j = 1:m
        V(end+1,:) = V(end,:) + E(j,:);
    end
    subplot(nrow,ncol,k)
    patch(V(:,1), V(:,2), cols(k,:), 'FaceAlpha',0.3, 'EdgeColor',cols(k,:), 'LineWidth',1.5);
    axis equal off
    title(sprintf('%d', k))
end
sgtitle(sprintf('All optimal convex %d-gons (num: %d)', m, s));

% ——————— 3. Visualization (optimal tight frames) ————————
plotVectorSetsTri(vectorSets, 'ColorMode','single','ShowIndex',true);

end


%% Function 0
function K = giveAnOptIntlist(m)
% Input:
%     m — an integer m ≥ 3 with an odd factor
% Output:
%     K — an integer list K = [k_1,...,k_r] corresponding to a certain 
%         optimal convex m-gon, where r denotes the smallest odd 
%         prime factor of m, and k_1 = ... = k_r = m / r


if m < 3 || floor(m) ~= m
    error('The input parameter m must be a positive integer with m ≥ 3!');
end

mm = m;
while mod(mm, 2) == 0
    mm = mm / 2;
end

if mm == 1
    error('The input parameter m must admit at least one odd factor!');
end

r = mm;
for k = 3:2:floor(sqrt(mm))
    if mod(mm, k) == 0
        r = k;
    end
end
K = ones(1,r)*m/r;
end


%% Function 1
function angInt = giveAngles(m,K)
% Input:
%     m — an integer m ≥ 3 with an odd factor
%     K — an integer list K = [k_1,...,k_r] corresponding to an 
%         optimal convex m-gon
% Output:
%     angInt — 1×m matrix, the counterclockwise angles of the 
%              optimal convex polygon associated with the list K,
%              from e_1 to e_i, i=1,...,m, normalized by (π/m)

r = length(K);
angInt = zeros(1,m);
angInt(1:K(1)) = 0:(-1+K(1));
angInt((K(1)+1):sum(K(1:2))) = (m+K(1)):(m-1+sum(K(1:2)));
for i = 2:(r-1)/2
    angInt((sum(K(1:(2*i-2)))+1):sum(K(1:(2*i-1)))) = ...
        sum(K(1:(2*i-2))):(-1+sum(K(1:(2*i-1))));
    angInt((sum(K(1:(2*i-1)))+1):sum(K(1:(2*i)))) = ...
        (m+sum(K(1:(2*i-1)))):(m-1+sum(K(1:(2*i))));
end
angInt((sum(K(1:(r-1)))+1):sum(K(1:r))) = ...
    sum(K(1:(r-1))):(-1+sum(K(1:r)));
% angInt = angInt*pi/m;
angInt = sort(angInt);
end


%% Function 2 (main)
function [angleSets,vectorSets] = allOptPolygons(m, tol)
% Find all distinct optimal convex m-gons obtained by "reversing"
% all edge vectors in a zero-sum edge vector subset of 
% a given optimal convex m-gon (see Lemma 4.2)
%
% Note: Polygons that differ only by translations, rotations, flips, 
%   or their combination are considered equivalent.
%
% Input:
%     m — an integer m ≥ 3 with an odd factor
%     tol — Tolerance for detecting a zero-sum vector set (optional, default 1e-8)
%
% Output:
%   angleSets  — the counterclockwise angles of all optimal polygons, 
%            s×m matrix. Each row corresponds to the set of 
%            the counterclockwise angles of one optimal polygon, 
%            from e_1 to e_i, i=1,...,m.
%   vectorSets — all optimal tight frames, 
%             2×m×s matrix. Each slice (a 2×m matrix) corresponds to 
%             one optimal tight frame.

if nargin < 2
    tol = 1e-8;
end

%% 0) Generate an initial optimal convex m-gon, normalize, add it to the result list
K = giveAnOptIntlist(m);
theta0Int = giveAngles(m,K);
theta0Int = mod(theta0Int - theta0Int(1), 2*m);
theta0Int = sort(theta0Int);
thetaIntSets = theta0Int;

%% 1) Construct the edge vectors E_j = [cosθ_j; sinθ_j]
E = [cos(theta0Int(:)*pi/m), sin(theta0Int(:)*pi/m)];  % m×2

%% 2) Enumerate all subsets I (excluding ∅ and the full set)
AI = 1:m;
for mask = 1:(2^m - 2)
    I = logical(bitget(mask, 1:m));
    % Verify if the sum of the subset equals the zero vector
    if norm(sum(E(I, :), 1)) < tol
        AII = AI(I);
        fprintf('[');fprintf('%d,', AII(1:end-1));
        fprintf('%d', AII(end));fprintf(']');fprintf('\n');
        % reversing all edge vectors in I
        theta = theta0Int;
        theta(I) = mod(theta0Int(I) + m, 2*m);
        theta = mod(theta - theta(1), 2*m);
        theta = sort(theta);
        % Remove polygons that are equivalent under 
        % translations, rotations, flips, or any combination thereof
        flag = 1;
        for j=1:size(thetaIntSets,1)
            for i = 1:m
                theta_C = mod(theta - theta(i), 2*m);
                theta_C = sort(theta_C);
                theta_Cflip = mod(-theta_C, 2*m);
                theta_Cflip = sort(theta_Cflip);
                if (norm(theta_C-thetaIntSets(j,:), 1) < tol)...
                        || (norm(theta_Cflip-thetaIntSets(j,:), 1) < tol)
                    flag = 0;
                    break;
                end
            end
            if flag == 0
                break;
            end
        end
        if flag == 1
            thetaIntSets = [thetaIntSets; theta];
        end
    end
end
%% 3) Show results
angleSets = thetaIntSets*pi/m;
s = size(angleSets,1);
% clc;
fprintf('Found %d distinct optimal convex %d-gons.\n', s, m);
fprintf('Found %d distinct optimal %d-vector sets (tight frames).\n', s, m);
vectorSets = zeros(2,m,s);
for i = 1:s
    vectorSets(:,:,i) = [cos(angleSets(i,:)/2);sin(angleSets(i,:)/2)];
end
end


%% Function 3
function plotVectorSetsTri(vectorSets, varargin)
% Plot a 2×m×s vector set with arrows represented 
% by line segments and small triangular heads
%
% plotVectorSetsTri(vectorSets)
% plotVectorSetsTri(vectorSets, 'AxisUnified', true, 'HeadLength', 0.12, ...)
%
% Input:
%   vectorSets - 2 x m x s array, where each 2×m matrix represents one set of vectors (x; y)
%
% Optional Name/Value Parameters:
%   'AxisUnified'     - logical, whether to unify axes across all subplots (default: true)
%   'Decimals'        - number of decimal places for coordinate display (default: 3)
%   'ShowIndex'       - logical, whether to display only indices instead of vector values (default: false)
%   'PlotScale'       - scale factor for vector lengths during plotting 
%                       (does not affect displayed labels) (default: 1)
%   'LineWidth'       - width of the line segments (default: 1.2)
%   'HeadLength'      - arrowhead length as a fraction of vector length (0-1) (default: 0.08)
%   'HeadWidth'       - arrowhead width as a fraction of vector length (0-1) (default: 0.06)
%   'FillHeads'       - logical, whether to fill the arrowhead triangles (default: true)
%   'MarkerOriginSize'- size of the origin marker (default: 6)
%   'ColorMode'       - 'byIndex' (each vector has a different color) or 
%                       'single' (all vectors same color) (default: 'byIndex')
%   'BaseColor'       - RGB triplet or MATLAB color used when ColorMode='single' 
%                       (default: [0 0.4470 0.7410])
%
% Example:
%   plotVectorSetsTri(vectorSets, 'HeadLength',0.10, 'HeadWidth',0.06, 'PlotScale',0.85);
%
% Notes:
%   - HeadLength and HeadWidth are relative 
%     to the vector length (shorter vectors have proportionally shorter heads)
%   - Very short vectors are automatically shortened 
%     or represented as points to avoid distorted triangles
%
% ---------------- Parameter Parsing ----------------
p = inputParser;
addParameter(p,'AxisUnified', true, @(x) islogical(x) || (x==0) || (x==1));
addParameter(p,'Decimals', 3, @(x) isnumeric(x) && x>=0 && mod(x,1)==0);
addParameter(p,'ShowIndex', false, @(x) islogical(x) || (x==0) || (x==1));
addParameter(p,'PlotScale', 1, @(x) isnumeric(x) && x>0);
addParameter(p,'LineWidth', 1.2, @isnumeric);
addParameter(p,'HeadLength', 0.08, @(x) isnumeric(x) && x>=0 && x<1);
addParameter(p,'HeadWidth', 0.06, @(x) isnumeric(x) && x>=0 && x<1);
addParameter(p,'FillHeads', true, @(x) islogical(x) || (x==0) || (x==1));
addParameter(p,'MarkerOriginSize', 6, @isnumeric);
addParameter(p,'ColorMode', 'byIndex', @(x) ischar(x) || isstring(x));
addParameter(p,'BaseColor', [0 0.4470 0.7410]);
parse(p, varargin{:});
opts = p.Results;

% ---------------- Input Validation ----------------
if ndims(vectorSets) < 2 || size(vectorSets,1) ~= 2
    error('vectorSets must be a 2 x m x s array.');
end
m = size(vectorSets,2);
s = size(vectorSets,3);

allX = reshape(vectorSets(1,:,:),1,[]);
allY = reshape(vectorSets(2,:,:),1,[]);
xmin = min([0, allX]);
xmax = max([0, allX]);
ymin = min([0, allY]);
ymax = max([0, allY]);
if xmin == xmax
    xmin = xmin - 1; xmax = xmax + 1;
end
if ymin == ymax
    ymin = ymin - 1; ymax = ymax + 1;
end
padX = 0.08*(xmax - xmin);
padY = 0.08*(ymax - ymin);
globalXLim = [xmin - padX, xmax + padX];
globalYLim = [ymin - padY, ymax + padY];

% Color settings
if strcmpi(opts.ColorMode,'byIndex')
    colors = lines(m);
else
    colors = repmat(reshape(opts.BaseColor,1,3), m, 1);
end

% Subplot layout
nRows = ceil(sqrt(s));
nCols = ceil(s / nRows);

figure('Name','All optimal tight frames','Color','w');
tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

for idx = 1:s
    nexttile;
    X = squeeze(vectorSets(1,:,idx));
    Y = squeeze(vectorSets(2,:,idx));
    Xp = X * opts.PlotScale;
    Yp = Y * opts.PlotScale;

    hold on;
    plot(0,0,'k+','MarkerSize',opts.MarkerOriginSize,'LineWidth',1);

    for j = 1:m
        xp = Xp(j); yp = Yp(j);
        col = colors(mod(j-1,size(colors,1))+1, :);
        L = hypot(xp, yp);
        if L < eps
            plot(0,0,'.','Color',col,'MarkerSize',8);
            endpt = [0;0];
            basept = [0;0];
        else
            uv = [xp; yp] / L;
            perp = [-uv(2); uv(1)];
            headLen = opts.HeadLength * L;
            headWid = opts.HeadWidth * L;

            if headLen > 0.8*L
                headLen = 0.8*L;
            end
            if headWid < 1e-6
                headWid = 1e-6;
            end

            endpt = [xp; yp];
            basept = endpt - headLen * uv;
            plot([0, basept(1)], [0, basept(2)], '-', 'LineWidth', opts.LineWidth, 'Color', col);
            left = basept + (headWid/2) * perp;
            right = basept - (headWid/2) * perp;

            if opts.FillHeads
                patch('XData',[endpt(1), left(1), right(1)], ...
                    'YData',[endpt(2), left(2), right(2)], ...
                    'FaceColor', col, 'EdgeColor', 'none');
            else
                patch('XData',[endpt(1), left(1), right(1)], ...
                    'YData',[endpt(2), left(2), right(2)], ...
                    'FaceColor', 'none', 'EdgeColor', col, 'LineWidth', max(0.5, opts.LineWidth/1.2));
            end
        end

        if opts.ShowIndex
            lab = sprintf('%d', j);
        else
            lab = sprintf('(%.*f, %.*f)', opts.Decimals, X(j), opts.Decimals, Y(j));
        end

        ax = gca;
        xrange = diff(ax.XLim);
        yrange = diff(ax.YLim);
        dx = 0.02 * xrange;
        dy = 0.02 * yrange;

        if hypot(endpt(1), endpt(2)) < eps
            txtpos = [dx, dy];
        else
            signx = sign(endpt(1)); if signx==0, signx = 1; end
            signy = sign(endpt(2)); if signy==0, signy = 1; end
            txtpos = [endpt(1) + signx*dx, endpt(2) + signy*dy];
        end
        text(txtpos(1), txtpos(2), lab, 'FontSize', 8, 'Color', col, 'Interpreter', 'none');
    end % j loop

    if opts.AxisUnified
        xlim(globalXLim);
        ylim(globalYLim);
    else
        xmin_i = min([0; X(:)]);
        xmax_i = max([0; X(:)]);
        ymin_i = min([0; Y(:)]);
        ymax_i = max([0; Y(:)]);
        if xmin_i == xmax_i
            xmin_i = xmin_i - 1; xmax_i = xmax_i + 1;
        end
        if ymin_i == ymax_i
            ymin_i = ymin_i - 1; ymax_i = ymax_i + 1;
        end
        padXi = 0.08*(xmax_i - xmin_i);
        padYi = 0.08*(ymax_i - ymin_i);
        xlim([xmin_i-padXi, xmax_i+padXi]);
        ylim([ymin_i-padYi, ymax_i+padYi]);
    end

    axis equal;
    grid on;
    title(sprintf('Opt tight frame %d', idx));
    hold off;
end

set(gcf,'Color','w');

end

