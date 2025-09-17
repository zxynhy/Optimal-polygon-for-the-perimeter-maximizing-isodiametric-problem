function optGon(m)
% Input:
%     m — an integer m ≥ 3 with an odd factor
% Output:
%     (1) all optimal tight frames (m-vector sets) in R^2
%         that achieve the minimal condition number
%         for the stability of phase retrieval;
%     (2) all optimal convex m-polygons in R^2
%         for the perimeter-maximizing isodiametric problem

clc

% ——————— 1. Computation of all results ————————
tic
[edgeVectorSets,frameVectorSets] = allOptPolygons(m);
toc

% ——————— 2. Visualization (optimal convex polygons) ————————
S = size(edgeVectorSets,3);
cols = lines(S);
ncol = ceil(sqrt(S));
nrow = ceil(S/ncol);

figure('Name','All optimal convex m-gons','Color','w');
for k = 1:S
    E = edgeVectorSets(:,:,k)';   % m×2
    V = [0,0];
    for j = 1:m
        V(end+1,:) = V(end,:) + E(j,:);
    end
    subplot(nrow,ncol,k)
    patch(V(:,1), V(:,2), cols(k,:), 'FaceAlpha',0.3, 'EdgeColor',cols(k,:), 'LineWidth',1.5);
    axis equal off
    title(sprintf('%d', k))
end
sgtitle(sprintf('All optimal convex %d-gons (num: %d)', m, S));

% ——————— 3. Visualization (optimal tight frames) ————————
plotVectorSetsTri(frameVectorSets, 'ColorMode','single','ShowIndex',true);

end



%% ---------- Auxiliary Function 1 ----------
function [edgeVectorSets,frameVectorSets] = allOptPolygons(m,fmtDecimals)
% Find all distinct optimal convex m-gons
% and all distinct optimal tight frames (m-vector sets)
% obtained by solving Problem \ref{pr}, the partition of roots of unity problem
%
% Note: Polygons that differ only by scaling, translation, rotation, flip,
%       or their combination are considered equivalent.
%
% Input:
%     m — an integer m ≥ 3 with an odd factor
%     fmtDecimals — Tolerance for detecting zero-sum (optional, default 1e-12)
%
% Output:
%   edgeVectorSets — edge vector sets of all optimal polygons, 2×m×s matrix.
%                    Each slice (2×m matrix) corresponds to the edge vectors
%                    arranged counterclockwise of one optimal m-gon.
%
%   frameVectorSets — all optimal tight frames, 2×m×s matrix.
%                     Each slice (2×m matrix) corresponds to
%                     one optimal tight frame in R^2 with m vectors

% Input checks
if ~(isscalar(m) && m>=3 && m==floor(m))
    error('The input parameter m must be a positive integer with m ≥ 3!');
end
mm = m;
while mod(mm, 2) == 0
    mm = mm / 2;
end
if mm == 1
    error('The input parameter m must have at least one odd factor!');
end

if nargin < 2 || isempty(fmtDecimals)
    fmtDecimals = 12;
end

% Find all solutions of the partition of roots of unity problem
partition_sols = find_eqPartition(m, fmtDecimals);

% Turn into opt gons and opt frames
S = size(partition_sols,1);
fprintf('Found %d distinct optimal convex %d-gons.\n', S, m);
fprintf('Found %d distinct optimal tight frames (%d-vector sets).\n', S, m);
edgeVectorSets = zeros(2,m,S);
frameVectorSets = zeros(2,m,S);
for k = 1:S
    eps = partition_sols(k,:);
    th = zeros(1,m);
    for j=1:m
        if eps(j) == 1
            th(j) = (j-1)*pi/m;
        else
            th(j) = (j-1)*pi/m + pi;
        end
    end
    th = sort(th);
    edgeVectorSets(:,:,k) = [cos(th);sin(th)];
    frameVectorSets(:,:,k) = [cos(th/2);sin(th/2)];
end

end



%% ---------- Auxiliary Function 2 ----------
function solutions = find_eqPartition(m, fmtDecimals)
% Input:
%     m — an integer m ≥ 3
%     fmtDecimals — Tolerance for detecting zero-sum (optional, default 1e-12)
% Output:
%     All vectors eps = (eps_0,...,eps_{m-1})^T with eps_j in {+1,-1}
%     such that g(eps) = sum_{j=0}^{m-1} eps_j * zeta_m^j = 0,
%     where zeta_m = exp(1i*pi/m)
%
% Notes:
% - The code supports two modes: for small m (default m<=22) it uses brute-force
%   enumeration (2^m). For larger m it uses a meet-in-the-middle (MITM)
%   strategy to reduce time complexity to roughly O(2^{m/2}).
% - Found solutions are merged according to the equivalence relations given in
%   the paper. The canonical representative of each equivalence class is chosen 
%   as the lexicographically smallest binary string obtained by encoding +1 -> '1', -1 -> '0'.
%
% Usage example:
% partition_sols = find_eqPartition(m, fmtDecimals);
%
% Caution:
% - For large m (e.g. m>30) the search space can still be very large and may
%   require substantial time and memory. Adjust your resources accordingly.

% default parameters
if nargin < 2 || isempty(fmtDecimals)
    fmtDecimals = 12;
end

% input checks
if ~(isscalar(m) && m>=3 && m==floor(m))
    error('The input parameter m must be a positive integer with m ≥ 3!');
end
mm = m;
while mod(mm, 2) == 0
    mm = mm / 2;
end
if mm == 1
    error('The input parameter m must have at least one odd factor!');
end

zeta = exp(1i*pi/m); % zeta_m
idx = 0:(m-1);
% method selection: brute force for small m, MITM otherwise
if m <= 22
    method = 'bruteforce';
else
    method = 'mitm';
end

switch method
    case 'bruteforce'
        N = 2^m;
        sols = zeros(0, m);
        cnt = 0;
        for mask = 0:(N-1)
            bits = bitget(mask, 1:m); % 1/0
            eps = 2*bits - 1; % {+1,-1}
            s = sum(eps .* (zeta .^ idx));
            if abs(real(s)) < 10^(-fmtDecimals) && abs(imag(s)) < 10^(-fmtDecimals)
                cnt = cnt + 1;
                sols(cnt, :) = eps;
            end
        end
    case 'mitm'
        n1 = floor(m/2);
        n2 = m - n1;
        idx1 = 0:(n1-1);
        idx2 = n1:(m-1);
        N1 = 2^n1;
        N2 = 2^n2;
        map = containers.Map();
        for mask = 0:(N2-1)
            bits = bitget(mask, 1:n2);
            eps2 = 2*bits - 1; % n2
            s2 = sum(eps2 .* (zeta .^ idx2));
            key = complex_key(s2, fmtDecimals);
            if isKey(map, key)
                list = map(key);
                list{end+1} = eps2; %#ok<AGROW>
                map(key) = list;
            else
                map(key) = {eps2};
            end
        end
        sols = zeros(0, m);
        cnt = 0;
        for mask = 0:(N1-1)
            bits = bitget(mask, 1:n1);
            eps1 = 2*bits - 1; % n1
            s1 = sum(eps1 .* (zeta .^ idx1));
            target = -s1;
            keyt = complex_key(target, fmtDecimals);
            if isKey(map, keyt)
                list = map(keyt);
                for k = 1:numel(list)
                    eps2 = list{k};
                    eps = [eps1, eps2];
                    cnt = cnt + 1;
                    sols(cnt, :) = eps;
                end
            end
        end
    otherwise
        error('Unknown method');
end

% fprintf('Found %d candidate sequences (before merging equivalence classes)\n', size(sols,1));

% Merge solutions according to the paper's equivalence relations and keep
% one canonical representative per equivalence class
canonical_map = containers.Map();
for i = 1:size(sols,1)
    eps = sols(i,:);
    crep = canonical_rep(eps);
    if ~isKey(canonical_map, crep)
        canonical_map(crep) = eps; % store one representative for this class
    end
end
keysList = keys(canonical_map);
S = numel(keysList);
solutions = zeros(S, m);
for i = 1:S
    eps_temp = canonical_map(keysList{i});
    if eps_temp(1) == -1
        eps_temp = -eps_temp;
    end
    solutions(i,:) = eps_temp;
end
% fprintf('After merging, %d equivalence classes remain\n', S);
end


function key = complex_key(c, decimals)
fmt = sprintf('%%.%df_%%.%df', decimals, decimals);
key = sprintf(fmt, real(c), imag(c));
end

function crep = canonical_rep(eps)
% Compute a canonical representative string for the equivalence class of
% eps under the operations: shift S_k (including k=0 as "no shift"),
% reversal R, and global sign flip alpha ∈ {−1,+1}.
% The canonical representative is the lexicographically smallest binary
% string obtained by encoding +1 -> '1' and -1 -> '0'.
m = numel(eps);
bestStr = '';
ks = [0, 1:(m-1)];
for k = ks
    for applyR = 0:1
        for alpha = [-1, 1]
            e = eps;
            if k > 0
                e = S_k(e, k);
            end
            if applyR == 1
                e = R_op(e);
            end
            e = alpha * e;
            s = vec_to_binstr(e);
            if isempty(bestStr)
                bestStr = s;
            else
                % lexicographic comparison
                C = sort({bestStr, s});
                if strcmp(C{1}, s)
                    bestStr = s;
                end
            end
        end
    end
end
crep = bestStr;
end

function s = vec_to_binstr(e)
s = repmat('0', 1, numel(e));
s(e==1) = '1';
end

function y = R_op(x)
y = x(end:-1:1);
end

function y = S_k(x, k)
m = numel(x);
if k < 0 || k >= m
    error('In S_k, k must satisfy 1 <= k <= m-1');
end
y = zeros(size(x));
y(1:(m-k)) = x((1+k):m);
y((m-k+1):m) = -x(1:k);
end



%% ---------- Auxiliary Function 3 ----------
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
