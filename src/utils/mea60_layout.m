function layout = mea60_layout()
%MEA60_LAYOUT Canonical 60-electrode MEA geometry as a struct of tables.
%
%   layout = MEA60_LAYOUT() returns a struct describing the Multi Channel
%   Systems 60MEA100/10iR 8x8 electrode grid (four corners absent). This
%   is the single source of truth for the channel-index to grid-position
%   mapping used by topographical_map.m, plot_network_on_mea.m, and any
%   future spatial analysis.
%
% INPUTS:
%   (none)
%
% OUTPUTS:
%   layout  -  Struct with fields:
%       .nCh               (scalar, 60)
%       .gridSize          ([8 8])
%       .channelList       (1 x 60) TDT channel indices 1..60
%       .gridRow           (1 x 60) grid row for each channel index
%       .gridCol           (1 x 60) grid column for each channel index
%       .xCoord            (1 x 60) Cartesian x coordinate (= gridCol)
%       .yCoord            (1 x 60) Cartesian y coordinate (= 9 - gridRow)
%                                   so that row 1 plots at the top
%       .chToRC            containers.Map: chIdx -> [row col]
%       .mcsLabels         (1 x 60) original MCS electrode labels
%                          (rowDigit*10 + colDigit)
%       .distanceMatrix    (60 x 60) pairwise Euclidean distance in
%                          electrode pitches (1 unit = 100 um)
%
% Notes:
%   - The electrode-label table was extracted from
%     src/archive/other_scripts/CartesianPlot_histo_Rasters_Dataset.m and
%     cross-checked against docs/60StandardMEA_Layout.pdf.
%   - Channel indices are 1-based (TDT convention). Grid rows/cols are
%     1..8 (MCS label convention: rowDigit is the tens digit, colDigit is
%     the units digit).

    chList = 1:60;
    elecLabels = [ ...
        47, 48, 46, 45, 38, 37, 28, 36, 27, 17, ...
        26, 16, 35, 25, 15, 14, 24, 34, 13, 23, ...
        12, 22, 33, 21, 32, 31, 44, 43, 41, 42, ...
        52, 51, 53, 54, 61, 62, 71, 63, 72, 82, ...
        73, 83, 64, 74, 84, 85, 75, 65, 86, 76, ...
        87, 77, 66, 78, 67, 68, 55, 56, 58, 57];

    rows = floor(elecLabels / 10);
    cols = mod(elecLabels, 10);

    layout = struct();
    layout.nCh         = numel(chList);
    layout.gridSize    = [8, 8];
    layout.channelList = chList;
    layout.gridRow     = rows;
    layout.gridCol     = cols;

    % Cartesian: x increases with column, y flipped so row 1 is at top.
    layout.xCoord = cols;
    layout.yCoord = 9 - rows;

    layout.chToRC = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for k = 1:numel(chList)
        layout.chToRC(int32(chList(k))) = [rows(k), cols(k)];
    end

    layout.mcsLabels = elecLabels;

    % Pairwise distances on the grid (units = electrode pitch).
    dx = layout.xCoord(:) - layout.xCoord(:).';
    dy = layout.yCoord(:) - layout.yCoord(:).';
    layout.distanceMatrix = sqrt(dx.^2 + dy.^2);
end
