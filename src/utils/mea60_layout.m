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
%       .channelList       (1 x 60) PZ5 channel indices (physical numbering,
%                          skipping ground channels 16,32,48,64)
%       .gridRow           (1 x 60) grid row for each channel
%       .gridCol           (1 x 60) grid column for each channel
%       .xCoord            (1 x 60) Cartesian x coordinate (= gridCol)
%       .yCoord            (1 x 60) Cartesian y coordinate (= 9 - gridRow)
%                                   so that row 1 plots at the top
%       .chToRC            containers.Map: PZ5 chIdx -> [row col]
%       .mcsLabels         (1 x 60) original MCS electrode labels
%                          (rowDigit*10 + colDigit)
%       .mcsPins           (1 x 60) MCS 1-dim pin numbers (1-60)
%       .distanceMatrix    (60 x 60) pairwise Euclidean distance in
%                          electrode pitches (1 unit = 100 um)
%
% Notes:
%   - The electrode-label table was extracted from
%     src/archive/other_scripts/CartesianPlot_histo_Rasters_Dataset.m and
%     cross-checked against docs/60StandardMEA_Layout.pdf.
%   - The TDT PZ5 digitizer uses 4 banks of 16 channels.  The 16th
%     channel of each bank is ground (PZ5 ch 16,32,48,64).  The 60
%     signal channels occupy PZ5 ch [1:15, 17:31, 33:47, 49:63],
%     which map sequentially to MCS 1-dim pins 1-60.
%   - MCS labels: tens digit = column, units digit = row.
%     (MCS convention: label = col*10 + row.  E.g. electrode 23 is
%      column 2, row 3.  See docs/60StandardMEA_Layout.pdf page 2.)
%   - On iR-model MEAs (60MEA100/10iR), MCS pin 15 (PZ5 ch 15) is the
%     internal reference electrode — a large reference pad, not a
%     recording electrode (59 recording + 1 ref).
%     See docs/60StandardMEA_Layout.pdf and cfg.channels.reference.

    % PZ5 physical channel numbers (60 signal channels, ground skipped).
    chList = [1:15, 17:31, 33:47, 49:63];

    % MCS 1-dim pin numbers (sequential 1-60).
    mcsPins = 1:60;

    % MCS electrode labels indexed by 1-dim pin number.
    % Pin 1 -> electrode 47, pin 2 -> electrode 48, etc.
    elecLabels = [ ...
        47, 48, 46, 45, 38, 37, 28, 36, 27, 17, ...
        26, 16, 35, 25, 15, 14, 24, 34, 13, 23, ...
        12, 22, 33, 21, 32, 31, 44, 43, 41, 42, ...
        52, 51, 53, 54, 61, 62, 71, 63, 72, 82, ...
        73, 83, 64, 74, 84, 85, 75, 65, 86, 76, ...
        87, 77, 66, 78, 67, 68, 55, 56, 58, 57];

    cols = floor(elecLabels / 10);   % tens digit = column (MCS convention)
    rows = mod(elecLabels, 10);     % ones digit = row    (MCS convention)

    layout = struct();
    layout.nCh         = numel(chList);
    layout.gridSize    = [8, 8];
    layout.channelList = chList;
    layout.gridRow     = rows;
    layout.gridCol     = cols;
    layout.mcsPins     = mcsPins;

    % Cartesian: x increases with column, y flipped so row 1 is at top.
    layout.xCoord = cols;
    layout.yCoord = 9 - rows;

    % Map from PZ5 channel index to [row, col].
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
