function validate_channel_mapping()
%VALIDATE_CHANNEL_MAPPING  Automated check of the 60MEA channel mapping.
%
%   VALIDATE_CHANNEL_MAPPING() runs a battery of checks against the
%   mapping defined in mea60_layout.m and the MCS 60StandardMEA datasheet
%   (docs/60StandardMEA_Layout.pdf, page 2).
%
%   Prints PASS/FAIL for each check.  Any failure indicates a mapping
%   inconsistency that must be resolved before spatial analyses are valid.
%
% See also: MEA60_LAYOUT, PROJECT_CONFIG.

    fprintf('\n=== MEA Channel Mapping Validation ===\n\n');
    nPass = 0;
    nFail = 0;

    layout = mea60_layout();

    % --- Check 1: 60 channels ---
    [nPass, nFail] = check(layout.nCh == 60, ...
        '60 signal channels', nPass, nFail);

    % --- Check 2: unique electrode labels ---
    [nPass, nFail] = check(numel(unique(layout.mcsLabels)) == 60, ...
        'All electrode labels are unique', nPass, nFail);

    % --- Check 3: PZ5 channel list ---
    expectedCh = [1:15, 17:31, 33:47, 49:63];
    [nPass, nFail] = check(isequal(layout.channelList, expectedCh), ...
        'PZ5 channel list matches [1:15, 17:31, 33:47, 49:63]', nPass, nFail);

    % --- Check 4: absent corners (11, 18, 81, 88) not in labels ---
    absentLabels = [11, 18, 81, 88];
    [nPass, nFail] = check(~any(ismember(absentLabels, layout.mcsLabels)), ...
        'Absent corners (11, 18, 81, 88) not in electrode labels', nPass, nFail);

    % --- Check 5: reference electrode (15) is present ---
    [nPass, nFail] = check(any(layout.mcsLabels == 15), ...
        'Reference electrode (label 15) is present', nPass, nFail);

    % --- Check 6: all labels have valid digits ---
    tens = floor(layout.mcsLabels / 10);
    ones = mod(layout.mcsLabels, 10);
    [nPass, nFail] = check(all(tens >= 1 & tens <= 8) && all(ones >= 1 & ones <= 8), ...
        'All labels have digits in range 1-8', nPass, nFail);

    % --- Check 7: verify pin-to-electrode mapping against MCS datasheet ---
    % Complete mapping from MCS 60StandardMEA_Layout.pdf page 2:
    %   Pin k -> electrode label mcsRef(k)
    mcsRef = [ ...
        47, 48, 46, 45, 38, 37, 28, 36, 27, 17, ...   % pins  1-10
        26, 16, 35, 25, 15, 14, 24, 34, 13, 23, ...   % pins 11-20
        12, 22, 33, 21, 32, 31, 44, 43, 41, 42, ...   % pins 21-30
        52, 51, 53, 54, 61, 62, 71, 63, 72, 82, ...   % pins 31-40
        73, 83, 64, 74, 84, 85, 75, 65, 86, 76, ...   % pins 41-50
        87, 77, 66, 78, 67, 68, 55, 56, 58, 57];      % pins 51-60
    [nPass, nFail] = check(isequal(layout.mcsLabels, mcsRef), ...
        'All 60 pin-to-electrode entries match MCS datasheet', nPass, nFail);

    % --- Check 8: distance matrix symmetry ---
    D = layout.distanceMatrix;
    [nPass, nFail] = check(max(abs(D - D.'), [], 'all') < 1e-12, ...
        'Distance matrix is symmetric', nPass, nFail);

    % --- Check 9: distance matrix diagonal is zero ---
    [nPass, nFail] = check(all(diag(D) == 0), ...
        'Distance matrix diagonal is zero', nPass, nFail);

    % --- Check 10: known distance spot-checks ---
    % Electrodes 44 and 45 are adjacent (same column, adjacent rows)
    % so their distance should be 1.0 electrode pitch.
    idx44 = find(layout.mcsLabels == 44);
    idx45 = find(layout.mcsLabels == 45);
    d_44_45 = D(idx44, idx45);
    [nPass, nFail] = check(abs(d_44_45 - 1.0) < 1e-10, ...
        sprintf('Distance(44, 45) = %.4f (expected 1.0)', d_44_45), nPass, nFail);

    % Electrodes 44 and 54 are adjacent (same row, adjacent columns)
    idx54 = find(layout.mcsLabels == 54);
    d_44_54 = D(idx44, idx54);
    [nPass, nFail] = check(abs(d_44_54 - 1.0) < 1e-10, ...
        sprintf('Distance(44, 54) = %.4f (expected 1.0)', d_44_54), nPass, nFail);

    % Electrodes 44 and 55 are diagonal neighbours -> distance = sqrt(2)
    idx55 = find(layout.mcsLabels == 55);
    d_44_55 = D(idx44, idx55);
    [nPass, nFail] = check(abs(d_44_55 - sqrt(2)) < 1e-10, ...
        sprintf('Distance(44, 55) = %.4f (expected %.4f)', d_44_55, sqrt(2)), nPass, nFail);

    % Electrodes 21 and 78 are at opposite corners of the active area
    % col diff = |7-2| = 5, row diff = |8-1| = 7, dist = sqrt(74)
    idx21 = find(layout.mcsLabels == 21);
    idx78 = find(layout.mcsLabels == 78);
    d_21_78 = D(idx21, idx78);
    expected_21_78 = sqrt(5^2 + 7^2);
    [nPass, nFail] = check(abs(d_21_78 - expected_21_78) < 1e-10, ...
        sprintf('Distance(21, 78) = %.4f (expected %.4f)', d_21_78, expected_21_78), nPass, nFail);

    % --- Check 11: cross-validate with project_config ---
    cfg = project_config();
    [nPass, nFail] = check(isequal(layout.channelList, cfg.channels.default), ...
        'layout.channelList == cfg.channels.default', nPass, nFail);
    [nPass, nFail] = check(cfg.channels.reference == 15, ...
        'cfg.channels.reference == 15', nPass, nFail);

    % --- Check 12: legacy mapping cross-check ---
    % Spot-check a few entries from CartesianPlot_histo_Rasters_Dataset.m
    legacy = containers.Map( ...
        {1, 15, 24, 37, 52, 60}, ...
        {47, 15, 21, 71, 77, 57});
    legacyKeys = keys(legacy);
    allMatch = true;
    for k = 1:numel(legacyKeys)
        pin = legacyKeys{k};
        if layout.mcsLabels(pin) ~= legacy(pin)
            allMatch = false;
        end
    end
    [nPass, nFail] = check(allMatch, ...
        'Spot-check entries match legacy CartesianPlot mapping', nPass, nFail);

    % --- Summary ---
    fprintf('\n--- Summary ---\n');
    fprintf('Passed: %d / %d\n', nPass, nPass + nFail);
    if nFail > 0
        fprintf(2, 'FAILED: %d checks — review mapping before running spatial analyses.\n', nFail);
    else
        fprintf('All checks passed. Channel mapping is consistent.\n');
    end
    fprintf('\n');
end


function [nPass, nFail] = check(condition, label, nPass, nFail)
    if condition
        fprintf('  [PASS] %s\n', label);
        nPass = nPass + 1;
    else
        fprintf(2, '  [FAIL] %s\n', label);
        nFail = nFail + 1;
    end
end
