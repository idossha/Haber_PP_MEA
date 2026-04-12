function colors = paired_plot_colors(study)
%PAIRED_PLOT_COLORS Standard color palette for paired baseline/treatment plots.
%
%   colors = PAIRED_PLOT_COLORS(study) returns a struct of RGB triplets and
%   transparency values used by the spike-rate / burst-rate figure scripts.
%   The 'doi' palette uses a blue treatment color; the 'ket' palette uses
%   the orange ketanserin color from the proof-of-concept figures.
%
% INPUTS:
%   study  -  'doi' or 'ket' (case-insensitive).
%
% OUTPUTS:
%   colors  -  Struct with fields:
%       baseline, treatment, increase, decrease     (paired-line colors)
%       baselineBox, treatmentBox                   (boxplot face colors)
%       baselineOutline, treatmentOutline           (boxplot edge colors)
%       alphaDots                                   (scalar in [0,1])

    colors.baseline         = [0.50, 0.50, 0.50];
    colors.increase         = [0.50, 0.72, 0.55];
    colors.decrease         = [0.82, 0.58, 0.58];
    colors.alphaDots        = 0.4;
    colors.baselineBox      = [0.88, 0.88, 0.88];
    colors.baselineOutline  = [0.65, 0.65, 0.65];

    switch lower(study)
        case 'doi'
            colors.treatment        = [0.25, 0.40, 0.60];
            colors.treatmentBox     = [0.573, 0.694, 0.843];
            colors.treatmentOutline = [0.35, 0.50, 0.70];
        case 'ket'
            colors.treatment        = [0.90, 0.50, 0.12];
            colors.treatmentBox     = [0.95, 0.72, 0.48];
            colors.treatmentOutline = [0.76, 0.48, 0.22];
        otherwise
            error('paired_plot_colors:UnknownStudy', ...
                'study must be one of {''doi'', ''ket''}.');
    end
end
