function step = nice_tick_step(span, targetTickCount)
%NICE_TICK_STEP Choose a 1/2/5 * 10^n tick step for readable axes.
%
%   step = NICE_TICK_STEP(span, targetTickCount) returns a "nice" axis tick
%   step from the {1, 2, 5} x 10^n family that produces approximately
%   targetTickCount ticks across the requested span. This is the same
%   helper used by every figure script in the proof-of-concept pipeline.
%
% INPUTS:
%   span             -  Total axis span (>0). For asymmetric % axes use
%                       (yMax - yMin), otherwise just yMax.
%   targetTickCount  -  Approximate number of major ticks desired (>=2).
%
% OUTPUTS:
%   step  -  Scalar tick step in data units.

    if span <= 0 || targetTickCount < 2
        step = 1;
        return;
    end
    rawStep = span / (targetTickCount - 1);
    mag  = 10^floor(log10(rawStep));
    frac = rawStep / mag;
    if frac <= 1
        niceFrac = 1;
    elseif frac <= 2
        niceFrac = 2;
    elseif frac <= 5
        niceFrac = 5;
    else
        niceFrac = 10;
    end
    step = niceFrac * mag;
end
