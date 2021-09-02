% newState = gmtOlsa_updateUserState(nCorrect, state, par)
%
% Update the internal state of the Olsa test given the number of number of 
% correct words in the current trial, the currect state and some static
% parameters of the test. If the test is adaptive, the change in SNR is:
%   dSNR = maxAbs( dSNR_raw, -deltaCorrect * par.minAdjust )
%     with       dSNR_raw = - deltaCorrect * par.fi_scale * par.fi_base^(-nReversals) / par.slope
%     and    deltaCorrect = (fractionCorrect - par.target)
%
% For par.minStep=0, the procedure simplifies to the one described in 
% Brand&Kollmeier (JASA, 2002). See paper for details.
% 
% INPUT:
%    nCorrect - # of correct words in the current trial
%    state - current state of the Olsa run; struct containing:
%               state.speechLvl - current speech level
%               state.noiseLvL  - current noise level
%               state.nReversals - current number of reversals
%               state.snrDirection - current direction of change in SNR (-1 / 0 / +1)
%    par - parameters controlling the adaptive behavior; struct containing:
%               par.adaptMode - adaptation mode: 'static' / 'speech' / 'noise'
%               par.target    - targeted fraction of correct responses;
%                               (0 < par.target < 1)
%               par.slope     - slope determining the rate of change [1/dB]
%               par.fi_scale  - scale for determining the rate of change 
%               par.fi_base   - base for determining the rate of change
%               par.minStep - min. step size for a deviation from target of 1.0
% OUTPUT:
%   newState - update state (same format as above)
%
% Change log:
%   10 Dec 2012, P.Hehrmann - created
%   02 Aug 2013, PH - added minimum step-size functionality (par.minAdjust), actually called minStep?
% 
function newState = matrix_updateUserState(nCorrect, state, par)
    
    % parAdapt.slope = 0.15;
    % parAdapt.fi_scale = 1.5;
    % parAdapt.fi_base = 1.4;

    % parameter checks
    assert(all(isfield(state,{'speechLevel_dB','noiseLevel_dB','nReversals','snrDirection'})),...
            'Illegal argument: check help for correct format of argument ''state''');
    assert(all(isfield(par,{'target','slope','fi_base','fi_scale'})), ...
            'Illegal argument: check help for correct format of argument ''par''');
    assert(par.target > 0 || par.target < 1, 'par.target must be in (0, 1).');
    assert(par.slope > 0, 'par.target must be positive.');
    assert(par.fi_base > 0, 'par.fi_base must be positive.');
    assert(par.fi_scale > 0, 'par.fi_scale must be positive.');        

    % if par.minStep is undefined -> set to 0  (for bw compatibility)
    if ~isfield(par, 'minStep')
        par.minStep = 0;
    end
    assert(par.minStep >= 0, 'par.minStep must be non-negative.');
       
    fracCorrect = nCorrect / 5;
    
    % update SNR direction (up/down)
    if fracCorrect > par.target
        snrDirection = -1; % SNR is decreasing
    elseif fracCorrect < par.target
        snrDirection = +1; % SNR is increasing
    else
        snrDirection = state.snrDirection; % on target: maintain previous direction
    end
    
    % update number of reversals
    if (snrDirection == -state.snrDirection)
        nReversals = state.nReversals + 1; % reversal has occurred
    else
        nReversals = state.nReversals;
    end
    
    % compute deltaSnr_raw
    fi = par.fi_scale*par.fi_base^(-nReversals);
    deltaCorrect = fracCorrect - par.target;
    deltaSnr_raw = - fi * deltaCorrect/par.slope; % f9
    
    % increase deltaSnr to chosen minimum if neccessary
    % choose between official rule and minimum stepsize
    deltaSnr = sign(deltaSnr_raw) * max(abs([deltaSnr_raw, par.minStep*deltaCorrect]));
%    [deltaSnr_raw, par.minStep*deltaCorrect]
    
    % compute new presentation levels
    switch par.adaptMode
        case{'off', 'static'}
            speechLevel_dB = state.speechLevel_dB;
            noiseLevel_dB = state.noiseLevel_dB;
        case{'speech'}
            speechLevel_dB = state.speechLevel_dB + deltaSnr;
            noiseLevel_dB = state.noiseLevel_dB;
        case{'noise'}
            speechLevel_dB = state.speechLevel_dB;
            noiseLevel_dB = state.noiseLevel_dB - deltaSnr;
    end
    
    % set new state
    newState.speechLevel_dB = speechLevel_dB;
    newState.noiseLevel_dB = noiseLevel_dB;
    newState.nReversals = nReversals;
    newState.snrDirection = snrDirection;
end