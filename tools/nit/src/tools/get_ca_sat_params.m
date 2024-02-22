function [dF_F_max, K_d] = get_ca_sat_params(indicator)
%%% *[dF_F_max, K_d] = get_ca_sat_params(indicator)*
%%%
%%% ### Description
%%% Returns the parameters used for saturating the calcium fluorescence
%%% according to: Helmchen, F. (2012). “Calcium imaging,” in Handbook of Neural Activity
%%% Measurement, eds R. Brette and A. Destexhe (Cambridge: Cambridge University
%%% Press), Eqn. (10.25). The values returned are taken from: Chen et al. (2013).
%%% "Ultrasensitive fluorescent proteins for imaging neuronal activity." Nature, 499(7458), 295-300.
%%%
%%% ### Inputs:
%%% - *indicator*: fluorescence indicator, must be one of `"GCaMP6f"` or `"GCaMP6s"`
%%%
%%% ### Outputs:
%%% - *dF_F_max*: maximum $\Delta F/F_0$ at saturation (from Chen et al.
%%% - *K_d*: affinity constant
%%% - *n_H*: Hill's exponent
%%%

switch indicator
    case {'GCaMP6s','gcamp6s'}
        dF_F_max = 16.8;     % 160 AP
        K_d = 144e-9;
    case {'GCaMP6f','gcamp6f'}
        dF_F_max = 13.1;     % 160 AP
        K_d = 290e-9;
end
end
