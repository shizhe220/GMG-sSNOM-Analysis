GMG#3_4x1um_GMGregion_rightedge_v1(07232026).csv -- column notes (2026-07-23, internal use)

graphene_4x1_manual, lp2 (right edge / interior), wn = 860-1000 cm-1. All q_p columns
are in raw cm^-1 (not this project's internal 1e5 cm^-1 convention). Sibling of lp1's
own GMG#3_4x1um_GMGregion_leftedge_v2(07222026).csv -- same columns/convention.

hankel_qp_cm-1, sqrtx_qp_cm-1
    Single-wave real-space fit (Hankel / 1-sqrt(x)). CAUTION: 5 of 15 wn
    (930,941,950,980,991) are bound-saturation artifacts of the fit's default
    lam_bound_factor=(0.7,1.2) -- NOT genuine convergence (verified via a wide
    lam_bound_factor search: no genuine local minimum exists near the other methods'
    consensus for these 5). Flagged in each row's `note`. Two more (911,960) were
    originally saturated too but were fixed by re-tuning lam0_guess (911) and lam0_guess
    + a widened lam_bound_factor just for that wn (960) -- both now genuine.

hankel_2wave_qp_cm-1
    Two-wave Hankel fit (edge-launched "q" term + tip-launched "2q" term). Filled for
    wn=900-1000 only -- use this instead of hankel_qp_cm-1 for that range (matches lp1's
    convention exactly).

fft_amp_qp_cm-1, fft_phase_qp_cm-1, fft_complex_qp_cm-1 (+ _fwhm_cm-1 each)
    FFT peak position (Lorentzian-fit) from each of the 3 channels. All 3 always
    computed and stored regardless of which is "recommended".

fft_recommended_channel, fft_recommended_qp_cm-1, fft_recommended_fwhm_cm-1
    "amp" for all 15 wn -- unlike lp1, no per-wn amp-vs-phase channel study has been
    done for lp2 yet, so this just defaults to the project's standard channel rather
    than guessing. (lp1's convention of using phase for 991/1000 does NOT carry over
    here -- checked, and lp2's per-channel FWHM doesn't follow the same pattern; e.g.
    1000cm-1's phase-channel FWHM is actually worse than amp's for lp2.)

peak2to3_spacing_qp_cm-1
    Real-space peak2-peak3 spacing / 2, independent of any fitting model. Filled for
    all 15 wn after several rounds of manual peak-position review (see
    log_and_conventions.md for the full correction history -- some wn needed the user's
    direct visual check to catch a mis-picked valley/shoulder). 930 and 1000cm-1 are
    genuine, unresolved outliers vs. hankel_2wave/fft (~40-50% lower) even after
    correction -- the real-space signal at those two wn has multiple close sub-peaks
    rather than one clean oscillation lobe, so this isn't a picking bug, just an open
    question about which sub-feature the spacing method should be reading. Flagged in
    their `note`.

note
    Short in-CSV flags for the caveats above -- kept brief since this file may be
    shared; see this README for the detailed version, for our own records.
