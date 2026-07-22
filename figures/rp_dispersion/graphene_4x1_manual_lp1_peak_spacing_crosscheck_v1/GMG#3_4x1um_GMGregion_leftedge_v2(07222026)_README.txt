GMG#3_4x1um_GMGregion_leftedge_v2(07222026).csv -- column notes (2026-07-22, internal use)

graphene_4x1_manual, lp1 (left edge), wn = 860-1000 cm-1. All q_p columns are in raw
cm^-1 (not this project's internal 1e5 cm^-1 convention).

hankel_qp_cm-1, sqrtx_qp_cm-1
    Single-wave real-space fit (Hankel / 1-sqrt(x)). Valid across the full 860-1000
    range but becomes unreliable above ~900cm-1 -- see hankel_2wave_qp_cm-1 below.

hankel_2wave_qp_cm-1
    Two-wave Hankel fit (edge-launched "q" term + tip-launched "2q" term, sharing one
    physical q). Filled for wn=900-1000 only. Above ~900cm-1 the edge-launched wave is
    no longer negligible (visible in the FFT's own 2-peak decomposition), so the
    single-wave columns above become unreliable there -- use this column instead of
    hankel_qp_cm-1 for wn>=900 (most visible at 1000cm-1: hankel_qp_cm-1=373999 is a
    clear outlier vs. hankel_2wave_qp_cm-1=242187, which agrees much better with every
    other method).

fft_amp_qp_cm-1, fft_phase_qp_cm-1, fft_complex_qp_cm-1 (+ _fwhm_cm-1 each)
    FFT peak position (Lorentzian-fit) from each of the 3 channels (amplitude, phase,
    complex = amp*exp(i*phase)). All 3 are always computed and stored.

fft_recommended_channel, fft_recommended_qp_cm-1, fft_recommended_fwhm_cm-1
    Which single channel to plot/use by default, pre-selected per row -- read this
    instead of guessing between the 3 fft_* columns above. "amp" (the project default)
    for wn=860-980. "phase" for wn=991,1000: at these two wn the amp-channel FFT peak
    is essentially unresolved (FWHM is 85-118% of the peak position itself, i.e.
    comparable to or larger than the peak position). Phase channel is a real
    improvement at 991 (FWHM/qp drops to 21%). At 1000, phase channel's FWHM/qp is
    still 57% (still broad) and its qp value (426864) is a visible outlier vs. every
    other method for that wn (hankel_2wave=242187, peak2to3_spacing=182200) -- treat
    the 1000cm-1 point as flagged/low-confidence regardless of which channel is used.

peak2to3_spacing_qp_cm-1
    Real-space peak2-peak3 spacing / 2 (independent of both Hankel and FFT). Filled for
    all 15 wn.

note
    Short in-CSV flag for the two items above (two-wave hankel above 900cm-1, phase
    channel at 991/1000) -- kept brief on purpose since this file is meant to be shared;
    this README is the detailed version, for our own records.
