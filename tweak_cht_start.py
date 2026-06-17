import nbformat as nbf

filepath = 'fitting_pipeline.ipynb'
with open(filepath, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

code_cht = r"""# === Tunable Parameters for CHT ===
x_start_cht = 0.0   # 🌟 Change this to >0 (e.g. 0.05) if you really want to cut the left edge!
lam0_guess_um = 0.420
L_cutoff = 0.9      # Your global cutoff
k_fit_range = (20, 35) # CHT k-space fit window
# ==================================

# Apply the optional left boundary cut
mask_cht = x_f >= x_start_cht
x_f_cht = x_f[mask_cht]
sig_f_cht = sig_f[mask_cht]

q_re_guess = 2 * np.pi / lam0_guess_um
p0_guess = [np.nanmax(np.abs(sig_f_cht)), q_re_guess, 0.5, 0.0]

fit_params_cht, k_arr, T_data, T_mod, mod_sig_cplx = nanoftir.fit_cht_peaks(
    x_f_cht, sig_f_cht, L=L_cutoff, k_fit_range=k_fit_range, p0=p0_guess, k_plot_range=(0.1, 40)
)

A_fit, q_re_fit, q_im_fit, phase_fit = fit_params_cht
lam_cht = (2 * np.pi / q_re_fit) * 1000
damping_cht = q_re_fit / q_im_fit
q_p = q_re_fit + 1j * q_im_fit

print(f"CHT Result: Wavelength = {lam_cht:.1f} nm, Damping = {damping_cht:.1f}")

# Reconstruct the real-space model
mod_sig_fit = np.real(mod_sig_cplx)

# Envelope calculation
x_safe = np.maximum(x_f_cht, 1e-5)
import scipy.special as sp
envelope = A_fit * np.abs(sp.hankel1(0, 2 * q_p * x_safe))

# ==========================================
# 📊 Plotting in Paper Style + k-space
# ==========================================
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), dpi=100)

c_data = '#555555'  
c_fit = '#b2182b'   
c_env_a = '#fddbc7' 
c_env_b_blue = '#92c5de' 
c_env_b_red = '#f4a582'  

# Panel (a)
ax = axes[0]
x_nm = x_f_cht * 1000 # Use the sliced x array for plotting
ax.fill_between(x_nm, envelope, -envelope, color=c_env_a, alpha=0.6)
ax.plot(x_nm, envelope, color='#d6604d', linestyle=':', lw=1.5)
ax.plot(x_nm, -envelope, color='#d6604d', linestyle=':', lw=1.5)
ax.plot(x_f * 1000, sig_f, marker='x', markersize=7, markeredgewidth=1.5, linestyle='None', color='lightgray', label='Discarded Data')
ax.plot(x_nm, sig_f_cht, marker='x', markersize=7, markeredgewidth=1.5, linestyle='None', color=c_data, label='Fitted Data')
ax.plot(x_nm, mod_sig_fit, color=c_fit, lw=2, label='Fit')
ax.set_xlabel('Distance from edge (nm)', fontweight='bold')
ax.set_ylabel(r'Re $\xi_{\mathbf{opt}}$ (a.u.)', fontweight='bold')
ax.set_xlim(-20, L_cutoff * 1000 + 20)
ax.set_yticks([])
txt_a = rf"$\lambda_p = {lam_cht:.1f}$ nm"
ax.text(0.95, 0.95, txt_a, transform=ax.transAxes, fontsize=12, va='top', ha='right')

# Panel (b)
ax = axes[1]
sqrt_x = np.sqrt(x_safe)
sig_fit_sqrt = sig_f_cht * sqrt_x
mod_sig_fit_sqrt = mod_sig_fit * sqrt_x
envelope_sqrt = envelope * sqrt_x

q_im_ideal = q_re_fit / 70.0
q_p_ideal = q_re_fit + 1j * q_im_ideal
envelope_ideal = A_fit * np.abs(sp.hankel1(0, 2 * q_p_ideal * x_safe))
envelope_sqrt_ideal = envelope_ideal * sqrt_x

ax.fill_between(x_nm, envelope_sqrt_ideal, -envelope_sqrt_ideal, color=c_env_b_blue, alpha=0.8)
ax.fill_between(x_nm, envelope_sqrt, -envelope_sqrt, color=c_env_b_red, alpha=0.8)
ax.plot(x_nm, envelope_sqrt_ideal, color='#053061', linestyle=':', lw=1.5)
ax.plot(x_nm, -envelope_sqrt_ideal, color='#053061', linestyle=':', lw=1.5)
ax.plot(x_nm, envelope_sqrt, color='#67001f', linestyle='--', lw=1.5)
ax.plot(x_nm, -envelope_sqrt, color='#67001f', linestyle='--', lw=1.5)
ax.plot(x_nm, sig_fit_sqrt, marker='x', markersize=7, markeredgewidth=1.5, linestyle='None', color=c_data)
ax.plot(x_nm, mod_sig_fit_sqrt, color=c_fit, lw=2)

ax.set_xlabel('Distance from edge (nm)', fontweight='bold')
ax.set_ylabel(r'Re $\xi_{\mathbf{opt}} \times \sqrt{\mathbf{x}}$ (a.u.)', fontweight='bold')
ax.set_xlim(-20, L_cutoff * 1000 + 20)
ax.set_yticks([])

max_env = np.max(envelope_sqrt)
ax.annotate(rf"$\gamma_p^{{-1}} = {damping_cht:.1f}$", 
            xy=(x_nm[-1]*0.7, envelope_sqrt[-1]*1.0),
            xytext=(x_nm[-1]*0.4, envelope_sqrt[-1]*1.0 + 0.5 * max_env),
            color='#67001f', fontsize=12, fontweight='bold',
            arrowprops=dict(arrowstyle="->", color='#67001f', lw=1.5, linestyle='--'))

ax.annotate(rf"$\gamma_p^{{-1}} = 70$", 
            xy=(x_nm[-1]*0.7, -envelope_sqrt_ideal[-1]*1.0),
            xytext=(x_nm[-1]*0.4, -envelope_sqrt_ideal[-1]*1.0 - 0.5 * max_env),
            color='#053061', fontsize=12, fontweight='bold',
            arrowprops=dict(arrowstyle="->", color='#053061', lw=1.5, linestyle=':'))

# Panel (c): k-space
ax = axes[2]
ax.plot(k_arr, np.abs(T_data), 'o', color=c_data, markersize=4, label='Data $|T(k)|$')
ax.plot(k_arr, np.abs(T_mod), '-', color=c_fit, lw=2, label='Fit $|T(k)|$')
ax.axvspan(k_fit_range[0], k_fit_range[1], color='gray', alpha=0.2, label='Fit Range')

k_peak = 2 * q_re_fit
ax.axvline(k_peak, color='blue', linestyle='--', lw=1, label=rf'Peak $\approx 2q_p$')

ax.set_xlabel(r'Momentum $k$ ($\mu$m$^{-1}$)', fontweight='bold')
ax.set_ylabel(r'$|T(k)|$ (a.u.)', fontweight='bold')
ax.set_xlim(0, 40)
ax.legend()
txt_c = rf"$q_p = {q_re_fit:.2f} + i{q_im_fit:.2f}$ $\mu$m$^{{-1}}$"
ax.text(0.95, 0.5, txt_c, transform=ax.transAxes, fontsize=12, va='center', ha='right')

plt.tight_layout()
fig
"""

count = 0
for cell in nb.cells:
    if cell.cell_type == 'code' and '=== Tunable Parameters for CHT ===' in cell.source:
        cell.source = code_cht
        count += 1

with open(filepath, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print(f"CHT cell successfully updated. Replaced {count} cells.")
