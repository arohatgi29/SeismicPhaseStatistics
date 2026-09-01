### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ b8a71a38-ee02-48d7-8ceb-1e9fa411cafe
using CSV, DataFrames, FFTW, Random, Plots

# ╔═╡ db63ff6f-86f4-441e-82df-4db7339c1cf3
using CircStats

# ╔═╡ 6b4e6ef5-776c-48f2-a805-5d0e71f7274d
using StatsBase

# ╔═╡ 954c55bc-bb1b-43e7-a364-b486ea070114
using SpecialFunctions

# ╔═╡ f5138299-f504-49ea-afa6-ec040596dbf0
using Distributions

# ╔═╡ 3b71e1f0-0c4a-11f1-0aa6-43a57fc47120
md"""
# **Phase variance as a seismic quality-control attribute**
"""

# ╔═╡ 67f35fa0-dfa6-424d-ba6b-7e9005936ff1
md"""
Akshika Rohatgi¹, Andrey Bakulin¹, and Sergey Fomel¹

¹ University of Texas at Austin, Bureau of Economic Geology
"""

# ╔═╡ f40c721d-f78c-41e8-97ed-efd0625d88c3
md"""
# **Abstract**
"""

# ╔═╡ a6d2d852-80d0-4074-83e2-53854ad1e5de
md"""
Seismic wavefields recorded on land are strongly distorted by near-surface heterogeneity, which introduces trace-specific, frequency-dependent phase perturbations that persist even after advanced time processing. Conventional processing relies primarily on surface-consistent deconvolution, which targets long- to mid-wavelength phase variability under the highly simplified assumption of surface consistency, equalizing large-scale trends and taming variability through overdetermination. However, this approximation is inherently unable to correct localized, non-surface-consistent phase distortions, and its effectiveness further degrades when such effects dominate, as is often the case for point-receiver data.

A separate and equally important limitation is that conventional workflows provide no direct, quantitative measure of phase reliability. Phase quality is therefore assessed only indirectly, typically through amplitude behavior or visual inspection, leaving residual phase disorder largely undiagnosed.

We introduce phase variance as a seismic quality-control attribute by treating seismic phases as circular random variables and analyzing local trace ensembles using circular statistics. This data-driven measure quantifies localized phase dispersion without phase unwrapping, enabling analysis of local phase trends and fluctuations without global assumptions or wavelet models. Phase variance is computed automatically and provides frequency-by-frequency classification of the data, ranging from coherent signal behavior to fully randomized, noise-dominated phase. Synthetic tests confirm that phase variance reliably captures imposed phase perturbations and their frequency dependence.

Application of phase variance analysis to field prestack land data shows that conventional processing reduces phase variability primarily in the low-to-intermediate frequency range and struggles within the noise cone, while the highest and lowest frequencies often show little improvement in phase coherence. Phase variance operates automatically over the full prestack volume, from shallow to deep, and frequency by frequency, providing a consistent, human-independent metric for defining effective bandwidth based on phase coherence and supporting phase-sensitive workflows such as AVO, migration, and full-waveform inversion.
"""

# ╔═╡ a400b618-d29d-46bd-87e0-09f4b92bc574
md"""
# **Introduction**
"""

# ╔═╡ cd96e0e5-e393-4f43-a7a2-ec0e15d95784
md"""
Seismic wavefields recorded at the surface are inevitably distorted by near-surface heterogeneity and small-scale elastic contrasts. These perturbations scatter energy and introduce complex, frequency-dependent distortions of amplitude and phase (Sato, Fehler, and Maeda, 2012). In land seismic processing, such effects have traditionally been mitigated through the extensive use of source and receiver arrays during acquisition (Meunier, 2011), followed by surface-consistent processing methods. Surface-consistent deconvolution and residual statics model the data using operators assigned separately to sources and receivers, implicitly assuming that near-surface distortions repeat in a predictable manner across traces (Taner, Koehler, and Alhilali, 1974; Taner and Koehler, 1981). Randomness is accommodated only in a limited sense, typically through trace-dependent static time shifts parameterized by one operator per source and one per receiver. These methods were deliberately designed to limit the number of degrees of freedom and overdetermine the estimation problem, enabling robust first-order corrections for statics and coupling, but without being grounded in a physical model of wave propagation through heterogeneous media.

Cary and Nagarajappa (2014a, 2014b) explicitly demonstrated that surface-consistent deconvolution can leave substantial residual phase errors and proposed additional surface-consistent phase corrections coupled with residual statics estimation. However, the underlying assumption remained that phase corrections are frequency-independent and surface-consistent. As a result, these approaches remain fundamentally unable to address the trace-to-trace, frequency-dependent phase variability that is abundant in modern single-sensor data. Several extensions have sought to relax strict surface consistency by introducing non-surface-consistent, time-varying statics that allow trace-dependent perturbations (Reilly et al., 2010). While these methods come closest to acknowledging the increased complexity of near-surface effects within conventional processing, they remain fundamentally constrained. By construction, they address only travel-time corrections and do not account for localized, frequency-dependent phase perturbations arising from small-scale scattering. Consequently, the available degrees of freedom remain insufficient to represent waveform distortions that vary continuously along the trace and differ for each source–receiver pair.

In practice, wave propagation through heterogeneous near-surface media is not surface consistent. Small-scale heterogeneity induces forward scattering and interference that generate non-surface-consistent, trace-specific phase perturbations (Bakulin et al., 2022a). These perturbations are not purely kinematic and cannot be corrected by time shifts alone. Instead, they modify the phase locally in time and frequency, even after advanced statics and deconvolution have been applied. As a result, residual phase disorder persists throughout prestack gathers, revealing a fundamental gap between classical processing assumptions and the physics of wave propagation in heterogeneous media (Stork et al., 2020; Bakulin et al., 2020).

Ensemble-based processing methods implicitly acknowledge the importance of phase coherence by exploiting redundancy in multi-channel seismic data. Historically, this redundancy was realized directly in acquisition through the use of source and receiver arrays, where local stacking in the time domain improved phase coherence but simultaneously conditioned phase and amplitude, a limitation perceived as detrimental to amplitude fidelity and high-frequency preservation (Newman and Mahoney, 1973). In seismology, phase-weighted stacks were introduced to explicitly prioritize phase coherence for signal detection, reflecting the recognized outsized importance of recovering correct phase from multi-channel recordings (Schimmel and Paulssen, 1997; Schimmel and Gallart, 2007). Phase-estimation and enhancement methods based on local attributes and blind deconvolution were later proposed to stabilize phase behavior (van der Baan, 2008; van der Baan and Fomel, 2009; Fomel and van der Baan, 2014; Holt and Lubrano, 2020). More recently, with the widespread adoption of point-receiver acquisition, ensemble-based phase conditioning has been reintroduced into reflection seismic processing through seismic time-frequency masking, where phase manipulation is decoupled from amplitude effects (Bakulin et al., 2023). While these methods effectively reduce incoherent phase locally, they treat phase variability primarily as a nuisance to be suppressed and do not explicitly analyze phase as a statistical quantity or quantify the underlying phase distributions that govern coherence loss.

Despite decades of advances in noise attenuation and imaging, and long-standing recognition of the critical role of phase coherence, seismic processing still lacks a routine, quantitative, and objective measure of phase reliability. In practice, phase quality is assessed indirectly through visual inspection, stack response, or amplitude-based attributes, relying on ad-hoc heuristics and human judgment rather than a standardized metric that isolates phase behavior or captures its frequency-dependent variability. This limitation is most evident in land environments with strong near-surface heterogeneity, where pronounced phase variability persists despite careful conventional processing (Bakulin et al., 2020; Stork, 2020). The problem becomes acute in desert environments such as the Middle East, where near-surface heterogeneity is particularly strong (Bakulin et al., 2022a; 2024), and reaches its extreme in highly scattering media such as basalts and crystalline rocks, where seismic imaging is often restricted to very low frequencies (Ziolkowski et al., 2003).

Recent work has reframed near-surface scattering as a multiplicative process that randomizes seismic phase, giving rise to speckle-like behavior in seismic wavefields (Bakulin et al., 2022a; 2023; Rohatgi et al., 2025). Analogous phenomena are well documented in optics and acoustics, where volumetric scatterers produce frequency-dependent phase fluctuations and coherence loss that are best described statistically (Goodman, 2007; Abbott and Thurstone, 1979). In seismic data, forward scattering on meter-scale heterogeneity generates multiple near-ballistic arrivals whose interference leads to rapid, trace-to-trace phase variations. These phase perturbations often dominate waveform degradation and increase with frequency, driving the wavefield toward a speckle regime in which individual traces appear chaotic, yet ensemble statistics remain stable and physically meaningful.

Modern land seismic acquisition is characterized by dense spatial sampling with point receivers rather than source-receiver arrays. While this increases spatial resolution and preserves high-frequency information, it also makes small-scale near-surface heterogeneity and trace-to-trace phase variability far more apparent within localized, densely sampled ensembles. In such settings, propagation through heterogeneous near-surface media renders the observed seismic phase effectively random from trace to trace.

This study adopts an explicitly stochastic viewpoint. Although the seismic source is deterministic, the resulting wavefield observed in prestack data reflects the cumulative effects of near-surface scattering and interference. Importantly, seismic phase is inherently circular rather than linear and must be treated accordingly. We introduce an ensemble-based framework that models seismic phase as a circular random variable and applies circular statistics to quantify phase variability within prestack gathers. We define a phase-variance attribute that summarizes frequency-dependent phase disorder. By shifting from single-trace phase estimates to ensemble phase distributions, the proposed framework provides a physically grounded measure of scattering-induced phase noise and establishes a transparent link among near-surface heterogeneity, phase disorder, and seismic coherence.
"""

# ╔═╡ a911602a-79e2-4b77-9be7-cdf2462dda18
md"""
# **PHASE BEHAVIOR AND STATISTICS IN SEISMIC ANALYSIS**
"""

# ╔═╡ dc709536-5ebb-4d6e-b196-1040262c5f7a
md"""
To motivate a statistical description of seismic phase, we start with controlled examples that show how phase coherence changes within local trace ensembles. Rather than interpreting phase trace by trace, we treat the phases at each frequency as a sample drawn from an underlying distribution and track how that distribution evolves as noise increases. We therefore analyze a sequence of cases spanning the full range of signal-to-noise ratios (SNR), anchored by two end members, a pure signal reference and a pure noise reference (Figure 1). These idealized limits provide intuition for the behavior encountered in real prestack gathers and establish a baseline for interpreting phase variability in more complex settings. Figure 1 organizes the sequence from fully coherent to fully randomized phase behavior. In the noise-free limit (SNR → ∞), phases are aligned across traces, and the phase distribution is maximally concentrated. As noise is introduced, phase alignment degrades and the distribution broadens. In the opposite limit (SNR → −∞), the phases are uniformly distributed on [−π, π], indicating complete loss of coherence. The following subsections describe these three regimes in detail and use them to motivate circular statistical measures of phase direction and dispersion.
"""

# ╔═╡ 4f9d39a7-a42f-4025-b4a4-23861b7474f8
T = 0.002

# ╔═╡ 3fd44e63-acf6-4be1-af3e-c80c5675dee8
fs = 1/T

# ╔═╡ e914ad8a-969a-4078-8786-7d0c7c686c23
num_wiggles=10000

# ╔═╡ b4515d5f-f5e6-4f2d-af72-b501b2c9902b
begin
	using HTTP
	file_id = "1suyqkO5jKfG3As5iub7iuWMDFXibwLH4"
	url = "https://drive.google.com/uc?export=download&id=$file_id"
	response = HTTP.get(url)
	klauder = DataFrame(CSV.File(IOBuffer(response.body))).signal[500:672]
	clean_traces = repeat(klauder, 1, num_wiggles)
end

# ╔═╡ 666ed7c7-b3bf-48f9-b7ec-13b8b22e6853
t = range(start=0, stop=length(klauder)*T, length=length(klauder))

# ╔═╡ 541d34f6-8957-42ba-8935-76bf644be6ed
function plot_wiggles(traces, t; step=200, scale=200, lw=1, color=:black, title="")
    p = plot(size=(800, 600), title=title, ylabel="Time (s)", xlabel="Trace",
             yflip=true, grid=false, legend=false)
    for i in 1:step:size(traces, 2)
        trc = traces[:, i]
        plot!(p, trc .* scale .+ i, t, linecolor=color, lw=lw)
    end
    return p
end

# ╔═╡ 22cd54ed-1fde-4fd3-a86f-23d5314792f6
function trace_spectrum(data, dt; extras=false)
	nt, ntr = size(data)
	freq = fftfreq(nt, 1/dt)
	pos_inds = findall(x -> x > 0, freq) # excludes DC
	nfreq = length(pos_inds)
	phase_vals = zeros(Float32, nfreq, ntr)
	amplitude_vals = zeros(Float32, nfreq, ntr)
	f_zero = zeros(ComplexF64, 1, ntr)
	for k in 1:ntr
	fft_trace = fft(data[:, k])
	f_zero[1, k] = fft_trace[1]
	phase_vals[:, k] = angle.(fft_trace[pos_inds])
	amplitude_vals[:, k] = abs.(fft_trace[pos_inds])
end
	return extras ? (phase_vals, amplitude_vals, freq[pos_inds], f_zero) :
	(phase_vals, amplitude_vals)
end

# ╔═╡ 55c749dd-532e-47ec-902c-55e050daa829
phase_clean, amplitude_clean, positive_freq, f_zero =
trace_spectrum(clean_traces, T, extras=true)

# ╔═╡ a5e0c66e-aa8a-4c0c-988b-0e27f2ef2cba
function circle(R; n=200)
    θ = range(0, 2π, length=n)
    return cos.(θ) .* R, sin.(θ) .* R
end

# ╔═╡ 081e2585-fae8-47ad-b055-6e58f9d4dc34
function sector(θ1, θ2, r; n=50)
    θ = range(θ1, θ2, length=n)
    x = vcat(0.0, r .* cos.(θ), 0.0)
    y = vcat(0.0, r .* sin.(θ), 0.0)
    return Shape(x, y)
end

# ╔═╡ a2205f1e-e936-4124-beb6-786dbc522fb2
function circ_hist_with_vm(phases, binwidth, μ=nothing, κ=nothing;
    title_str="",
    fillcol="#6fa76f",
    linecol=:red,
    lw=3.5,
    show_vm=true
)
    n = max(1, round(Int, 2π/binwidth))
    binwidth = 2π/n
    bins = range(0, stop=2π, length=n+1)
    data = mod.(phases, 2π)
    data .+= 10*eps(float(eltype(phases)))
    hist = fit(Histogram, data, bins, closed=:left)
    edges = hist.edges[1]
    R = maximum(hist.weights)

    p = plot(aspect_ratio=:equal, xlim=(-R, R), ylim=(-R, R), legend=false, grid=false, title=title_str)
    plot!(p, circle(R); lc=:black, fill=nothing)
    plot!(p, sector.(edges[1:end-1], edges[2:end], hist.weights), fill=fillcol, lw=1, showaxis=false)

    # only draw von Mises curve if requested and parameters are provided
    if show_vm && !isnothing(μ) && !isnothing(κ)
        θ = range(0, 2π, length=2000)
        r = @. 1 / (2π * besseli(0, κ)) * exp(κ * cos(θ - μ))
        rmax = maximum(r)
        r = rmax == 0 ? r : (R * r / rmax)
        x = r .* cos.(θ)
        y = r .* sin.(θ)
        plot!(p, x, y, color=linecol, lw=lw)
    end

    return p
end

# ╔═╡ 8cccd756-0e8e-46f8-bdb9-698b9e90b62c
function add_noise(signal, num_wiggles, desired_SNR_dB)
	noise_power = var(signal) / 10^(desired_SNR_dB / 10)
	noisy_traces = zeros(Float64, length(signal), num_wiggles)
		for i in 1:num_wiggles
			noisy_traces[:,i] = signal .+ sqrt(noise_power) .* randn(size(signal))
		end
		return noisy_traces
end

# ╔═╡ f9bf018d-b213-4673-ad05-10c1384b0f45
noise_only = add_noise(klauder, num_wiggles, -100);

# ╔═╡ 2715160b-7b81-4ebc-9d36-853be73bf925
phase_noise_only, amplitude_noise_only, positive_freq_noise_only, f_zero_noise_only =
trace_spectrum(noise_only, T, extras=true)

# ╔═╡ 38b2cb6b-a303-49dc-b60c-7085f48b9648
noise_only_plot = noise_only*0.00002;

# ╔═╡ 02d68919-6983-47f0-bb86-2707c2485de5
add_noise_traces5dB = add_noise(klauder, num_wiggles, -5);

# ╔═╡ e18a3bea-f46e-49e5-bac7-5da9a532ffee
phase_add_noise, amplitude_add_noise, positive_freq_add_noise, f_zero_add_noise =
trace_spectrum(add_noise_traces5dB, T, extras=true)

# ╔═╡ fc174366-77a5-4e98-855e-e38439a8507b
function kappa_from_R(R)
    if R < 0.53
        return 2R + R^3 + (5R^5)/6
    elseif R < 0.85
        return -0.4 + 1.39R + 0.43/(1 - R)
    else
        return 1/(3R - 4R^2 + R^3)
    end
end

# ╔═╡ 216596c5-bfa8-4abf-aebb-cc7ce3ad0d4c
function multiplicative_noise(clean_traces, phase_clean, amp_clean, f_zero, positive_freq, fs; tau=0.004)
    # ── compute frequency-dependent kappa from von Mises model 
    sig = 2π .* positive_freq .* tau
    num_samples = length(sig)
    kappa_in = zeros(length(sig))
    for (i, M) in enumerate(sig)
        phi = randn(num_samples) .* M
        z = exp.(1.0im .* phi)
        phase = angle.(z)
        r = circ_r(phase)
        kappa_in[i] = circ_kappa(r)
    end
    # trim to match phase_clean rows
    nf = size(phase_clean, 1)
    kappa_in = kappa_in[1:nf]
    # ── apply multiplicative phase perturbations 
    perturbed_phases = similar(phase_clean)
    for i in 1:size(phase_clean, 1)
        pp = phase_clean[i, :]
        kappat = kappa_in[i]
        von_mises_dist = VonMises(0.0, kappat)
        random = rand(von_mises_dist, size(phase_clean, 2))
        perturbed_phases[i, :] .= pp .+ random
    end
    # ── reconstruct perturbed signal via inverse FFT 
    new_f = amp_clean .* exp.(1.0im .* perturbed_phases)
    perturbed_signal_half = zeros(ComplexF64, size(new_f, 1), size(new_f, 2))
    for i in 1:size(new_f, 2)
        perturbed_signal_half[:, i] = reverse(conj(new_f[:, i]))
    end
    perturbed_signal = vcat(f_zero, new_f, perturbed_signal_half)
    ifft_perturbed_signal = zeros(Float64, size(perturbed_signal, 1), size(perturbed_signal, 2))
    for i in 1:size(perturbed_signal, 2)
        ifft_perturbed_signal[:, i] = real(ifft(perturbed_signal[:, i]))
    end
    return ifft_perturbed_signal
end

# ╔═╡ 54c4cabf-b788-4512-94ec-8d32c4c9a5b7
mul_noise_traces = multiplicative_noise(
    clean_traces, phase_clean, amplitude_clean, f_zero, positive_freq, fs; tau=0.004)

# ╔═╡ 68879a15-a986-4985-8c05-73ccb092c990
phase_mul_noise, amplitude_mul_noise, positive_freq_mul_noise, f_zero_mul_noise =
    trace_spectrum(mul_noise_traces, T, extras=true)

# ╔═╡ 5b0b748d-126e-4615-b909-40adbff7f685
md"""
## **Signal only (SNR = ∞)**
"""

# ╔═╡ aa88a1e0-4586-443e-9d0f-1c1a29eb330b
md"""
We begin with the noise-free reference case. Here the time-domain traces are perfectly repeatable within the local ensemble (Figure 1a), and the frequency-domain phases are coherently aligned across traces at all frequencies (Figure 1b). To describe this behavior statistically, we examine the distribution of phase values across the ensemble at each frequency, rather than individual phase spectra. For a multichannel dataset with ``N`` traces, define the set of phases at frequency ``\omega`` as
"""

# ╔═╡ 95c4fc77-9389-4576-a22e-ff2ada4c0ca6
# Cell 2 - the equation alone
md"""
```math
\Theta(\omega) = \{\theta_1(\omega),\, \theta_2(\omega),\, \ldots,\, \theta_N(\omega)\} \qquad (1)
```
"""

# ╔═╡ cc90f86b-3d00-4ce2-abf6-4eb5f7490c9c
md"""
In the absence of noise, all phases are identical, ``\theta_1(\omega) = \theta_2(\omega) = \theta_N(\omega) = \theta_0(\omega)``, so the phase distribution collapses to a single value, which defines the lower bound of phase variability. This collapse is illustrated at 17 Hz and 52 Hz in Figures 1d and 1f, respectively, where phases are tightly concentrated, and the circular variance is zero (``V = 0.00``). This limiting case provides a baseline against which phase perturbations from additive noise and scattering-induced multiplicative distortions can be quantified.
"""

# ╔═╡ 08a4cdc5-ff4a-44c8-96d1-4f3701e32314
begin
    V1d = round(max(0.0, 1 - circ_r(phase_clean[6, :])),  digits=2)
    V1f = round(max(0.0, 1 - circ_r(phase_clean[18, :])), digits=2)

    p1a = plot_wiggles(clean_traces, t; title="Clean traces (SNR = ∞)")

    p1b = heatmap(1:200:10000, positive_freq[1:25], phase_clean[1:25, 1:200:10000],
                  yflip=true, c=:coolwarm, colorbar=false, clims=(-π, π),
                  ylabel="Frequency (Hz)", title="Clean phases")

    p1c = circ_hist_with_vm(phase_clean[6, :],  2π/1000; title_str="17 Hz", show_vm=false)
    p1e = circ_hist_with_vm(phase_clean[18, :], 2π/1000; title_str="52 Hz", show_vm=false)

    p1d = plot(xlims=(-π, π), ylims=(0, 1.1), title="V = $V1d",
               legend=false, grid=false, ylabel="Normalized Density")
    vline!(p1d, [mean(phase_clean[6, :])],  color=:black, lw=3)

    p1f = plot(xlims=(-π, π), ylims=(0, 1.1), title="V = $V1f",
               legend=false, grid=false)
    vline!(p1f, [mean(phase_clean[18, :])], color=:black, lw=3)

    l = @layout [
        [a{0.5h}; b{0.5h}]  [c e; d f]
    ]

    plot(p1a, p1b, p1c, p1e, p1d, p1f,
         layout=l,
         size=(1400, 700),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         top_margin=4Plots.mm)
end

# ╔═╡ 6a3bf061-1449-4d9f-b258-0f7318c3b837
md"""
## **Noise only (SNR = -∞)**
"""

# ╔═╡ 79b745b1-8ed8-4ccc-850e-658fd52145e3
md"""
At the opposite extreme, we consider the noise-only limit, which defines the upper bound of phase variability (Figure 1g). With no coherent signal present, the frequency-domain phases are random and show no alignment across traces (Figure 1h). For a local ensemble of ``N`` traces, the phase set at frequency ``\omega`` is again given by equation (1), but the phases are now mutually uncorrelated across channels. As illustrated in Figures 1i and 1k, the phases form a uniform circular distribution over ``[-\pi, \pi]`` at all frequencies,
"""

# ╔═╡ d66e5b68-0fd1-4e96-855b-8dd2b42977fd
md"""
```math
\theta_k \sim U(-\pi, \pi), \qquad k = 1, 2, \ldots, N. \qquad (2)
```
"""

# ╔═╡ 21990862-4208-46a0-84cb-f2478cbeccdf
md"""
This case corresponds to a complete loss of phase coherence, with circular variance equal to one. Goodman (2007) showed that when a wavefield is modeled as a superposition of many independent scatterers with random phases, the resulting phase of the complex field is uniformly distributed on the circle. In the seismic context, Bakulin et al. (2022a) observed the same limiting behavior within a random multiplicative (speckle) framework, confirming that the noise-only end member exhibits uniform phase statistics.
"""

# ╔═╡ 60a4d061-ca42-422b-b3f8-952b4c89453a
begin
    # circular variance values
    V1j = round(max(0.0, 1 - circ_r(phase_noise_only[6, :])),  digits=2)
    V1l = round(max(0.0, 1 - circ_r(phase_noise_only[18, :])), digits=2)

    # circular mean and kappa for noise only
    mean_noise_only_17  = circ_mean(phase_noise_only[6, :]).μ
    kappa_noise_only_17 = kappa_from_R(circ_r(phase_noise_only[6, :]))
    mean_noise_only_52  = circ_mean(phase_noise_only[18, :]).μ
    kappa_noise_only_52 = kappa_from_R(circ_r(phase_noise_only[18, :]))

    θ_range_n = range(-π, π, length=500)

    # 1g — wiggles
    p1g = plot_wiggles(noise_only_plot, t; title="Noise only traces (SNR = -∞)")

    # 1h — phase heatmap
    p1h = heatmap(1:200:10000, positive_freq_noise_only[1:25], phase_noise_only[1:25, 1:200:10000],
                  yflip=true, c=:coolwarm, colorbar=false, clims=(-π, π),
                  ylabel="Frequency (Hz)", title="Noise-only phases",
                  xlabel="Trace Index")

    # 1i — circular histogram 17 Hz
    p1i = circ_hist_with_vm(phase_noise_only[6, :], 2π/50,
                             mean_noise_only_17, kappa_noise_only_17;
                             title_str="17 Hz", show_vm=true)

    # 1k — circular histogram 52 Hz
    p1k = circ_hist_with_vm(phase_noise_only[18, :], 2π/50,
                             mean_noise_only_52, kappa_noise_only_52;
                             title_str="52 Hz", show_vm=true)

    # 1j — linear phase distribution 17 Hz
    p1j = plot(xlims=(-π, π), ylims=(0, 1), title="V = $V1j",
               legend=false, grid=false, ylabel="Normalized Density")
    histogram!(p1j, phase_noise_only[6, :], normalized=:pdf, bins=50,
               fillcolor="#6fa76f", linecolor=:black, lw=0.5)
    pdf_vm_n1 = @. exp(kappa_noise_only_17 * cos(θ_range_n - mean_noise_only_17)) /
                   (2π * besseli(0, kappa_noise_only_17))
    plot!(p1j, θ_range_n, pdf_vm_n1, color=:red, lw=2)

    # 1l — linear phase distribution 52 Hz
    p1l = plot(xlims=(-π, π), ylims=(0, 1), title="V = $V1l",
               legend=false, grid=false)
    histogram!(p1l, phase_noise_only[18, :], normalized=:pdf, bins=50,
               fillcolor="#6fa76f", linecolor=:black, lw=0.5)
    pdf_vm_n2 = @. exp(kappa_noise_only_52 * cos(θ_range_n - mean_noise_only_52)) /
                   (2π * besseli(0, kappa_noise_only_52))
    plot!(p1l, θ_range_n, pdf_vm_n2, color=:red, lw=2)

    # ── layout 
    l_noise = @layout [
        [a{0.5h}; b{0.5h}]  [c e; d f]
    ]

    plot(p1g, p1h, p1i, p1k, p1j, p1l,
         layout=l_noise,
         size=(1400, 700),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         top_margin=4Plots.mm)
end

# ╔═╡ 62f2f1f8-284a-471c-b348-e3a10242476a
md"""
## **Noisy signal (-∞ < SNR < ∞)**
"""

# ╔═╡ 2ac31116-3f26-428b-87d3-80c566af2c19
md"""
Most seismic data lie between the two end members and therefore contain both coherent reflections and noise. In this intermediate regime, reflection events are contaminated by either additive noise (Figure 1m) or multiplicative distortions associated with scattering (Figure 1s), so phase coherence is partially preserved but systematically perturbed, leading to broadened phase distributions within local ensembles.

As a conceptual model, we write the phase of the ``k``-th trace at frequency ``\omega`` as
"""

# ╔═╡ c883fcfb-1b21-4836-8ef4-2ce3321d5369
md"""
```math
\theta_k(\omega) = \theta_0(\omega) + \delta_k(\omega), \qquad k = 1, 2, \ldots, N, \qquad (3)
```
"""

# ╔═╡ f6e2cacb-c0ca-44ed-827c-4b7fcaa2b0b0
md"""
where ``\theta_0(\omega)`` denotes the unperturbed signal phase and ``\delta_k(\omega)`` is a symmetric random perturbation. As noise increases, the dispersion of ``\{\theta_k(\omega)\}`` around ``\theta_0(\omega)`` grows and the distribution broadens. To reproduce phase behavior observed in challenging land data, we adopt the multiplicative noise model of Bakulin et al. (2022a, 2023), which captures the key feature that scattering-induced phase perturbations increase with frequency; accordingly, phase distributions remain compact at low frequencies but broaden substantially at higher frequencies (Figures 1v and 1x). A practical complication is that seismic phase is wrapped to ``[-\pi, \pi]``, which makes linear summaries and linear histograms misleading, for example, through apparent bimodality caused purely by wrapping (Figures 1r and 1x). Although phase unwrapping is common in radar, InSAR, and optical interferometry (Zebker, 1998), unwrapping becomes unstable in noisy or incoherent regions (Feigl et al., 2007) and is problematic for seismic traces (Shatilo, 1992). These issues motivate treating seismic phase as circular data and analyzing it with circular statistics (Mardia and Jupp, 2000; Mohammad et al., 2021; Rohatgi et al., 2025). Circular statistics operate directly on wrapped phases, avoid boundary artifacts at ``-\pi`` and ``\pi``, and provide robust measures of mean phase direction and dispersion, forming the basis for the phase-variance attributes introduced next (Bakulin et al., 2025; Rohatgi et al., 2025).
"""

# ╔═╡ 1cbc72ad-b5e2-4709-8006-cfd5cb8eacbf
begin
    # circular variance values
    V1p = round(max(0.0, 1 - circ_r(phase_add_noise[6, :])),  digits=2)
    V1r = round(max(0.0, 1 - circ_r(phase_add_noise[18, :])), digits=2)

    # circular mean and kappa for additive noise
    mean_add_noise_17  = circ_mean(phase_add_noise[6, :]).μ
    kappa_add_noise_17 = kappa_from_R(circ_r(phase_add_noise[6, :]))
    mean_add_noise_52  = circ_mean(phase_add_noise[18, :]).μ
    kappa_add_noise_52 = kappa_from_R(circ_r(phase_add_noise[18, :]))

    θ_range = range(-π, π, length=500)

    # 1m — wiggles
    p1m = plot_wiggles(add_noise_traces5dB, t; title="Additive noise (SNR = -5 dB)")

    # 1n — phase heatmap
    p1n = heatmap(1:200:10000, positive_freq_add_noise[1:25], phase_add_noise[1:25, 1:200:10000],
                  yflip=true, c=:coolwarm, colorbar=false, clims=(-π, π),
                  ylabel="Frequency (Hz)", title="Noisy phases with additive noise",
                  xlabel="Trace Index")

    # 1o — circular histogram 17 Hz
    p1o = circ_hist_with_vm(phase_add_noise[6, :], 2π/50,
                             mean_add_noise_17, kappa_add_noise_17;
                             title_str="17 Hz", show_vm=true)

    # 1q — circular histogram 52 Hz
    p1q = circ_hist_with_vm(phase_add_noise[18, :], 2π/50,
                             mean_add_noise_52, kappa_add_noise_52;
                             title_str="52 Hz", show_vm=true)

    # 1p — linear phase distribution 17 Hz
    p1p = plot(xlims=(-π, π), ylims=(0, 1), title="V = $V1p",
               legend=false, grid=false, ylabel="Normalized Density")
    histogram!(p1p, phase_add_noise[6, :], normalized=:pdf, bins=50,
               fillcolor="#6fa76f", linecolor=:black, lw=0.5)
    pdf_vm1 = @. exp(kappa_add_noise_17 * cos(θ_range - mean_add_noise_17)) /
                 (2π * besseli(0, kappa_add_noise_17))
    plot!(p1p, θ_range, pdf_vm1, color=:red, lw=2)

    # 1r — linear phase distribution 52 Hz
    p1r = plot(xlims=(-π, π), ylims=(0, 1), title="V = $V1r",
               legend=false, grid=false)
    histogram!(p1r, phase_add_noise[18, :], normalized=:pdf, bins=50,
               fillcolor="#6fa76f", linecolor=:black, lw=0.5)
    pdf_vm2 = @. exp(kappa_add_noise_52 * cos(θ_range - mean_add_noise_52)) /
                 (2π * besseli(0, kappa_add_noise_52))
    plot!(p1r, θ_range, pdf_vm2, color=:red, lw=2)

    # ── layout 
    l_add = @layout [
        [a{0.5h}; b{0.5h}]  [c e; d f]
    ]

    plot(p1m, p1n, p1o, p1q, p1p, p1r,
         layout=l_add,
         size=(1400, 700),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         top_margin=4Plots.mm)
end

# ╔═╡ 8315a8e9-bbe4-480e-9b76-16d9de08401f
begin
    # circular variance values
    V1v = round(max(0.0, 1 - circ_r(phase_mul_noise[6, :])),  digits=2)
    V1x = round(max(0.0, 1 - circ_r(phase_mul_noise[18, :])), digits=2)

    # circular mean and kappa for multiplicative noise
    mean_mul_noise_17  = circ_mean(phase_mul_noise[6, :]).μ
    kappa_mul_noise_17 = kappa_from_R(circ_r(phase_mul_noise[6, :]))
    mean_mul_noise_52  = circ_mean(phase_mul_noise[18, :]).μ
    kappa_mul_noise_52 = kappa_from_R(circ_r(phase_mul_noise[18, :]))

    θ_range_m = range(-π, π, length=500)

    # 1s — wiggles
    p1s = plot_wiggles(mul_noise_traces, t; title="Multiplicative noise")

    # 1t — phase heatmap
    p1t = heatmap(1:200:10000, positive_freq_mul_noise[1:25], phase_mul_noise[1:25, 1:200:10000],
                  yflip=true, c=:coolwarm, colorbar=false, clims=(-π, π),
                  ylabel="Frequency (Hz)", title="Noisy phases with multiplicative noise",
                  xlabel="Trace Index")

    # 1u — circular histogram 17 Hz
    p1u = circ_hist_with_vm(phase_mul_noise[6, :], 2π/50,
                             mean_mul_noise_17, kappa_mul_noise_17;
                             title_str="17 Hz", show_vm=true)

    # 1w — circular histogram 52 Hz
    p1w = circ_hist_with_vm(phase_mul_noise[18, :], 2π/50,
                             mean_mul_noise_52, kappa_mul_noise_52;
                             title_str="52 Hz", show_vm=true)

    # 1v — linear phase distribution 17 Hz
    p1v = plot(xlims=(-π, π), ylims=(0, 1), title="V = $V1v",
               legend=false, grid=false, ylabel="Normalized Density")
    histogram!(p1v, phase_mul_noise[6, :], normalized=:pdf, bins=50,
               fillcolor="#6fa76f", linecolor=:black, lw=0.5)
    pdf_vm_m1 = @. exp(kappa_mul_noise_17 * cos(θ_range_m - mean_mul_noise_17)) /
                   (2π * besseli(0, kappa_mul_noise_17))
    plot!(p1v, θ_range_m, pdf_vm_m1, color=:red, lw=2)

    # 1x — linear phase distribution 52 Hz
    p1x = plot(xlims=(-π, π), ylims=(0, 1), title="V = $V1x",
               legend=false, grid=false)
    histogram!(p1x, phase_mul_noise[18, :], normalized=:pdf, bins=50,
               fillcolor="#6fa76f", linecolor=:black, lw=0.5)
    pdf_vm_m2 = @. exp(kappa_mul_noise_52 * cos(θ_range_m - mean_mul_noise_52)) /
                   (2π * besseli(0, kappa_mul_noise_52))
    plot!(p1x, θ_range_m, pdf_vm_m2, color=:red, lw=2)

    # ── layout 
    l_mul = @layout [
        [a{0.5h}; b{0.5h}]  [c e; d f]
    ]

    plot(p1s, p1t, p1u, p1w, p1v, p1x,
         layout=l_mul,
         size=(1400, 700),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         top_margin=4Plots.mm)
end

# ╔═╡ 1d99aa32-f8b7-481d-885c-3b46d1bf48fd
md"""
# **CIRCULAR STATISTICS FOR SEISMIC PHASE ANALYSIS**
"""

# ╔═╡ c28e595e-3857-46c6-929b-44a3ba6806e9
md"""
## **Circular Mean : Measure of mean phase direction**
"""

# ╔═╡ 5468f447-9a9a-4fdc-a6a4-7d080434a841
md"""
For an ensemble of seismic traces, the complex spectra can be normalized to unit magnitude so that only phase information is retained, independent of amplitude. This normalization ensures that phase statistics are not biased by variations in signal strength across traces. Figure 2 illustrates this procedure: noisy complex spectra (Figure 2a) are represented in the complex plane and then projected onto the unit circle (Figure 2b), preserving phase while removing amplitude effects.
"""

# ╔═╡ b00bf530-e245-45cb-8f87-1c33c394ef7f
md"""
At a given frequency ``\omega``, the phase of each trace is represented as a unit phasor
"""

# ╔═╡ 0f37ae56-928d-4a73-bf0a-81262d1ea982
md"""
```math
e^{i\theta_1(\omega)},\, e^{i\theta_2(\omega)},\, \ldots,\, e^{i\theta_N(\omega)}, \qquad (4)
```
"""

# ╔═╡ 7cb4d395-b266-46c4-8e3f-787c18b075d5
md"""
where for a complex signal ``X_k(\omega)``, ``\theta_k(\omega) = \arg[X_k(\omega)]`` is the phase of the ``k``-th trace and ``N`` is the number of traces in the local ensemble.

The preferred phase direction is characterized by the circular mean ``\bar{\theta}(\omega)``, defined as the direction of the resultant vector formed by summing the unit phasors. At a fixed frequency, the cosine and sine components are averaged as:
"""

# ╔═╡ 4022a458-2e6b-4b79-92cf-575ac2fd63e9
md"""
```math
\bar{C}(\omega) = \frac{1}{N}\sum_{k=1}^{N}\cos\theta_k(\omega), \qquad
\bar{S}(\omega) = \frac{1}{N}\sum_{k=1}^{N}\sin\theta_k(\omega). \qquad (5)
```

The length of the mean resultant vector is
```math
\bar{R}(\omega) = \sqrt{\bar{C}^2(\omega) + \bar{S}^2(\omega)}, \qquad (6)
```

and the circular mean direction is given by
```math
\bar{\theta}(\omega) = \operatorname{atan2}\!\left(\bar{S}(\omega),\, \bar{C}(\omega)\right). \qquad (7)
```

The mean direction is undefined when ``\bar{R}(\omega) = 0``, corresponding to a uniform phase distribution with no preferred direction.

It is important to emphasize that the circular mean ``\bar{\theta}`` is not equivalent to the arithmetic average ``(\theta_1 + \theta_2 + \cdots + \theta_N)/N``. Instead, it represents the direction of the geometric vector sum of the unit phasors and correctly accounts for the periodic nature of phase. This distinction is critical when analyzing wrapped seismic phases, particularly in the presence of noise.
"""

# ╔═╡ 87a0577d-0a4d-4e6e-b5b7-ecbee0f7cf51
begin
    # generate additive noise at -10 dB
    add_noise_10dB = add_noise(klauder, num_wiggles, -10)
    phase_add_noise_10dB, amplitude_add_noise_10dB, positive_freq_add_noise_10dB, f_zero_add_noise_10dB =
        trace_spectrum(add_noise_10dB, T, extras=true)

    freq_idx = 18  # 52 Hz

    signal_vector_10 = amplitude_add_noise_10dB[freq_idx, :] .* exp.(1.0im .* phase_add_noise_10dB[freq_idx, :])
    phase_vector_10  = exp.(1.0im .* phase_add_noise_10dB[freq_idx, :])

    xs_phase_10  = real.(phase_vector_10)
    ys_phase_10  = imag.(phase_vector_10)
    xs_signal_10 = real.(signal_vector_10)
    ys_signal_10 = imag.(signal_vector_10)

    step = 10
    idx_10 = 1:step:length(xs_phase_10)

    xs_phase_sub_10  = xs_phase_10[idx_10]
    ys_phase_sub_10  = ys_phase_10[idx_10]
    xs_signal_sub_10 = xs_signal_10[idx_10]
    ys_signal_sub_10 = ys_signal_10[idx_10]

    mean_phase_10 = circ_mean(phase_add_noise_10dB[freq_idx, :]).μ
    mean_r_10     = circ_r(phase_add_noise_10dB[freq_idx, :])

    # 2a — noisy complex spectra
    p2a = quiver(zeros(length(xs_signal_sub_10)), zeros(length(ys_signal_sub_10));
                 quiver=(xs_signal_sub_10, ys_signal_sub_10),
                 aspect_ratio=1,
                 xlabel="Re", ylabel="Im",
                 title="Phase vectors with noisy amplitudes",
                 legend=false, linecolor="#6fa76f", arrow=true, lw=0.5, grid=false)

    # 2b — unit phasors
    p2b = quiver(zeros(length(xs_phase_sub_10)), zeros(length(ys_phase_sub_10));
                 quiver=(xs_phase_sub_10, ys_phase_sub_10),
                 aspect_ratio=1,
                 xlabel="Re", ylabel="Im",
                 title="Phase vectors with normalized amplitudes",
                 legend=false, linecolor="#6fa76f", arrow=true, lw=0.5, grid=false)
    quiver!(p2b, [0.0], [0.0];
            quiver=([mean_r_10 * cos(mean_phase_10)], [mean_r_10 * sin(mean_phase_10)]),
            linecolor=:red, arrow=true, lw=3)
    annotate!(p2b, mean_r_10 * cos(mean_phase_10) + 0.05, mean_r_10 * sin(mean_phase_10) + 0.1,
              text("Mean phase\ndirection", :red, 9, :left))

    plot(p2a, p2b,
         layout=(1, 2),
         size=(1200, 500),
         fontsize=12, labelfontsize=12,
         legendfontsize=12, tickfontsize=12,
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         dpi=300)
end

# ╔═╡ 248b59c9-e468-41c6-a6ab-ca582d94c48e
md"""
## **Circular Variance: Measure of phase variability**
"""

# ╔═╡ bf0ed95f-cc64-418d-9fde-851657638a2f
md"""
While the mean direction characterizes the central tendency of seismic phases within an ensemble, it is equally important to quantify how tightly the phases cluster around that direction. This variability is measured by the mean resultant length ``\bar{R}`` (equation 6), which takes values between 0 and 1. Values of ``\bar{R}`` close to 1 indicate strong phase coherence, whereas values near 0 indicate highly dispersed, incoherent phases,
```math
0 \le \bar{R}(\omega) \le 1. \qquad (8)
```

Although the mean resultant length is the fundamental quantity in circular statistics, it is often more intuitive, particularly by analogy with linear statistics, to express phase variability in terms of circular variance, defined as
```math
V(\omega) = 1 - \bar{R}(\omega), \qquad (9)
```

where ``0 \le V(\omega) \le 1``.

In the ideal signal-only case, all phase vectors are aligned and point in the same direction, yielding ``\bar{R}(\omega) = 1`` and ``V(\omega) = 0``, corresponding to perfect phase coherence (Figure 3a). In contrast, when the data consist solely of noise, phase vectors are uniformly distributed around the unit circle, resulting in ``\bar{R}(\omega) \approx 0`` and ``V(\omega) \approx 1``, indicating maximal phase variability (Figure 3b). Most realistic seismic situations lie between these two extremes. For example, with moderate additive noise at ``-5`` dB, the phase vectors retain a preferred direction but exhibit noticeable dispersion, producing ``\bar{R}(\omega) \approx 0.73`` and ``V(\omega) \approx 0.27`` (Figure 3c). Under stronger or frequency-dependent phase perturbations, such as those associated with multiplicative noise, the phase vectors spread further and the circular variance increases accordingly, reaching values around ``V(\omega) \approx 0.5`` (Figure 3d). Although circular statistics provide a general framework for analyzing wrapped seismic phases, no specific parametric distribution is required to compute the mean direction or circular variance.
"""

# ╔═╡ b0b90f69-c265-4c54-bd86-5d62e97b04f7
md"""
In this study, the von Mises distribution is introduced as a convenient conceptual model, as it is the circular analogue of the normal distribution, providing an intuitive interpretation of phase clustering, dispersion, and a well-defined center of gravity when a single dominant signal event is present within the analysis window.

However, our analysis does not depend on the von Mises assumption. Instead, we primarily use the circular variance ``V(\omega)``, which is defined directly from the mean resultant length and can be estimated nonparametrically from phase ensembles. Circular variance is a robust, model-free observable that remains stable under strong noise contamination, avoids phase-unwrapping ambiguities, and can be computed consistently across frequencies, offsets, and processing stages, making it well suited as a practical QC and diagnostic metric for seismic phase behavior.
"""

# ╔═╡ f6523673-ad48-450f-8e95-4a8f4750dfa3
begin
    θcirc_fig3 = range(-π, π; length=500)
    cx_fig3, cy_fig3 = cos.(θcirc_fig3), sin.(θcirc_fig3)

    function make_phasor(phases, ttl; step=20, lw=1)
        pv = exp.(1.0im .* phases)
        x  = real.(pv)[1:step:end]
        y  = imag.(pv)[1:step:end]
        p  = quiver(zeros(length(x)), zeros(length(y));
                    quiver=(x, y),
                    aspect_ratio=1,
                    xlabel="Re", ylabel="Im",
                    legend=false, linecolor="#6fa76f",
                    arrow=true, lw=lw, grid=false,
                    title=ttl)
        plot!(p, cx_fig3, cy_fig3, color=:black, lw=2)
        # resultant vector R
        R  = circ_r(phases)
        θm = circ_mean(phases).μ
        quiver!(p, [0.0], [0.0];
                quiver=([R * cos(θm)], [R * sin(θm)]),
                linecolor=:black, arrow=true, lw=3)
        xlims!(p, -1.2, 1.2); ylims!(p, -1.2, 1.2)
        return p
    end

    p3a = make_phasor(phase_clean[18, :],      "Signal Only";          step=80, lw=4)
    p3b = make_phasor(phase_noise_only[18, :], "Noise Only";           step=20, lw=1)
    p3c = make_phasor(phase_add_noise[18, :],  "Additive Noise";       step=20, lw=1)
    p3d = make_phasor(phase_mul_noise[18, :],  "Multiplicative Noise"; step=20, lw=1)

    plot(p3a, p3b, p3c, p3d,
         layout=(2, 2),
         size=(1100, 1100),
         margin=5Plots.mm,
         fontsize=12, labelfontsize=12,
         legendfontsize=12, tickfontsize=12,
         xguidefontsize=12, yguidefontsize=12,
         dpi=300)
end

# ╔═╡ 433b4239-faab-4750-8228-fb36d7074c40
md"""
## **The von Mises Distribution**
"""

# ╔═╡ 8b382333-ff48-4174-b895-1cb8a5c28fff
md"""
The von Mises distribution, introduced by von Mises (1918) as the circular analogue of the normal distribution, has been widely used to model directional or wrapped data. It has found broad application in acoustics, optical interferometry, and InSAR, where phase is inherently periodic and must be analyzed without unwrapping (Feigl et al., 2009; Jiang et al., 2019; Takamichi et al., 2018; Nugraha et al., 2019). In these fields, von Mises models are commonly used to represent phase variability, reconstruct phase information, and analyze noisy or wrapped phase measurements in a probabilistically consistent manner. In seismology, Gosselin (2022) applied von Mises distributions to model uncertainty in circular phase spectra of surface-wave dispersion, while Rohatgi et al. (2025) used circular statistics to quantify phase coherence within local seismic ensembles.

This distribution is parameterized by the mean direction ``\bar{\theta}`` and the concentration parameter ``\kappa``. Its probability density function is given as:
```math
f(\theta;\, \bar{\theta},\, \kappa) = \frac{1}{2\pi I_0(\kappa)}\, e^{\kappa \cos(\theta - \bar{\theta})}, \qquad (10)
```

where ``I_0`` denotes the modified Bessel function of the first kind and order zero, defined as
```math
I_0(\kappa) = \frac{1}{2\pi}\int_0^{2\pi} e^{\kappa \cos\theta}\, d\theta. \qquad (11)
```

The function ``I_0`` also admits the power series expansion
```math
I_0(\kappa) = \sum_{m=0}^{\infty} \frac{1}{(m!)^2}\left(\frac{\kappa}{2}\right)^{2m}. \qquad (12)
```

The von Mises distribution is unimodal and symmetric about its mean direction. When ``\kappa = 0``, it reduces to the uniform distribution on ``[-\pi, \pi]``, corresponding to completely random phases. As ``\kappa`` increases, the distribution becomes increasingly concentrated around ``\bar{\theta}`` and approaches a normal distribution in shape, while preserving phase periodicity. Figure 4a illustrates this behavior for a fixed mean direction of 2.2 radians and increasing ``\kappa``, showing the transition from nearly uniform to strongly clustered phase distributions.

The concentration parameter ``\kappa`` can be estimated from the mean resultant length ``\bar{R}`` using standard approximations (Fisher, 1993),
```math
\kappa \approx \begin{cases} 2\bar{R} + \bar{R}^3 + \dfrac{5\bar{R}^5}{6}, & \bar{R} < 0.53 \\[8pt] -0.4 + 1.39\bar{R} + \dfrac{0.43}{1 - \bar{R}}, & 0.53 \le \bar{R} < 0.85 \\[8pt] \dfrac{1}{3\bar{R} - 4\bar{R}^2 + \bar{R}^3}, & \bar{R} \ge 0.85 \end{cases} \qquad (13)
```
"""

# ╔═╡ be3e9677-d5f8-44dd-94fd-f6de2d7a7a4a
begin
    function wrap_phase(unwrapped_phase)
        return mod.(unwrapped_phase .+ π, 2π) .- π
    end

    mean_value = 2.2

    kappa1 = 0.1
    x_vm1  = wrap_phase(rand(VonMises(mean_value, kappa1), 1000000))
    V_vm1  = round(1 - circ_r(x_vm1), digits=2)

    kappa2 = 2
    x_vm2  = wrap_phase(rand(VonMises(mean_value, kappa2), 1000000))
    V_vm2  = round(1 - circ_r(x_vm2), digits=2)

    kappa3 = 10
    x_vm3  = wrap_phase(rand(VonMises(mean_value, kappa3), 1000000))
    V_vm3  = round(1 - circ_r(x_vm3), digits=2)

    p4a = histogram(x_vm1, bins=100, normalize=true,
                    linecolor=:blue, linewidth=2, fillalpha=0,
                    label="θ=2.2, κ=$kappa1, V=$V_vm1",
                    xlabel="Phase (radians)", ylabel="Normalized Density",
                    title="von Mises Distribution",
                    grid=false, legend=:topleft,
                    fontsize=12, labelfontsize=12,
                    legendfontsize=12, tickfontsize=12,
                    xguidefontsize=12, yguidefontsize=12,
                    dpi=300)

    histogram!(p4a, x_vm2, bins=100, normalize=true,
               linecolor=:red, linewidth=2, fillalpha=0,
               label="θ=2.2, κ=$kappa2, V=$V_vm2")

    histogram!(p4a, x_vm3, bins=100, normalize=true,
               linecolor=:green, linewidth=2, fillalpha=0,
               label="θ=2.2, κ=$kappa3, V=$V_vm3")

    p4a
end

# ╔═╡ cf4f405c-f54f-4b6d-ac5a-b9f7c36c27f8
md"""
Circular variance ``V`` and the von Mises concentration parameter ``\kappa`` both describe the dispersion of phase vectors around the mean direction, but from complementary perspectives. Circular variance is defined directly from the mean resultant length ``V = 1 - \bar{R}`` and provides a bounded, geometrically intuitive measure of phase variability, where values near zero indicate strong phase coherence and values near one indicate incoherent or noise-dominated behavior. In contrast, ``\kappa`` is specific to the von Mises model and quantifies the strength of clustering in an unbounded manner, ranging from zero to infinity. As a result, ``\kappa`` can be more sensitive to noise and may fluctuate strongly between neighboring ensembles, particularly at low coherence. When the von Mises model is assumed, however, ``V`` and ``\kappa`` are uniquely related through ``\bar{R}`` and can be converted analytically or numerically, as illustrated in Figure 4b. For clarity, Table 1 summarizes the key circular statistical quantities used in this study and their physical interpretation.
"""

# ╔═╡ 48e8a718-f908-4541-af03-af0302dc83eb
md"""
| Attribute | Symbol | Range | Interpretation |
|:----------|:------:|:-----:|:---------------|
| Circular mean | ``\bar{\theta}`` | ``[-\pi, \pi]`` | Mean (central) phase direction |
| Mean resultant length | ``\bar{R}`` | ``[0, 1]`` | Degree of phase alignment (clustering strength) |
| Circular variance | ``V = 1 - \bar{R}`` | ``[0, 1]`` | Phase dispersion or incoherence |
| Concentration (von Mises) | ``\kappa`` | ``[0, \infty)`` | Phase clustering strength (inverse dispersion) |

*Table 1: Summary of key circular statistical attributes used for ensemble-based seismic phase analysis.*
"""

# ╔═╡ 2d0e8b6a-c822-4447-9a52-761e253f5a29
begin
    kappa_values = 0.001:0.1:20
    circ_var_values = map(κ -> 1.0 - circ_r(wrap_phase.(rand(VonMises(0.0, κ), 100_001))), kappa_values)

    p4b = plot(kappa_values, circ_var_values,
               xlabel="κ", ylabel="Circular Variance",
               title="Circular Variance vs κ",
               legend=false, linewidth=4, color=:black, grid=false,
               tickfontsize=12, xguidefontsize=14, yguidefontsize=14,
               titlefontsize=16, dpi=300, size=(700, 450))

    p4b
end

# ╔═╡ 470bee5b-4ea5-4a02-ba24-a508968cbf0c
md"""
## **QUANTIFYING VARIABILITY IN SEISMIC PHASE**
"""

# ╔═╡ ef5433c7-c551-452c-8d55-545629f1d26a
md"""
Although phase statistics are estimated directly from the data, without fitting any specific parametric distribution, the von Mises distribution remains an important interpretive reference. When a single coherent reflection dominates a local ensemble, the phase distribution is expected to be unimodal, with the mean direction representing the signal phase and the dispersion reflecting noise. This provides a simple and intuitive signal–plus–noise picture in circular phase space.

Real seismic data are more complex. Noise can itself be partially coherent, such as ground roll, guided waves, or other organized wave modes, and multiple signals may interfere or cross within the same analysis window. In such cases, phase distributions may broaden asymmetrically or become multimodal. This behavior does not merely indicate stronger noise, but rather the absence of a single dominant signal phase. Importantly, the apparent phase distribution may depend on the number of traces in the ensemble: small local ensembles can appear unstable or noisy, whereas larger ensembles tend to reveal more stable phase structure. Even under these conditions, the unimodal von Mises case remains a useful reference state, providing a baseline for distinguishing signal-dominated behavior from mixed or incoherent phase regimes.

Building on this conceptual framework, we demonstrate how circular variance can be used as a quantitative quality-control attribute for seismic phase. Circular variance measures the spread of phase distributions within a local ensemble. Because seismic phase variability depends on both offset and frequency, these statistics must be computed locally using sliding windows over ensembles of traces. This approach enables objective tracking of phase variability across a seismic gather and provides a means to assess the impact of processing steps.

The workflow for computing circular variance is summarized as follows. Each local time-space window is analyzed in the frequency domain using the following procedure. Each trace ``x_k(t)`` within the window is transformed to the frequency domain,
"""

# ╔═╡ d58cee5e-b7f0-4daa-b17d-3e5b498b7ef8
md"""
```math
X_k(\omega) = \mathcal{F}\{x_k(t)\} = |X_k(\omega)|\, e^{i\theta_k(\omega)},
```

where ``|X_k(\omega)|`` is the spectral amplitude and ``\theta_k(\omega)`` the spectral phase.

- To isolate phase information, each spectrum is normalized to unit magnitude,
```math
\tilde{X}_k(\omega) = \frac{X_k(\omega)}{|X_k(\omega)|} = e^{i\theta_k(\omega)}.
```

- The set ``\{\tilde{X}_k(\omega)\}_{k=1}^{N}`` thus forms a phase ensemble at frequency ``\omega``.

- Phase dispersion within the local ensemble is quantified using the circular variance,
```math
V(\omega) = 1 - \left|\frac{1}{N}\sum_{k=1}^{N} e^{i\theta_k(\omega)}\right|,
```

where ``V(\omega) \in [0, 1]``.

- The phase variance ``V(\omega)`` is assigned to the center of the analysis window and mapped across time and frequency.

The ensemble size ``N`` is chosen to ensure statistically stable and interpretable estimates, as discussed in Appendix A. To assess the reliability of this metric, we first apply it to a controlled synthetic example. This allows us to demonstrate how circular variance responds to known phase perturbations, providing a benchmark before extending the method to more complex field data.
"""

# ╔═╡ d36129ee-22ae-463f-b133-f2e596544edf
md"""
## **Synthetic example with phase perturbations**
"""

# ╔═╡ 78e5e4f0-2422-44c7-9ae3-1421aedbbbd7
md"""
The choice of ensemble size is dictated by the role of the analysis and is kept consistent between synthetic tests and field applications. In high-noise conditions, larger ensembles are required to see through noise and obtain statistically stable estimates, whereas in low-noise settings smaller ensembles may suffice. Accordingly, large ensembles of 10,000 traces analyzed using 2,000-trace sliding windows are employed to obtain statistically stable estimates of circular variance when it is used as a diagnostic and QC attribute. This ensures that the measured phase variability reflects the imposed perturbations rather than sampling noise.

To investigate the variability of seismic phase within an ensemble, we designed a controlled synthetic experiment. Phase perturbations were introduced by drawing random realizations from a von Mises distribution, consistent with the conceptual model described earlier. For each trace, a circular variance ``V`` was prescribed and systematically varied across the ensemble to impose laterally varying phase noise. Figure 5a illustrates the imposed variation in phase perturbations as a function of trace number. The variance decreases from left to right, mimicking a realistic acquisition scenario in which near-offset traces fall within a noise cone characterized by stronger phase perturbations, while far-offset traces exhibit increasingly coherent reflection behavior. The resulting synthetic seismic section is shown in Figure 5b.
"""

# ╔═╡ 79eec12e-ecaa-42ae-87fe-d80f636c0a1c
function offset_dependent_perturbation(phase_clean, amplitude_clean, f_zero, t, num_wiggles;
                                        kappa_start=0.2, kappa_stop=2.0)
    Random.seed!(2025)
    kappa_values = range(start=kappa_start, stop=kappa_stop, length=num_wiggles)
    var_values   = [1 - besseli(1, k)/besseli(0, k) for k in kappa_values]
    phase_perturbed = similar(phase_clean)
    for i in 1:size(phase_clean, 2)
        pp     = phase_clean[:, i]
        kappat = kappa_values[i]
        random = rand(VonMises(0.0, kappat), size(phase_clean, 1))
        phase_perturbed[:, i] .= pp .+ random
    end
    phase_perturbed = mod.(phase_perturbed .+ π, 2π) .- π
    new_f = amplitude_clean .* exp.(1.0im .* phase_perturbed)
    perturbed_signal_half = zeros(ComplexF64, size(new_f, 1), num_wiggles)
    for i in 1:size(new_f, 2)
        perturbed_signal_half[:, i] = reverse(conj(new_f[:, i]))
    end
    perturbed_signal = vcat(f_zero, new_f, perturbed_signal_half)
    ifft_perturbed = zeros(Float64, size(perturbed_signal, 1), num_wiggles)
    for i in 1:size(perturbed_signal, 2)
        ifft_perturbed[:, i] = real(ifft(perturbed_signal[:, i]))
    end
    return ifft_perturbed, phase_perturbed, var_values
end

# ╔═╡ 54c9b2ec-0c9d-4f23-b615-0da0d4023b9c
offset_traces, phase_offset, var_imposed = offset_dependent_perturbation(
    phase_clean, amplitude_clean, f_zero, t, num_wiggles;
    kappa_start=0.2, kappa_stop=2.0)

# ╔═╡ 988b7168-7c2a-4cae-a21c-0487b1477073
begin
    p5b = plot_wiggles(offset_traces, t; title="Noisy traces")

    p5c = heatmap(1:10000, positive_freq[1:25], phase_offset[1:25, 1:10000],
                  yflip=true, c=:coolwarm, colorbar=false, clims=(-π, π),
                  ylabel="Frequency (Hz)", title="Offset-dependent perturbed phases",
                  xlabel="Trace Index")

    l_5bc = @layout [a{0.5h}; b{0.5h}]

    plot(p5b, p5c,
         layout=l_5bc,
         size=(1000, 800),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         top_margin=4Plots.mm,
         dpi=300)
end

# ╔═╡ 2b0dfbbc-33de-4aed-bbfc-c49f66055c33
md"""

The baseline signal was generated using a Klauder wavelet, defined as the autocorrelation of a vibroseis sweep commonly employed in land seismic surveys. This wavelet served as a clean reference, from which an ensemble of 10,000 identical traces was constructed. In the frequency domain, the unperturbed signal at channel ``k`` and frequency ``\omega`` is given by
```math
X_k(\omega) = A_0(\omega)\, e^{i\theta_0(\omega)}, \qquad (14)
```

where ``A_0(\omega)`` denotes the true amplitude spectrum and ``\theta_0(\omega)`` the corresponding true phase. Phase distortions were then introduced by applying spatially dependent perturbations, such that the perturbed signal becomes
```math
X_k(\omega) = A_0(\omega)\, e^{i\{\theta_0(\omega) + \theta_k^{\delta}(\omega)\}}, \qquad (15)
```

where ``\theta_k^{\delta}(\omega)`` represents the imposed perturbation.

Within each window, amplitudes were normalized to isolate phase effects, ensuring that the computed statistics reflect phase variability alone. The resulting phase ensemble in the frequency domain is shown in Figure 5c. Sliding the window across the ensemble allows phase distributions to be evaluated as a function of offset and frequency. Representative phase distributions extracted from near-offset and far-offset windows are shown in Figures 5e and 5f, respectively. These distributions illustrate the transition from highly dispersed phases at near offsets to more concentrated phase behavior at far offsets. For example, at 17 Hz the near-offset window yields a circular variance of 0.86, indicating strong phase variability, whereas the far-offset window produces a lower value of 0.32, reflecting greater phase stability.
"""

# ╔═╡ 2be746d9-d942-4e29-933d-c28d0dec843e
begin
    # circular stats for near and far offset at 17 Hz
    mean_c_near_offset  = circ_mean(phase_offset[6, 1:2000]).μ
    kappa_c_near_offset = kappa_from_R(circ_r(phase_offset[6, 1:2000]))
    mean_c_far_offset   = circ_mean(phase_offset[6, 8000:10000]).μ
    kappa_c_far_offset  = kappa_from_R(circ_r(phase_offset[6, 8000:10000]))

    function add_circular_mean_arrow!(p, μ; scale=1.0, color=:red, lw=3)
        μw = mod(μ, 2π)
        xlims_p = Plots.xlims(p)
        R = maximum(abs.(xlims_p)) * scale
        quiver!(p, [0.0], [0.0];
                quiver=([R * cos(μw)], [R * sin(μw)]),
                color=color, lw=lw, arrow=true, label=false)
    end

    function add_circle_annotations!(p)
        xlims_p = Plots.xlims(p)
        R = maximum(abs.(xlims_p))
        # -π,π on the left
        annotate!(p, -R * 1.05, 0.0,  text("-π,\nπ",  :black, 8, :right))
        # 0 on the right
        annotate!(p,  R * 1.05, 0.0,  text("0",       :black, 8, :left))
    end

    # 5e — near offset circular histogram
    p5e = circ_hist_with_vm(phase_offset[6, 1:2000], 2π/50,
                             mean_c_near_offset, kappa_c_near_offset;
                             title_str="17 Hz (near offset)", show_vm=true)
    add_circular_mean_arrow!(p5e, mean_c_near_offset)
    add_circle_annotations!(p5e)

    # 5f — far offset circular histogram
    p5f = circ_hist_with_vm(phase_offset[6, 8000:10000], 2π/50,
                             mean_c_far_offset, kappa_c_far_offset;
                             title_str="17 Hz (far offset)", show_vm=true)
    add_circular_mean_arrow!(p5f, mean_c_far_offset)
    add_circle_annotations!(p5f)

    plot(p5e, p5f,
         layout=(1, 2),
         size=(800, 400),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         dpi=300)
end

# ╔═╡ 1b31cff4-9f96-490c-97ec-45004449ae7e
md"""
The circular variance computed within each window is assigned to the window center before advancing to the next position. Repeating this procedure produces a new seismic attribute, referred to here as phase variance. The resulting phase-variance map (Figure 5d) captures the spatial and spectral evolution of phase variability across the synthetic gather. As expected, phase variance is highest at near offsets and decreases steadily toward far offsets, consistent with the imposed perturbations. Because the perturbations in this test were identical across frequencies, the frequency dependence of the phase-variance map is uniform. The framework naturally extends to cases where phase noise varies with frequency, although the present example focuses on spatial variability for clarity.

The average phase variance across all traces (Figure 5g) shows close agreement between the imposed and recovered values, confirming that circular variance reliably quantifies phase perturbations in a controlled setting. These synthetic results demonstrate that circular variance provides a robust and interpretable measure of seismic phase variability. Having validated that ``V(\omega)`` accurately recovers known phase perturbations in synthetic data, we now apply it as a QC observable to field seismic gathers to diagnose where phase becomes unreliable.
"""

# ╔═╡ bb3bc3d2-6f4a-4021-997a-2c0613ca814d
function local_phase_variance(phases; batch_size=2000)
    n_freqs, n_traces = size(phases)
    half_batch = div(batch_size, 2)
    variance = Matrix{Float64}(undef, n_freqs, n_traces)

    for j in 1:n_traces
        start_idx = max(1, j - half_batch)
        end_idx   = min(n_traces, j + half_batch - 1)

        for i in 1:n_freqs
            window_phases = phases[i, start_idx:end_idx]
            r = circ_r(window_phases)
            variance[i, j] = 1 - r
        end
    end

    return variance
end

# ╔═╡ bbedd657-28fa-42b0-9418-296bf80feead
variance_offset = local_phase_variance(phase_offset; batch_size=2000);

# ╔═╡ 43de8089-ac37-4b06-ab88-750d152e367a
p5d = heatmap(1:10000, positive_freq[1:size(variance_offset, 1)], variance_offset,
              yflip=false, c=cgrad(:jet, rev=true), colorbar=true, clims=(0, 1),
              ylabel="Frequency (Hz)", title="Local Phase Variance",
              xlabel="Trace Index", size=(1000,400),
              left_margin=8Plots.mm, bottom_margin=8Plots.mm, dpi=300)

# ╔═╡ fd33d1bf-3ef2-474c-ae4e-df30b18e8157
begin
    p5g = plot(var_imposed, label="Variance in", lw=3, grid=false)
    plot!(p5g, variance_offset[6, :], label="Variance out", lw=3,
          linestyle=:dash, color=:black,
          ylabel="Phase Variance", xlabel="Trace Index",
          title="True vs Estimated Phase Variance",
          fontsize=12, labelfontsize=12,
          legendfontsize=12, tickfontsize=12,
          xguidefontsize=12, yguidefontsize=12,
          size=(600, 500), grid=false, dpi=300)

    p5g
end

# ╔═╡ efb01d36-8fe3-413b-9edb-d13c574acd28
md"""
## **Monitoring phase variance in field data**
"""

# ╔═╡ 9763d963-8e50-42e6-b9f5-eedfaa60571e
md"""
Monitoring phase quality is critical in seismic processing because phase integrity directly controls imaging, inversion, and attribute analysis. Despite its importance, conventional workflows provide few tools for quantifying phase variability or tracking phase noise in a systematic way. As a result, phase distortions caused by near-surface effects, scattering, or processing artifacts may persist and propagate through the workflow without being explicitly diagnosed. Circular variance provides a practical phase quality-control attribute that allows phase behavior to be evaluated objectively as a function of offset and frequency.

To illustrate this capability, we analyze two prestack seismic gathers extracted from a 3D CMP supergather, one before and one after conventional time processing. A 300-ms time window containing a deep reflection event was selected. Phase statistics were computed using a sliding-window approach with local ensembles of 2,000 traces across a total of 10,000 traces. The gathers shown in Figures 6a and 6b are displayed after normal moveout correction and statics. The conventional time-processing sequence included linear noise attenuation, refraction statics, random noise removal, two passes of surface-consistent deconvolution, post-deconvolution noise attenuation targeting linear, random, and burst noise, merge-phase and static matching, surface-consistent scaling, velocity analysis, and residual statics correction.

The effect of processing on seismic phase can be examined by inspecting phase distributions within local time–frequency windows. At 16 Hz in a mid-offset local window, the pre-processed data exhibit a broad and weakly organized phase distribution, indicative of substantial phase variability (Figure 6c). After conventional processing, the same window shows noticeably tighter phase clustering (Figure 6d), with circular variance decreasing from 0.81 to 0.59. This reduction reflects improved phase coherence at that frequency and offset and demonstrates how phase variance provides a quantitative measure of processing impact on seismic phase behavior in field data.
"""

# ╔═╡ f434f74b-4b9e-4d95-b283-228452094247
function read_bin_seismic(file_id, rows, cols, ::Type{T}) where {T}
    url = "https://drive.google.com/uc?export=download&id=$file_id"
    response = HTTP.get(url)
    data = Matrix{T}(undef, rows, cols)
    read!(IOBuffer(response.body), data)
    return data
end

# ╔═╡ eac93f35-2256-499d-9f43-d92e88900fef
function trace_normalize(data::Matrix{T}) where {T <: AbstractFloat}
    normalized = similar(data)
    for i in 1:size(data, 2)
        max_val = maximum(abs.(data[:, i]))
        normalized[:, i] .= max_val != 0 ? data[:, i] ./ max_val : data[:, i]
    end
    return normalized
end

# ╔═╡ 61079dd0-875f-43ea-aad6-9ebb1e45826c
raw_supergather = read_bin_seismic("1NZ1hvEVtE7rhHpaMAmbeKTeiAUZuj47y", 149, 10000,
Float32);

# ╔═╡ c69d6959-c396-42b6-b08c-d5b8fade2cb3
processed_supergather = read_bin_seismic("1bPv6pk3eYzkHKU4tIiXmbGr9Uc26bIm4", 149,
10000, Float32);

# ╔═╡ d5a723cf-72d0-414d-842b-867e2d366289
t_seismic =range(start=0, step=0.002, length=149)

# ╔═╡ d0121087-e4e9-4a9b-8e29-32cc4ac0ea97
phase_raw, amplitude_raw, positive_freq_raw, f_zero_raw =
trace_spectrum(raw_supergather, T, extras=true)

# ╔═╡ 9edc3c0a-c93b-458e-a5b6-dd89ca189faf
variance_raw = local_phase_variance(phase_raw; batch_size=2000);

# ╔═╡ 6ac9ae73-269f-47ac-b80f-605b890b5b85
phase_processed, amplitude_processed, positive_freq_processed, f_zero_processed =
trace_spectrum(processed_supergather, T, extras=true)

# ╔═╡ 62bb18af-58f0-440b-8f61-80315f4d59f9
variance_processed = local_phase_variance(phase_processed; batch_size=2000);

# ╔═╡ 40e2b933-9409-49ed-9e4e-47151650aa95
begin
    n_traces_field = size(raw_supergather, 2)
    offset_axis = range(0, stop=22, length=n_traces_field)
    raw_norm = trace_normalize(raw_supergather)
    processed_norm = trace_normalize(processed_supergather)

    # ── 6a — raw seismic heatmap 
    p6a = heatmap(1:size(raw_norm, 2), t_seismic, raw_norm,
                  yflip=true, c=:grays, colorbar=false,
                  ylabel="Time (s)", xlabel="Offset (kft)",
                  title="Unprocessed Land Seismic Field Data")

    # ── 6b — processed seismic heatmap
    p6b = heatmap(1:size(processed_norm, 2), t_seismic, processed_norm,
                  yflip=true, c=:grays, colorbar=true,
                  ylabel="Time (s)", xlabel="Offset (kft)",
                  title="Processed Land Seismic Field Data")

    # ── 6c — raw phase heatmap 
    p6c = heatmap(offset_axis, positive_freq_raw[1:20],
                  phase_raw[1:20, 1:n_traces_field],
                  yflip=true, c=:coolwarm, colorbar=false, clims=(-π, π),
                  ylabel="Frequency (Hz)", xlabel="Offset (kft)",
                  title="Raw - Phase Ensemble")

    # ── 6d — processed phase heatmap 
    p6d = heatmap(offset_axis, positive_freq_processed[1:20],
                  phase_processed[1:20, 1:n_traces_field],
                  yflip=true, c=:coolwarm, colorbar=true, clims=(-π, π),
                  ylabel="Frequency (Hz)", xlabel="Offset (kft)",
                  title="Processed - Phase Ensemble")

    # ── 6g — raw phase variance heatmap 
    p6g = heatmap(offset_axis, positive_freq_raw[1:20],
                  variance_raw[1:20, :],
                  yflip=true, c=cgrad(:jet, rev=true), colorbar=false, clims=(0, 1),
                  ylabel="Frequency (Hz)", xlabel="Offset (kft)",
                  title="Raw - Circular Variance")

    # ── 6h — processed phase variance heatmap 
    p6h = heatmap(offset_axis, positive_freq_processed[1:20],
                  variance_processed[1:20, :],
                  yflip=true, c=cgrad(:jet, rev=true), colorbar=true, clims=(0, 1),
                  ylabel="Frequency (Hz)", xlabel="Offset (kft)",
                  title="Processed - Circular Variance")

    # ── layout 
    l_fig6 = @layout [
        a b
        c d
        e f
    ]

    plot(p6a, p6b, p6c, p6d, p6g, p6h,
         layout=l_fig6,
         size=(1600, 1200),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         top_margin=4Plots.mm,
         dpi=300)
end

# ╔═╡ b23715bb-fe6a-419e-8015-012109e8774c
md"""
To generalize these observations, we computed circular variance across all frequencies and offsets using overlapping sliding time windows (Figures 6g and 6h). Although reflectors appear visually clearer after conventional processing, the phase-variance maps reveal a more nuanced picture. At higher frequencies, circular variance remains close to one, indicating that phase behavior is still dominated by noise despite apparent amplitude enhancement. The spatial distribution of variance highlights systematic trends. Near offsets, which fall within the noise cone, exhibit consistently high phase variance, while mid to far offsets show progressively lower variance. For example, at 16 Hz, the unprocessed near-offset data yield a variance of 0.93, characteristic of noise-dominated phase behavior, whereas mid to far offsets decrease to about 0.71, indicating partial coherence. After processing, these same offsets show a marked reduction in variance to approximately 0.40, quantitatively confirming that conventional processing improves phase coherence primarily at lower frequencies and away from the noise cone (Figures 6j and 6l).

"""

# ╔═╡ e72b071e-be58-45a0-8c5d-53250d577d23
begin
    # circular stats for raw and processed at 16 Hz (freq index 5)
    mean_c_16_raw       = circ_mean(phase_raw[5, 2500:7500]).μ
    kappa_c_16_raw      = kappa_from_R(circ_r(phase_raw[5, 2500:7500]))
    mean_c_16_processed = circ_mean(phase_processed[5, 2500:7500]).μ
    kappa_c_16_processed = kappa_from_R(circ_r(phase_processed[5, 2500:7500]))

    # 6c inset — raw circular histogram
    p6c_circ = circ_hist_with_vm(phase_raw[5, 2500:7500], 2π/50,
                                  mean_c_16_raw, kappa_c_16_raw;
                                  title_str="16 Hz (raw)", show_vm=true)
    add_circular_mean_arrow!(p6c_circ, mean_c_16_raw)
    add_circle_annotations!(p6c_circ)

    # 6d inset — processed circular histogram
    p6d_circ = circ_hist_with_vm(phase_processed[5, 2500:7500], 2π/50,
                                  mean_c_16_processed, kappa_c_16_processed;
                                  title_str="16 Hz (processed)", show_vm=true)
    add_circular_mean_arrow!(p6d_circ, mean_c_16_processed)
    add_circle_annotations!(p6d_circ)

    plot(p6c_circ, p6d_circ,
         layout=(1, 2),
         size=(800, 400),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         dpi=300)
end

# ╔═╡ a9b554b5-4fe6-455f-aa9a-fb001827619e
md"""
This trend is consistent with stack-based signal-to-noise ratio estimates (Figure 7a), computed following Bakulin et al. (2022b). Prestack SNR improves from roughly ``-22`` dB to ``-12`` dB after processing, reflecting substantial suppression of random noise. However, SNR averages energy over the full frequency band and therefore masks important broadband behavior. In contrast, circular variance isolates phase variability directly and reveals that improvements are confined to a limited frequency range.
"""

# ╔═╡ c2450094-75fe-4cc8-9a3c-028e30ec2f36
begin
    # ── raw stats 
    mean_var_near_raw  = round(mean(vec(variance_raw[5, 1:2000])),    digits=2)
    mean_var_far_raw   = round(mean(vec(variance_raw[5, 6000:8000])), digits=2)

    p6i = histogram(vec(variance_raw[5, 1:3000]),
                    label="Near offsets", title="V̄ = $mean_var_near_raw",
                    fill="#FF0000", bins=10, grid=false, normalize=true,
                    legend=:topleft,
                    xtickfontsize=10, ytickfontsize=10)
    xlims!(p6i, 0.0, 1.0); ylims!(p6i, 0, 48)

    p6j = histogram(vec(variance_raw[5, 6000:8000]),
                    label="Far offsets", title="V̄ = $mean_var_far_raw",
                    fill="#4EA72E", bins=10, grid=false, normalize=true,
                    xtickfontsize=10, ytickfontsize=10)
    xlims!(p6j, 0.0, 1.0); ylims!(p6j, 0, 48)

    # ── processed stats 
    mean_var_near_processed = round(mean(vec(variance_processed[5, 1:3000])),    digits=2)
    mean_var_far_processed  = round(mean(vec(variance_processed[5, 6000:8000])), digits=2)

    p6k = histogram(vec(variance_processed[5, 1:3000]),
                    label="Near offsets", title="V̄ = $mean_var_near_processed",
                    fill="#FF0000", bins=10, grid=false, normalize=true,
                    legend=:topleft,
                    xtickfontsize=10, ytickfontsize=10)
    xlims!(p6k, 0.0, 1.0); ylims!(p6k, 0, 52)

    p6l = histogram(vec(variance_processed[5, 6000:8000]),
                    label="Far offsets", title="V̄ = $mean_var_far_processed",
                    fill="#4EA72E", bins=10, grid=false, normalize=true,
                    xtickfontsize=10, ytickfontsize=10)
    xlims!(p6l, 0.0, 1.0); ylims!(p6l, 0, 52)

    plot(p6i, p6j, p6k, p6l,
         layout=(1, 4),
         size=(1400, 350),
         left_margin=6Plots.mm,
         bottom_margin=6Plots.mm,
         dpi=300)
end

# ╔═╡ bb960156-be5a-4a0e-b43a-c66b541c8850
md"""
Amplitude spectra before and after processing (Figure 7b) show that conventional workflows often boost high frequencies through repeated deconvolution and filtering, giving the appearance of increased bandwidth. Amplitude-based measures alone, however, do not capture phase stability and can therefore be misleading. Phase-variance spectra (Figure 7c) provide a complementary and more physically meaningful perspective. By examining circular variance as a function of frequency, averaged over offsets, it becomes possible to define an effective frequency band in terms of phase coherence rather than amplitude alone. In practice, one can select a threshold variance corresponding to acceptable phase stability for a given application. Frequencies below this threshold retain coherent phase and are suitable for phase-sensitive workflows such as full-waveform inversion and AVO analysis, whereas frequencies above it contributes primarily incoherent energy.
"""

# ╔═╡ b70ac738-fa0b-4606-8070-d86666cd66f8
function calculate_semblance(window)
    N, M = size(window)
    numerator   = sum((sum(window[i, j] for j in 1:M))^2 for i in 1:N)
    denominator = M * sum(sum(window[i, j]^2 for j in 1:M) for i in 1:N)
    return numerator / denominator
end

# ╔═╡ bae896e2-30e7-4580-a360-f275aaa22b5c
function compute_snr(data; step_size=1, xdims=2000, ydims=149)
    rows = size(data, 1)
    cols = size(data, 2)
    snr_values = Float64[]
    for i in 1:step_size:rows-ydims+1
        for j in 1:step_size:cols-xdims+1
            window = data[i:i+ydims-1, j:j+xdims-1]
            sem    = calculate_semblance(window)
            snr    = sem / (1 - sem)
            snr_db = 10 * log10(snr)
            push!(snr_values, snr_db)
        end
    end
	 return snr_values
end

# ╔═╡ 5e22327e-e153-4f3f-8932-dee5bb764cf6
begin
    SNR_value_raw       = compute_snr(raw_supergather)
    SNR_value_processed = compute_snr(processed_supergather)

    p7a = plot(SNR_value_raw, lw=3, label="Unprocessed", color=:black,
               ylabel="SNR (dB)", xlabel="Trace Index", grid=false,
               fontsize=12, xtickfontsize=12, ytickfontsize=12,
               guidefontsize=12, legendfontsize=12, size=(1000, 400))
    plot!(p7a, SNR_value_processed, lw=3, label="Processed",
          color=:black, linestyle=:dash)

    p7a
end

# ╔═╡ 588f8c8c-c3a7-430d-b4e8-65ed743c1e42
begin
    # ── amplitude spectra ─────────────────────────────────────────────────────
    amplitude_av_raw       = sum(amplitude_raw, dims=2) / size(amplitude_raw, 2)
    amplitude_av_db_raw    = 20 .* log10.(amplitude_av_raw)
    amplitude_av_db_raw    = amplitude_av_db_raw .- maximum(amplitude_av_db_raw)

    amplitude_av_processed    = sum(amplitude_processed, dims=2) / size(amplitude_processed, 2)
    amplitude_av_db_processed = 20 .* log10.(amplitude_av_processed)
    amplitude_av_db_processed = amplitude_av_db_processed .- maximum(amplitude_av_db_processed)

    p7b = plot(positive_freq_raw, amplitude_av_db_raw,
               label="Unprocessed", linewidth=3, color=:black, grid=false,
               xlabel="Frequency (Hz)", ylabel="Amplitude (dB)",
               fontsize=12, xtickfontsize=12, ytickfontsize=12,
               guidefontsize=12, legendfontsize=12)
    plot!(p7b, positive_freq_raw, amplitude_av_db_processed,
          label="Processed", linewidth=3, color=:black, linestyle=:dash)
    xlims!(p7b, 0, 80)
    ylims!(p7b, minimum(amplitude_av_db_processed), 0)

    # ── phase variance spectra ────────────────────────────────────────────────
    r_raw_5000       = [1 - circ_r(phase_raw[i, 2500:7500])       for i in 1:size(phase_raw, 1)]
    r_processed_5000 = [1 - circ_r(phase_processed[i, 2500:7500]) for i in 1:size(phase_processed, 1)]

    p7c = plot(positive_freq_raw, r_raw_5000,
               label=false, lw=3, color=:black, grid=false,
               xlabel="Frequency (Hz)", ylabel="Circular (Phase) Variance",
               fontsize=12, xtickfontsize=12, ytickfontsize=12,
               guidefontsize=12, legendfontsize=12)
    plot!(p7c, positive_freq_processed, r_processed_5000,
          label=false, lw=3, linestyle=:dash, color=:black)
    xlims!(p7c, 0, 80)
    ylims!(p7c, -0.1, 1)

    plot(p7b, p7c,
         layout=(1, 2),
         size=(1000, 400),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         dpi=300)
end

# ╔═╡ 5b0e55a7-200a-4d1a-a45e-d51be98c8fbc
md"""
These observations suggest that seismic processing should be evaluated not only by amplitude and SNR metrics but also by its ability to systematically reduce phase variance. Monitoring phase variance provides a direct means to assess whether a given processing step improves or degrades phase integrity and offers a quantitative guide for parameter selection. Although we examine only two stages here, raw data and fully time-processed data, the same framework can be applied step by step throughout the processing sequence.

Identifying frequency bands with acceptable phase variance naturally leads to strategies for mitigating phase variability within those bands. In practice, seismic processing aims not only to diagnose phase distortions but also to compensate for them. Bakulin et al. (2023) showed that when the underlying signal is consistent across traces, the phase of the clean signal can be recovered using phase masking based on local ensembles, because the phase spectrum of the ensemble expectation matches that of the signal. In this case, the circular mean phase provides a convenient, amplitude-independent estimate of the signal phase.
"""

# ╔═╡ 249b65b6-44e6-41b2-8efe-03add61dea65
md"""
## **DISCUSSION**
"""

# ╔═╡ 7f58875c-f87b-4fd0-893b-b5557e2b0365
md"""
Treating seismic phase as a circular variable reframes phase variability from a difficult, trace-level attribute into a measurable ensemble property that can be evaluated systematically across frequency and offset. Circular variance, here called phase variance ``V(\omega)``, operates directly on wrapped phases and avoids the ambiguities of phase unwrapping. Low ``V(\omega)`` corresponds to tightly clustered phases and coherent signal behavior, whereas high ``V(\omega)`` indicates strong phase dispersion caused by noise, scattering, or other incoherent effects. In addition, the circular mean phase provides a useful estimate of the underlying signal phase within an ensemble. Although not the focus of this study, this estimate may support additional quality control of phase trends and serve as a key ingredient in phase-enhancement approaches such as seismic time-frequency masking (Bakulin et al., 2023) and phase-based filtering.

Synthetic tests confirm that ``V(\omega)`` reliably recovers imposed phase variability as a function of frequency and offset. Application to field data shows that conventional time-domain processing primarily reduces phase variance at low frequencies and mid-to-far offsets. At higher frequencies, phase variance often remains elevated, even when amplitude spectra suggest apparent bandwidth extension. This discrepancy highlights a key limitation of amplitude-based diagnostics and motivates defining an effective bandwidth in terms of phase coherence rather than amplitude alone.

Phase variance maps and spectra provide a direct way to identify where phase information becomes reliable after processing. By selecting an application-dependent threshold on ``V(\omega)``, practitioners can objectively delineate the frequency range suitable for phase-sensitive workflows, such as AVO analysis, full-waveform inversion, and migration. Unlike visual inspection or SNR metrics, ``V(\omega)`` isolates phase behavior and reveals frequency-dependent limitations that would otherwise remain hidden. The statistical stability of ``V(\omega)`` and its dependence on ensemble size are examined in Appendix A.
"""

# ╔═╡ 0ad0323c-6c8f-4a0a-984f-47d574a974fc


# ╔═╡ e5bedb7e-3cb2-48a5-b27c-196102bfcfbb
md"""
# **Conclusions**
"""

# ╔═╡ 7b6fca1e-e3ae-49b9-be66-75f9dcae238d
md"""
We have introduced circular statistics as a new and necessary framework for analyzing seismic phases. By treating phase as a circular variable rather than a trace-level artifact, phase variability becomes a measurable ensemble property that can be quantified, tracked, and acted upon systematically. Circular variance is not merely another quality-control attribute but a physically grounded descriptor of phase coherence that directly reflects the reliability of phase information across frequency and offset.

Our results show that amplitude-based diagnostics alone can substantially overestimate usable bandwidth. Conventional processing often amplifies high-frequency amplitudes while leaving phase behavior noise-dominated. Phase variance directly exposes this limitation, enabling the definition of an effective bandwidth based on phase coherence rather than spectral amplitude. This represents a fundamental shift in how seismic data quality is assessed for phase-sensitive applications.

We propose using phase variance as a standard prestack diagnostic, evaluated alongside amplitude spectra and SNR. Tracking ``V(\omega)`` through a processing sequence provides an objective way to identify where phase integrity improves or deteriorates, to guide parameter selection, and to prevent the propagation of phase distortions that would otherwise remain invisible. Importantly, this diagnostic view naturally enables phase-conditioning and masking strategies driven by explicit objectives, such as achieving a prescribed level of phase fidelity for imaging or inversion, rather than by visual QC or amplitude-based criteria.

More broadly, this work establishes circular statistics as a practical and extensible foundation for seismic phase analysis and provides the statistical basis for future phase-enhancement and masking approaches that operate directly on ensemble phase behavior. The framework applies equally to QC, processing design, and phase-sensitive imaging and inversion. We anticipate that phase-based diagnostics and conditioning methods built on circular statistics will become increasingly important as seismic workflows move toward higher frequencies, denser sampling, and greater sensitivity to subtle phase effects.
"""

# ╔═╡ 423e4c92-8645-4174-bec6-a568921b0d72
md"""
# **Acknowledgments**
"""

# ╔═╡ 89bbd044-b74b-4245-ba0b-20b5673d15f5
md"""
We  acknowledge  Fairfield  Geotechnologies  for  granting  permission to use the data presented in this study.
"""

# ╔═╡ 237f7c93-a97f-4e35-8071-1964bc78c207
md"""
# **References**
"""

# ╔═╡ 22f6de88-8754-4513-9c6d-09be0234d81f
md"""
Abbott, J. G., and T. H. Thurstone, 1979, Acoustic speckle: Theory and experimental analysis: Ultrasonic Imaging, 1, 303–324, doi: 10.1016/0161-7346(79)90021-2.

Bakulin, A., I. Silvestrov, I. Dmitriev, D. Neklyudov, M. Protasov, K. Gadylshin, and others, 2020, Nonlinear beamforming for enhancement of 3D prestack land seismic data: Geophysics, 85, 1MJ–Z13, doi: 10.1190/geo2019-0631.1.

Bakulin, A., D. Neklyudov, and I. Silvestrov, 2022a, Multiplicative seismic noise caused by small-scale near-surface scattering and its transformation during stacking: Geophysics, 87(5), V419–V435, doi: 10.1190/geo2021-0632.1.

Bakulin, A., I. Silvestrov, and M. Protasov, 2022b, Signal-to-noise ratio computation for challenging land single-sensor seismic data: Geophysical Prospecting, 70(3), 629–638, doi: 10.1111/1365-2478.13197.

Bakulin, A., D. Neklyudov, and I. Silvestrov, 2023, Seismic time-frequency masking for suppression of seismic speckle noise: Geophysics, 88(5), V371–V385, doi: 10.1190/geo2022-0701.1.

Bakulin, A., D. Neklyudov, and I. Silvestrov, 2024, The impact of receiver arrays on suppressing seismic speckle scattering noise caused by meter-scale near-surface heterogeneity: Geophysics, 89, V551–V565, doi: 10.1190/geo2023-0489.1.

Bakulin, A., A. Rohatgi, and S. Fomel, 2025, Statistical analysis of seismic phase variability in dense data: 86th EAGE Annual Conference & Exhibition, EAGE Technical Program Expanded Abstracts, doi: 10.3997/2214-4609.2025101166.

Cary, P., and N. Nagarajappa, 2014a, Robust surface-consistent residual statics and phase correction — Part 1: GeoConvention 2014, FOCUS: Adapt, Refine, Sustain, Calgary, Alberta, Canada.

Cary, P., and N. Nagarajappa, 2014b, Robust surface-consistent residual statics and phase correction — Part 2: GeoConvention 2014, FOCUS, Calgary, Alberta, Canada.

Feigl, K. L., R. M. Smith, F. D. Behn, and A. F. Sheehan, 2009, A method for modeling radar interferograms without phase unwrapping: Geophysical Journal International, 176, 491–497, doi: 10.1111/j.1365-246X.2008.03931.x.

Fisher, N. I., 1993, Statistical analysis of circular data: Cambridge University Press.

Fomel, S., and van der Baan, M., 2014, Local skewness attribute as a seismic phase detector: Interpretation, 2, SA49–SA56, doi: 10.1190/INT-2013-0080.1.

Goodman, J. W., 2007, Speckle phenomena in optics: Theory and applications: Roberts & Company.

Gosselin, J. M., 2022, Probabilistic inversion of circular phase spectra: Application to two-station phase-velocity dispersion estimation: Geophysical Journal International, 229, 1492–1517, doi: 10.1093/gji/ggac021.

Holt, R., and Lubrano, A., 2020, Stabilizing the phase of onshore 3D seismic data: Geophysics, 85(6), V473–V479, doi: 10.1190/geo2019-0695.1.

Jiang, Y., J. Hu, Z. Li, and Q. Liu, 2019, Bayesian inversion of wrapped satellite interferometric phase: Geophysical Journal International, 219, 1500–1518, doi: 10.1093/gji/ggz379.

Mardia, K. V., and P. E. Jupp, 2000, Directional statistics: John Wiley & Sons.

Meunier, J., 2011, Seismic acquisition from yesterday to tomorrow: 2011 Distinguished Instructor Course, Distinguished Instructor Series, No. 14, SEG and EAGE.

Mohammad, H. H., S. Z. Satari, and W. N. S. W. Yusoff, 2021, Review on circular–linear regression models: Journal of Physics: Conference Series, 1988, 012108, doi: 10.1088/1742-6596/1988/1/012108.

Newman, P., and J. T. Mahoney, 1973, Patterns—With a pinch of salt: Geophysical Prospecting, 21(2), 197–219, doi: 10.1111/j.1365-2478.1973.tb00023.x.

Nugraha, A. A., K. Kobayashi, T. Toda, S. Sakti, and S. Nakamura, 2019, A deep generative model of speech complex spectrograms: IEEE/ACM Transactions on Audio, Speech, and Language Processing, 27, 63–76, doi: 10.1109/TASLP.2018.2879810.

Reilly, J. M., A. P. Shatilo, and Z. J. Shevchek, 2010, The case for separate sensor processing: Meeting the imaging challenge in a producing carbonate field in the Middle East: The Leading Edge, 29, 1240–1249, doi: 10.1190/1.3496914.

Rohatgi, A., A. Bakulin, and S. Fomel, 2024, Analyzing the impact of additive and multiplicative noise on seismic data analysis: Fourth International Meeting for Applied Geoscience & Energy, SEG Technical Program Expanded Abstracts, doi: 10.1190/image2024-4086176.1.

Rohatgi, A., A. Bakulin, and S. Fomel, 2025, Data-driven analysis of seismic phase using circular statistics: The Leading Edge, 44, 683–691, doi: 10.1190/tle44090683.1.

Sato, H., M. C. Fehler, and T. Maeda, 2012, Seismic wave propagation and scattering in the heterogeneous Earth, 2nd ed.: Springer, doi: 10.1007/978-3-642-23029-5.

Schimmel, M., and J. Gallart, 2007, Frequency-dependent phase coherence for noise suppression in seismic array data: Journal of Geophysical Research: Solid Earth, 112, B04303, doi: 10.1029/2006JB004680.

Schimmel, M., and H. Paulssen, 1997, Noise reduction and detection of weak, coherent signals through phase-weighted stacks: Geophysical Journal International, 130, 497–505, doi: 10.1111/j.1365-246X.1997.tb05664.x.

Shatilo, A. P., 1992, Seismic phase unwrapping: Methods, results, problems: Geophysical Prospecting, 40, 211–225, doi: 10.1111/j.1365-2478.1992.tb00080.x.

Stork, C., 2020, How does the thin near surface of the Earth produce 10–100 times more noise on land seismic data than on marine data?: First Break, 38, 67–75, doi: 10.3997/1365-2397.fb2020062.

Takamichi, S., K. Kobayashi, T. Toda, and H. Saruwatari, 2018, Phase reconstruction from amplitude spectrograms based on a von Mises-distribution deep neural network: Interspeech 2018 Proceedings, 2027–2031.

Taner, M. T., F. Koehler, and K. A. Alhilali, 1974, Estimation and correction of near-surface time anomalies: Geophysics, 39, 441–463, doi: 10.1190/1.1440441.

Taner, M. T., and F. Koehler, 1981, Surface-consistent corrections: Geophysics, 46, 17–22, doi: 10.1190/1.1441133.

van der Baan, M., 2008, Time-varying wavelet estimation and deconvolution by kurtosis maximization: Geophysics, 73, V11–V18, doi: 10.1190/1.2831936.

van der Baan, M., and Fomel, S., 2009, Nonstationary phase estimation using regularized local kurtosis maximization: Geophysics, 74, A75–A80, doi: 10.1190/1.3213533.

von Mises, R., 1918, Über die "Ganzzahligkeit" der Atomgewichte und verwandte Fragen: Physikalische Zeitschrift, 19, 490–500.

Zebker, H. A., and Y. Lu, 1998, Phase unwrapping algorithms for radar interferometry: A systematic comparison: Journal of the Optical Society of America A, 15, 586–600, doi: 10.1364/JOSAA.15.000586.

Ziolkowski, A., P. Hanssen, R. Gatliff, X.-Y. Li, and H. Jakubowicz, 2003, Use of low frequencies for sub-basalt imaging: Geophysical Prospecting, 51, 169–182, doi: 10.1046/j.1365-2478.2003.00363.x.
"""

# ╔═╡ f85c0466-490c-4f2a-b189-277c74657dfa
md"""
# **APPENDIX A: EFFECT OF ENSEMBLE SIZE ON THE STABILITY OF PHASE VARIANCE ESTIMATE**
"""

# ╔═╡ de8a7d62-7429-47c4-8bb0-94ba01803536
md"""
Circular variance is a statistical quantity estimated from a finite ensemble of phase samples. As such, the variance itself is a random variable whose stability depends on the number of traces included in the ensemble. Importantly, small ensemble sizes do not render phase-variance estimates uninterpretable; rather, they increase estimator variability and make interpretation more sensitive to local fluctuations. Interpretable phase-variance estimates therefore require an ensemble size appropriate to the objective of the analysis. To assess the effect of ensemble size, we computed circular variance using sliding windows containing 10, 100, and 1000 traces (Figure A1). For small ensembles, the estimated variance exhibits strong local fluctuations due to limited sampling, revealing fine-scale detail but with increased statistical scatter. As the number of traces increases, these fluctuations are progressively suppressed, and the estimated variance converges toward a smooth and stable trend that reflects the underlying phase dispersion rather than estimator noise. Importantly, the large-scale trend in phase variability is preserved across all window sizes, while increasing the ensemble size improves statistical stability.

This behavior is expected and reflects the nested nature of the problem: individual phase values are random, and circular variance is itself an estimate derived from those random samples. Consequently, the ensemble size must be chosen based on the analysis objective. When the goal is to resolve only the underlying trend in phase variability, large ensembles are required, particularly when circular variance is used as a diagnostic or quality-control attribute where stability and robustness are critical. Smaller ensembles may be useful when emphasizing localized structure, but at the expense of increased variability and interpretational uncertainty.

This distinction underpins the consistent use of different ensemble sizes throughout this study and ensures that circular variance is applied appropriately to each task.
"""

# ╔═╡ 82900dcf-62a0-41b6-8b47-33d5dab19744
begin
	variance_perturbed_10   = local_phase_variance(phase_offset; batch_size=10)
    variance_perturbed_100  = local_phase_variance(phase_offset; batch_size=100)
    variance_perturbed_1000 = local_phase_variance(phase_offset; batch_size=1000)

    variance_out_10   = variance_perturbed_10[6, :]
    variance_out_100  = variance_perturbed_100[6, :]
    variance_out_1000 = variance_perturbed_1000[6, :]

    p_appendix = plot(variance_out_10, label="N = 10", lw=0.5,
                      linestyle=:dash, color=:black, alpha=0.5,
                      ylabel="Phase Variance", xlabel="Trace Index",
                      grid=false, fontsize=12, labelfontsize=12,
                      legendfontsize=12, tickfontsize=12,
                      xguidefontsize=12, yguidefontsize=12,
                      size=(1000, 300))

    plot!(p_appendix, variance_out_100, label="N = 100",
          lw=1, linestyle=:dash, color=:black, alpha=0.7)

    plot!(p_appendix, variance_out_1000, label="N = 1000",
          lw=3, color=:black, grid=false)

    plot!(p_appendix, var_imposed, label="Imposed",
          lw=3, color="#499af2", grid=false)

    p_appendix
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CircStats = "2f6764a1-d620-4564-9394-76eb7c776766"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.10.15"
CircStats = "~1.0.4"
DataFrames = "~1.8.1"
Distributions = "~0.25.123"
FFTW = "~1.10.0"
HTTP = "~1.11.0"
Plots = "~1.41.5"
SpecialFunctions = "~2.7.2"
StatsBase = "~0.34.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "fa6dc16daadd6fe0ad553667673f0e965b26ee13"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "856ecd7cebb68e5fc87abecd2326ad59f0f911f3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.43"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a21c5464519504e41e0cbc91f0188e8ca23d7440"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+1"

[[deps.CircStats]]
deps = ["Distributions", "HypothesisTests", "LinearAlgebra", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "ecfe2e9a260c4723026b4a71460cf0420def9e40"
uuid = "2f6764a1-d620-4564-9394-76eb7c776766"
version = "1.0.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "c761b00e7755700f9cdf5b02039939d1359330e1"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.1.0"

[[deps.CommonSolve]]
git-tree-sha1 = "78ea4ddbcf9c241827e7035c3a03e2e456711470"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.6"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Printf", "Random", "Roots", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "b67985bd11331ccef26109a6269dbaae01474a72"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.11.6"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c9cbeda6aceffc52d8a0017e71db27c7a7c0beaf"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.5+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "1cc8ad0762e59e713ee3ef28f9b78b2c9f4ca078"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.5"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "c5a07210bd060d6a8491b0ccdee2fa0235fc00bf"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.1.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "10a488dbecb88a9679c8f357d383d7d83dcc748d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.13"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2700b235561b0335d5bef7097a111dc513b8655e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.2"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a3c1536470bf8c5e02096ad4853606d7c8f62721"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.2"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "28145feabf717c5d65c1d5e09747ee7b1ff3ed13"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.3"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e015f211ebb898c8180887012b938f3851e719ac"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.55+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─3b71e1f0-0c4a-11f1-0aa6-43a57fc47120
# ╟─67f35fa0-dfa6-424d-ba6b-7e9005936ff1
# ╟─f40c721d-f78c-41e8-97ed-efd0625d88c3
# ╟─a6d2d852-80d0-4074-83e2-53854ad1e5de
# ╟─a400b618-d29d-46bd-87e0-09f4b92bc574
# ╟─cd96e0e5-e393-4f43-a7a2-ec0e15d95784
# ╟─a911602a-79e2-4b77-9be7-cdf2462dda18
# ╟─dc709536-5ebb-4d6e-b196-1040262c5f7a
# ╠═b8a71a38-ee02-48d7-8ceb-1e9fa411cafe
# ╠═4f9d39a7-a42f-4025-b4a4-23861b7474f8
# ╠═3fd44e63-acf6-4be1-af3e-c80c5675dee8
# ╠═e914ad8a-969a-4078-8786-7d0c7c686c23
# ╠═b4515d5f-f5e6-4f2d-af72-b501b2c9902b
# ╠═666ed7c7-b3bf-48f9-b7ec-13b8b22e6853
# ╠═541d34f6-8957-42ba-8935-76bf644be6ed
# ╠═22cd54ed-1fde-4fd3-a86f-23d5314792f6
# ╠═55c749dd-532e-47ec-902c-55e050daa829
# ╠═db63ff6f-86f4-441e-82df-4db7339c1cf3
# ╠═6b4e6ef5-776c-48f2-a805-5d0e71f7274d
# ╠═a5e0c66e-aa8a-4c0c-988b-0e27f2ef2cba
# ╠═081e2585-fae8-47ad-b055-6e58f9d4dc34
# ╠═a2205f1e-e936-4124-beb6-786dbc522fb2
# ╠═8cccd756-0e8e-46f8-bdb9-698b9e90b62c
# ╠═f9bf018d-b213-4673-ad05-10c1384b0f45
# ╠═2715160b-7b81-4ebc-9d36-853be73bf925
# ╠═38b2cb6b-a303-49dc-b60c-7085f48b9648
# ╠═02d68919-6983-47f0-bb86-2707c2485de5
# ╠═e18a3bea-f46e-49e5-bac7-5da9a532ffee
# ╠═954c55bc-bb1b-43e7-a364-b486ea070114
# ╠═fc174366-77a5-4e98-855e-e38439a8507b
# ╠═f5138299-f504-49ea-afa6-ec040596dbf0
# ╠═216596c5-bfa8-4abf-aebb-cc7ce3ad0d4c
# ╠═54c4cabf-b788-4512-94ec-8d32c4c9a5b7
# ╠═68879a15-a986-4985-8c05-73ccb092c990
# ╟─5b0b748d-126e-4615-b909-40adbff7f685
# ╟─aa88a1e0-4586-443e-9d0f-1c1a29eb330b
# ╟─95c4fc77-9389-4576-a22e-ff2ada4c0ca6
# ╟─cc90f86b-3d00-4ce2-abf6-4eb5f7490c9c
# ╠═08a4cdc5-ff4a-44c8-96d1-4f3701e32314
# ╟─6a3bf061-1449-4d9f-b258-0f7318c3b837
# ╟─79b745b1-8ed8-4ccc-850e-658fd52145e3
# ╟─d66e5b68-0fd1-4e96-855b-8dd2b42977fd
# ╟─21990862-4208-46a0-84cb-f2478cbeccdf
# ╠═60a4d061-ca42-422b-b3f8-952b4c89453a
# ╟─62f2f1f8-284a-471c-b348-e3a10242476a
# ╟─2ac31116-3f26-428b-87d3-80c566af2c19
# ╟─c883fcfb-1b21-4836-8ef4-2ce3321d5369
# ╟─f6e2cacb-c0ca-44ed-827c-4b7fcaa2b0b0
# ╠═1cbc72ad-b5e2-4709-8006-cfd5cb8eacbf
# ╠═8315a8e9-bbe4-480e-9b76-16d9de08401f
# ╟─1d99aa32-f8b7-481d-885c-3b46d1bf48fd
# ╟─c28e595e-3857-46c6-929b-44a3ba6806e9
# ╟─5468f447-9a9a-4fdc-a6a4-7d080434a841
# ╟─b00bf530-e245-45cb-8f87-1c33c394ef7f
# ╟─0f37ae56-928d-4a73-bf0a-81262d1ea982
# ╟─7cb4d395-b266-46c4-8e3f-787c18b075d5
# ╟─4022a458-2e6b-4b79-92cf-575ac2fd63e9
# ╠═87a0577d-0a4d-4e6e-b5b7-ecbee0f7cf51
# ╟─248b59c9-e468-41c6-a6ab-ca582d94c48e
# ╟─bf0ed95f-cc64-418d-9fde-851657638a2f
# ╟─b0b90f69-c265-4c54-bd86-5d62e97b04f7
# ╠═f6523673-ad48-450f-8e95-4a8f4750dfa3
# ╟─433b4239-faab-4750-8228-fb36d7074c40
# ╟─8b382333-ff48-4174-b895-1cb8a5c28fff
# ╠═be3e9677-d5f8-44dd-94fd-f6de2d7a7a4a
# ╟─cf4f405c-f54f-4b6d-ac5a-b9f7c36c27f8
# ╟─48e8a718-f908-4541-af03-af0302dc83eb
# ╠═2d0e8b6a-c822-4447-9a52-761e253f5a29
# ╟─470bee5b-4ea5-4a02-ba24-a508968cbf0c
# ╟─ef5433c7-c551-452c-8d55-545629f1d26a
# ╟─d58cee5e-b7f0-4daa-b17d-3e5b498b7ef8
# ╟─d36129ee-22ae-463f-b133-f2e596544edf
# ╟─78e5e4f0-2422-44c7-9ae3-1421aedbbbd7
# ╠═79eec12e-ecaa-42ae-87fe-d80f636c0a1c
# ╠═54c9b2ec-0c9d-4f23-b615-0da0d4023b9c
# ╠═988b7168-7c2a-4cae-a21c-0487b1477073
# ╟─2b0dfbbc-33de-4aed-bbfc-c49f66055c33
# ╠═2be746d9-d942-4e29-933d-c28d0dec843e
# ╟─1b31cff4-9f96-490c-97ec-45004449ae7e
# ╠═bb3bc3d2-6f4a-4021-997a-2c0613ca814d
# ╠═bbedd657-28fa-42b0-9418-296bf80feead
# ╠═43de8089-ac37-4b06-ab88-750d152e367a
# ╠═fd33d1bf-3ef2-474c-ae4e-df30b18e8157
# ╟─efb01d36-8fe3-413b-9edb-d13c574acd28
# ╟─9763d963-8e50-42e6-b9f5-eedfaa60571e
# ╠═f434f74b-4b9e-4d95-b283-228452094247
# ╠═eac93f35-2256-499d-9f43-d92e88900fef
# ╠═61079dd0-875f-43ea-aad6-9ebb1e45826c
# ╠═c69d6959-c396-42b6-b08c-d5b8fade2cb3
# ╠═d5a723cf-72d0-414d-842b-867e2d366289
# ╠═d0121087-e4e9-4a9b-8e29-32cc4ac0ea97
# ╠═9edc3c0a-c93b-458e-a5b6-dd89ca189faf
# ╠═6ac9ae73-269f-47ac-b80f-605b890b5b85
# ╠═62bb18af-58f0-440b-8f61-80315f4d59f9
# ╠═40e2b933-9409-49ed-9e4e-47151650aa95
# ╟─b23715bb-fe6a-419e-8015-012109e8774c
# ╠═e72b071e-be58-45a0-8c5d-53250d577d23
# ╟─a9b554b5-4fe6-455f-aa9a-fb001827619e
# ╠═c2450094-75fe-4cc8-9a3c-028e30ec2f36
# ╟─bb960156-be5a-4a0e-b43a-c66b541c8850
# ╠═bae896e2-30e7-4580-a360-f275aaa22b5c
# ╠═b70ac738-fa0b-4606-8070-d86666cd66f8
# ╠═5e22327e-e153-4f3f-8932-dee5bb764cf6
# ╠═588f8c8c-c3a7-430d-b4e8-65ed743c1e42
# ╟─5b0e55a7-200a-4d1a-a45e-d51be98c8fbc
# ╟─249b65b6-44e6-41b2-8efe-03add61dea65
# ╟─7f58875c-f87b-4fd0-893b-b5557e2b0365
# ╟─0ad0323c-6c8f-4a0a-984f-47d574a974fc
# ╟─e5bedb7e-3cb2-48a5-b27c-196102bfcfbb
# ╟─7b6fca1e-e3ae-49b9-be66-75f9dcae238d
# ╟─423e4c92-8645-4174-bec6-a568921b0d72
# ╟─89bbd044-b74b-4245-ba0b-20b5673d15f5
# ╟─237f7c93-a97f-4e35-8071-1964bc78c207
# ╟─22f6de88-8754-4513-9c6d-09be0234d81f
# ╟─f85c0466-490c-4f2a-b189-277c74657dfa
# ╟─de8a7d62-7429-47c4-8bb0-94ba01803536
# ╠═82900dcf-62a0-41b6-8b47-33d5dab19744
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
