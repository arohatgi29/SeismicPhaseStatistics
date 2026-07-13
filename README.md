<div align="center">
  <img src="logo.png" alt="SeismicPhaseStatistics Logo" width="800"/>
</div>

Code and data accompanying publications on Data-Driven Seismic Phase Analysis using Circular Statistics.

Each subdirectory corresponds to a published or submitted paper and is self-contained: you can reproduce a paper's results by working entirely within its folder.

## Background: Phase Wrapping and Unwrapping

Instantaneous phase, by construction, is only defined modulo $2\pi$: the arctangent used to compute it from a complex (analytic) signal always returns a value in $(-\pi, \pi]$. As phase accumulates over time, it "wraps" around this interval every time it crosses a $\pm\pi$ boundary, producing the sawtooth pattern below. The true, continuously accumulating phase (unwrapped) is recovered by adding or subtracting multiples of $2\pi$ at each discontinuity so the sequence increases smoothly.

<div align="center">
  <img src="assets/Phase_wrap_unwrap.jpg" alt="Phase wrapping and unwrapping illustration" width="700"/>
  <p><em>Figure 3 from Rohatgi et al., 2025, <a href="https://doi.org/10.1190/tle44090683.1">The Leading Edge</a>.</em></p>
</div>

This distinction matters throughout the circular-statistics framework used in this repository: raw (wrapped) phase differences between traces or gathers can appear artificially large or discontinuous near the $\pm\pi$ boundary, which is exactly the kind of artifact that circular statistics (e.g., circular variance, mean resultant length) are designed to handle correctly — unlike linear statistics, which implicitly assume phase has already been unwrapped or lives on a line rather than a circle.

## Papers in this repository

- **Data-driven analysis of seismic phase using circular statistics**
  A. Rohatgi, A. Bakulin, and S. Fomel, *The Leading Edge*, 44(9), 683–691, 2025.
  
  [![DOI](https://img.shields.io/badge/DOI-10.1190%2Ftle44090683.1-2ea44f?style=flat-square)](https://doi.org/10.1190/tle44090683.1)
  [![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.20055837-1682D4?style=flat-square&logo=zenodo&logoColor=white)](https://doi.org/10.5281/zenodo.20055837)
  Folder: [`Rohatgi-et-al-2025-TLE/`](./Rohatgi-et-al-2025-TLE)

- **Amplitude-invariant phase masking for coherence recovery in scattered wavefields**
  A. Rohatgi, A. Bakulin, and S. Fomel — *JASA Express Letters*, 2026 (under revision).
  
  [![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.20059964-1682D4?style=flat-square&logo=zenodo&logoColor=white)](https://doi.org/10.5281/zenodo.20059964)
  Folder: [`Rohatgi-et-al-2026-JASA/`](./Rohatgi-et-al-2026-JASA)

<!-- Add future papers here as new entries -->

## How to reproduce a paper's results

1. Clone the repo:
```bash
   git clone https://github.com/arohatgi29/SeismicPhaseStatistics.git
   cd SeismicPhaseStatistics
```
2. Navigate to the relevant paper folder.
3. Follow the instructions in that folder's `README.md`.

Each paper folder lists its own dependencies and entry points. Software versions used for the published results are pinned in each folder.

## License

Released under the [MIT License](./LICENSE) unless noted otherwise within a paper folder.

## Contact

Akshika Rohatgi — [GitHub @arohatgi29](https://github.com/arohatgi29)

Issues and questions are welcome via the [Issues tab](https://github.com/arohatgi29/SeismicPhaseStatistics/issues).
