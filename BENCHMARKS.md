# Benchmarks

Wall-clock measurements collected while developing this code, kept as a reference for anyone planning runs on similar hardware. All runs use WENO advection, 3rd-order Runge-Kutta time-stepping, and Float32 (Oceananigans' default on GPU), at CFL ≈ 0.8.

## Hardware tested

| GPU | Architecture | VRAM | Memory | Bandwidth | FP32 |
|---|---|---|---|---|---|
| RTX 2000 Ada | Ada Lovelace | 16 GB | GDDR6 | ~256 GB/s | ~12 TFLOPS |
| RTX 4080 SUPER | Ada Lovelace | 16 GB | GDDR6X | ~717 GB/s | ~49 TFLOPS |
| A100-40G SXM4 | Ampere | 40 GB | HBM2e | ~1555 GB/s | ~19.5 TFLOPS |

The A100 measurements were taken on a shared HPC cluster; the consumer cards on local workstations.

## Timings

Each cell shows time per iteration and total wall time to reach `t = 60`. VRAM in parentheses where measured.

| Grid | RTX 2000 Ada (16 GB) | RTX 4080 SUPER (16 GB) | A100-40G |
|---|---|---|---|
| 2000² | 0.30 s/iter | 0.11 s/iter — 28 min (~4 GB) | — |
| ~4000² | — | 0.38 s/iter — 3 hr (~4 GB) | 0.09 s/iter — 36 min |
| 8000² | — | ~1.2 s/iter — crashed at t ≈ 34 after 10 hr | 0.24 s/iter — 3.6 hr (~12 GB) |
| 12000² | — | — | 0.52 s/iter — 11.6 hr (~31 GB) |

The 12000² A100 run is the case shown in the paper.

## Memory budget

For a 2D Float32 simulation, each `Nx × Nz` field occupies `4·Nx·Nz` bytes — about 64 MiB at 4096². Observed VRAM at runtime is ~4 GiB, i.e. roughly 64 fields, far more than the four prognostic fields `(u, w, b, p)` you might expect.

The extras come from halos, "old" and "new" copies for time-stepping, RK3 substage tendencies, the pressure solver's residual and search-direction vectors, and per-direction temporary fluxes. As a rough rule of thumb, each prognostic field expands to ~15–20 backing arrays, so four prognostic fields give ~60–80 arrays. Switching to a less memory-hungry time-stepper (Heun, AB2) reduces this somewhat at a small accuracy cost.

This sets the upper grid size on each card: ~8000² is the practical limit on the 16 GB consumer cards (and is fragile, as the crashed 4080 run shows), 12000² is comfortable on the 40 GB A100.

## Memory bandwidth, not FLOPS, sets performance

Two observations stand out:

- The A100 (~19.5 TFLOPS FP32) is ~5× faster than the RTX 4080 SUPER (~49 TFLOPS FP32), despite less than half the peak compute throughput.
- The 4080 vs RTX 2000 measured ratio (~2.7×) sits on the **bandwidth ratio** (~2.8×), not the FLOPS ratio (~4.1×).

Both are consistent with this DNS being **memory-bandwidth limited**. Stencil derivatives, halo updates, and the pressure-Poisson solve have low arithmetic intensity (1–5 FLOPs per byte moved), so the GPU spends most of each timestep waiting on VRAM rather than computing. The relevant figure of merit is bandwidth: A100 HBM2e ≈ 1.55 TB/s vs 4080 GDDR6X ≈ 0.72 TB/s — a 2.2× advantage that accounts for most of the 5× wall-clock speedup. The remaining factor comes from the A100's better latency hiding (more concurrent warps and a larger register file feeding the cores while memory fetches complete), the wider 5120-bit HBM bus and physically stacked memory giving lower access latency than GDDR6X's 256-bit bus, and the SXM4 form factor sustaining full power and clocks indefinitely while the 4080 thermally throttles under continuous load.

A useful diagnostic on your own runs: if `nvidia-smi` shows the GPU pinned at ~100% utilisation but the iteration time scales with the card's memory bandwidth rather than its TFLOPS, you're memory-bound. NVIDIA Nsight Compute can confirm this by reporting a low arithmetic-intensity / high memory-stall fraction for the dominant kernels.

**Practical implication.** When choosing GPU hardware for this kind of 2D DNS, prioritise memory bandwidth and total VRAM over peak FP32 throughput. Data-centre cards with HBM memory (A100, H100, V100) are usually a much better fit than gaming cards of similar list price. A 32 GB V100 (HBM2, ~900 GB/s) would comfortably handle grids up to ~11500² and is a sensible fallback if A100/H100 access isn't available.
