/*
 * ca_cuda_sinc.h — C-linkage declarations for CUDA kernel wrappers
 *
 * Include from integrated.h when USE_CUDA is defined.
 */

#ifndef CA_CUDA_SINC_H_
#define CA_CUDA_SINC_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Allocate device grids and zero d_and_count */
void cuda_alloc_grids(void);

/* Upload host grids to device (call once after init) */
void cuda_upload_grid(void *h_grid, void *h_grid_next);

/* Run one pulse step on GPU (BFS + activation flags) */
void cuda_pulse_step(int tick);

/* Run one sinc step on GPU (wave eq + Bresenham + TTL + AND + damping copy-back) */
void cuda_sinc_step(int tick);

/* Copy grid from device to host for rendering / convergence check */
void cuda_download_grid(void *h_grid);

/* Copy and_count from device to host */
void cuda_download_and_count(int *h_and_count, int n);

/* Download AND2/AND3 per-tick counters from device */
void cuda_download_and_counters(int *h_and2, int *h_and3);

/* Reset per-tick AND counters on device to zero */
void cuda_reset_and_counters(void);

/* Upload full grid back to device (after convergence detection on host) */
void cuda_upload_grid_full(void *h_grid);

/* Free device memory */
void cuda_free(void);

#ifdef __cplusplus
}
#endif

#endif /* CA_CUDA_SINC_H_ */
