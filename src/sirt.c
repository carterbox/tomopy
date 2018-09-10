// Copyright (c) 2015, UChicago Argonne, LLC. All rights reserved.

// Copyright 2015. UChicago Argonne, LLC. This software was produced
// under U.S. Government contract DE-AC02-06CH11357 for Argonne National
// Laboratory (ANL), which is operated by UChicago Argonne, LLC for the
// U.S. Department of Energy. The U.S. Government has rights to use,
// reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
// UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
// ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
// modified to produce derivative works, such modified software should
// be clearly marked, so as not to confuse it with the version available
// from ANL.

// Additionally, redistribution and use in source and binary forms, with
// or without modification, are permitted provided that the following
// conditions are met:

//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.

//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in
//       the documentation and/or other materials provided with the
//       distribution.

//     * Neither the name of UChicago Argonne, LLC, Argonne National
//       Laboratory, ANL, the U.S. Government, nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago
// Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "utils.h"


void
sirt(
    const float *data, int dy, int dt, int dx,
	const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter)
{
    // int i, s, p, d, n; // preferred loop order
    // For each slice
    for (int s=0; s<dy; s++)
    {
        float *recon_slice = recon + s*ngridx*ngridy;
        float *gridx = (float *)malloc((ngridx+1)*sizeof(float));
        float *gridy = (float *)malloc((ngridy+1)*sizeof(float));
        assert(gridx != NULL && gridy != NULL);
        float mov;
        preprocessing(ngridx, ngridy, dx, center[s],
            &mov, gridx, gridy);
            // Outputs: mov, gridx, gridy
        float *all_dist, *all_sum_dist2;
        int *all_indi, *ray_start, *ray_stride;
        compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov,
            ngridx, ngridy, &ray_start, &ray_stride, &all_indi, &all_dist,
            &all_sum_dist2);
            // Outputs: ray_start, ray_stride, all_indi, all_dist
        free(gridx);
        free(gridy);
        // For each iteration
        for (int i=0; i<num_iter; i++)
        {
            float *simdata = calloc(dt*dx, sizeof *simdata);
            float *sum_dist = calloc(ngridx*ngridy, sizeof *sum_dist);
            float *update = calloc(ngridx*ngridy, sizeof *update);
            assert(simdata != NULL && sum_dist != NULL && update != NULL);
            // For each projection angle
            for (int p=0; p<dt; p++)
            {
                // For each detector pixel
                for (int d=0; d<dx; d++)
                {
                    int ray = d + dx*p;
                    float *dist = all_dist + ray_start[ray];
                    int *indi = all_indi + ray_start[ray];
                    float sum_dist2 = all_sum_dist2[ray];
                    if (sum_dist2 != 0.0)
                    {
                        // Calculate simdata
                        calc_simdata(0, p, d, ngridx, ngridy, dt, dx,
                            ray_stride[ray]+1, indi, dist, recon_slice,
                            simdata); // Output: simdata
                        // Update
                        int ind_data = d + dx*(p + dt*s);
                        int ind_sim = d + dx*p;
                        float upd = (data[ind_data]-simdata[ind_sim])/sum_dist2;
                        for (int n=0; n<ray_stride[ray]; n++)
                        {
                            update[indi[n]] += upd*dist[n];
                            sum_dist[indi[n]] += dist[n];
                        }
                    }
                }
            }
            for (int n = 0; n < ngridx*ngridy; n++) {
                if (sum_dist[n] > 0) {
                    recon_slice[n] += update[n]/sum_dist[n];
                }
            }
            free(simdata);
            free(sum_dist);
            free(update);
        }
        free(all_dist);
        free(all_sum_dist2);
        free(all_indi);
        free(ray_start);
        free(ray_stride);
    }
}


void
sirt_fly_rotation(
    const float *data, int dy, int dt, int dx,
    const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter, int nmask, int *mask)
{
    sirt_convolve(data, dy, dt, dx,
        center, theta,
        recon, ngridx, ngridy, num_iter,
        nmask, mask);
}


void
sirt_convolve(
    const float *data, int dy, int dt, int dx,
    const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter,
    int nmask, bool *mask)
{
    int step = 1;
    assert(step > 0 && "Step must be positive or else infinite loop.");
    // int i, s, p, b, d, n; // preferred loop order
    // For each slice
    for (int s=0; s<dy; s++)
    {
        float *recon_slice = recon + s*ngridx*ngridy;
        float *gridx = (float *)malloc((ngridx+1)*sizeof(float));
        float *gridy = (float *)malloc((ngridy+1)*sizeof(float));
        assert(gridx != NULL && gridy != NULL);
        float mov;
        preprocessing(ngridx, ngridy, dx, center[s],
            &mov, gridx, gridy);
            // Outputs: mov, gridx, gridy
        float *all_dist, *all_sum_dist2;
        int *all_indi, *ray_start, *ray_stride;
        compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov,
            ngridx, ngridy, &ray_start, &ray_stride, &all_indi, &all_dist,
            &all_sum_dist2);
            // Outputs: ray_start, ray_stride, all_indi, all_dist
        free(gridx);
        free(gridy);
        // For each iteration
        for (int i=0; i<num_iter; i++)
        {
            float *simdata = calloc(dt*dx, sizeof *simdata);
            assert(simdata != NULL);
            float *update = calloc(ngridx * ngridy, sizeof *update);
            float *nupdate = calloc(ngridx * ngridy, sizeof *nupdate);
            assert(update != NULL && nupdate != NULL);
            // For each projection angle, simulate data acquisition
            for (int p=0; p<dt; p++)
            {
                // For each detector pixel
                for (int d=0; d<dx; d++)
                {
                    int ray = d + dx*p;
                    float *dist = all_dist + ray_start[ray];
                    int *indi = all_indi + ray_start[ray];
                    if (all_sum_dist2[ray] != 0.0)
                    {
                        // Calculate simdata
                        calc_simdata(0, p, d, ngridx, ngridy, dt, dx,
                            ray_stride[ray]+1, indi, dist, recon_slice,
                            simdata); // Output: simdata
                    }
                }
            }
            // For each projection angle, pool data and compute updates
            for (int p=0; p<dt; p+=step)
            {
                // Initialize buffers to zero
                float *pool_sim = calloc(dx, sizeof *pool_sim);
                float *pool_sum_dist2 = calloc(dx, sizeof *pool_sum_dist2);
                assert(pool_sim != NULL && pool_sum_dist2 != NULL);
                // For each code element
                for (int b=0; b<nmask; b++)
                {
                    if (mask[b])
                    {
                        int p1 = (p+b) % dt;
                        // For each detector pixel
                        for (int d=0; d<dx; d++)
                        {
                            int ray = d + dx*(p1);
                            if (all_sum_dist2[ray] != 0.0)
                            {
                                // Calculate pool data
                                pool_sum_dist2[d] += all_sum_dist2[ray];
                                int ind_sim = d + dx*p1;
                                pool_sim[d] += simdata[ind_sim];
                            }
                        }
                    }
                }
                // For each code element
                for (int b=0; b<nmask; b++)
                {
                    if (mask[b])
                    {
                        int p1 = (p+b) % dt;
                        // For each detector pixel
                        for (int d=0; d<dx; d++)
                        {
                            if (pool_sum_dist2[d] > 0)
                            {
                                // Compute update
                                int ind_data = d + dx*(p1 + dt*s);
                                float pool_upd = (data[ind_data] - pool_sim[d])
                                                  / pool_sum_dist2[d];
                                // Update
                                int ray = d + dx*(p1);
                                float *dist = all_dist + ray_start[ray];
                                int *indi = all_indi + ray_start[ray];
                                for (int n=0; n<ray_stride[ray]; n++)
                                {
                                    update[indi[n]] += pool_upd*dist[n];
                                    nupdate[indi[n]] += dist[n];
                                }
                            }
                        }
                    }
                }
                free(pool_sim);
                free(pool_sum_dist2);
            }
            for (int n=0; n<(ngridx*ngridy); n++){
                if (nupdate[n] > 0) {
                    recon_slice[n] += update[n] / nupdate[n];
                }
            }
            free(simdata);
            free(update);
            free(nupdate);
        }
        free(ray_start);
        free(ray_stride);
        free(all_indi);
        free(all_dist);
        free(all_sum_dist2);
    }
}
