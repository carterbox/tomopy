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
art(
    const float *data, int dy, int dt, int dx,
    const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter)
{
    float *gridx = (float *)malloc((ngridx+1)*sizeof(float));
    float *gridy = (float *)malloc((ngridy+1)*sizeof(float));
    assert(gridx != NULL && gridy != NULL);
    float mov;
    preprocessing(ngridx, ngridy, dx, center[0],
        &mov, gridx, gridy); // Outputs: mov, gridx, gridy

    float *all_dist, *all_sum_dist2;
    int *all_indi, *ray_start, *ray_stride;
    compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov,
        ngridx, ngridy, &ray_start, &ray_stride, &all_indi, &all_dist,
        &all_sum_dist2);
        // Outputs: ray_start, ray_stride, all_indi, all_dist

    free(gridx);
    free(gridy);

    int i, s, p, d, n; // preferred loop order
    for (i=0; i<num_iter; i++)
    {
        printf("art: iteration %d\n", i);
        // initialize simdata to zero
        float *simdata = calloc((dt*dy*dx), sizeof *simdata);
        assert(simdata != NULL);
        // For each slice
        for (s=0; s<dy; s++)
        {
            // For each projection angle
            for (p=0; p<dt; p++)
            {
                // For each detector pixel
                for (d=0; d<dx; d++)
                {
                    int ray = d + dx*p;
                    float *dist = all_dist + ray_start[ray];
                    int *indi = all_indi + ray_start[ray];
                    float sum_dist2 = all_sum_dist2[ray];
                    if (sum_dist2 != 0.0)
                    {
                        // Calculate simdata
                        calc_simdata(s, p, d, ngridx, ngridy, dt, dx,
                            ray_stride[ray]+1, indi, dist, recon,
                            simdata); // Output: simdata
                        // Update
                        int ind_data = d + dx*(p + dt*s);
                        int ind_recon = s*ngridx*ngridy;
                        float upd = (data[ind_data]-simdata[ind_data])/sum_dist2;
                        for (n=0; n<ray_stride[ray]; n++)
                        {
                            recon[indi[n]+ind_recon] += upd*dist[n];
                        }
                    }
                }
            }
        }
        free(simdata);
    }
    free(ray_start);
    free(ray_stride);
    free(all_indi);
    free(all_dist);
    free(all_sum_dist2);
}


/**
 * @param data Measurements collected at each position. The size of data is dy, dt, dx
 * @param bin The number of angles to be grouped together.
 * @param mask The weights of each of the angles, i.e. the convolution kernel.
 * @param theta The angles at which measurements were collected the size is dt

 Given a series of measurements, [data], collected at angles, [theta]. Pool
 adjacent measurements together using a convolutional [mask] of size [bin].

 If the size of [data] is 5 and [bin] is 4, then the [mask] is evaluated at 2
 positions: indices 0 and 1.
 */
void
art_fly_rotation(
    const float *data, int dy, int dt, int dx,
    const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter, int bin, int *mask)
{
    float *gridx = (float *)malloc((ngridx+1)*sizeof(float));
    float *gridy = (float *)malloc((ngridy+1)*sizeof(float));
    assert(gridx != NULL && gridy != NULL);
    float* simdata;

    float mov;

    preprocessing(ngridx, ngridy, dx, center[0],
        &mov, gridx, gridy); // Outputs: mov, gridx, gridy

    float *all_dist, *all_sum_dist2;
    int *all_indi, *ray_start, *ray_stride;
    compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov,
        ngridx, ngridy, &ray_start, &ray_stride, &all_indi, &all_dist,
        &all_sum_dist2);
        // Outputs: ray_start, ray_stride, all_indi, all_dist

    free(gridx);
    free(gridy);

    float *update;
    int *nupdate;

    float *dist;
    int *indi;
    float *sum_dist2;
    int ray, ind_data, ind_recon;
    int s, p, d, i, n, b;
    float pool_sim, pool_data, pool_sum_dist2, pool_upd;

    for (i=0; i<num_iter; i++)
    {
        // initialize simdata to zero
        simdata = calloc(dy*dt*dx, sizeof(float));
        sum_dist2 = calloc(dt * dx, sizeof *sum_dist2);
        assert(simdata != NULL);
        assert(sum_dist2 != NULL);
        // For each slice
        for (s=0; s<dy; s++)
        {
            // For each projection angle
            for (p=0; p<dt; p++)
            {
                // For each detector pixel
                for (d=0; d<dx; d++)
                {
                    ray = d + dx*p;
                    dist = all_dist + ray_start[ray];
                    indi = all_indi + ray_start[ray];
                    // Calculate dist*dist
                    for (n=0; n<ray_stride[ray]; n++)
                    {
                        sum_dist2[ray] += dist[n]*dist[n];
                    }
                    if (sum_dist2[ray] != 0.0)
                    {
                        // Calculate simdata
                        calc_simdata(s, p, d, ngridx, ngridy, dt, dx,
                            ray_stride[ray]+1, indi, dist, recon,
                            simdata); // Output: simdata
                    }
                }
                if ((p+1) % bin == 0)
                {
                    update = calloc(ngridx * ngridy * dy, sizeof *update);
                    nupdate = calloc(ngridx * ngridy * dy, sizeof *nupdate);
                    assert(update != NULL && nupdate != NULL);
                    // For each detector pixel
                    for (d=0; d<dx; d++)
                    {
                        // Simulate pooled data
                        pool_sim = 0; pool_data = 0; pool_sum_dist2 = 0;
                        for (b=0; b<bin; b++)
                        {
                            if (mask[b] > 0) {
                                ray = d + dx*(p-b);
                                pool_sum_dist2 += sum_dist2[ray];
                            }
                        }
                        if (pool_sum_dist2 > 0)
                        {
                            for (b=0; b<bin; b++)
                            {
                                if (mask[b] > 0)
                                {
                                    int p1 = p-b;
                                    ind_data = d+p1*dx+s*dt*dx;
                                    pool_sim += simdata[ind_data];
                                    pool_data += data[ind_data];
                                }
                            }
                            // Compute update
                            pool_upd = (pool_data - pool_sim) / pool_sum_dist2;

                            // Update
                            for (b=0; b<bin; b++)
                            {
                                if (mask[b] > 0) {
                                    ray = d + dx*(p-b);
                                    dist = all_dist + ray_start[ray];
                                    indi = all_indi + ray_start[ray];
                                    ind_recon = s*ngridx*ngridy;
                                    for (n=0; n<ray_stride[ray]; n++)
                                    {
                                        update[indi[n]+ind_recon] += pool_upd*dist[n];
                                        nupdate[indi[n]+ind_recon] += 1;
                                    }
                                }
                            }
                        }
                    }
                    for (n=0; n<(ngridx*ngridy*dy); n++){
                        if (nupdate[n] > 0) {
                            recon[n] += update[n] / nupdate[n];
                        }
                    }
                    free(update);
                    free(nupdate);
                }
            }
        }
        free(simdata);
        free(sum_dist2);
    }
    free(ray_start);
    free(ray_stride);
    free(all_indi);
    free(all_dist);
    free(all_sum_dist2);
}

void
art_convolve(
    const float *data, int dy, int dt, int dx,
    const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter, int bin, int *mask)
{
    float *gridx = (float *)malloc((ngridx+1)*sizeof(float));
    float *gridy = (float *)malloc((ngridy+1)*sizeof(float));
    assert(gridx != NULL && gridy != NULL);
    float mov;
    preprocessing(ngridx, ngridy, dx, center[0],
        &mov, gridx, gridy); // Outputs: mov, gridx, gridy

    float *all_dist, *all_sum_dist2;
    int *all_indi, *ray_start, *ray_stride;
    compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov,
        ngridx, ngridy, &ray_start, &ray_stride, &all_indi, &all_dist,
        &all_sum_dist2);
        // Outputs: ray_start, ray_stride, all_indi, all_dist

    free(gridx);
    free(gridy);

    int i, s, p, b, d, n; // preferred loop order
    for (i=0; i<num_iter; i++)
    {
        // For each slice
        for (s=0; s<dy; s++)
        {
            // For each projection angle
            for (p=bin-1; p<dt; p++)
            {
                // initialize simdata to zero
                float *simdata = calloc(dy*dt*dx, sizeof *simdata);
                assert(simdata != NULL);
                float *update = calloc(ngridx * ngridy * dy, sizeof *update);
                int *nupdate = calloc(ngridx * ngridy * dy, sizeof *nupdate);
                assert(update != NULL && nupdate != NULL);
                // For each code element
                for (b=0; b<bin; b++)
                {
                    if (mask[b] > 0)
                    {
                        // For each detector pixel
                        for (d=0; d<dx; d++)
                        {
                            int ray = d + dx*(p-b);
                            float *dist = all_dist + ray_start[ray];
                            int *indi = all_indi + ray_start[ray];
                            if (all_sum_dist2[ray] != 0.0)
                            {
                                // Calculate simdata
                                calc_simdata(s, (p-b), d, ngridx, ngridy, dt, dx,
                                    ray_stride[ray]+1, indi, dist, recon,
                                    simdata); // Output: simdata
                            }
                        }
                    }
                }
                // For each detector pixel
                for (d=0; d<dx; d++)
                {
                    // Simulate pooled data
                    float pool_sim = 0;
                    float pool_data = 0;
                    float pool_sum_dist2 = 0;
                    // For each code element
                    for (b=0; b<bin; b++)
                    {
                        if (mask[b] > 0) {
                            int ray = d + dx*(p-b);
                            pool_sum_dist2 += all_sum_dist2[ray];
                        }
                    }
                    if (pool_sum_dist2 > 0)
                    {
                        for (b=0; b<bin; b++)
                        {
                            if (mask[b] > 0)
                            {
                                int p1 = p-b;
                                int ind_data = d + dx*(p1 + dt*s);
                                pool_sim += simdata[ind_data];
                                pool_data += data[ind_data];
                            }
                        }
                        // Compute update
                        float pool_upd = (pool_data - pool_sim) / pool_sum_dist2;
                        // Update
                        for (b=0; b<bin; b++)
                        {
                            if (mask[b] > 0) {
                                int ray = d + dx*(p-b);
                                float* dist = all_dist + ray_start[ray];
                                int *indi = all_indi + ray_start[ray];
                                int ind_recon = s*ngridx*ngridy;
                                for (n=0; n<ray_stride[ray]; n++)
                                {
                                    update[indi[n]+ind_recon] += pool_upd*dist[n];
                                    nupdate[indi[n]+ind_recon] += 1;
                                }
                            }
                        }
                    }
                }
                for (n=0; n<(ngridx*ngridy*dy); n++){
                    if (nupdate[n] > 0) {
                        recon[n] += update[n] / nupdate[n];
                    }
                }
                free(simdata);
                free(update);
                free(nupdate);
            }
        }
    }
    free(ray_start);
    free(ray_stride);
    free(all_indi);
    free(all_dist);
    free(all_sum_dist2);
}
