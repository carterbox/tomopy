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
    float *gridx = (float *)malloc((ngridx+1)*sizeof(float));
    float *gridy = (float *)malloc((ngridy+1)*sizeof(float));
    float *coordx = (float *)malloc((ngridy+1)*sizeof(float));
    float *coordy = (float *)malloc((ngridx+1)*sizeof(float));
    float *ax = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *ay = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *bx = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *by = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *coorx = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *coory = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *dist = (float *)malloc((ngridx+ngridy)*sizeof(float));
    int *indi = (int *)malloc((ngridx+ngridy)*sizeof(int));

    assert(coordx != NULL && coordy != NULL &&
        ax != NULL && ay != NULL && by != NULL && bx != NULL &&
        coorx != NULL && coory != NULL && dist != NULL && indi != NULL);

    int s, p, d, i, n;
    int quadrant;
    float theta_p, sin_p, cos_p;
    float mov, xi, yi;
    int asize, bsize, csize;
    float *simdata;
    float upd;
    int ind_data, ind_recon;
    float *sum_dist;
    float sum_dist2;
    float *update;

    for (i=0; i<num_iter; i++)
    {
        simdata = (float *)calloc((dt*dy*dx), sizeof(float));

        // For each slice
        for (s=0; s<dy; s++)
        {
            preprocessing(ngridx, ngridy, dx, center[s],
                &mov, gridx, gridy); // Outputs: mov, gridx, gridy

            sum_dist = (float *)calloc((ngridx*ngridy), sizeof(float));
            update = (float *)calloc((ngridx*ngridy), sizeof(float));

            // For each projection angle
            for (p=0; p<dt; p++)
            {
                // Calculate the sin and cos values
                // of the projection angle and find
                // at which quadrant on the cartesian grid.
                theta_p = fmod(theta[p], 2*M_PI);
                quadrant = calc_quadrant(theta_p);
                sin_p = sinf(theta_p);
                cos_p = cosf(theta_p);

                // For each detector pixel
                for (d=0; d<dx; d++)
                {
                    // Calculate coordinates
                    xi = -ngridx-ngridy;
                    yi = (1-dx)/2.0+d+mov;
                    calc_coords(
                        ngridx, ngridy, xi, yi, sin_p, cos_p, gridx, gridy,
                        coordx, coordy);

                    // Merge the (coordx, gridy) and (gridx, coordy)
                    trim_coords(
                        ngridx, ngridy, coordx, coordy, gridx, gridy,
                        &asize, ax, ay, &bsize, bx, by);

                    // Sort the array of intersection points (ax, ay) and
                    // (bx, by). The new sorted intersection points are
                    // stored in (coorx, coory). Total number of points
                    // are csize.
                    sort_intersections(
                        quadrant, asize, ax, ay, bsize, bx, by,
                        &csize, coorx, coory);

                    // Calculate the distances (dist) between the
                    // intersection points (coorx, coory). Find the
                    // indices of the pixels on the reconstruction grid.
                    calc_dist(
                        ngridx, ngridy, csize, coorx, coory,
                        indi, dist);

                    // Calculate simdata
                    calc_simdata(s, p, d, ngridx, ngridy, dt, dx,
                        csize, indi, dist, recon,
                        simdata); // Output: simdata


                    // Calculate dist*dist
                    sum_dist2 = 0.0;
                    for (n=0; n<csize-1; n++)
                    {
                        sum_dist2 += dist[n]*dist[n];
                        sum_dist[indi[n]] += dist[n];
                    }

                    // Update
                    if (sum_dist2 != 0.0)
                    {
                        ind_data = d+p*dx+s*dt*dx;
                        upd = (data[ind_data]-simdata[ind_data])/sum_dist2;
                        for (n=0; n<csize-1; n++)
                        {
                            update[indi[n]] += upd*dist[n];
                        }
                    }
                }
            }

            for (n = 0; n < ngridx*ngridy; n++) {
                if (sum_dist[n] != 0.0) {
                    ind_recon = s*ngridx*ngridy;
                    recon[n+ind_recon] += update[n]/sum_dist[n];
                }
            }

            free(sum_dist);
            free(update);
        }

        free(simdata);
    }

    free(gridx);
    free(gridy);
    free(coordx);
    free(coordy);
    free(ax);
    free(ay);
    free(bx);
    free(by);
    free(coorx);
    free(coory);
    free(dist);
    free(indi);
}


void
sirt_fly_rotation(
    const float *data, int dy, int dt, int dx,
    const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter, int bin, int *mask)
{
    printf("SIRT fly rotation\n");
    float *gridx = (float *)malloc((ngridx+1)*sizeof(float));
    float *gridy = (float *)malloc((ngridy+1)*sizeof(float));
    assert(gridx != NULL && gridy != NULL);
    float* simdata = (float *)malloc((dy*dt*dx)*sizeof(float));
    assert(simdata != NULL);

    float mov;

    preprocessing(ngridx, ngridy, dx, center[0],
        &mov, gridx, gridy); // Outputs: mov, gridx, gridy

    float *all_dist;
    int *all_indi, *ray_start, *ray_stride;
    compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov,
        ngridx, ngridy, &ray_start, &ray_stride, &all_indi, &all_dist);
        // Outputs: ray_start, ray_stride, all_indi, all_dist

    free(gridx);
    free(gridy);

    float *update = malloc(ngridx * ngridy * dy * sizeof *update);
    int *nupdate = malloc(ngridx * ngridy * dy * sizeof *nupdate);
    assert(update != NULL && nupdate != NULL);

    float *dist;
    int *indi;
    float *sum_dist2 = malloc(sizeof *sum_dist2 * dt * dx);
    assert(sum_dist2 != NULL);
    int ray, ind_data, ind_recon;
    int s, p, d, i, n, b;
    float pool_sim, pool_data, pool_sum_dist2, pool_upd;

    for (i=0; i<num_iter; i++)
    {
        // initialize simdata to zero
        memset(update, 0, ngridx * ngridy * dy * sizeof *update);
        memset(nupdate, 0, ngridx * ngridy * dy * sizeof *nupdate);
        memset(simdata, 0, dy*dt*dx*sizeof(float));
        memset(sum_dist2, 0, sizeof *sum_dist2 * dt * dx);
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
                    // For each slice
                    for (s=0; s<dy; s++)
                    {
                        // Calculate simdata
                        calc_simdata(s, p, d, ngridx, ngridy, dt, dx,
                            ray_stride[ray]+1, indi, dist, recon,
                            simdata); // Output: simdata
                    }
                }
            }
            if ((p+1) % bin == 0)
            {
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
                        // For each slice
                        for (s=0; s<dy; s++)
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
                }
            }
        }
        for (n=0; n<(ngridx*ngridy*dy); n++){
            if (nupdate[n] > 0) {
                recon[n] += update[n] / nupdate[n];
            }
        }
    }
    free(simdata);
    free(ray_start);
    free(ray_stride);
    free(all_indi);
    free(all_dist);
    free(update);
    free(nupdate);
    free(sum_dist2);
}


void
sirt_convolve(
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

    float *all_dist;
    int *all_indi, *ray_start, *ray_stride;
    compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov,
        ngridx, ngridy, &ray_start, &ray_stride, &all_indi, &all_dist);
        // Outputs: ray_start, ray_stride, all_indi, all_dist

    free(gridx);
    free(gridy);

    float *update = malloc(ngridx * ngridy * dy * sizeof *update);
    int *nupdate = malloc(ngridx * ngridy * dy * sizeof *nupdate);
    assert(update != NULL && nupdate != NULL);
    float* simdata = (float *)malloc((dy*dt*dx)*sizeof(float));
    assert(simdata != NULL);
    float *sum_dist2 = malloc(sizeof *sum_dist2 * dt * dx);
    assert(sum_dist2 != NULL);

    #pragma omp parallel
    {
        float *dist;
        int *indi;
        int ray, ind_data, ind_recon;
        int s, p, d, i, n, b;
        float pool_sim, pool_data, pool_sum_dist2, pool_upd;

        for (i=0; i<num_iter; i++)
        {
            #pragma omp single
            {
                printf("sirt_convolve: iteration %d\n", i);
                // initialize simdata to zero
                memset(update, 0, ngridx * ngridy * dy * sizeof *update);
                memset(nupdate, 0, ngridx * ngridy * dy * sizeof *nupdate);
                memset(simdata, 0, dy*dt*dx*sizeof(float));
                memset(sum_dist2, 0, sizeof *sum_dist2 * dt * dx);
            }
            // For each projection angle
            #pragma omp for
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
                        // For each slice
                        for (s=0; s<dy; s++)
                        {
                            // Calculate simdata
                            calc_simdata(s, p, d, ngridx, ngridy, dt, dx,
                                ray_stride[ray]+1, indi, dist, recon,
                                simdata); // Output: simdata
                        }
                    }
                }
            }
            #pragma omp single
            for (p=0; p<dt; p++)
            {
                if (p >= bin - 1)
                {
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
                            // For each slice
                            for (s=0; s<dy; s++)
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
                    }
                }
            }
            #pragma omp single
            for (n=0; n<(ngridx*ngridy*dy); n++){
                if (nupdate[n] > 0) {
                    recon[n] += update[n] / nupdate[n];
                }
            }
        }
    }
    free(all_dist);
    free(all_indi);
    free(ray_start);
    free(ray_stride);
    free(update);
    free(nupdate);
    free(simdata);
    free(sum_dist2);
}
