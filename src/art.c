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


void compute_indices_and_lengths(
const float * const theta, const int dt, const int dx,
const float *gridx, const float *gridy, const float mov,
const int ngridx, const int ngridy,
int ** const ray_start, int ** const ray_stride,
int ** const indices, float ** const distances)
{
    int ** const indi_list = malloc(sizeof *indi_list * dx*dt);
    float ** const dist_list = malloc(sizeof *dist_list * dx*dt);
    *ray_start = malloc(sizeof **ray_start * dx*dt);
    *ray_stride = malloc(sizeof **ray_stride * dx*dt);
    assert(indi_list != NULL && dist_list != NULL &&
        ray_start != NULL && ray_stride != NULL);

    #pragma omp parallel
    {
        int quadrant;
        float theta_p, sin_p, cos_p;
        float xi, yi;
        int asize, bsize, csize;

        float *coordx = (float *)malloc((ngridy+1)*sizeof(float));
        float *coordy = (float *)malloc((ngridx+1)*sizeof(float));
        float *ax = (float *)malloc((ngridx+ngridy)*sizeof(float));
        float *ay = (float *)malloc((ngridx+ngridy)*sizeof(float));
        float *bx = (float *)malloc((ngridx+ngridy)*sizeof(float));
        float *by = (float *)malloc((ngridx+ngridy)*sizeof(float));
        float *coorx = (float *)malloc((ngridx+ngridy)*sizeof(float));
        float *coory = (float *)malloc((ngridx+ngridy)*sizeof(float));
        assert(coordx != NULL && coordy != NULL &&
            ax != NULL && ay != NULL && by != NULL && bx != NULL &&
            coorx != NULL && coory != NULL);

        // For each projection angle
        #pragma omp for
        for (int p=0; p<dt; p++)
        {
            // For each detector pixel
            for (int d=0; d<dx; d++)
            {
                float *dist = (float *)malloc((ngridx+ngridy)*sizeof(float));
                int *indi = (int *)malloc((ngridx+ngridy)*sizeof(int));
                assert(dist != NULL && indi != NULL);
                // Calculate the sin and cos values
                // of the projection angle and find
                // at which quadrant on the cartesian grid.
                theta_p = fmod(theta[p], 2*M_PI);
                quadrant = calc_quadrant(theta_p);
                sin_p = sinf(theta_p);
                cos_p = cosf(theta_p);
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
                // Save the intersections and lengths from this ray
                int ray = d + p*dx;
                indi_list[ray] = indi;
                dist_list[ray] = dist;
                (*ray_stride)[ray] = csize - 1;
            }
        }
        free(coordx);
        free(coordy);
        free(ax);
        free(ay);
        free(bx);
        free(by);
        free(coorx);
        free(coory);
        // Compute the total length of the combined arrays and the start of each
        // block
        #pragma omp single
        {
            int sum_ray_stride = 0;
            for (int ray=0; ray<dx*dt; ray++)
            {
                (*ray_start)[ray] = sum_ray_stride;
                sum_ray_stride += (*ray_stride)[ray];
            }
            // Copy all of the intersections and distances into one array each
            *indices = malloc(sizeof **indices * sum_ray_stride);
            *distances = malloc(sizeof **distances * sum_ray_stride);
        }
        #pragma omp barrier

        #pragma omp for nowait
        for (int ray=0; ray<dx*dt; ray++)
        {
            int j = (*ray_start)[ray];
            memcpy(&(*indices)[j], indi_list[ray],
                sizeof **indices * (*ray_stride)[ray]);
            free(indi_list[ray]);
            memcpy(&(*distances)[j], dist_list[ray],
                sizeof **distances * (*ray_stride)[ray]);
            free(dist_list[ray]);
        }
    }
    free(indi_list);
    free(dist_list);
}

void
art(
    const float *data, int dy, int dt, int dx,
    const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter)
{
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

    float *dist;
    int *indi;
    float sum_dist2, upd;
    int ray, ind_data, ind_recon;
    int s, p, d, i, n;

    for (i=0; i<num_iter; i++)
    {
        // initialize simdata to zero
        memset(simdata, 0, dy*dt*dx*sizeof(float));
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
                sum_dist2 = 0.0;
                for (n=0; n<ray_stride[ray]; n++)
                {
                    sum_dist2 += dist[n]*dist[n];
                }
                if (sum_dist2 != 0.0)
                {
                    // For each slice
                    for (s=0; s<dy; s++)
                    {
                        // Calculate simdata
                        calc_simdata(s, p, d, ngridx, ngridy, dt, dx,
                            ray_stride[ray]+1, indi, dist, recon,
                            simdata); // Output: simdata
                        // Update
                        ind_data = d+p*dx+s*dt*dx;
                        ind_recon = s*ngridx*ngridy;
                        upd = (data[ind_data]-simdata[ind_data])/sum_dist2;
                        for (n=0; n<ray_stride[ray]; n++)
                        {
                        	recon[indi[n]+ind_recon] += upd*dist[n];
                        }
                    }
                }
            }
        }
    }
    free(simdata);
    free(ray_start);
    free(ray_stride);
    free(all_indi);
    free(all_dist);
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
                memset(update, 0, ngridx * ngridy * dy * sizeof *update);
                memset(nupdate, 0, ngridx * ngridy * dy * sizeof *nupdate);
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
art_convolve(
    const float *data, int dy, int dt, int dx,
    const float *center, const float *theta,
    float *recon, int ngridx, int ngridy, int num_iter, int bin, int *mask)
{
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
        // For each projection angle
        for (p=bin-1; p<dt; p++)
        {
            // initialize simdata to zero
            memset(simdata, 0, dy*dt*dx*sizeof(float));
            memset(sum_dist2, 0, sizeof *sum_dist2 * dt * dx);
            memset(update, 0, ngridx * ngridy * dy * sizeof *update);
            memset(nupdate, 0, ngridx * ngridy * dy * sizeof *nupdate);
            for (b=0; b<bin; b++)
            {
                // For each detector pixel
                for (d=0; d<dx; d++)
                {
                    ray = d + dx*(p-b);
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
