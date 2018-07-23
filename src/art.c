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
    float* simdata = (float *)malloc((dy*dt*dx)*sizeof(float));

    assert(coordx != NULL && coordy != NULL &&
        ax != NULL && ay != NULL && by != NULL && bx != NULL &&
        coorx != NULL && coory != NULL && dist != NULL &&
        indi != NULL && simdata != NULL);

    int s, p, d, i, n;
    int quadrant;
    float theta_p, sin_p, cos_p;
    float mov, xi, yi;
    int asize, bsize, csize;
    float upd;
    int ind_data, ind_recon;

    for (i=0; i<num_iter; i++)
    {
        // initialize simdata to zero
        memset(simdata, 0, dy*dt*dx*sizeof(float));

        preprocessing(ngridx, ngridy, dx, center[0],
            &mov, gridx, gridy); // Outputs: mov, gridx, gridy

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


                // Calculate dist*dist
                float sum_dist2 = 0.0;
                for (n=0; n<csize-1; n++)
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
                            csize, indi, dist, recon,
                            simdata); // Output: simdata

                        // Update
                        ind_data = d+p*dx+s*dt*dx;
                        ind_recon = s*ngridx*ngridy;
                        upd = (data[ind_data]-simdata[ind_data])/sum_dist2;
                        for (n=0; n<csize-1; n++)
                        {
                        	recon[indi[n]+ind_recon] += upd*dist[n];
                        }
                    }
                }
            }
        }
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
    free(simdata);
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
    float *coordx = (float *)malloc((ngridy+1)*sizeof(float));
    float *coordy = (float *)malloc((ngridx+1)*sizeof(float));
    float *ax = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *ay = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *bx = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *by = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *coorx = (float *)malloc((ngridx+ngridy)*sizeof(float));
    float *coory = (float *)malloc((ngridx+ngridy)*sizeof(float));
    assert(gridx != NULL && gridy != NULL && coordx != NULL &&
           coordy != NULL && ax != NULL && ay != NULL && by != NULL &&
           bx != NULL && coorx != NULL && coory != NULL);

    float *simdata = malloc(dy * dt * dx * sizeof *simdata);
    assert(simdata != NULL);


    int pool_buffer_size = bin * dx;
    float* sum_dist2 = malloc(pool_buffer_size * sizeof *sum_dist2);
    float *dist = malloc(pool_buffer_size * (ngridx + ngridy) * sizeof *dist);
    float *dist_b;
    int *indi   = malloc(pool_buffer_size * (ngridx + ngridy) * sizeof *indi);
    int *indi_b;
    int *csize  = malloc(pool_buffer_size * sizeof *csize);
    assert(sum_dist2 != NULL);
    assert(dist != NULL && indi != NULL && csize != NULL);

    float *update = malloc(ngridx * ngridy * sizeof *update);
    int *nupdate = malloc(ngridx * ngridy * sizeof *nupdate);
    assert(update != NULL && nupdate != NULL);

    int s, p, d, i, n, b;
    int quadrant;
    float theta_p, sin_p, cos_p;
    float mov, xi, yi;
    int asize, bsize;
    int ind_data, ind_recon;

    for (i=0; i<num_iter; i++)
    {
        printf("num_iter=%i\n", i);

        preprocessing(ngridx, ngridy, dx, center[0],
            &mov, gridx, gridy); // Outputs: mov, gridx, gridy

        memset(simdata, 0, dy * dt * dx * sizeof *simdata);

        // For each slice
        for (s=0; s<dy; s++)
        {
            // For each projection angle
            for (p=0; p<dt-bin+1; p+=bin)
            {
                memset(update, 0, ngridx * ngridy * sizeof *update);
                memset(nupdate, 0, ngridx * ngridy * sizeof *nupdate);
                memset(sum_dist2, 0, pool_buffer_size * sizeof *sum_dist2);
                // For each detector pixel
                for (d=0; d<dx; d++)
                {
                    // For binned of projection angle
                    for (b=0; b<bin; b++)
                    {
                        // Choose where to store indices and lengths for this
                        // projection.
                        int ind_buffer = b + (d * bin);
                        indi_b = indi + ind_buffer * (ngridx + ngridy);
                        dist_b = dist + ind_buffer * (ngridx + ngridy);

                        // calculate lengths and intersections ---->
                        // Calculate the sin and cos values
                        // of the projection angle and find
                        // at which quadrant on the cartesian grid.
                        theta_p = fmod(theta[p+b], 2*M_PI);
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
                        &csize[ind_buffer], coorx, coory);

                        // Calculate the distances (dist) between the
                        // intersection points (coorx, coory). Find the
                        // indices of the pixels on the reconstruction grid.
                        calc_dist(
                        ngridx, ngridy, csize[ind_buffer], coorx, coory,
                        indi_b, dist_b);
                        // <---- calculate lengths and intersections

                        // Calculate the dot product of the intersection lengths
                        for (n=0; n<csize[ind_buffer]-1; n++)
                        {
                            sum_dist2[ind_buffer] += dist_b[n] * dist_b[n];
                        }

                        if (sum_dist2[ind_buffer] != 0.0)
                        {
                            calc_simdata(
                                s, p+b, d, ngridx, ngridy, dt, dx,
                                csize[ind_buffer], indi_b, dist_b, recon,
                                simdata); // Output: simdata
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
                    float pool_upd;
                    for (b=0; b<bin; b++)
                    {
                        int ind_buffer = b + (d * bin);
                        if (mask[b] > 0) {
                            pool_sum_dist2 += sum_dist2[ind_buffer];
                        }
                    }
                    if (pool_sum_dist2 > 0) {
                        for (b=0; b<bin; b++)
                        {
                            int p1 = p+b;
                            if (mask[b] > 0) {
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
                            int ind_buffer = b + (d * bin);
                            indi_b = indi + ind_buffer * (ngridx + ngridy);
                            dist_b = dist + ind_buffer * (ngridx + ngridy);
                            for (n=0; n<csize[ind_buffer]-1; n++)
                            {
                                float upd = pool_upd*dist_b[n];
                                // printf("update %d -> %f\n", n, upd);
                                update[indi_b[n]] += upd;
                                nupdate[indi_b[n]] += 1;
                            }
                        }
                    }
                }
                ind_recon = s*ngridx*ngridy;
                for (n=0; n<(ngridx*ngridy); n++){
                    if (nupdate[n] > 0) {
                        recon[n+ind_recon] += update[n] / nupdate[n];
                    }
                }
            }
        }
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
    free(simdata);
    free(sum_dist2);
    free(indi);
    free(dist);
    free(csize);
    // free(update);
}
