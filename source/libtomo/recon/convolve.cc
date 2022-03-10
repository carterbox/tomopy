
#include "recon.h"
#include "utils.h"

#include <stdint.h>

void
compute_indices_and_lengths(const float* const theta, const int dt, const int dx,
                            const float* gridx, const float* gridy, const float mov,
                            const int ngridx, const int ngridy, int** const ray_start,
                            int** const ray_stride, int** const indices,
                            float** const distances, float** const sum_distances2)
{
    int** const   indi_list = malloc(sizeof *indi_list * dx * dt);
    float** const dist_list = malloc(sizeof *dist_list * dx * dt);
    *ray_start              = malloc(sizeof **ray_start * dx * dt);
    *ray_stride             = malloc(sizeof **ray_stride * dx * dt);
    assert(indi_list != NULL && dist_list != NULL && *ray_start != NULL &&
           *ray_stride != NULL);
    *sum_distances2 = calloc(dx * dt, sizeof **sum_distances2);
    assert(*sum_distances2 != NULL);

#pragma omp parallel
    {
        int   quadrant;
        float theta_p, sin_p, cos_p;
        float xi, yi;
        int   asize, bsize, csize;

        float* coordx = (float*) malloc((ngridy + 1) * sizeof(float));
        float* coordy = (float*) malloc((ngridx + 1) * sizeof(float));
        float* ax     = (float*) malloc((ngridx + ngridy) * sizeof(float));
        float* ay     = (float*) malloc((ngridx + ngridy) * sizeof(float));
        float* bx     = (float*) malloc((ngridx + ngridy) * sizeof(float));
        float* by     = (float*) malloc((ngridx + ngridy) * sizeof(float));
        float* coorx  = (float*) malloc((ngridx + ngridy) * sizeof(float));
        float* coory  = (float*) malloc((ngridx + ngridy) * sizeof(float));
        assert(coordx != NULL && coordy != NULL && ax != NULL && ay != NULL &&
               by != NULL && bx != NULL && coorx != NULL && coory != NULL);

// For each projection angle
#pragma omp for
        for(int p = 0; p < dt; p++)
        {
            // For each detector pixel
            for(int d = 0; d < dx; d++)
            {
                float* dist = (float*) malloc((ngridx + ngridy) * sizeof(float));
                int*   indi = (int*) malloc((ngridx + ngridy) * sizeof(int));
                assert(dist != NULL && indi != NULL);
                // Calculate the sin and cos values
                // of the projection angle and find
                // at which quadrant on the cartesian grid.
                theta_p  = fmod(theta[p], 2 * M_PI);
                quadrant = calc_quadrant(theta_p);
                sin_p    = sinf(theta_p);
                cos_p    = cosf(theta_p);
                // Calculate coordinates
                xi = -ngridx - ngridy;
                yi = (1 - dx) * 0.5f + d + mov;
                calc_coords(ngridx, ngridy, xi, yi, sin_p, cos_p, gridx, gridy, coordx,
                            coordy);
                // Merge the (coordx, gridy) and (gridx, coordy)
                trim_coords(ngridx, ngridy, coordx, coordy, gridx, gridy, &asize, ax, ay,
                            &bsize, bx, by);
                // Sort the array of intersection points (ax, ay) and
                // (bx, by). The new sorted intersection points are
                // stored in (coorx, coory). Total number of points
                // are csize.
                sort_intersections(quadrant, asize, ax, ay, bsize, bx, by, &csize, coorx,
                                   coory);
                // Calculate the distances (dist) between the
                // intersection points (coorx, coory). Find the
                // indices of the pixels on the reconstruction grid.
                calc_dist(ngridx, ngridy, csize, coorx, coory, indi, dist);
                // Save the intersections and lengths from this ray
                int ray            = d + p * dx;
                indi_list[ray]     = indi;
                dist_list[ray]     = dist;
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
            for(int ray = 0; ray < dx * dt; ray++)
            {
                (*ray_start)[ray] = sum_ray_stride;
                sum_ray_stride += (*ray_stride)[ray];
            }
            *indices   = malloc(sizeof **indices * sum_ray_stride);
            *distances = malloc(sizeof **distances * sum_ray_stride);
        }
#pragma omp barrier

#pragma omp for nowait
        for(int ray = 0; ray < dx * dt; ray++)
        {
            // Compute the squared sum of the distances for each ray
            for(int n = 0; n < (*ray_stride)[ray]; n++)
            {
                (*sum_distances2)[ray] += dist_list[ray][n] * dist_list[ray][n];
            }
            // Copy all of the intersections and distances into one array each
            int j = (*ray_start)[ray];
            memcpy(&(*indices)[j], indi_list[ray], sizeof **indices * (*ray_stride)[ray]);
            free(indi_list[ray]);
            memcpy(&(*distances)[j], dist_list[ray],
                   sizeof **distances * (*ray_stride)[ray]);
            free(dist_list[ray]);
        }
    }
    free(indi_list);
    free(dist_list);
}

/**
Given a series of convolved measurements, [data], collected at angles,
[theta]. Pool adjacent angles together using a convolutional [mask] of
size [nmask].

Example
-------
[raw_data] = [3.0, 1.0, 2.0, 0.0, 2.0, 0.0, 0.0]
[mask]     = [true, true]
[nmask]    = 2
[data]     = [4.0, 3.0, 2.0, 2.0, 2.0, 0.0, 3.0]

Note about the convolution mask:
    1. It is left aligned.
    2. It wraps around from the left to the right edge.

Then to reconstruct, we compare [data] with mask * simdata([theta]).

@param data Measurements collected at each position. The size of data is
    dy, dt, dx. The data is preconvolved.
@param nmask The number of angles to be grouped together.
@param mask A intean mask for the angles, i.e. the convolution kernel.
@param theta The angles at which measurements were collected the size is dt
 */
void
art_convolve(const float* data, int dy, int dt, int dx, const float* center,
             const float* theta, float* recon, int ngridx, int ngridy, int num_iter,
             int nmask, int* mask, int* ind_block)
{
    // int s, i, p, b, d, n; // preferred loop order
    // For each slice
    for(int s = 0; s < dy; s++)
    {
        float* recon_slice = recon + s * ngridx * ngridy;
        float* gridx       = (float*) malloc((ngridx + 1) * sizeof(float));
        float* gridy       = (float*) malloc((ngridy + 1) * sizeof(float));
        assert(gridx != NULL && gridy != NULL);
        float mov;
        preprocessing(ngridx, ngridy, dx, center[s], &mov, gridx, gridy);
        // Outputs: mov, gridx, gridy
        float *all_dist, *all_sum_dist2;
        int *  all_indi, *ray_start, *ray_stride;
        compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov, ngridx, ngridy,
                                    &ray_start, &ray_stride, &all_indi, &all_dist,
                                    &all_sum_dist2);
        // Outputs: ray_start, ray_stride, all_indi, all_dist
        free(gridx);
        free(gridy);
        // For each iteration
        for(int i = 0; i < num_iter; i++)
        {
            // For each projection angle
            for(int k = 0; k < dt; k++)
            {
                int p = ind_block[k];
                // Initialize buffers to zero
                float* simdata = calloc(dt * dx, sizeof *simdata);
                assert(simdata != NULL);
                float* update  = calloc(ngridx * ngridy, sizeof *update);
                float* nupdate = calloc(ngridx * ngridy, sizeof *nupdate);
                assert(update != NULL && nupdate != NULL);
                float* pool_sim       = calloc(dx, sizeof *pool_sim);
                float* pool_sum_dist2 = calloc(dx, sizeof *pool_sum_dist2);
                assert(pool_sim != NULL && pool_sum_dist2 != NULL);
                // For each code element
                for(int b = 0; b < nmask; b++)
                {
                    if(mask[b])
                    {
                        int p1 = (p + b) % dt;
                        // For each detector pixel
                        for(int d = 0; d < dx; d++)
                        {
                            int    ray  = d + dx * (p1);
                            float* dist = all_dist + ray_start[ray];
                            int*   indi = all_indi + ray_start[ray];
                            if(all_sum_dist2[ray] != 0.0)
                            {
                                // Calculate simdata
                                calc_simdata(0, p1, d, ngridx, ngridy, dt, dx,
                                             ray_stride[ray] + 1, indi, dist, recon_slice,
                                             simdata);  // Output: simdata
                                // Calculate pool data
                                pool_sum_dist2[d] += all_sum_dist2[ray];
                                int ind_sim = d + dx * p1;
                                pool_sim[d] += simdata[ind_sim];
                            }
                        }
                    }
                }
                // For each detector pixel
                for(int d = 0; d < dx; d++)
                {
                    if(pool_sum_dist2[d] > 0)
                    {
                        int   ind_data = d + dx * (p + dt * s);
                        float pool_upd =
                            (data[ind_data] - pool_sim[d]) / pool_sum_dist2[d];
                        // For each code element
                        for(int b = 0; b < nmask; b++)
                        {
                            if(mask[b])
                            {
                                int    ray  = d + dx * ((p + b) % dt);
                                float* dist = all_dist + ray_start[ray];
                                int*   indi = all_indi + ray_start[ray];
                                for(int n = 0; n < ray_stride[ray]; n++)
                                {
                                    update[indi[n]] += pool_upd * dist[n];
                                    nupdate[indi[n]] += dist[n];
                                }
                            }
                        }
                    }
                }
                for(int n = 0; n < (ngridx * ngridy); n++)
                {
                    if(nupdate[n] > 0)
                    {
                        recon_slice[n] += update[n] / nupdate[n];
                    }
                }
                free(simdata);
                free(update);
                free(nupdate);
                free(pool_sim);
                free(pool_sum_dist2);
            }
        }
        free(ray_start);
        free(ray_stride);
        free(all_indi);
        free(all_dist);
        free(all_sum_dist2);
    }
}

void
art_fly_rotation(const float* data, int dy, int dt, int dx, const float* center,
                 const float* theta, float* recon, int ngridx, int ngridy, int num_iter,
                 int nmask, int* mask)
{
    return;
    int*  something;
    int* bmask;
    art_convolve(data, dy, dt, dx, center, theta, recon, ngridx, ngridy, num_iter, nmask,
                 bmask, something);
}

void
sirt_fly_rotation(const float* data, int dy, int dt, int dx, const float* center,
                  const float* theta, float* recon, int ngridx, int ngridy, int num_iter,
                  int nmask, int* mask)
{
    sirt_convolve(data, dy, dt, dx, center, theta, recon, ngridx, ngridy, num_iter, nmask,
                  mask);
}

void
sirt_convolve(const float* data, int dy, int dt, int dx, const float* center,
              const float* theta, float* recon, int ngridx, int ngridy, int num_iter,
              int nmask, int* mask)
{
    int step = 1;
    assert(step > 0 && "Step must be positive or else infinite loop.");
    // int i, s, p, b, d, n; // preferred loop order
    // For each slice
    for(int s = 0; s < dy; s++)
    {
        float* recon_slice = recon + s * ngridx * ngridy;
        float* gridx       = (float*) malloc((ngridx + 1) * sizeof(float));
        float* gridy       = (float*) malloc((ngridy + 1) * sizeof(float));
        assert(gridx != NULL && gridy != NULL);
        float mov;
        preprocessing(ngridx, ngridy, dx, center[s], &mov, gridx, gridy);
        // Outputs: mov, gridx, gridy
        float *all_dist, *all_sum_dist2;
        int *  all_indi, *ray_start, *ray_stride;
        compute_indices_and_lengths(theta, dt, dx, gridx, gridy, mov, ngridx, ngridy,
                                    &ray_start, &ray_stride, &all_indi, &all_dist,
                                    &all_sum_dist2);
        // Outputs: ray_start, ray_stride, all_indi, all_dist
        free(gridx);
        free(gridy);
        // For each iteration
        for(int i = 0; i < num_iter; i++)
        {
            float* simdata = calloc(dt * dx, sizeof *simdata);
            assert(simdata != NULL);
            float* update  = calloc(ngridx * ngridy, sizeof *update);
            float* nupdate = calloc(ngridx * ngridy, sizeof *nupdate);
            assert(update != NULL && nupdate != NULL);
            // For each projection angle, simulate data acquisition
            for(int p = 0; p < dt; p++)
            {
                // For each detector pixel
                for(int d = 0; d < dx; d++)
                {
                    int    ray  = d + dx * p;
                    float* dist = all_dist + ray_start[ray];
                    int*   indi = all_indi + ray_start[ray];
                    if(all_sum_dist2[ray] != 0.0)
                    {
                        // Calculate simdata
                        calc_simdata(0, p, d, ngridx, ngridy, dt, dx, ray_stride[ray] + 1,
                                     indi, dist, recon_slice,
                                     simdata);  // Output: simdata
                    }
                }
            }
            // For each projection angle, pool data and compute updates
            for(int p = 0; p < dt; p += step)
            {
                // Initialize buffers to zero
                float* pool_sim       = calloc(dx, sizeof *pool_sim);
                float* pool_sum_dist2 = calloc(dx, sizeof *pool_sum_dist2);
                assert(pool_sim != NULL && pool_sum_dist2 != NULL);
                // For each code element
                for(int b = 0; b < nmask; b++)
                {
                    if(mask[b])
                    {
                        int p1 = (p + b) % dt;
                        // For each detector pixel
                        for(int d = 0; d < dx; d++)
                        {
                            int ray = d + dx * (p1);
                            if(all_sum_dist2[ray] != 0.0)
                            {
                                // Calculate pool data
                                pool_sum_dist2[d] += all_sum_dist2[ray];
                                int ind_sim = d + dx * p1;
                                pool_sim[d] += simdata[ind_sim];
                            }
                        }
                    }
                }
                // For each code element
                for(int b = 0; b < nmask; b++)
                {
                    if(mask[b])
                    {
                        int p1 = (p + b) % dt;
                        // For each detector pixel
                        for(int d = 0; d < dx; d++)
                        {
                            if(pool_sum_dist2[d] > 0)
                            {
                                // Compute update
                                int   ind_data = d + dx * (p1 + dt * s);
                                float pool_upd =
                                    (data[ind_data] - pool_sim[d]) / pool_sum_dist2[d];
                                // Update
                                int    ray  = d + dx * (p1);
                                float* dist = all_dist + ray_start[ray];
                                int*   indi = all_indi + ray_start[ray];
                                for(int n = 0; n < ray_stride[ray]; n++)
                                {
                                    update[indi[n]] += pool_upd * dist[n];
                                    nupdate[indi[n]] += dist[n];
                                }
                            }
                        }
                    }
                }
                free(pool_sim);
                free(pool_sum_dist2);
            }
            for(int n = 0; n < (ngridx * ngridy); n++)
            {
                if(nupdate[n] > 0)
                {
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
