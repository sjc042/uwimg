#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for (i = 0; i < n; ++i)
    {
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i % im.w;
    d.p.y = i / im.w;
    d.data = calloc(w * w * im.c, sizeof(float));
    d.n = w * w * im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for (c = 0; c < im.c; ++c)
    {
        float cval = im.data[c * im.w * im.h + i];
        for (dx = -w / 2; dx < (w + 1) / 2; ++dx)
        {
            for (dy = -w / 2; dy < (w + 1) / 2; ++dy)
            {
                float val = get_pixel(im, c, i / im.w + dy, i % im.w + dx);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for (i = -9; i < 10; ++i)
    {
        set_pixel(im, 0, y, x + i, 1);
        set_pixel(im, 0, y + i, x, 1);
        set_pixel(im, 1, y, x + i, 0);
        set_pixel(im, 1, y + i, x, 0);
        set_pixel(im, 2, y, x + i, 1);
        set_pixel(im, 2, y + i, x, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for (i = 0; i < n; ++i)
    {
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    int w = ((int)(6 * sigma)) | 1; // Comment: i.e. 1|1 = 1, 2|1 = 3, 3|1 = 3, 4|1 = 5
    image f = make_image(1, 1, w);
    // Method 1
    float G;
    for (int i = w / -2; i < w / 2 + 1; ++i)
    {
        G = 1 / (sqrt(TWOPI) * sigma) * exp(-1 * (i * i) / (2 * sigma * sigma));
        set_pixel(f, 0, 0, i + w / 2, G);
    }
    // // Method 1
    // for(int i = 0; i < w; ++i){
    //     float x = w / 2. - i - .5;
    //     float val = 1./(TWOPI*sigma*sigma) * exp(-(x*x)/(2*sigma*sigma));
    //     set_pixel(f, 0, 0, i, val);
    // }

    l1_normalize(f);
    return f;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{

    if (0)
    {
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    }
    else
    {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.

        // Method 1
        // image gauss_row = make_1d_gaussian(sigma);
        // image gauss_col = make_image(1, gauss_row.w, 1);
        // float pixVal;
        // for (int i = 0; i < gauss_row.w; ++i)    // Comment: transpose gauss_row
        // {
        //     gauss_col.data[i] = gauss_row.data[i];
        // }
        // im = convolve_image(im, gauss_row, 1);
        // im = convolve_image(im, gauss_col, 1);
        // return copy_image(im);

        // Method 2
        image g = make_1d_gaussian(sigma);
        image s1 = convolve_image(im, g, 1);
        g.h = g.w;
        g.w = 1;
        image s2 = convolve_image(s1, g, 1);
        free_image(g);
        free_image(s1);
        return s2;
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(3, im.h, im.w);
    // TODO: calculate structure matrix for im.
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image Ix = convolve_image(im, gx_filter, 0); // !!! May need to change to non-preserve, 1-->0
    image Iy = convolve_image(im, gy_filter, 0); // !!! May need to change to non-preserve, 1-->0
    image Ix2 = make_image(Ix.c, Ix.h, Ix.w);
    image Iy2 = make_image(Iy.c, Iy.h, Iy.w);
    image IxIy = make_image(Ix.c, Ix.h, Ix.w);
    for (int i = 0; i < Ix.c * Ix.h * Ix.w; ++i)
    {
        Ix2.data[i] = Ix.data[i] * Ix.data[i];
        Iy2.data[i] = Iy.data[i] * Iy.data[i];
        IxIy.data[i] = Ix.data[i] * Iy.data[i];
    }
    Ix2 = smooth_image(Ix2, sigma);
    Iy2 = smooth_image(Iy2, sigma);
    IxIy = smooth_image(IxIy, sigma);
    for (int h = 0; h < S.h; ++h)
    {
        for (int w = 0; w < S.w; ++w)
        {
            float Ix2_val = get_pixel(Ix2, 0, h, w);
            float Iy2_val = get_pixel(Iy2, 0, h, w);
            float IxIy_val = get_pixel(IxIy, 0, h, w);
            set_pixel(S, 0, h, w, Ix2_val);
            set_pixel(S, 1, h, w, Iy2_val);
            set_pixel(S, 2, h, w, IxIy_val);
        }
    }
    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(1, S.h, S.w);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    // det(S) = Ix_2(channel 1 of S) * Iy_2(channel 2 of S) - IxIy*IxIy(channel 3 of S)
    // trace(S) = Ix_2 + Iy_2
    // image det_S = make_image(1, S.h, S.w);
    // image trace_S = make_image(1, S.h, S.w);
    float alpha = 0.06;
    float det_S;
    float trace_S;
    float Ix2_val;
    float Iy2_val;
    float IxIy_val;
    float r;
    for (int h = 0; h < S.h; ++h)
    {
        for (int w = 0; w < S.w; ++w)
        {
            Ix2_val = get_pixel(S, 0, h, w);
            Iy2_val = get_pixel(S, 1, h, w);
            IxIy_val = get_pixel(S, 2, h, w);
            det_S = Ix2_val * Iy2_val - IxIy_val * IxIy_val;
            trace_S = Ix2_val + Iy2_val;
            r = det_S - alpha * trace_S * trace_S;
            set_pixel(R, 0, h, w, r);
        }
    }
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    image r = make_image(im.c, im.h, im.w);
    for (int c = 0; c < im.c; ++c)
    {
        for (int col = 0; col < im.h; ++col)
        {
            for (int row = 0; row < im.w; ++row)
            {
                int found = 0;
                float pixVal = get_pixel(im, c, col, row);
                for (int y = 0; y < (2 * w + 1); y++)
                {
                    for (int x = 0; x < (2 * w + 1); x++)
                    {
                        float neighbor_val = get_pixel(im, c, col + y - w, row + x - w);
                        if (neighbor_val > pixVal)
                        {
                            set_pixel(r, c, col, row, -999999);
                            found = 1;
                            break;
                        }
                    }
                    if (found)
                    {
                        break;
                    }
                }
                if (!found)
                {
                    set_pixel(r, c, col, row, pixVal);
                }
            }
        }
    }

    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    // Count number of responses over threshold
    int count = 0;
    for (int i = 0; i < Rnms.c * Rnms.h * Rnms.w; ++i)
    {
        if (Rnms.data[i] > thresh)
            ++count;
    }
    *n = count; // <- set *n equal to number of corners in image.

    descriptor *d = calloc(count, sizeof(descriptor));
    // Fill in array *d with descriptors of corners, use describe_index()
    int index = 0;
    for (int i = 0; i < Rnms.c * Rnms.h * Rnms.w; ++i)
    {
        if (Rnms.data[i] > thresh)
        {
            d[index++] = describe_index(im, i);
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
