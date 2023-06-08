#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// Draws a line on an image with color corresponding to the direction of line
// image im: image to draw line on
// float x, y: starting point of line
// float dx, dy: vector corresponding to line angle and magnitude
void draw_line(image im, float y, float x, float dy, float dx)
{
    assert(im.c == 3);
    float angle = 6 * (atan2(dy, dx) / TWOPI + .5);
    int index = floor(angle);
    float f = angle - index;
    float r, g, b;
    if (index == 0)
    {
        r = 1;
        g = f;
        b = 0;
    }
    else if (index == 1)
    {
        r = 1 - f;
        g = 1;
        b = 0;
    }
    else if (index == 2)
    {
        r = 0;
        g = 1;
        b = f;
    }
    else if (index == 3)
    {
        r = 0;
        g = 1 - f;
        b = 1;
    }
    else if (index == 4)
    {
        r = f;
        g = 0;
        b = 1;
    }
    else
    {
        r = 1;
        g = 0;
        b = 1 - f;
    }
    float i;
    float d = sqrt(dx * dx + dy * dy);
    for (i = 0; i < d; i += 1)
    {
        int xi = x + dx * i / d;
        int yi = y + dy * i / d;
        set_pixel(im, 0, yi, xi, r);
        set_pixel(im, 1, yi, xi, g);
        set_pixel(im, 2, yi, xi, b);
    }
}

// Make an integral image or summed area table from an image
// image im: image to process
// returns: image I such that I[x,y] = sum{i<=x, j<=y}(im[i,j])
image make_integral_image(image im)
{
    image integ = make_image(im.c, im.h, im.w);
    // TODO: fill in the integral image
    for (int c = 0; c < im.c; ++c)
    {
        for (int h = 0; h < im.h; ++h)
        {
            for (int w = 0; w < im.w; ++w)
            {
                float i_xy = get_pixel(im, c, h, w);
                float top_left_pixVal = ((w - 1) < 0 || (h - 1) < 0) ? 0 : get_pixel(integ, c, h - 1, w - 1);
                float left_pixVal = ((w - 1) < 0) ? 0 : get_pixel(integ, c, h, w - 1);
                float top_pixVal = ((h - 1) < 0) ? 0 : get_pixel(integ, c, h - 1, w);
                float I_xy = i_xy + top_pixVal + left_pixVal - top_left_pixVal;
                set_pixel(integ, c, h, w, I_xy);
            }
        }
    }
    return integ;
}

// Apply a box filter to an image using an integral image for speed
// image im: image to smooth
// int s: window size for box filter
// returns: smoothed image
image box_filter_image(image im, int s)
{
    int i, j, k;
    image integ = make_integral_image(im);
    image S = make_image(im.c, im.h, im.w);
    // TODO: fill in S using the integral image.
    for (int i = 0; i < im.c; ++i) // iterate through channels of S
    {
        for (int y = 0; y < im.h; ++y) // iterate through rows of S
        {
            for (int x = 0; x < im.w; ++x) // iterate through columns of S
            {
                int left = (x < s / 2) ? 0 : (x - s / 2);
                int right = ((x + s / 2) < (im.w - 1)) ? (x + s / 2) : (im.w - 1);
                int top = (y < s / 2) ? 0 : (y - s / 2);
                int bottom = ((y + s / 2) < (im.h - 1)) ? (y + s / 2) : (im.h - 1);
                int width = right - left + 1;
                int height = bottom - top + 1;
                float D = get_pixel(integ, i, bottom, right);
                float A = get_pixel(integ, i, top - 1, left - 1);
                float B = get_pixel(integ, i, top - 1, right);
                float C = get_pixel(integ, i, bottom, left - 1);
                float sum;
                if (left == 0 && top == 0)
                {
                    sum = D;
                }
                else if (left > 0 && top == 0)
                {
                    sum = D - C;
                }
                else if (left == 0 && top > 0)
                {
                    sum = D - B;
                }
                else if (left > 0 && top > 0)
                {
                    sum = D + A - B - C;
                }
                set_pixel(S, i, y, x, sum / (width * height));
            }
        }
    }

    free_image(integ);
    return S;
}

// Calculate the time-structure matrix of an image pair.
// image im: the input image.
// image prev: the previous image in sequence.
// int s: window size for smoothing.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          3rd channel is IxIy, 4th channel is IxIt, 5th channel is IyIt.
image time_structure_matrix(image im, image prev, int s)
{
    int i;
    int converted = 0;
    if (im.c == 3)
    {
        converted = 1;
        im = rgb_to_grayscale(im);
        prev = rgb_to_grayscale(prev);
    }

    // TODO: calculate gradients, structure components, and smooth them
    image S = make_image(5, im.h, im.w);
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image Ix = convolve_image(im, gx_filter, 0);
    image Iy = convolve_image(im, gy_filter, 0);
    image Ix2 = make_image(Ix.c, Ix.h, Ix.w);
    image Iy2 = make_image(Iy.c, Iy.h, Iy.w);
    image IxIy = make_image(Ix.c, Ix.h, Ix.w);
    image IxIt = make_image(Ix.c, Ix.h, Ix.w);
    image IyIt = make_image(Iy.c, Iy.h, Iy.w);
    for (int i = 0; i < Ix.c * Ix.h * Ix.w; ++i)
    {
        Ix2.data[i] = Ix.data[i] * Ix.data[i];
        Iy2.data[i] = Iy.data[i] * Iy.data[i];
        IxIy.data[i] = Ix.data[i] * Iy.data[i];
        IxIt.data[i] = Ix.data[i] * (im.data[i] - prev.data[i]);
        IyIt.data[i] = Iy.data[i] * (im.data[i] - prev.data[i]);
    }
    Ix2 = box_filter_image(Ix2, s);
    Iy2 = box_filter_image(Iy2, s);
    IxIy = box_filter_image(IxIy, s);
    IxIt = box_filter_image(IxIt, s);
    IyIt = box_filter_image(IyIt, s);

    for (int h = 0; h < S.h; ++h)
    {
        for (int w = 0; w < S.w; ++w)
        {
            float Ix2_val = get_pixel(Ix2, 0, h, w);
            float Iy2_val = get_pixel(Iy2, 0, h, w);
            float IxIy_val = get_pixel(IxIy, 0, h, w);
            float IxIt_val = get_pixel(IxIt, 0, h, w);
            float IyIt_val = get_pixel(IyIt, 0, h, w);
            set_pixel(S, 0, h, w, Ix2_val);
            set_pixel(S, 1, h, w, Iy2_val);
            set_pixel(S, 2, h, w, IxIy_val);
            set_pixel(S, 3, h, w, IxIt_val);
            set_pixel(S, 4, h, w, IyIt_val);
        }
    }

    if (converted)
    {
        free_image(im);
        free_image(prev);
    }
    return S;
}

// Calculate the velocity given a structure image
// image S: time-structure image
// int stride:
image velocity_image(image S, int stride)
{
    image v = make_image(3, S.h / stride, S.w / stride);
    int i, j;
    matrix M = make_matrix(2, 2);
    matrix xyT = make_matrix(2, 1);
    matrix Vxy;
    for (j = (stride - 1) / 2; j < S.h; j += stride)
    {
        for (i = (stride - 1) / 2; i < S.w; i += stride)
        {
            float Ixx = S.data[i + S.w * j + 0 * S.w * S.h];
            float Iyy = S.data[i + S.w * j + 1 * S.w * S.h];
            float Ixy = S.data[i + S.w * j + 2 * S.w * S.h];
            float Ixt = S.data[i + S.w * j + 3 * S.w * S.h];
            float Iyt = S.data[i + S.w * j + 4 * S.w * S.h];
            float vx, vy;
            // TODO: calculate vx and vy using the flow equation
            if ((Ixx * Iyy - Ixy * Ixy) > 0 || (Ixx * Iyy - Ixy * Ixy) < 0)
            {
                M.data[0][0] = Ixx;
                M.data[0][1] = Ixy;
                M.data[1][0] = Ixy;
                M.data[1][1] = Iyy;
                matrix M_inv = matrix_invert(M);
                xyT.data[0][0] = -Ixt;
                xyT.data[1][0] = -Iyt;
                Vxy = matrix_mult_matrix(M_inv, xyT);
                vx = Vxy.data[0][0];
                vy = Vxy.data[1][0];
            }
            else
            {
                vx = 0;
                vy = 0;
            }

            set_pixel(v, 0, j / stride, i / stride, vx);
            set_pixel(v, 1, j / stride, i / stride, vy);
        }
    }
    free_matrix(M);
    return v;
}

// Draw lines on an image given the velocity
// image im: image to draw on
// image v: velocity of each pixel
// float scale: scalar to multiply velocity by for drawing
void draw_flow(image im, image v, float scale)
{
    int stride = im.w / v.w;
    int i, j;
    for (j = (stride - 1) / 2; j < im.h; j += stride)
    {
        for (i = (stride - 1) / 2; i < im.w; i += stride)
        {
            float dx = scale * get_pixel(v, 0, j / stride, i / stride);
            float dy = scale * get_pixel(v, 1, j / stride, i / stride);
            if (fabs(dx) > im.w)
                dx = 0;
            if (fabs(dy) > im.h)
                dy = 0;
            draw_line(im, j, i, dy, dx);
        }
    }
}

// Constrain the absolute value of each image pixel
// image im: image to constrain
// float v: each pixel will be in range [-v, v]
void constrain_image(image im, float v)
{
    int i;
    for (i = 0; i < im.w * im.h * im.c; ++i)
    {
        if (im.data[i] < -v)
            im.data[i] = -v;
        if (im.data[i] > v)
            im.data[i] = v;
    }
}

// Calculate the optical flow between two images
// image im: current image
// image prev: previous image
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// returns: velocity matrix
image optical_flow_images(image im, image prev, int smooth, int stride)
{
    image S = time_structure_matrix(im, prev, smooth);
    image v = velocity_image(S, stride);
    constrain_image(v, 6);
    image vs = smooth_image(v, 2);
    free_image(v);
    free_image(S);
    return vs;
}

// Run optical flow demo on webcam
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// int div: downsampling factor for images from webcam
void optical_flow_webcam(int smooth, int stride, int div)
{
#ifdef OPENCV
    void *cap;
    // What video stream you open
    // cap = open_video_stream(0, 1, 0, 0, 0);
    cap = open_video_stream(0, 0, 0, 0, 30);
    printf("%ld\n", cap);
    if (!cap)
    {
        fprintf(stderr, "couldn't open\n");
        exit(0);
    }
    image prev = get_image_from_stream(cap);
    printf("%d %d\n", prev.w, prev.h);
    image prev_c = nn_resize(prev, prev.h / div, prev.w / div);
    image im = get_image_from_stream(cap);
    image im_c = nn_resize(im, im.h / div, im.w / div);
    while (im.data)
    {
        image copy = copy_image(im);
        image v = optical_flow_images(im_c, prev_c, smooth, stride);
        draw_flow(copy, v, smooth * div * 2);
        int key = show_image(copy, "flow", 5);
        free_image(v);
        free_image(copy);
        free_image(prev);
        free_image(prev_c);
        prev = im;
        prev_c = im_c;
        if (key != -1)
        {
            key = key % 256;
            printf("%d\n", key);
            if (key == 27)
                break;
        }
        im = get_image_from_stream(cap);
        im_c = nn_resize(im, im.h / div, im.w / div);
    }
#else
    fprintf(stderr, "Must compile with OpenCV\n");
#endif
}
