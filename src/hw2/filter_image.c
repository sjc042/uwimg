#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    float sum = 0;
    for (int c = 0; c < im.c; ++c)
    {
        for (int h = 0; h < im.h; ++h)
        {
            for (int w = 0; w < im.w; ++w)
            {
                sum += get_pixel(im, c, h, w);
            }
        }
    }
    for (int c = 0; c < im.c; ++c)
    {
        for (int h = 0; h < im.h; ++h)
        {
            for (int w = 0; w < im.w; ++w)
            {
                float pixVal = get_pixel(im, c, h, w);
                set_pixel(im, c, h, w, pixVal / sum);
            }
        }
    }
}

image make_box_filter(int w)
{
    // // O(c*w*w)
    image filter = make_image(1, w, w);
    for (int y = 0; y < w; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            set_pixel(filter, 0, y, x, (float)1 / (w * w));
        }
    }

    // O(2*c*w*w)
    // image filter = make_image(1, w, w);
    // for(int z = 0; z < c; ++z) {
    //     for(int y = 0; y < w; ++y) {
    //         for(int x = 0; x < w; ++x) {
    //             set_pixel(filter, z, y, x, 1);
    //         }
    //     }
    // }
    // l1_normalize(filter);
    return filter;
}

image convolve_image(image im, image filter, int preserve)
{
    // assert(filter.c == 1 || filter.c == im.c);
    // image img_conv_multiC;
    // image img_conv_1c;
    // image filter_conv = filter;

    // if (filter_conv.c != im.c) // make a multi-channel filter to preserve original image channels
    // {
    //     filter_conv = make_image(im.c, filter.h, filter.w);
    //     for (int h = 0; h < filter.h; ++h)
    //     {
    //         for (int w = 0; w < filter.w; ++w)
    //         {
    //             for (int c = 0; c < im.c; ++c)
    //             {
    //                 float pixVal = get_pixel(filter, 0, h, w); // Question: is it ok to set all channels of filter to original value ???
    //                 set_pixel(filter_conv, c, h, w, pixVal);
    //             }
    //         }
    //     }
    // }
    // // at this poit, filter_conv.c == im.c
    // int ker_len = filter_conv.w;
    // img_conv_multiC = make_image(im.c, im.h, im.w);
    // for (int c = 0; c < im.c; ++c) // go through all pixels of original image and do convolution
    // {
    //     for (int h = 0; h < im.h; ++h)
    //     {
    //         for (int w = 0; w < im.w; ++w)
    //         {
    //             float new_pixVal = 0;
    //             // do convolution on current pixel with filter_w*filter_w for-loop
    //             for (int row = 0; row < ker_len; ++row)      // assuming a NxN filter, does not handle 1d decomposed filter
    //             {
    //                 for (int col = 0; col < ker_len; ++col)  // assuming a NxN filter, does not handle 1d decomposed filter
    //                 {
    //                     float og_pixVal = get_pixel(im, c, h + row - ker_len / 2, w + col - ker_len / 2);
    //                     float coef = get_pixel(filter_conv, c, row, col);
    //                     new_pixVal += og_pixVal * coef;
    //                 }
    //             }
    //             set_pixel(img_conv_multiC, c, h, w, new_pixVal);
    //         }
    //     }
    // }
    // // img_conv_multiC is done
    // if (preserve == 1) // return img_conv_multiC, with same number of channels as original image
    // {
    //     return img_conv_multiC;
    // }
    // else // produce and return img_conv_1c, a one-channel image
    // {
    //     {
    //         // collapse channels into one by summing over spatial and channel dimensions
    //         img_conv_1c = make_image(1, im.h, im.w);

    //         for (int h = 0; h < img_conv_multiC.h; ++h)
    //         {
    //             for (int w = 0; w < img_conv_multiC.w; ++w)
    //             {
    //                 float pixVal_sum = 0;
    //                 for (int c = 0; c < img_conv_multiC.c; ++c)
    //                 {
    //                     pixVal_sum += get_pixel(img_conv_multiC, c, h, w);
    //                 }
    //                 set_pixel(img_conv_1c, 0, h, w, pixVal_sum);
    //             }
    //         }
    //         return img_conv_1c;
    //     }
    // }
    image out = make_image(preserve ? im.c : 1, im.h, im.w);
    assert(im.c == filter.c || filter.c == 1);
    int i, j, k, dx, dy;
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){      // looping through all pixels in image
                for(dy = 0; dy < filter.h; ++dy){
                    for(dx = 0; dx < filter.w; ++dx){
                        float weight = get_pixel(filter, (filter.c == 1) ? 0 : k, dy, dx);
                        float val = get_pixel(im, k, j+dy-filter.h/2, i+dx-filter.w/2);
                        if(preserve) out.data[i + im.w*(j+im.h*k)] += val*weight;
                        else         out.data[i + im.w*j] += val*weight;
                    }
                }
            }
        }
    }
    return out;

}

image make_highpass_filter()
{
    image highPass_filter = make_image(1, 3, 3);
    float values[] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
    for (int h = 0; h < 3; ++h)
    {
        for (int w = 0; w < 3; ++w)
        {
            set_pixel(highPass_filter, 0, h, w, values[w + h * 3]);
        }
    }
    return highPass_filter;
}

image make_sharpen_filter()
{
    image sharpen_filter = make_image(1, 3, 3);
    float values[] = {0, -1, 0, -1, 5, -1, 0, -1, 0};
    for (int h = 0; h < 3; ++h)
    {
        for (int w = 0; w < 3; ++w)
        {
            set_pixel(sharpen_filter, 0, h, w, values[w + h * 3]);
        }
    }
    return sharpen_filter;
}

image make_emboss_filter()
{
    // TODO
    image emboss_filter = make_image(1, 3, 3);
    float values[] = {-2, -1, 0, -1, 1, 1, 0, 1, 2};
    for (int h = 0; h < 3; ++h)
    {
        for (int w = 0; w < 3; ++w)
        {
            set_pixel(emboss_filter, 0, h, w, values[w + h * 3]);
        }
    }
    return emboss_filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: We want to preserve number of channels on sharpen and emboss filters so convolved images look enhanced comparing to original
//         Might not need to preserve channels for highpass filter since it only gives frequency information to extract features.

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Highpass filter might need post-processing to eliminate noise from the high frequency features.
//         Similar with Sharpen and Emboss filters, filter values add up to one. But the values for those are also large and noisy.

image make_gaussian_filter(float sigma)
{
    // TODO
    int ker_len = (int)sigma * 6 + 1;
    image gaussian_filter = make_image(1, ker_len, ker_len);
    float G;
    for (int h = ker_len / -2; h < ker_len / 2 + 1; ++h)
    {
        for (int w = ker_len / -2; w < ker_len / 2 + 1; ++w)
        {
            G = 1 / (TWOPI * sigma * sigma) * exp(-1 * (w * w + h * h) / (2 * sigma * sigma));
            set_pixel(gaussian_filter, 0, h + ker_len / 2, w + ker_len / 2, G);
        }
    }
    l1_normalize(gaussian_filter);
    return gaussian_filter;
}

image add_image(image a, image b)
{
    // Checking images' sizes
    assert(a.h == b.h && a.w == b.w);
    int high_c = a.c > b.c ? a.c : b.c;
    image im_add = make_image(high_c, a.h, a.w);

    for (int h = 0; h < a.h; ++h)
    {
        for (int w = 0; w < a.w; ++w)
        {
            float pixVal;
            for (int c = 0; c < high_c; ++c)
            {
                if (a.c < b.c)
                {
                    pixVal = get_pixel(a, 0, h, w) + get_pixel(b, c, h, w);
                }
                else if (a.c > b.c)
                {
                    pixVal = get_pixel(a, c, h, w) + get_pixel(b, 0, h, w);
                }
                else
                {
                    pixVal = get_pixel(a, c, h, w) + get_pixel(b, c, h, w);
                }
                set_pixel(im_add, c, h, w, pixVal);
            }
        }
    }
    return im_add;
}

image sub_image(image a, image b)
{
    assert(a.h == b.h && a.w == b.w);
    // if image don't have the same channels, upsampling image with fewer channels to match that of more
    int high_c = a.c > b.c ? a.c : b.c;
    // make upsampled to the same number of channels
    image im_sub = make_image(high_c, a.h, a.w);

    for (int h = 0; h < a.h; ++h)
    {
        for (int w = 0; w < a.w; ++w)
        {
            float pixVal;
            for (int c = 0; c < high_c; ++c)
            {
                if (a.c < b.c)
                {
                    pixVal = get_pixel(a, 0, h, w) - get_pixel(b, c, h, w);
                }
                else if (a.c > b.c)
                {
                    pixVal = get_pixel(a, c, h, w) - get_pixel(b, 0, h, w);
                }
                else
                {
                    pixVal = get_pixel(a, c, h, w) - get_pixel(b, c, h, w);
                }
                set_pixel(im_sub, c, h, w, pixVal);
            }
        }
    }
    return im_sub;
}

image make_gx_filter()
{
    image gx_filter = make_image(1, 3, 3);
    float values[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    for (int h = 0; h < 3; ++h)
    {
        for (int w = 0; w < 3; ++w)
        {
            set_pixel(gx_filter, 0, h, w, values[w + h * 3]);
        }
    }
    return gx_filter;
}

image make_gy_filter()
{
    image gy_filter = make_image(1, 3, 3);
    float values[] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    for (int h = 0; h < 3; ++h)
    {
        for (int w = 0; w < 3; ++w)
        {
            set_pixel(gy_filter, 0, h, w, values[w + h * 3]);
        }
    }
    return gy_filter;
}

void feature_normalize(image im)
{
    // to normalize and visualize
    // loop through all pixels in im
    //  get minimum pixel value, maximum pixel value, range=max-min
    // if range == 0, set whole image to zero
    // else normalize all pixel to real value [0,1]
    float min = 256; // setting to 1 is enough, setting to 256 to be safe
    float max = 0;
    float range = max - min;
    float pixVal;
    for (int c = 0; c < im.c; ++c)
    {
        for (int h = 0; h < im.h; ++h)
        {
            for (int w = 0; w < im.w; ++w)
            {
                pixVal = get_pixel(im, c, h, w);
                min = pixVal < min ? pixVal : min;
                max = pixVal > max ? pixVal : max;
                range = max - min;
            }
        }
    }
    // loop through all pixels in im and normalize pixel value to [0,1]
    for (int c = 0; c < im.c; ++c)
    {
        for (int h = 0; h < im.h; ++h)
        {
            for (int w = 0; w < im.w; ++w)
            {
                pixVal = get_pixel(im, c, h, w);
                set_pixel(im, c, h, w, (pixVal - min) / range);
            }
        }
    }
}

image *sobel_image(image im)
{
    // get grad_x image
    // get grad_y image
    // loop through all pixels in each synchronously
    //      calculate magnitude and direction
    image grad_x_filter = make_gx_filter();
    image grad_y_filter = make_gy_filter();
    image im_grad_x = convolve_image(im, grad_x_filter, 0);
    image im_grad_y = convolve_image(im, grad_y_filter, 0);

    float grad_x_val;
    float grad_y_val;
    float magnitude;
    float angle;

    image im_magnitude = make_image(1, im.h, im.w);
    image im_angle = make_image(1, im.h, im.w);

    int c = 0;
    for (int h = 0; h < im_grad_x.h; ++h)
    {
        for (int w = 0; w < im_grad_x.w; ++w)
        {
            grad_x_val = get_pixel(im_grad_x, c, h, w);
            grad_y_val = get_pixel(im_grad_y, c, h, w);
            magnitude = sqrtf(grad_x_val * grad_x_val + grad_y_val * grad_y_val);
            angle = atan2f(grad_y_val, grad_x_val);
            set_pixel(im_magnitude, c, h, w, magnitude);
            set_pixel(im_angle, c, h, w, angle);
        }
    }
    image *result = calloc(2, sizeof(image));
    result[0] = im_magnitude;
    result[1] = im_angle;
    return result;
}

image colorize_sobel(image im)
{
    // preprocess image by smoothing with gaussian filter
    image gauss_f = make_gaussian_filter(2);
    im = convolve_image(im, gauss_f, 1);
    image *res = sobel_image(im);
    image mag = res[0];
    image theta = res[1];
    feature_normalize(mag);
    feature_normalize(theta);
    // make a HSV image with magnitude being saturation and value, theta being hue
    image im_return = make_image(3, im.h, im.w);
    for(int h = 0; h < im.h; ++h) {
        for(int w = 0; w < im.w; ++w) {
            float pixVal_mag = get_pixel(mag, 0, h, w);
            float pixVal_theta = get_pixel(theta, 0, h, w);

            set_pixel(im_return, 0, h, w, pixVal_theta);
            set_pixel(im_return, 1, h, w, pixVal_mag);
            set_pixel(im_return, 2, h, w, pixVal_mag);
        }
    }
    // convert HSV image to RGV image
    hsv_to_rgb(im_return);
    return im_return;
}