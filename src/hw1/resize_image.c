#include <math.h>
#include "image.h"

float nn_interpolate(image im, int c, float h, float w)
{
    float pixelVal = get_pixel(im, c, round(h), round(w));
    return pixelVal;
}

image nn_resize(image im, int h, int w)
{
    image im_nnResized = make_image(im.c, h, w);
    float h_og = im.h;
    float w_og = im.w;
    
    float a_h = h_og / h;
    float b_h = 0.5 * a_h - 0.5;
    float h_resize2og;
    
    float a_w = w_og / w;
    float b_w = 0.5 * a_w - 0.5;
    float w_resize2og;

    float pixelVal_og;
    for(int c = 0; c < im.c; ++c) {
        for(int row = 0; row < h; ++row) {
            h_resize2og = a_h * row + b_h;
            for(int col = 0; col < w; ++col) {
                w_resize2og = a_w * col + b_w;
                pixelVal_og = nn_interpolate(im, c, h_resize2og, w_resize2og);
                set_pixel(im_nnResized, c, row, col, pixelVal_og);
            }
        }
    }
    return im_nnResized;
}

float bilinear_interpolate(image im, int c, float h, float w)
{
    int h_top = h;
    int h_bottom = (h + 1);
    int w_left = w;
    int w_right = (w + 1);
    float topLeft_pixVal = get_pixel(im, c, h_top, w_left);
    float topRight_pixVal = get_pixel(im, c, h_top, w_right);
    float bottomLeft_pixVal = get_pixel(im, c, h_bottom, w_left);
    float bottomRight_pixVal = get_pixel(im, c, h_bottom, w_right);

    float d_top = h - h_top;
    float d_bottom = 1 - d_top;
    float d_left = w - w_left;
    float d_right = 1 - d_left;

    float q1 = topLeft_pixVal * d_bottom + bottomLeft_pixVal * d_top;
    float q2 = topRight_pixVal * d_bottom + bottomRight_pixVal * d_top;

    float pixVal = q1 * d_right + q2 * d_left;

    return pixVal;
}

image bilinear_resize(image im, int h, int w)
{
    image im_biResized = make_image(im.c, h, w);
    float h_og = im.h;
    float w_og = im.w;
    
    float a_h = h_og / h;
    float b_h = 0.5 * a_h - 0.5;
    float h_resize2og;
    
    float a_w = w_og / w;
    float b_w = 0.5 * a_w - 0.5;
    float w_resize2og;

    float pixelVal_og;
    for(int c = 0; c < im.c; ++c) {
        for(int row = 0; row < h; ++row) {
            h_resize2og = a_h * row + b_h;
            for(int col = 0; col < w; ++col) {
                w_resize2og = a_w * col + b_w;
                pixelVal_og = bilinear_interpolate(im, c, h_resize2og, w_resize2og);
                set_pixel(im_biResized, c, row, col, pixelVal_og);
            }
        }
    }
    return im_biResized;
}

//// !! Compiler reports invalid initializer error when helper function is called. !!
// image resize_image(image im, int h, int w, char resize_method)
// {
//     //// for each channel, loop through all rows and columns
//     //// for each pixel location, performs bilinear interpolation
//     //      Bilinear interpolation:
//     image im_resized = make_image(im.c, h, w);
//     float h_og = im.h;
//     float w_og = im.w;
    
//     float a_h = h_og / h;
//     float b_h = 0.5 * a_h - 0.5;
//     float h_resize2og;
    
//     float a_w = w_og / w;
//     float b_w = 0.5 * a_w - 0.5;
//     float w_resize2og;

//     float pixelVal_og;

//     for(int c = 0; c < im.c; ++c) {
//         for(int row = 0; row < h; ++row) {
//             h_resize2og = a_h * row + b_h;
//             for(int col = 0; col < w; ++col) {
//                 w_resize2og = a_w * col + b_w;
//                 if(!strcmp(resize_method, "bilinear interpolation")) {
//                     pixelVal_og = bilinear_interpolate(im, c, h_resize2og, w_resize2og);
//                 } else {
//                     pixelVal_og = nn_interpolate(im, c, h_resize2og, w_resize2og);
//                 }
//                 set_pixel(im_resized, c, row, col, pixelVal_og);
//             }
//         }
//     }
//     return im_resized;
// }
