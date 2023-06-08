#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int c, int h, int w)
{
    int im_dim[3] = {im.c, im.h, im.w};
    int get_coord[3] = {c, h, w};
    // Apply clamp padding
    for (int i = 0; i < 3; ++i) {
        if(get_coord[i] < 0) {
            get_coord[i] = 0;
        } else if(get_coord[i] > (im_dim[i]-1)) {
            get_coord[i] = im_dim[i] - 1;
        }
    }
    int C = get_coord[0];
    int H = get_coord[1];
    int W = get_coord[2];
    float pixel_value = im.data[W + H*im.w + C*im.w*im.h];
    return pixel_value;
}

void set_pixel(image im, int c, int h, int w, float v)
{
    int set_coord[3] = {c, h, w};
    int im_dim[3] = {im.c, im.h, im.w};
    for (int i = 0; i < 3; ++i) {
        if(set_coord[i] < 0 || set_coord[i] > im_dim[i]-1) {
            return;
        }
    }
    im.data[w + h*im.w + c*im.w*im.h] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.c, im.h, im.w);
    for(int c = 0; c < im.c; ++c) {
        for(int h = 0; h < im.h; ++h) {
            for(int w = 0; w < im.w; ++w) {
                float pixelVal = get_pixel(im, c, h, w);
                set_pixel(copy, c, h, w, pixelVal);
            }
        }
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(1, im.h, im.w);
    for(int h = 0; h < im.h; ++h) {
        for(int w = 0; w < im.w; ++w) {
            float RGB[3];
            for(int c = 0; c < 3; ++c) {
                RGB[c] = get_pixel(im, c, h, w);
            }
            float R = RGB[0];
            float G = RGB[1];
            float B = RGB[2];

            float Y_luma = 0.299 * R + 0.587 * G + 0.114 * B;
            set_pixel(gray, 0, h, w, Y_luma);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    for(int h = 0; h < im.h; ++h) {
        for(int w = 0; w < im.w; ++w) {
            float pixelVal = get_pixel(im, c, h, w) + v;
            set_pixel(im, c, h, w, pixelVal);
        }
    }
}

void clamp_image(image im)
{
    for(int c = 0; c < im.c; ++c) {
        for(int h = 0; h < im.h; ++h) {
            for(int w = 0; w < im.w; ++w) {
                float pixelVal = get_pixel(im, c, h, w);
                if(pixelVal < 0) {
                    set_pixel(im, c, h, w, 0);
                } else if(pixelVal > 1) {
                    set_pixel(im, c, h, w, 1);
                }
            }
        }
    }
}

// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    for(int h = 0; h < im.h; ++h) {
        for(int w = 0; w < im.w; ++w) {
            float RGB[3];
            for(int c = 0; c < 3; ++c) {
                RGB[c] = get_pixel(im, c, h, w);
            }
            float R = RGB[0];
            float G = RGB[1];
            float B = RGB[2];
            // getting S, V
            float M = three_way_max(R, G, B);
            float V = M;
            float m = three_way_min(R, G, B);
            float C = V - m;
            float S = (V > 0) ? (C / V) : 0;
            // getting H
            float H_1;
            if((C-0.002)<0 && 0<(C+0.002)) {               // C==0
                H_1 = 0;
            } else if((V-0.002)<R && R<(V+0.002)) {        // V==R
                H_1 = (G - B) / C;
            } else if((V-0.002)<G && G<(V+0.002)) {        // V==G
                H_1 = (B - R) / C + 2;
            } else if((V-0.002)<B && B<(V+0.002)) {        // V==B
                H_1 = (R - G) / C + 4;
            }
            float H = (H_1 < 0) ? (H_1 / 6 + 1) : (H_1 / 6);
            float HSV[3] = {H, S, V};
            for (int c = 0; c < 3; ++c) {
                set_pixel(im, c, h, w, HSV[c]);
            }
        }
    }
}

void hsv_to_rgb(image im)
{
    for(int h = 0; h < im.h; ++h) {
        for(int w = 0; w < im.w; ++w) {
            float HSV[3];
            float RGB[3];
            for(int c = 0; c < 3; ++c) {
                HSV[c] = get_pixel(im, c, h, w);
            }

            float H = HSV[0];
            float S = HSV[1];
            float V = HSV[2];
            float C = V * S;
            float m = V - C;
            float H_1 = H * 6;
            float X = C * (1 - fabs(fmodf(H_1, 2) - 1));

            if(H_1 >= 0 && H_1 < 1) {
                RGB[0] = C;
                RGB[1] = X;
                RGB[2] = 0;
            } else if(H_1 >= 1 && H_1 < 2) {
                RGB[0] = X;
                RGB[1] = C;
                RGB[2] = 0;
            } else if(H_1 >= 2 && H_1 < 3) {
                RGB[0] = 0;
                RGB[1] = C;
                RGB[2] = X;
            } else if(H_1 >= 3 && H_1 < 4) {
                RGB[0] = 0;
                RGB[1] = X;
                RGB[2] = C;
            } else if(H_1 >= 4 && H_1 < 5) {
                RGB[0] = X;
                RGB[1] = 0;
                RGB[2] = C;
            } else {        // H >= 5 && H < 6
                RGB[0] = C;
                RGB[1] = 0;
                RGB[2] = X;
            }
            for(int c = 0; c < 3; ++c) {
                RGB[c] += m;
                set_pixel(im, c, h, w, RGB[c]);
            }
        }
    }
}