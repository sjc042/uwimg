#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// Comparator for matches
// const void *a, *b: pointers to the matches to compare.
// returns: result of comparison, 0 if same, 1 if a > b, -1 if a < b.
int match_compare(const void *a, const void *b)
{
    match *ra = (match *)a;
    match *rb = (match *)b;
    if (ra->distance < rb->distance)
        return -1;
    else if (ra->distance > rb->distance)
        return 1;
    else
        return 0;
}

// Helper function to create 2d points.
// float x, y: coordinates of point.
// returns: the point.
point make_point(float x, float y)
{
    point p;
    p.x = x;
    p.y = y;
    return p;
}

// Place two images side by side on canvas, for drawing matching pixels.
// image a, b: images to place.
// returns: image with both a and b side-by-side.
image both_images(image a, image b)
{
    image both = make_image(a.c > b.c ? a.c : b.c, a.h > b.h ? a.h : b.h, a.w + b.w);
    int i, j, k;
    for (k = 0; k < a.c; ++k)
    {
        for (j = 0; j < a.h; ++j)
        {
            for (i = 0; i < a.w; ++i)
            {
                set_pixel(both, k, j, i, get_pixel(a, k, j, i));
            }
        }
    }
    for (k = 0; k < b.c; ++k)
    {
        for (j = 0; j < b.h; ++j)
        {
            for (i = 0; i < b.w; ++i)
            {
                set_pixel(both, k, j, i + a.w, get_pixel(b, k, j, i));
            }
        }
    }
    return both;
}

// Draws lines between matching pixels in two images.
// image a, b: two images that have matches.
// match *matches: array of matches between a and b.
// int n: number of matches.
// int inliers: number of inliers at beginning of matches, drawn in green.
// returns: image with matches drawn between a and b on same canvas.
image draw_matches(image a, image b, match *matches, int n, int inliers)
{
    image both = both_images(a, b);
    int i, j;
    for (i = 0; i < n; ++i)
    {
        int bx = matches[i].p.x;
        int ex = matches[i].q.x;
        int by = matches[i].p.y;
        int ey = matches[i].q.y;
        for (j = bx; j < ex + a.w; ++j)
        {
            int r = (float)(j - bx) / (ex + a.w - bx) * (ey - by) + by;
            set_pixel(both, 0, r, j, i < inliers ? 0 : 1);
            set_pixel(both, 1, r, j, i < inliers ? 1 : 0);
            set_pixel(both, 2, r, j, 0);
        }
    }
    return both;
}

// Draw the matches with inliers in green between two images.
// image a, b: two images to match.
// matches *
image draw_inliers(image a, image b, matrix H, match *m, int n, float thresh)
{
    int inliers = model_inliers(H, m, n, thresh);
    image lines = draw_matches(a, b, m, n, inliers);
    return lines;
}

// Find corners, match them, and draw them between two images.
// image a, b: images to match.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
image find_and_draw_matches(image a, image b, float sigma, float thresh, int nms)
{
    int an = 0;
    int bn = 0;
    int mn = 0;
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    mark_corners(a, ad, an);
    mark_corners(b, bd, bn);
    image lines = draw_matches(a, b, m, mn, 0);

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);
    return lines;
}

// Calculates L1 distance between two floating point arrays.
// float *a, *b: arrays to compare.
// int n: number of values in each array.
// returns: l1 distance between arrays (sum of absolute differences).
float l1_distance(float *a, float *b, int n)
{
    // TODO: return the correct number.
    float sum = 0;
    for (int i = 0; i < n; ++i)
    {
        sum += fabs(a[i] - b[i]);
    }
    return sum;
}

// Finds best matches between descriptors of two images.
// descriptor *a, *b: array of descriptors for pixels in two images.
// int an, bn: number of descriptors in arrays a and b.
// int *mn: pointer to number of matches found, to be filled in by function.
// returns: best matches found. each descriptor in a should match with at most
//          one other descriptor in b.
match *match_descriptors(descriptor *a, int an, descriptor *b, int bn, int *mn)
{
    int i, j;

    // We will have at most an matches.
    *mn = an;
    match *m = calloc(an, sizeof(match));
    for (j = 0; j < an; ++j)
    {
        // TODO: for every descriptor in a, find best match in b.
        float least_dist = 999999;
        int bind = -1;
        for (i = 0; i < bn; ++i)
        {
            float l1_dist = l1_distance(a[j].data, b[i].data, a[j].n);
            if (least_dist > l1_dist)
            {
                least_dist = l1_dist;
                bind = i;
            }
        }
        // record ai as the index in *a and bi as the index in *b.]
        m[j].ai = j;
        m[j].bi = bind; // <- should be index in b.
        m[j].p = a[j].p;
        m[j].q = b[bind].p;
        m[j].distance = least_dist; // <- should be the smallest L1 distance!
    }

    int count = 1;
    int *seen = calloc(bn, sizeof(int));
    for (int y = 0; y < bn; ++y) // initialize all elements to -1 instead of 0
    {
        seen[y] = -1;
    }

    // TODO: we want matches to be injective (one-to-one).
    // Sort matches based on distance using match_compare and qsort.
    // Then throw out matches to the same element in b. Use seen to keep track.
    // Each point should only be a part of one match.
    // Some points will not be in a match.
    // In practice just bring good matches to front of list, set *mn.
    qsort(m, an, sizeof(match), match_compare);
    seen[0] = m[0].bi;
    for (int z = 1; z < an; z++) // iterate through m
    {
        // iterate through each match in sorted m based on distance low to high
        // if the match's b index is seen before,
        //      drop the current match
        int found_duplicate = 0;
        for (int y = 0; y < bn; ++y) // iterate through seen
        {
            if (m[z].bi == seen[y])
            {
                for (int x = z; x < an - 1; ++x)
                {
                    m[x] = m[x + 1];
                }
                --z;
                --an;
                found_duplicate = 1;
                break;
            }
        }
        if (!found_duplicate)
        {
            seen[z] = m[z].bi;
            count++;
        }
    }
    *mn = count;
    free(seen);
    return m;
}

// Apply a projective transformation to a point.
// matrix H: homography to project point.
// point p: point to project.
// returns: point projected using the homography.
point project_point(matrix H, point p)
{
    matrix c = make_matrix(3, 1);
    // TODO: project point p with homography H.
    // Remember that homogeneous coordinates are equivalent up to scalar.
    // Have to divide by.... something...
    c.data[0][0] = p.x;
    c.data[1][0] = p.y;
    c.data[2][0] = 1; // c=[x, y, w] while arbituarily setting w to 1
    matrix P = matrix_mult_matrix(H, c);
    double x = P.data[0][0];
    double y = P.data[1][0];
    double w = P.data[2][0];
    point q = make_point(x / w, y / w); // Normoalize on new w
    return q;
}

// Calculate L2 distance between two points.
// point p, q: points.
// returns: L2 distance between them.
float point_distance(point p, point q)
{
    // TODO: should be a quick one.
    float dist = sqrtf((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y));
    return dist;
}

// Count number of inliers in a set of matches. Should also bring inliers
// to the front of the array.
// matrix H: homography between coordinate systems.
// match *m: matches to compute inlier/outlier.
// int n: number of matches in m.
// float thresh: threshold to be an inlier.
// returns: number of inliers whose projected point falls within thresh of
//          their match in the other image. Should also rearrange matches
//          so that the inliers are first in the array. For drawing.
int model_inliers(matrix H, match *m, int n, float thresh)
{
    int i;
    int count = 0;
    // TODO: count number of matches that are inliers
    // i.e. distance(H*p, q) < thresh
    // Also, sort the matches m so the inliers are the first 'count' elements.
    for (i = 0; i < n; ++i)
    {
        match curr = m[i];
        point p = curr.p;
        point q = curr.q;
        point q_proj = project_point(H, p);
        if (!(point_distance(q_proj, q) < thresh))
        {
            m[i] = m[n - 1];
            m[n - 1] = curr;
            n--;
            i--;
        }
    }
    count = n;
    return count;
}

// Randomly shuffle matches for RANSAC.
// match *m: matches to shuffle in place.
// int n: number of elements in matches.
void randomize_matches(match *m, int n)
{
    // TODO: implement Fisher-Yates to shuffle the array.
    int i;
    int j;

    for (i + 0; i < n - 2; ++i)
    {
        j = rand() % (n - i);
        match temp = m[i];
        m[i] = m[j];
        m[j] = temp;
    }
}

// Computes homography between two images given matching pixels.
// match *matches: matching points between images.
// int n: number of matches to use in calculating homography.
// returns: matrix representing homography H that maps image a to image b.
matrix compute_homography(match *matches, int n)
{
    matrix M = make_matrix(n * 2, 8);
    matrix b = make_matrix(n * 2, 1);
    int i;
    for (i = 0; i < M.rows; i += 2)
    {
        double x = matches[i / 2].p.x;
        double xp = matches[i / 2].q.x;
        double y = matches[i / 2].p.y;
        double yp = matches[i / 2].q.y;
        // TODO: fill in the matrices M and b.
        b.data[i][0] = xp;
        b.data[i + 1][0] = yp;
        double row[2][8] = {{x, y, 1, 0, 0, 0, -(x * xp), -(y * xp)}, {0, 0, 0, x, y, 1, -(x * yp), -(y * yp)}};
        for (int j = 0; j < 8; ++j)
        {
            M.data[i][j] = row[0][j];
            M.data[i + 1][j] = row[1][j];
        }
    }
    matrix a = solve_system(M, b);
    free_matrix(M);
    free_matrix(b);

    // If a solution can't be found, return empty matrix;
    matrix none = {0};
    if (!a.data)
        return none;

    matrix H = make_matrix(3, 3);
    // TODO: fill in the homography H based on the result in a.
    for (i = 0; i < H.rows - 1; ++i)
    {
        for (int j = 0; j < H.rows; ++j)
        {
            H.data[i][j] = a.data[i * 3 + j][0];
        }
    }
    H.data[2][0] = a.data[6][0];
    H.data[2][1] = a.data[7][0];
    H.data[2][2] = 1;
    free_matrix(a);

    return H;
}

// Perform RANdom SAmple Consensus to calculate homography for noisy matches.
// match *m: set of matches.
// int n: number of matches.
// float thresh: inlier/outlier distance threshold.
// int k: number of iterations to run.
// int cutoff: inlier cutoff to exit early.
// returns: matrix representing most common homography between matches.
matrix RANSAC(match *m, int n, float thresh, int k, int cutoff)
{
    int e;
    int best = 0;
    matrix Hb = make_translation_homography(256, 0);
    //         if it's better than the cutoff:
    //             return it immediately
    // if we get to the end return the best homography
    int pairs = 4;                  // need at least four matches to solve 8 degree of freedom matrix
    for (e = 0; e < k; ++e)
    {
        randomize_matches(m, n);
        matrix H = compute_homography(m, pairs);

        int count = model_inliers(H, m, n, thresh);
        if (count > best)
        {
            Hb = compute_homography(m, count);
            best = model_inliers(Hb, m, n, thresh);
            if (best > cutoff)
                return Hb;
        }
    }
    return Hb;
}

// Stitches two images together using a projective transformation.
// image a, b: images to stitch.
// matrix H: homography from image a coordinates to image b coordinates.
// returns: combined image stitched together.
image combine_images(image a, image b, matrix H)
{

    matrix Hinv = matrix_invert(H);
    // Project the corners of image b into image a coordinates.
    point c1 = project_point(Hinv, make_point(0, 0));
    point c2 = project_point(Hinv, make_point(b.w - 1, 0));
    point c3 = project_point(Hinv, make_point(0, b.h - 1));
    point c4 = project_point(Hinv, make_point(b.w - 1, b.h - 1));
    // Find top left and bottom right corners of image b warped into image a.
    point topleft, botright;
    botright.x = MAX(c1.x, MAX(c2.x, MAX(c3.x, c4.x)));
    botright.y = MAX(c1.y, MAX(c2.y, MAX(c3.y, c4.y)));
    topleft.x = MIN(c1.x, MIN(c2.x, MIN(c3.x, c4.x)));
    topleft.y = MIN(c1.y, MIN(c2.y, MIN(c3.y, c4.y)));

    int y_top, y_bot, x_left, x_right;
    y_top = topleft.y;
    y_bot = botright.y;
    x_left = topleft.x;
    x_right = botright.x;

    // Find how big our new image should be and the offsets from image a.
    int dx = MIN(0, topleft.x);
    int dy = MIN(0, topleft.y);
    int w = MAX(a.w, botright.x) - dx;
    int h = MAX(a.h, botright.y) - dy;

    // // Can disable this if you are making very big panoramas.
    // // Usually this means there was an error in calculating H.
    // if (w > 7000 || h > 7000)
    // {
    //     fprintf(stderr, "output too big, stopping\n");
    //     // return copy_image(a);
    // }

    int i, j, k;
    image c = make_image(a.c, h, w);

    // Paste image a into the new image offset by dx and dy.
    for (k = 0; k < a.c; ++k)
    {
        for (j = 0; j < a.h; ++j)
        {
            for (i = 0; i < a.w; ++i)
            {
                // TODO: fill in.
                float aVal = get_pixel(a, k, j, i);
                int h = j - dy;
                int w = i - dx;
                set_pixel(c, k, h, w, aVal);
            }
        }
    }

    // TODO: Paste in image b as well.
    // You should loop over some points in the new image (which? all?)
    // and see if their projection from a coordinates to b coordinates falls
    // inside of the bounds of image b. If so, use bilinear interpolation to
    // estimate the value of b at that projection, then fill in image c.

    // from top left to bottom right of image b projected in c coordinates,
    // project each point by H, see projected point is in b
    // if point is in b, get pixel value in b using bl interpolation
    for (int z = 0; z < c.c; ++z)
    {
        for (float y = y_top; y < y_bot; ++y)
        {
            for (float x = x_left; x < x_right; ++x)
            {

                point c_point = make_point(x, y);
                point c2b = project_point(H, c_point);
                if ((c2b.x < b.w && (c2b.x > 0 || (c2b.x < 0.002 && c2b.x > -0.002))) && (c2b.y < b.h && (c2b.y > 0 || (c2b.y < 0.002 && c2b.y > -0.002))))
                {
                    // get pixel value of b at [c2b.x, c2b.y, z] using bilinear interpolation
                    float bVal = bilinear_interpolate(b, z, c2b.y, c2b.x);
                    // set pixel value of c at [x, y, z]
                    set_pixel(c, z, y - dy, x - dx, bVal);
                }
            }
        }
    }

    return c;
}

// Create a panoramam between two images.
// image a, b: images to stitch together.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
// float inlier_thresh: threshold for RANSAC inliers. Typical: 2-5
// int iters: number of RANSAC iterations. Typical: 1,000-50,000
// int cutoff: RANSAC inlier cutoff. Typical: 10-100
image panorama_image(image a, image b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff)
{
    srand(10);
    int an = 0;
    int bn = 0;
    int mn = 0;

    // Calculate corners and descriptors
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);
    // Find matches
    match *m = match_descriptors(ad, an, bd, bn, &mn);
    // Run RANSAC to find the homography
    matrix H = RANSAC(m, mn, inlier_thresh, iters, cutoff);
    if (0)
    {
        // Mark corners and matches between images, turn this off!!
        mark_corners(a, ad, an);
        mark_corners(b, bd, bn);
        image inlier_matches = draw_inliers(a, b, H, m, mn, inlier_thresh);
        save_image(inlier_matches, "inliers");
    }

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);
    // Stitch the images together with the homography
    image comb = combine_images(a, b, H);
    return comb;
}

// Project an image onto a cylinder.
// image im: image to project.
// float f: focal length used to take image (in pixels).
// returns: image projected onto cylinder, then flattened.
image cylindrical_project(image im, float f)
{
    // TODO: project image onto a cylinder
    float theta;
    float h;
    float X;
    float Y;
    float Z;

    image im_copy = copy_image(im);
    for (int c = 0; c < im.c; ++c)
    {
        for (int h = 0; h < im.h; ++h)
        {
            for (int w = 0; w < im.w; ++w)
            {
                float theta = (w - w); // !!! Complete later !!! xc is 255
            }
        }
    }
    return im_copy;
}
