#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Rounding Functions
 */

double C441(double f) {
    return ceil(f - 0.00001);
}

double F441(double f) {
    return floor(f + 0.00001);
}

/*
 * Rasterization Classes
 */

typedef struct Triangle {
   double X[3];
   double Y[3];
   double Z[3];
   double color[3][3]; // color[2][0] is for V2, red channel
} Triangle;

typedef struct TriangleList {
   int numTriangles;
   Triangle *triangles;
} TriangleList;

double get_x_intercept(double ypos, double x1, double y1, double x2, double y2) {
    /*
     * Given a line between two points, calculate the X intercept
     * that crosses the line at a given y position.
     */
    if (x2 == x1) {
        // This line has no slope, x intercept is clear.
        return x1;
    }

    double slope = (y2 - y1) / (x2 - x1);
    double b = y1 - (slope * x1);
    return (ypos - b) / (slope);
}

double lerp(double a, double b, double t) {
    /*
     * Lerps a value between a and b.
     */
    return (t * a) + ((1 - t) * b);
}

typedef struct FragmentLerpData {
    // Data struct for holding fragment lerp data.
    double z[2];
    double color[3][2];
} FragmentLerpData;

struct FragmentLerpData *get_fragment_lerp(struct Triangle *tri, int v1, int v2, int v3,
                                           double l_xpos, double r_xpos, double ypos) {
    /*
     * Given a triangle, ordered vertices, left/right intercept
     * positions and a target ypos,
     * return the useful left/right data we get from lerping those vertices.
     */
    // First: We need to calculate the lerp amounts.
    // We'll do this by calculating a general theta
    // for us to lerp our values by.
    // We need to get some useful constants first.
    double x1 = tri->X[v1];
    double y1 = tri->Y[v1];
    double x2 = tri->X[v2];
    double y2 = tri->Y[v2];
    double x3 = tri->X[v3];
    double y3 = tri->Y[v3];

    // First, calculate the distance between the endpoints.
    double l_distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    double r_distance = sqrt(pow(x3 - x2, 2) + pow(y3 - y2, 2));

    // Then calculate the distance between one point and the checked point.
    double pl_distance = sqrt(pow(x2 - l_xpos, 2) + pow(y2 - ypos, 2));
    double pr_distance = sqrt(pow(x3 - r_xpos, 2) + pow(y3 - ypos, 2));

    // Calculate the 'theta' to lerp by.
    // But sanity check first just in case...
    if (l_distance == 0 || r_distance == 0) {
        // No distance, avoid dividing by zero.
        printf("get_fragment_lerp received bad distance\n");
        abort();
    }

    double l_theta = pl_distance / l_distance;
    double r_theta = pr_distance / r_distance;

    // Now that we've got the theta,
    // build the fragment lerp data to use.
    struct FragmentLerpData *fld = (struct FragmentLerpData *)malloc(sizeof(struct FragmentLerpData));

    // And populate it.
    fld->z[0] = lerp(tri->Z[v1], tri->Z[v2], l_theta);
    fld->z[1] = lerp(tri->Z[v2], tri->Z[v3], r_theta);
    for (int i = 0; i < 3; i++) {
        fld->color[i][0] = lerp(tri->color[v1][i], tri->color[v2][i], l_theta);
        fld->color[i][1] = lerp(tri->color[v2][i], tri->color[v3][i], r_theta);
    }

    // Return the FLD.
    // They are responsible for freeing it.
    return fld;
}

/*
 * Project 4 specs
 */

char *
ReadTuple3(char *tmp, double *v1, double *v2, double *v3)
{
    tmp++; /* left paren */
    *v1 = atof(tmp);
    while (*tmp != ',')
       tmp++;
    tmp += 2; // comma+space
    *v2 = atof(tmp);
    while (*tmp != ',')
       tmp++;
    tmp += 2; // comma+space
    *v3 = atof(tmp);
    while (*tmp != ')')
       tmp++;
    tmp++; /* right paren */
    return tmp;
}

TriangleList *get_triangles() {
    /*
     * Reads triangles (lovely, lovely triangles).
     */
     FILE *f = fopen("tris_w_r_rgb.txt", "r");
     if (f == NULL) {
        fprintf(stderr, "You must place the tris_w_r_rgb.txt file in the current directory.\n");
        exit(EXIT_FAILURE);
    }
    fseek(f, 0, SEEK_END);
    int numBytes = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (numBytes != 13488634) {
        fprintf(stderr, "Your tris_w_r_rgb.txt file is corrupted.  It should be 13488634 bytes, but you have %d.\n", numBytes);
        exit(EXIT_FAILURE);
    }

    char *buffer = (char *) malloc(numBytes);
    if (buffer == NULL) {
        fprintf(stderr, "Unable to allocate enough memory to load file.\n");
        exit(EXIT_FAILURE);
    }

    fread(buffer, sizeof(char), numBytes, f);

    char *tmp = buffer;
    int numTriangles = atoi(tmp);
    while (*tmp != '\n')
        tmp++;
    tmp++;

    if (numTriangles != 42281) {
        fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
        exit(EXIT_FAILURE);
    }

    TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
    tl->numTriangles = numTriangles;
    tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

    for (int i = 0 ; i < tl->numTriangles ; i++) {
        double x1, y1, z1, x2, y2, z2, x3, y3, z3;
        double r[3], g[3], b[3];

        tmp = ReadTuple3(tmp, &x1, &y1, &z1);
        tmp += 3; /* space+equal+space */
        tmp = ReadTuple3(tmp, r+0, g+0, b+0);
        tmp += 2; /* comma+space */
        tmp = ReadTuple3(tmp, &x2, &y2, &z2);
        tmp += 3; /* space+equal+space */
        tmp = ReadTuple3(tmp, r+1, g+1, b+1);
        tmp += 2; /* comma+space */
        tmp = ReadTuple3(tmp, &x3, &y3, &z3);
        tmp += 3; /* space+equal+space */
        tmp = ReadTuple3(tmp, r+2, g+2, b+2);
        tmp++;    /* newline */

        tl->triangles[i].X[0] = x1;
        tl->triangles[i].X[1] = x2;
        tl->triangles[i].X[2] = x3;
        tl->triangles[i].Y[0] = y1;
        tl->triangles[i].Y[1] = y2;
        tl->triangles[i].Y[2] = y3;
        tl->triangles[i].Z[0] = z1;
        tl->triangles[i].Z[1] = z2;
        tl->triangles[i].Z[2] = z3;
        tl->triangles[i].color[0][0] = r[0];
        tl->triangles[i].color[0][1] = g[0];
        tl->triangles[i].color[0][2] = b[0];
        tl->triangles[i].color[1][0] = r[1];
        tl->triangles[i].color[1][1] = g[1];
        tl->triangles[i].color[1][2] = b[1];
        tl->triangles[i].color[2][0] = r[2];
        tl->triangles[i].color[2][1] = g[2];
        tl->triangles[i].color[2][2] = b[2];
        //printf("Read triangle (%f, %f, %f) / (%f, %f, %f), (%f, %f, %f) / (%f, %f, %f), (%f, %f, %f) / (%f, %f, %f)\n", x1, y1, z1, r[0], g[0], b[0], x2, y2, z2, r[1], g[1], b[1], x3, y3, z3, r[2], g[2], b[2]);
    }

    free(buffer);
    return tl;
}

/*
 * Pixel Struct
 */

typedef struct Pixel {
    // Stores color data.
    unsigned char r;
    unsigned char g;
    unsigned char b;
    double z;
} Pixel;

Pixel *make_pixel(unsigned char r, unsigned char g, unsigned char b) {
    /*
     * Creates a pixel struct of a given RGB.
     */
    Pixel *pixel = (Pixel *)malloc(sizeof(Pixel));
    pixel->r = r;
    pixel->g = g;
    pixel->b = b;
    pixel->z = -1;
    return pixel;
}

void cleanup_pixel(Pixel *pixel) {
    /*
     * Cleans up the pixel struct.
     */
    free(pixel);
}

/*
 * Image Struct
 */

typedef struct Image {
    unsigned int width;
    unsigned int height;
    Pixel ***data;  // i'm seeing stars...
} Image;

Image *make_image(unsigned int width, unsigned int height) {
    /*
     * Creates a blank canvas.
     */
    Image *image = (Image *)malloc(sizeof(Image));
    image->width = width;
    image->height = height;

    // Make a pixel array for the columns.
    Pixel ***columns = (Pixel ***)malloc(width * sizeof(Pixel **));
    image->data = columns;

    // Populate the columns array.
    for (unsigned int i = 0; i < width; i++) {
        // Make a pixel array for the content of this column.
        Pixel **column = (Pixel **)malloc(height * sizeof(Pixel *));
        columns[i] = column;

        // Populate this column with pixels.
        for (unsigned int j = 0; j < height; j++) {
            Pixel *pixel = make_pixel(0, 0, 0);
            column[j] = pixel;
        }
    }

    // We are very done.
    return image;
}

void cleanup_image(Image *image) {
    /*
     * Cleans up an image struct.
     */
    // Iterate over the columns array.
    for (unsigned int i = 0; i < image->width; i++) {
        Pixel **column = image->data[i];

        // Cleanup each pixel in the column.
        for (unsigned int j = 0; j < image->height; j++)
            cleanup_pixel(column[j]);

        // Cleanup the column.
        free(column);
    }

    // Cleanup the rest of the image.
    free(image->data);
    free(image);
}

void write_image(Image *image, char *fname) {
    /*
     * Writes an image out to a file.
     */
    FILE *fp;
    fp = fopen(fname, "w+");

    // Print the file header.
    fprintf(fp, "P6\n%d %d\n255\n", image->width, image->height);

    // Iterate over each pixel, starting at the top left.
    for (unsigned int j = 0; j < image->height; j++) {
        for (unsigned int i = 0; i < image->width; i++) {
            // Get the pixel to print.
            // (don't ask why I had to do this weird math op)
            Pixel *pixel = (image->data)[i][(image->height) - j - 1];

            // Little memory optimization: just write the struct directly
            // (the rgb values are adjacent to eachother)
            fwrite(pixel, 1, 3, fp);
        }
    }

    // Close the file.
    fclose(fp);
}

/*
 * Draw Functions
 */

void draw_triangle(Image *image, Triangle *tri) {
    /*
     * Draws a triangle.
     */

    // Figure out the top and bottom of the triangle.
    double top    = fmax(fmax((tri->Y)[0], (tri->Y)[1]), (tri->Y)[2]);
    double bottom = fmin(fmin((tri->Y)[0], (tri->Y)[1]), (tri->Y)[2]);

    // Check for the vertices of the bottom and top triangles.
    int v_bot = (tri->Y[0] == bottom)            ? (0) : ((tri->Y[1] == bottom)            ? (1) : (2));
    int v_top = (tri->Y[0] == top && v_bot != 0) ? (0) : ((tri->Y[1] == top && v_bot != 1) ? (1) : (2));

    // The middle vertex will be whichever one the others aren't.
    int v_mid = (v_bot != 0 && v_top != 0) ? (0) : ((v_bot != 1 && v_top != 1) ? (1) : (2));

    // Good opportunity to get the y coordinate of the middle coordinate as well.
    double middle = (tri->Y)[v_mid];

    // OK, now we iterate across each scanline.
    for (double y = F441(bottom); y <= C441(top); y++) {
        // Make sure this scanline is onscreen.
        if (y < 0) continue;
        if (y >= (double)(image->height)) continue;

        // It's time to figure out our intercepts.
        double left = 0;
        double right = 0;

        double c_y = fmax(bottom, fmin(y, top)); // constrained y

        struct FragmentLerpData *fld = NULL;

        // Depending on where our Y coordinate is relative to the
        // middle vertex, we will have to calculate the intercepts differently.
        if (y < middle) {
            // We are looking below the middle coordinate.
            // Calculate two intercepts.
            double intercept_a = get_x_intercept(c_y,
                                                 tri->X[v_mid], tri->Y[v_mid],
                                                 tri->X[v_bot], tri->Y[v_bot]);
            double intercept_b = get_x_intercept(c_y,
                                                 tri->X[v_bot], tri->Y[v_bot],
                                                 tri->X[v_top], tri->Y[v_top]);

            // Determine left and right relative to the intercepts.
            left  = fmin(intercept_a, intercept_b);
            right = fmax(intercept_a, intercept_b);

            if (left == intercept_a)
                fld = get_fragment_lerp(tri, v_mid, v_bot, v_top, left, right, c_y);
            else
                fld = get_fragment_lerp(tri, v_top, v_bot, v_mid, left, right, c_y);
        } else if (y > middle) {
            // We are looking above the middle coordinate.
            // Calculate two intercepts.
            double intercept_a = get_x_intercept(c_y,
                                                 tri->X[v_top], tri->Y[v_top],
                                                 tri->X[v_bot], tri->Y[v_bot]);
            double intercept_b = get_x_intercept(c_y,
                                                 tri->X[v_top], tri->Y[v_top],
                                                 tri->X[v_mid], tri->Y[v_mid]);

            // Determine left and right relative to the intercepts.
            left  = fmin(intercept_a, intercept_b);
            right = fmax(intercept_a, intercept_b);

            if (left == intercept_a)
                fld = get_fragment_lerp(tri, v_bot, v_top, v_mid, left, right, c_y);
            else
                fld = get_fragment_lerp(tri, v_mid, v_top, v_bot, left, right, c_y);
        } else {
            // Our intercept lies right on the middle coordinate.
            // Does this triangle have a flat base?
            if (middle == bottom) {
                // The triangle has a flat base.
                // Our intercepts will be across the base of the triangle.
                left  = fmin(tri->X[v_mid], tri->X[v_bot]);
                right = fmax(tri->X[v_mid], tri->X[v_bot]);
                if (left == tri->X[v_mid])
                    fld = get_fragment_lerp(tri, v_top, v_mid, v_bot, left, right, c_y);
                else
                    fld = get_fragment_lerp(tri, v_top, v_bot, v_mid, left, right, c_y);
            } else if (middle == top) {
                // The triangle has a flat base at the top.
                // Our intercepts will be across the base of the triangle.
                left  = fmin(tri->X[v_mid], tri->X[v_top]);
                right = fmax(tri->X[v_mid], tri->X[v_top]);
                if (left == tri->X[v_mid])
                    fld = get_fragment_lerp(tri, v_mid, v_top, v_bot, left, right, c_y);
                else
                    fld = get_fragment_lerp(tri, v_top, v_mid, v_bot, left, right, c_y);
            } else {
                // The triangle does not have a flat base.
                // Calculate the intercept.
                double intercept = get_x_intercept(c_y,
                                                   tri->X[v_top], tri->Y[v_top],
                                                   tri->X[v_bot], tri->Y[v_bot]);

                // Determine left and right relative to this intercept.
                left  = fmin(tri->X[v_mid], intercept);
                right = fmax(tri->X[v_mid], intercept);

                if (left == intercept)
                    fld = get_fragment_lerp(tri, v_bot, v_top, v_mid, left, right, c_y);
                else
                    fld = get_fragment_lerp(tri, v_mid, v_top, v_bot, left, right, c_y);
            }
        }

        // if they're too close just ignore this point
        if (((right - left) < 0.00001) && (right - left) > -0.00001) continue;

        // we need valid lerp data
        if (fld == NULL) continue;

        // Draw the scanline now.
        for (double x = C441(left); x <= F441(right); x++) {
            // Make sure these coordinates are onscreen.
            if (x < 0) continue;
            if (x >= (double)(image->width)) continue;

            // Calculate the depth at this pixel.
            double theta = fmax(0, fmin((x - right) / (left - right), 1));
            double z = lerp(fld->z[0], fld->z[1], theta);

            // Get the pixel to write over.
            Pixel *pixel = (image->data)[(int)x][(int)y];

            // Are we deep enough?
            if (z <= pixel->z) continue;

            // Write over it.
            pixel->r = C441(lerp(fld->color[0][0], fld->color[0][1], theta) * 255);
            pixel->g = C441(lerp(fld->color[1][0], fld->color[1][1], theta) * 255);
            pixel->b = C441(lerp(fld->color[2][0], fld->color[2][1], theta) * 255);
            pixel->z = z;
        }

        // Cleanup.
        free(fld);
    }
}

/*
 * this is where the magic happens
 */

int main(int argc, char *argv[]) {
    /*
     * Create an image canvas and manipulate it with silly little colors.
     */
    Image *image = make_image(1000, 1000);

    // Draw over this image.
    TriangleList *tl = get_triangles();
    for (int i = 0; i < tl->numTriangles; i++)
        draw_triangle(image, &((tl->triangles)[i]));

    // Print this image.
    write_image(image, "proj1D_out.pnm");

    // Cleanup.
    cleanup_image(image);

    // We are done here.
    return EXIT_SUCCESS;
}
