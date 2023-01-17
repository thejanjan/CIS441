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
   unsigned char color[3];
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

/*
 * Project 3 specs
 */

TriangleList *get_triangles(int small_read) {
    /*
     * Reads triangles (lovely, lovely triangles).
     */
    FILE *f = fopen("tris.txt", "r");
    if (f == NULL) {
        fprintf(stderr, "You must place the tris.txt file in the current directory.\n");
        exit(EXIT_FAILURE);
    }
    fseek(f, 0, SEEK_END);
    int numBytes = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (numBytes != 241511792) {
        fprintf(stderr, "Your tris.txt file is corrupted.  It should be 241511792 bytes, but you only have %d.\n", numBytes);
        exit(EXIT_FAILURE);
    }

    if (small_read == 1) {
        numBytes = 10000;
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

    if (numTriangles != 2566541) {
        fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
        exit(EXIT_FAILURE);
    }

    if (small_read == 1)
        numTriangles = 100;

    TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
    tl->numTriangles = numTriangles;
    tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

    for (int i = 0 ; i < tl->numTriangles ; i++) {
        double x1, y1, x2, y2, x3, y3;
        int    r, g, b;
        /*
         * Weird: sscanf has a terrible implementation for large strings.
         * When I did the code below, it did not finish after 45 minutes.
         * Reading up on the topic, it sounds like it is a known issue that
         * sscanf fails here.  Stunningly, fscanf would have been faster.
         *     sscanf(tmp, "(%lf, %lf), (%lf, %lf), (%lf, %lf) = (%d, %d, %d)\n%n",
         *              &x1, &y1, &x2, &y2, &x3, &y3, &r, &g, &b, &numRead);
         *
         *  So, instead, do it all with atof/atoi and advancing through the buffer manually...
         */
        tmp++,
        x1 = atof(tmp);
        while (*tmp != ',')
            tmp++;
        tmp += 2; // comma+space
        y1 = atof(tmp);
        while (*tmp != ')')
            tmp++;
        tmp += 4; // right-paren+comma+space+left-paren
        x2 = atof(tmp);
        while (*tmp != ',')
            tmp++;
        tmp += 2; // comma+space
        y2 = atof(tmp);
        while (*tmp != ')')
            tmp++;
        tmp += 4; // right-paren+comma+space+left-paren
        x3 = atof(tmp);
        while (*tmp != ',')
            tmp++;
        tmp += 2; // comma+space
        y3 = atof(tmp);
        while (*tmp != ')')
            tmp++;
        tmp += 5; // right-paren+space+equal+space+left-paren
        r = atoi(tmp);
        while (*tmp != ',')
            tmp++;
        tmp += 2; // comma+space
        g = atoi(tmp);
        while (*tmp != ',')
            tmp++;
        tmp += 2; // comma+space
        b = atoi(tmp);
        while (*tmp != '\n')
            tmp++;
        tmp++; // onto next line

        tl->triangles[i].X[0] = x1;
        tl->triangles[i].X[1] = x2;
        tl->triangles[i].X[2] = x3;
        tl->triangles[i].Y[0] = y1;
        tl->triangles[i].Y[1] = y2;
        tl->triangles[i].Y[2] = y3;
        tl->triangles[i].color[0] = r;
        tl->triangles[i].color[1] = g;
        tl->triangles[i].color[2] = b;
        //printf("Read triangle %f, %f, %f, %f, %f, %f, %d, %d, %d\n", x1, y1, x2, y2, x3, y3, r, g, b);
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
} Pixel;

Pixel *make_pixel(unsigned char r, unsigned char g, unsigned char b) {
    /*
     * Creates a pixel struct of a given RGB.
     */
    Pixel *pixel = (Pixel *)malloc(sizeof(Pixel));
    pixel->r = r;
    pixel->g = g;
    pixel->b = b;
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

        // Depending on where our Y coordinate is relative to the
        // middle vertex, we will have to calculate the intercepts differently.
        if (y < middle) {
            // We are looking below the middle coordinate.
            // Calculate two intercepts.
            double intercept_a = get_x_intercept(fmax(y, bottom),
                                                 tri->X[v_mid], tri->Y[v_mid],
                                                 tri->X[v_bot], tri->Y[v_bot]);
            double intercept_b = get_x_intercept(fmax(y, bottom),
                                                 tri->X[v_bot], tri->Y[v_bot],
                                                 tri->X[v_top], tri->Y[v_top]);

            // Determine left and right relative to the intercepts.
            left  = fmin(intercept_a, intercept_b);
            right = fmax(intercept_a, intercept_b);
        } else if (y > middle) {
            // We are looking above the middle coordinate.
            // Calculate two intercepts.
            double intercept_a = get_x_intercept(fmin(y, top),
                                                 tri->X[v_top], tri->Y[v_top],
                                                 tri->X[v_bot], tri->Y[v_bot]);
            double intercept_b = get_x_intercept(fmin(y, top),
                                                 tri->X[v_top], tri->Y[v_top],
                                                 tri->X[v_mid], tri->Y[v_mid]);

            // Determine left and right relative to the intercepts.
            left  = fmin(intercept_a, intercept_b);
            right = fmax(intercept_a, intercept_b);
        } else {
            // Our intercept lies right on the middle coordinate.
            // Does this triangle have a flat base?
            if (middle == bottom) {
                // The triangle has a flat base.
                // Our intercepts will be across the base of the triangle.
                left  = fmin(tri->X[v_mid], tri->X[v_bot]);
                right = fmax(tri->X[v_mid], tri->X[v_bot]);
            } else if (middle == top) {
                // The triangle has a flat base at the top.
                // Our intercepts will be across the base of the triangle.
                left  = fmin(tri->X[v_mid], tri->X[v_top]);
                right = fmax(tri->X[v_mid], tri->X[v_top]);
            } else {
                // The triangle does not have a flat base.
                // Calculate the intercept.
                double intercept = get_x_intercept(y,
                                                   tri->X[v_top], tri->Y[v_top],
                                                   tri->X[v_bot], tri->Y[v_bot]);

                // Determine left and right relative to this intercept.
                left  = fmin(tri->X[v_mid], intercept);
                right = fmax(tri->X[v_mid], intercept);
            }
        }

        // if they're too close just ignore this point
        if (((right - left) < 0.00001) && (right - left) > -0.00001) continue;

        // Draw the scanline now.
        for (double x = C441(left); x <= F441(right); x++) {
            // Make sure these coordinates are onscreen.
            if (x < 0) continue;
            if (x >= (double)(image->width)) continue;

            // Get the pixel to write over.
            Pixel *pixel = (image->data)[(int)x][(int)y];

            // Write over it.
            pixel->r = tri->color[0];
            pixel->g = tri->color[1];
            pixel->b = tri->color[2];
        }
    }
}

/*
 * this is where the magic happens
 */

int main(int argc, char *argv[]) {
    /*
     * Create an image canvas and manipulate it with silly little colors.
     */
    Image *image = make_image(1786, 1344);

    // Draw over this image.
    TriangleList *tl = get_triangles(0);
    for (int i = 0; i < tl->numTriangles; i++)
        draw_triangle(image, &((tl->triangles)[i]));

    // Print this image.
    write_image(image, "proj1C_out.pnm");

    // Cleanup.
    cleanup_image(image);

    // We are done here.
    return EXIT_SUCCESS;
}
