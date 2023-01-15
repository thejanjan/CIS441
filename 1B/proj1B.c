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
   double x[3];
   double y[3];
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
 * Project 2 specs
 */

TriangleList *get_triangles(void) {
    /*
     * Generates a list of triangles to use.
     */
    TriangleList *tl = (TriangleList *)malloc(sizeof(TriangleList));
    tl->numTriangles = 100;
    tl->triangles = (Triangle *)malloc(sizeof(Triangle)*(tl->numTriangles));

    // Define the colors to populate the triangles with.
    unsigned char colors[6][3] = {
        {255,128,0}, {255, 0, 127}, {0,204,204},
        {76,153,0},  {255, 204, 204}, {204, 204, 0}
    };

    // A little math to define each triangle.
    for (int i = 0 ; i < (tl->numTriangles); i++) {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       tl->triangles[i].x[firstPt] = posI;
       if (i == 50)
           tl->triangles[i].x[firstPt] = -10;
       tl->triangles[i].y[firstPt] = posJ+10*(idxJ+1);
       tl->triangles[i].x[(firstPt+1)%3] = posI+105;
       tl->triangles[i].y[(firstPt+1)%3] = posJ;
       tl->triangles[i].x[(firstPt+2)%3] = posI+i;
       tl->triangles[i].y[(firstPt+2)%3] = posJ;
       if (i == 95)
          tl->triangles[i].y[firstPt] = 1050;
       tl->triangles[i].color[0] = colors[i%6][0];
       tl->triangles[i].color[1] = colors[i%6][1];
       tl->triangles[i].color[2] = colors[i%6][2];
   }

   // Return the list.
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
    double top    = fmax(fmax((tri->y)[0], (tri->y)[1]), (tri->y)[2]);
    double bottom = fmin(fmin((tri->y)[0], (tri->y)[1]), (tri->y)[2]);

    // Check for the vertices of the bottom and top triangles.
    int v_bot = (tri->y[0] == bottom) ? (0) : ((tri->y[1] == bottom) ? (1) : (2));
    int v_top = (tri->y[0] == top)    ? (0) : ((tri->y[1] == top)    ? (1) : (2));

    // SANITY CHECK: they better not be the same
    if (v_bot == v_top) {
        fprintf(stderr, "v_bot and v_top are equal\n");
        exit(EXIT_FAILURE);
    }

    // The middle vertex will be whichever one the others aren't.
    int v_mid = (v_bot != 0 && v_top != 0) ? (0) : ((v_bot != 1 && v_top != 1) ? (1) : (2));

    // Good opportunity to get the y coordinate of the middle coordinate as well.
    double middle = (tri->y)[v_mid];

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
            double intercept_a = get_x_intercept(y,
                                                 tri->x[v_mid], tri->y[v_mid],
                                                 tri->x[v_bot], tri->y[v_bot]);
            double intercept_b = get_x_intercept(y,
                                                 tri->x[v_mid], tri->y[v_mid],
                                                 tri->x[v_top], tri->y[v_top]);

            // Determine left and right relative to the intercepts.
            left  = fmin(intercept_a, intercept_b);
            right = fmax(intercept_a, intercept_b);
        } else if (y > middle) {
            // We are looking above the middle coordinate.
            // Calculate two intercepts.
            double intercept_a = get_x_intercept(y,
                                                 tri->x[v_top], tri->y[v_top],
                                                 tri->x[v_bot], tri->y[v_bot]);
            double intercept_b = get_x_intercept(y,
                                                 tri->x[v_top], tri->y[v_top],
                                                 tri->x[v_mid], tri->y[v_mid]);

            // Determine left and right relative to the intercepts.
            left  = fmin(intercept_a, intercept_b);
            right = fmax(intercept_a, intercept_b);
        } else {
            // Our intercept lies right on the middle coordinate.
            // Does this triangle have a flat base?
            if (middle == bottom) {
                // The triangle has a flat base.
                // Our intercepts will be across the base of the triangle.
                left  = fmin(tri->x[v_mid], tri->x[v_bot]);
                right = fmax(tri->x[v_mid], tri->x[v_bot]);
            } else {
                // The triangle does not have a flat base.
                // Calculate the intercept.
                double intercept = get_x_intercept(y,
                                                   tri->x[v_top], tri->y[v_top],
                                                   tri->x[v_bot], tri->y[v_bot]);

                // Determine left and right relative to this intercept.
                left  = fmin(tri->x[v_mid], intercept);
                right = fmax(tri->x[v_mid], intercept);
            }
        }

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
    Image *image = make_image(1000, 1000);

    // Draw over this image.
    TriangleList *tl = get_triangles();
    for (int i = 0; i < tl->numTriangles; i++)
        draw_triangle(image, &((tl->triangles)[i]));

    // Print this image.
    write_image(image, "proj1B_out.pnm");

    // Cleanup.
    cleanup_image(image);

    // We are done here.
    return EXIT_SUCCESS;
}
