#include <stdio.h>
#include <stdlib.h>

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

void draw_square(Image *image,
                 unsigned char r, unsigned char g, unsigned char b,
                 unsigned int xstart, unsigned int ystart,
                 unsigned int width, unsigned int height) {
    /*
     * Draws a square at given coordinates.
     * Uses a pixel struct for color reference.
     */
    for (unsigned int i = xstart; i < (xstart + width); i++) {
        for (unsigned int j = ystart; j < (ystart +height); j++) {
            // Get the pixel to write over.
            Pixel *pixel = (image->data)[i][j];

            // Write over it.
            pixel->r = r;
            pixel->g = g;
            pixel->b = b;
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
    Image *image = make_image(300, 300);

    // Draw over this image.
    draw_square(image,   0,   0,   0,   0, 200, 100, 100);
    draw_square(image, 128, 128, 128, 100, 200, 100, 100);
    draw_square(image, 255, 255, 255, 200, 200, 100, 100);

    draw_square(image, 255,   0,   0,   0, 100, 100, 100);
    draw_square(image,   0, 255,   0, 100, 100, 100, 100);
    draw_square(image,   0,   0, 255, 200, 100, 100, 100);

    draw_square(image, 255,   0, 255,   0,   0, 100, 100);
    draw_square(image,   0, 255, 255, 100,   0, 100, 100);
    draw_square(image, 255, 255,   0, 200,   0, 100, 100);

    // Print this image. (Sorry ISO C++)
    write_image(image, "proj1A_out.pnm");

    // Cleanup.
    cleanup_image(image);

    // We are done here.
    return EXIT_SUCCESS;
}
