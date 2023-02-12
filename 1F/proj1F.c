#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NORMALS

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

typedef struct {
    double X[3];
    double Y[3];
    double Z[3];
    double color[3][3]; // color[2][0] is for V2, red channel
#ifdef NORMALS
    double normals[3][3]; // normals[2][0] is for V2, x-component
    double shading[3];
#endif
} Triangle;

typedef struct {
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

typedef struct
{
    double          A[4][4];     // A[i][j] means row i, column j
} Matrix;


void
PrintMatrix(Matrix m)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        printf("(%.7f %.7f %.7f %.7f)\n", m.A[i][0], m.A[i][1], m.A[i][2], m.A[i][3]);
    }
}

Matrix *ComposeMatrices(Matrix M1, Matrix M2)
{
    Matrix *m_out = (Matrix *)malloc(sizeof (Matrix));
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            (m_out->A)[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                (m_out->A)[i][j] += M1.A[i][k]*M2.A[k][j];
        }
    return m_out;
}

void 
TransformPoint(Matrix m, const double *ptIn, double *ptOut)
{  
    ptOut[0] = ptIn[0]*m.A[0][0]
             + ptIn[1]*m.A[1][0]
             + ptIn[2]*m.A[2][0]
             + ptIn[3]*m.A[3][0];
    ptOut[1] = ptIn[0]*m.A[0][1]
             + ptIn[1]*m.A[1][1]
             + ptIn[2]*m.A[2][1]
             + ptIn[3]*m.A[3][1];
    ptOut[2] = ptIn[0]*m.A[0][2]
             + ptIn[1]*m.A[1][2]
             + ptIn[2]*m.A[2][2]
             + ptIn[3]*m.A[3][2];
    ptOut[3] = ptIn[0]*m.A[0][3]
             + ptIn[1]*m.A[1][3]
             + ptIn[2]*m.A[2][3]
             + ptIn[3]*m.A[3][3];
    ptOut[0] /= ptOut[3];
    ptOut[1] /= ptOut[3];
    ptOut[2] /= ptOut[3];
    ptOut[3] /= ptOut[3];
}

typedef struct
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
} Camera;

double SineParameterize(int curFrame, int nFrames, int ramp)
{  
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {        
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }        
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
} 

Camera       
GetCamera(int frame, int nframes)
{            
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0; 
    c.focus[1] = 0; 
    c.focus[2] = 0;
    c.up[0] = 0;    
    c.up[1] = 1;    
    c.up[2] = 0;    
    return c;       
}


typedef struct 
{
    double lightDir[3];  // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
} LightingParameters;


LightingParameters 
GetLighting(Camera c)
{
    LightingParameters lp;
    lp.Ka = 0.3;
    lp.Kd = 0.7;
    lp.Ks = 2.8;
    lp.alpha = 50.5;
    lp.lightDir[0] = c.position[0]-c.focus[0];
    lp.lightDir[1] = c.position[1]-c.focus[1];
    lp.lightDir[2] = c.position[2]-c.focus[2];
    double mag = sqrt(lp.lightDir[0]*lp.lightDir[0]
                    + lp.lightDir[1]*lp.lightDir[1]
                    + lp.lightDir[2]*lp.lightDir[2]);
    if (mag > 0)
    {
        lp.lightDir[0] /= mag;
        lp.lightDir[1] /= mag;
        lp.lightDir[2] /= mag;
    }

    return lp;
}

char *
Read3Numbers(char *tmp, double *v1, double *v2, double *v3)
{
    *v1 = atof(tmp);
    while (*tmp != ' ')
       tmp++;
    tmp++; /* space */
    *v2 = atof(tmp);
    while (*tmp != ' ')
       tmp++;
    tmp++; /* space */
    *v3 = atof(tmp);
    while (*tmp != ' ' && *tmp != '\n')
       tmp++;
    return tmp;
}

TriangleList *
Get3DTriangles()
{
   FILE *f = fopen("ws_tris.txt", "r");
   if (f == NULL)
   {
       fprintf(stderr, "You must place the ws_tris.txt file in the current directory.\n");
       exit(EXIT_FAILURE);
   }
   fseek(f, 0, SEEK_END);
   int numBytes = ftell(f);
   fseek(f, 0, SEEK_SET);
   if (numBytes != 3892295)
   {
       fprintf(stderr, "Your ws_tris.txt file is corrupted.  It should be 3892295 bytes, but you have %d.\n", numBytes);
       exit(EXIT_FAILURE);
   }

   char *buffer = (char *) malloc(numBytes);
   if (buffer == NULL)
   {
       fprintf(stderr, "Unable to allocate enough memory to load file.\n");
       exit(EXIT_FAILURE);
   }
   
   fread(buffer, sizeof(char), numBytes, f);

   char *tmp = buffer;
   int numTriangles = atoi(tmp);
   while (*tmp != '\n')
       tmp++;
   tmp++;
 
   if (numTriangles != 14702)
   {
       fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
       exit(EXIT_FAILURE);
   }

   TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
   tl->numTriangles = numTriangles;
   tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

   for (int i = 0 ; i < tl->numTriangles ; i++)
   {
       for (int j = 0 ; j < 3 ; j++)
       {
           double x, y, z;
           double r, g, b;
           double normals[3];
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
           tmp = Read3Numbers(tmp, &x, &y, &z);
           tmp += 3; /* space+slash+space */
           tmp = Read3Numbers(tmp, &r, &g, &b);
           tmp += 3; /* space+slash+space */
           tmp = Read3Numbers(tmp, normals+0, normals+1, normals+2);
           tmp++;    /* newline */

           tl->triangles[i].X[j] = x;
           tl->triangles[i].Y[j] = y;
           tl->triangles[i].Z[j] = z;
           tl->triangles[i].color[j][0] = r;
           tl->triangles[i].color[j][1] = g;
           tl->triangles[i].color[j][2] = b;
#ifdef NORMALS
           tl->triangles[i].normals[j][0] = normals[0];
           tl->triangles[i].normals[j][1] = normals[1];
           tl->triangles[i].normals[j][2] = normals[2];
#endif
       }
   }

   free(buffer);
   return tl;
}

/*
 * Math Ops
 */

void cross(double *a, double *b, double *out) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

double dot(double *a, double *b) {
    return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

void norm(double *a) {
    double size = sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
    a[0] /= size;
    a[1] /= size;
    a[2] /= size;
}

/*
 * Triangle Transformations
 */

Triangle *transform_triangle(Matrix M, Triangle *ref) {
    /*
     * Applies a given transformation onto a triangle,
     * and returns a new triangle.
     */
    // Make a copy of the triangle.
    Triangle *tri = (Triangle *)malloc(sizeof(Triangle));
    memcpy(tri, ref, sizeof(Triangle));

    // Apply matrix transformations.
    for (int i = 0; i < 3; i++) {
        double ptIn[4] = {
            ref->X[i], ref->Y[i], ref->Z[i], 1
        };
        double ptOut[4];
        TransformPoint(M, ptIn, ptOut);
        tri->X[i] = ptOut[0];
        tri->Y[i] = ptOut[1];
        tri->Z[i] = ptOut[2];
    }

    // Return it.
    return tri;
}

/*
 * Fragment Interpolation 
 */

double calculate_phong_shading(Triangle *tri, LightingParameters lp, Camera c, int vertex) {
    /*
    typedef struct {
        double lightDir[3];  // The direction of the light source
        double Ka;           // The coefficient for ambient lighting.
        double Kd;           // The coefficient for diffuse lighting.
        double Ks;           // The coefficient for specular lighting.
        double alpha;        // The exponent term for specular lighting.
    } LightingParameters;

    typedef struct {
        double          near, far;
        double          angle;
        double          position[3];
        double          focus[3];
        double          up[3];
    } Camera;
    */
    // Get the normal.
    double *normal = tri->normals[vertex];

    // Calculate the diffuse amount.
    double diffuse = fmax(0, dot(lp.lightDir, normal));

    // Calculate R and V vectors for specular.
    double v[3] = {
        c.position[0] - tri->X[vertex],
        c.position[1] - tri->Y[vertex],
        c.position[2] - tri->Z[vertex]
    };
    norm(v);

    // R will be 2*(L dot N)*N - L.
    double r[3] = {0, 0, 0};
    double r_mult = 2 * dot(normal, lp.lightDir);
    r[0] = (normal[0] * r_mult) - (lp.lightDir[0]);
    r[1] = (normal[1] * r_mult) - (lp.lightDir[1]);
    r[2] = (normal[2] * r_mult) - (lp.lightDir[2]);
    norm(r);

    // Now calculate specular.
    double specular = 1 * pow(fmax(0, dot(r, v)), lp.alpha);
    if (specular < 0) specular = -specular;

    return lp.Ka + (diffuse * lp.Kd) + (specular * lp.Ks);
}

void apply_shading_to_triangle(Triangle *tri, LightingParameters lp, Camera c) {
    tri->shading[0] = calculate_phong_shading(tri, lp, c, 0);
    tri->shading[1] = calculate_phong_shading(tri, lp, c, 1);
    tri->shading[2] = calculate_phong_shading(tri, lp, c, 2);
}

typedef struct FragmentLerpData {
    // Data struct for holding fragment lerp data.
    double z[2];
    double color[3][2];
    double shading[2];
} FragmentLerpData;

struct FragmentLerpData *get_fragment_lerp(Triangle *tri, int v1, int v2, int v3,
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
    fld->shading[0] = lerp(tri->shading[v1], tri->shading[v2], l_theta);
    fld->shading[1] = lerp(tri->shading[v2], tri->shading[v3], r_theta);

    // Return the FLD.
    // They are responsible for freeing it.
    return fld;
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

void initialize_screen(Image *image) {
    /*
     * Flushes the screen buffer.
     */
    for (unsigned int i = 0; i < image->width; i++) {
        Pixel **column = image->data[i];
        for (unsigned int j = 0; j < image->height; j++) {
            Pixel *pixel = column[j];

            pixel->r = 0;
            pixel->g = 0;
            pixel->b = 0;
            pixel->z = -1;
        }
    }
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
            double s = lerp(fld->shading[0], fld->shading[1], theta);;

            // Get the pixel to write over.
            Pixel *pixel = (image->data)[(int)x][(int)y];

            // Are we deep enough?
            if (z <= pixel->z) continue;

            // Write over it.
            pixel->r = fmin(255, C441(lerp(fld->color[0][0], fld->color[0][1], theta) * 255 * s));
            pixel->g = fmin(255, C441(lerp(fld->color[1][0], fld->color[1][1], theta) * 255 * s));
            pixel->b = fmin(255, C441(lerp(fld->color[2][0], fld->color[2][1], theta) * 255 * s));
            pixel->z = z;
        }

        // Cleanup.
        free(fld);
    }
}

void transform_render_triangles(Camera c, TriangleList *tl, Image *img) {
    // Valculate camera transform.
    double cf_o[3] = {
        c.position[0], c.position[1], c.position[2]
    };

    double cf_w[3] = {
        cf_o[0] - c.focus[0], cf_o[1] - c.focus[1], cf_o[2] - c.focus[2]
    };
    norm(cf_w);

    double cf_u[3] = {0, 0, 0};
    cross(c.up, cf_w, cf_u);
    norm(cf_u);

    double cf_v[3] = {0, 0, 0};
    cross(cf_w, cf_u, cf_v);
    norm(cf_v);

    double cf_t[3] = {
        -cf_o[0], -cf_o[1], -cf_o[2]
    };

    double cf_z[3] = {
        dot(cf_u, cf_t), dot(cf_v, cf_t), dot(cf_w, cf_t)
    };

    Matrix camera_transform = {{
        {cf_u[0], cf_v[0], cf_w[0], 0},
        {cf_u[1], cf_v[1], cf_w[1], 0},
        {cf_u[2], cf_v[2], cf_w[2], 0},
        {cf_z[0], cf_z[1], cf_z[2], 1}
    }};

    // Calculate view transform.
    double vt_cot = 1 / tan(c.angle / 2);
    double far = c.far;
    double near = c.near;

    Matrix view_transform = {{
        {vt_cot, 0, 0, 0},
        {0, vt_cot, 0, 0},
        {0, 0, (far + near) / (far - near), -1},
        {0, 0, (2 * far * near) / (far - near), 0}
    }};

    // calculate device transform
    int n = img->width;
    int m = img->height;
    Matrix device_transform = {{
        {n / 2, 0, 0, 0},
        {0, m / 2, 0, 0},
        {0, 0,     1, 0},
        {n/2, m/2, 0, 1}
    }};

    // Create composed matrix M.
    Matrix *CameraViewMatrix = ComposeMatrices(camera_transform, view_transform);
    Matrix *M = ComposeMatrices(*CameraViewMatrix, device_transform);
    free(CameraViewMatrix);

    #ifdef DEBUG_TRT
        printf("Camera Frame: U = %f, %f, %f\n", cf_u[0], cf_u[1], cf_u[2]);
        printf("Camera Frame: V = %f, %f, %f\n", cf_v[0], cf_v[1], cf_v[2]);
        printf("Camera Frame: W = %f, %f, %f\n", cf_w[0], cf_w[1], cf_w[2]);
        printf("Camera Frame: O = %f, %f, %f\n", cf_o[0], cf_o[1], cf_o[2]);
        printf("\nCamera Transform\n");
        PrintMatrix(camera_transform);
        printf("\nView Transform\n");
        PrintMatrix(view_transform);
        printf("\nTotal Transform\n");
        PrintMatrix(*M);

        printf("\n\nMapping (0, 36.46, 36.46)\n");
        double ptIn[4] = {0, 36.464466, 36.464466, 1};
        double ptOut[4];
        TransformPoint(*M, ptIn, ptOut);
        printf("Result: %f, %f, %f, %f\n", ptOut[0], ptOut[1], ptOut[2], ptOut[3]);
    #endif

    // Apply lighting onto all of the triangles.
    LightingParameters lp = GetLighting(c);

    // process each triangle over M and render it
    for (int i = 0; i < tl->numTriangles; i++) {
        // Get the reference triangle.
        Triangle *ref = &((tl->triangles)[i]);

        // Figure out the lighting for the triangle.
        apply_shading_to_triangle(ref, lp, c);

        // Get its tranformation, draw it and clean it up.
        Triangle *tri = transform_triangle(*M, ref);
        draw_triangle(img, tri);
        free(tri);
    }

    free(M);
}

/*
 * this is where the magic happens
 */

int main(int argc, char *argv[]) {
    /*
     * Create an image canvas and manipulate it with silly little colors.
     */

    Image *image = make_image(1000, 1000);
    TriangleList *tl = Get3DTriangles();

    for (int i = 0 ; i < 1000 ; i++) {
        initialize_screen(image);   
        Camera c = GetCamera(i, 1000);
        transform_render_triangles(c, tl, image);

        // Write to file.
        char fname[21];
        sprintf(fname, "proj1F_frame%04d.pnm", i);
        write_image(image, fname);

        // disable for movie mode!
        break;
    }

    // Cleanup.
    cleanup_image(image);
    free(tl->triangles);
    free(tl);

    // We are done here.
    return EXIT_SUCCESS;
}
