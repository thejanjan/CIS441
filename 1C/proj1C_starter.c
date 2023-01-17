
typedef struct
{
   double         X[3];
   double         Y[3];
   unsigned char color[3];
} Triangle;

typedef struct
{
   int numTriangles;
   Triangle *triangles;
} TriangleList;


TriangleList *
GetTriangles(int small_read)
{
   FILE *f = fopen("tris.txt", "r");
   if (f == NULL)
   {
       fprintf(stderr, "You must place the tris.txt file in the current directory.\n");
       exit(EXIT_FAILURE);
   }
   fseek(f, 0, SEEK_END);
   int numBytes = ftell(f);
   fseek(f, 0, SEEK_SET);
   if (numBytes != 241511792)
   {
       fprintf(stderr, "Your tris.txt file is corrupted.  It should be 241511792 bytes, but you only have %d.\n", numBytes);
       exit(EXIT_FAILURE);
   }

   if (small_read == 1)
   {
       numBytes = 10000;
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
 
   if (numTriangles != 2566541)
   {
       fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
       exit(EXIT_FAILURE);
   }

   if (small_read == 1)
       numTriangles = 100;

   TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
   tl->numTriangles = numTriangles;
   tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

   for (int i = 0 ; i < tl->numTriangles ; i++)
   {
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

    // code to read all of the triangles: 
    //    TriangleList *tl = GetTriangles(0);

    // code to read just the first 100 triangles: 
    //    TriangleList *tl = GetTriangles(0);

