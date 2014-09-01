/*
   Copyright (C) 2014 Richard Plangger <rich [at] pasra [dot] at>

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include "normalmap.h"


NormalmapVals nmapvals =
{
   .filter = FILTER_NONE,
   .minz = 0.0,
   .scale = 2.0,
   .wrap = 1,
   .height_source = 0,
   .alpha = ALPHA_NONE,
   .conversion = CONVERT_NONE,
   .dudv = DUDV_NONE,
   .xinvert = 0,
   .yinvert = 0,
   .swapRGB = 0,
   .contrast = 0.0,
   .alphamap_id = 0
};

static char * progname;
void usage() {
    (void) fprintf(stdout, "usage: %s [options] <image> <normal_map>\n"
""
"  Options:\n"
"    -s scale. default 1"
"    -f filter. default FILTER_NONE\n"
"      values: " EACH_FILTER_TYPE(STR_EACH_SEP_WS) "\n"
"",
progname); 
}

int main(int argc,char **argv) {
    int bflag, ch, fd;
    char * endptr;
    progname = argv[0];
    char found = 0;

    bflag = 0;
    // the following is not buffer overflow safe.
    // it is not safe to use it in a remote service
    while ((ch = getopt(argc, argv, "f:s:h")) != -1) {
         switch (ch) {
         case 's':
             nmapvals.scale = strtod(optarg, &endptr);
             break;
         case 'f':
             for (int i = 0; i < FILTER_TYPE_COUNT; i++) {
                 if (strcmp(FILTER_TYPE_NAMES[i], optarg) == 0) {
                     nmapvals.filter = i;
                     found = 1;
                     break;
                 }
             }
             if (!found) {
                 printf("unkown filter: %s\n", optarg);
             }
             break;
         case 'h':
         default:
                 usage();
                 exit(1);
         }
    }
    argc -= optind;
    argv += optind;

    if (argc != 2) {
        usage();
        exit(1);
    }
    ExceptionInfo *exception;
    Image *image, *images, *resize_image;
    ImageInfo *image_info;


    MagickCoreGenesis(*argv,MagickTrue);
    exception=AcquireExceptionInfo();
    image_info=CloneImageInfo((ImageInfo *) NULL);
    (void) strcpy(image_info->filename,argv[0]);
    images = ReadImage(image_info,exception);
    if (exception->severity != UndefinedException) CatchException(exception);
    if (images == (Image *) NULL)
        exit(1);

    while ((image=RemoveFirstImageFromList(&images)) != (Image *) NULL) {

        MagickPixelPacket pixel_packet;
        ImageInfo * nm_info = CloneImageInfo((ImageInfo *) NULL);
        Image * nm_image = NewMagickImage(nm_info, image->columns, image->rows, &pixel_packet);

        normalmap(image, nm_image, nmapvals);

        (void) strcpy(nm_image->filename, argv[1]);
        WriteImage(image_info, nm_image);
        (void)DestroyImageList(nm_image);
        (void)DestroyImageInfo(nm_info);
    }

    image_info=DestroyImageInfo(image_info);
    exception=DestroyExceptionInfo(exception);
    MagickCoreTerminus();
    return(0); 
}


