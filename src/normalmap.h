/*
   Copyright (C) 2002-2008 Shawn Kirst <skirst@insightbb.com>

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
#ifndef NORMAL_MAP
#define NORMAL_MAP

#include <magick/MagickCore.h>

#define ENUM_EACH(PARAM) PARAM,
#define STR_EACH(PARAM) #PARAM,
#define STR_EACH_SEP_WS(PARAM) #PARAM " "

#define EACH_FILTER_TYPE(DO) \
    DO(FILTER_NONE) \
    DO(FILTER_SOBEL_3x3) \
    DO(FILTER_SOBEL_5x5) \
    DO(FILTER_PREWITT_3x3) \
    DO(FILTER_PREWITT_5x5) \
    DO(FILTER_3x3) \
    DO(FILTER_5x5) \
    DO(FILTER_7x7) \
    DO(FILTER_9x9) \
    DO(MAX_FILTER_TYPE)
#define FILTER_TYPE_COUNT (10)

enum FILTER_TYPE { EACH_FILTER_TYPE(ENUM_EACH) };

static char * FILTER_TYPE_NAMES[] = {
    EACH_FILTER_TYPE(STR_EACH)
};

enum ALPHA_TYPE
{
   ALPHA_NONE, ALPHA_HEIGHT, ALPHA_INVERSE_HEIGHT, ALPHA_ZERO, ALPHA_ONE,
   ALPHA_INVERT, ALPHA_MAP, MAX_ALPHA_TYPE
};

enum CONVERSION_TYPE
{
   CONVERT_NONE, CONVERT_BIASED_RGB, CONVERT_RED, CONVERT_GREEN, 
   CONVERT_BLUE, CONVERT_MAX_RGB, CONVERT_MIN_RGB, CONVERT_COLORSPACE,
   CONVERT_NORMALIZE_ONLY, CONVERT_DUDV_TO_NORMAL, CONVERT_HEIGHTMAP,
   MAX_CONVERSION_TYPE
};

enum DUDV_TYPE
{
   DUDV_NONE, DUDV_8BIT_SIGNED, DUDV_8BIT_UNSIGNED, DUDV_16BIT_SIGNED,
   DUDV_16BIT_UNSIGNED,
   MAX_DUDV_TYPE
};

typedef struct
{
   int filter;
   double minz;
   double scale;
   int wrap;
   int height_source;
   int alpha;
   int conversion;
   int dudv;
   int xinvert;
   int yinvert;
   int swapRGB;
   double contrast;
   int32_t alphamap_id;
} NormalmapVals;

int32_t normalmap(Image * img, Image * nm_img, NormalmapVals nmapvals);


#endif
