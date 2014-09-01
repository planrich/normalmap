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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <unistd.h>
#include <magick/MagickCore.h>
#include "normalmap.h"

#define PREVIEW_SIZE 150

static
int sample_alpha_map(unsigned char *pixels, int x, int y, int w, int h, int sw, int sh);

static inline int icerp(int a, int b, int c, int d, int x)
{
   int p = (d - c) - (a - b);
   int q = (a - b) - p;
   int r = c - a;
   return((x * (x * (x * p + (q << 7)) + (r << 14)) + (b << 21)) >> 21);
}

void scale_pixels(unsigned char *dst, int dw, int dh, unsigned char *src, int sw, int sh, int bpp)
{
   int x, y, n, ix, iy, wx, wy, v;
   int a, b, c, d;
   int dstride = dw * bpp;
   unsigned char *s;
   
   for(y = 0; y < dh; ++y)
   {
      if(dh > 1)
      {
         iy = (((sh - 1) * y) << 7) / (dh - 1);
         if(y == dh - 1) --iy;
         wy = iy & 0x7f;
         iy >>= 7;
      }
      else
         iy = wy = 0;
      
      for(x = 0; x < dw; ++x)
      {
         if(dw > 1)
         {
            ix = (((sw - 1) * x) << 7) / (dw - 1);
            if(x == dw - 1) --ix;
            wx = ix & 0x7f;
            ix >>= 7;
         }
         else
            ix = wx = 0;
         
         s = src + ((iy - 1) * sw + (ix - 1)) * bpp;
         
         for(n = 0; n < bpp; ++n)
         {
            b = icerp(s[(sw + 0) * bpp],
                      s[(sw + 1) * bpp],
                      s[(sw + 2) * bpp],
                      s[(sw + 3) * bpp], wx);
            if(iy > 0)
            {
               a = icerp(s[      0],
                         s[    bpp],
                         s[2 * bpp],
                         s[3 * bpp], wx);
            }
            else
               a = b;
            
            c = icerp(s[(2 * sw + 0) * bpp],
                      s[(2 * sw + 1) * bpp],
                      s[(2 * sw + 2) * bpp],
                      s[(2 * sw + 3) * bpp], wx);
            if(iy < dh - 1)
            {
               d = icerp(s[(3 * sw + 0) * bpp],
                         s[(3 * sw + 1) * bpp],
                         s[(3 * sw + 2) * bpp],
                         s[(3 * sw + 3) * bpp], wx);
            }
            else
               d = c;
            
            v = icerp(a, b, c, d, wy);
            if(v < 0) v = 0;
            if(v > 255) v = 255;
            dst[(y * dstride) + (x * bpp) + n] = v;
            ++s;
         }
      }
   }
}

static const float oneover255 = 1.0f / 255.0f;

#ifndef min
#define min(a,b)  ((a)<(b) ? (a) : (b))
#endif
#ifndef max
#define max(a,b)  ((a)>(b) ? (a) : (b))
#endif

#define SQR(x)      ((x) * (x))
#define LERP(a,b,c) ((a) + ((b) - (a)) * (c))

static inline void NORMALIZE(float *v)
{
   float len = sqrtf(SQR(v[0]) + SQR(v[1]) + SQR(v[2]));
   
   if(len > 1e-04f)
   {
      len = 1.0f / len;
      v[0] *= len;
      v[1] *= len;
      v[2] *= len;
   }
   else
      v[0] = v[1] = v[2] = 0;
}

typedef struct
{
   int x,y;
   float w;
} kernel_element;

static void make_kernel(kernel_element *k, float *weights, int size)
{
   int x, y, idx;
   
   for(y = 0; y < size; ++y)
   {
      for(x = 0; x < size; ++x)
      {
         idx = x + y * size;
         k[idx].x = x - (size / 2);
         k[idx].y = (size / 2) - y;
         k[idx].w = weights[idx];
      }
   }
}

static void rotate_array(float *dst, float *src, int size)
{
   int x, y, newx, newy;
   
   for(y = 0; y < size; ++y)
   {
      for(x = 0; x < size; ++x)
      {
         newy = size - x - 1;
         newx = y;
         dst[newx + newy * size] = src[x + y * size];
      }
   }
}

int sample_alpha_map(unsigned char *pixels, int x, int y,
                            int w, int h, int sw, int sh)
{
   int ix, iy, wx, wy, v;
   int a, b, c, d;
   unsigned char *s;
   
   if(sh > 1)
   {
      iy = (((h - 1) * y) << 7) / (sh - 1);
      if(y == sh - 1) --iy;
      wy = iy & 0x7f;
      iy >>= 7;
   }
   else
      iy = wy = 0;
   
   if(sw > 1)
   {
      ix = (((w - 1) * x) << 7) / (sw - 1);
      if(x == sw - 1) --ix;
      wx = ix & 0x7f;
      ix >>= 7;
   }
   else
      ix = wx = 0;
   
   s = pixels + ((iy - 1) * w + (ix - 1));
   
   b = icerp(s[w + 0],
             s[w + 1],
             s[w + 2],
             s[w + 3], wx);
   if(iy > 0)
   {
      a = icerp(s[0],
                s[1],
                s[2],
                s[3], wx);
   }
   else
      a = b;
   
   c = icerp(s[2 * w + 0],
             s[2 * w + 1],
             s[2 * w + 2],
             s[2 * w + 3], wx);
   if(iy < sh - 1)
   {
      d = icerp(s[3 * w + 0],
                s[3 * w + 1],
                s[3 * w + 2],
                s[3 * w + 3], wx);
   }
   else
      d = c;
   
   v = icerp(a, b, c, d, wy);
            
   if(v <   0) v = 0;
   if(v > 255) v = 255;
   
   return((unsigned char)v);
}

static void make_heightmap(unsigned char *image, int w, int h, int bpp, NormalmapVals nmapvals)
{
   unsigned int i, num_pixels = w * h;
   int x, y;
   float v, hmin, hmax;
   float *s, *r;
   
   s = (float*)malloc(w * h * 3 * sizeof(float));
   if(s == 0)
   {
      perror("Memory allocation error!");
      return;
   }
   r = (float*)malloc(w * h * 4 * sizeof(float));
   if(r == 0)
   {
      free(s);
      perror("Memory allocation error!");
      return;
   }
   
   /* scale into 0 to 1 range, make signed -1 to 1 */
   for(i = 0; i < num_pixels; ++i)
   {
      s[3 * i + 0] = (((float)image[bpp * i + 0] / 255.0f) - 0.5) * 2.0f;
      s[3 * i + 1] = (((float)image[bpp * i + 1] / 255.0f) - 0.5) * 2.0f;
      s[3 * i + 2] = (((float)image[bpp * i + 2] / 255.0f) - 0.5) * 2.0f;
   }
   
   memset(r, 0, w * h * 4 * sizeof(float));

#define S(x, y, n) s[(y) * (w * 3) + ((x) * 3) + (n)]
#define R(x, y, n) r[(y) * (w * 4) + ((x) * 4) + (n)]
   
   /* top-left to bottom-right */
   for(x = 1; x < w; ++x)
      R(x, 0, 0) = R(x - 1, 0, 0) + S(x - 1, 0, 0);
   for(y = 1; y < h; ++y)
      R(0, y, 0) = R(0, y - 1, 0) + S(0, y - 1, 1);
   for(y = 1; y < h; ++y)
   {
      for(x = 1; x < w; ++x)
      {
         R(x, y, 0) = (R(x, y - 1, 0) + R(x - 1, y, 0) +
                       S(x - 1, y, 0) + S(x, y - 1, 1)) * 0.5f;
      }
   }

   /* top-right to bottom-left */
   for(x = w - 2; x >= 0; --x)
      R(x, 0, 1) = R(x + 1, 0, 1) - S(x + 1, 0, 0);
   for(y = 1; y < h; ++y)
      R(0, y, 1) = R(0, y - 1, 1) + S(0, y - 1, 1);
   for(y = 1; y < h; ++y)
   {
      for(x = w - 2; x >= 0; --x)
      {
         R(x, y, 1) = (R(x, y - 1, 1) + R(x + 1, y, 1) -
                       S(x + 1, y, 0) + S(x, y - 1, 1)) * 0.5f;
      }
   }

   /* bottom-left to top-right */
   for(x = 1; x < w; ++x)
      R(x, 0, 2) = R(x - 1, 0, 2) + S(x - 1, 0, 0);
   for(y = h - 2; y >= 0; --y)
      R(0, y, 2) = R(0, y + 1, 2) - S(0, y + 1, 1);
   for(y = h - 2; y >= 0; --y)
   {
      for(x = 1; x < w; ++x)
      {
         R(x, y, 2) = (R(x, y + 1, 2) + R(x - 1, y, 2) +
                       S(x - 1, y, 0) - S(x, y + 1, 1)) * 0.5f;
      }
   }

   /* bottom-right to top-left */
   for(x = w - 2; x >= 0; --x)
      R(x, 0, 3) = R(x + 1, 0, 3) - S(x + 1, 0, 0);
   for(y = h - 2; y >= 0; --y)
      R(0, y, 3) = R(0, y + 1, 3) - S(0, y + 1, 1);
   for(y = h - 2; y >= 0; --y)
   {
      for(x = w - 2; x >= 0; --x)
      {
         R(x, y, 3) = (R(x, y + 1, 3) + R(x + 1, y, 3) -
                       S(x + 1, y, 0) - S(x, y + 1, 1)) * 0.5f;
      }
   }

#undef S
#undef R

   /* accumulate, find min/max */
   hmin =  1e10f;
   hmax = -1e10f;
   for(i = 0; i < num_pixels; ++i)
   {
      r[4 * i] += r[4 * i + 1] + r[4 * i + 2] + r[4 * i + 3];
      if(r[4 * i] < hmin) hmin = r[4 * i];
      if(r[4 * i] > hmax) hmax = r[4 * i];
   }

   /* scale into 0 - 1 range */
   for(i = 0; i < num_pixels; ++i)
   {   
      v = (r[4 * i] - hmin) / (hmax - hmin);
      /* adjust contrast */
      v = (v - 0.5f) * nmapvals.contrast + v;
      if(v < 0) v = 0;
      if(v > 1) v = 1;
      r[4 * i] = v;
   }

   /* write out results */
   for(i = 0; i < num_pixels; ++i)
   {
      v = r[4 * i] * 255.0f;
      image[bpp * i + 0] = (unsigned char)v;
      image[bpp * i + 1] = (unsigned char)v;
      image[bpp * i + 2] = (unsigned char)v;
   }

   free(s);
   free(r);
}

int32_t normalmap(Image * img, Image * nm_img, NormalmapVals nmapvals)
{
   int x, y;
   int width, height, bpp, rowbytes, pw, ph, amap_w = 0, amap_h = 0;
   unsigned char *d, *dst, *s, *src, *tmp, *amap = 0;
   float *heights;
   float val, du, dv, n[3], weight;
   float rgb_bias[3];
   int i, num_elements = 0;
   kernel_element *kernel_du = 0;
   kernel_element *kernel_dv = 0;

   ExceptionInfo *exception;

   bpp = 4; // TODO drawable->bpp;
   
   if(nmapvals.filter < 0 || nmapvals.filter >= MAX_FILTER_TYPE)
      nmapvals.filter = FILTER_NONE;
   if(bpp != 4) nmapvals.height_source = 0;
   if(bpp != 4 && (nmapvals.dudv == DUDV_16BIT_SIGNED || nmapvals.dudv == DUDV_16BIT_UNSIGNED))
      nmapvals.dudv = DUDV_NONE;
   
   width = img->columns;
   height = img->rows;
   rowbytes = width * bpp;

   dst = malloc(width * height * bpp);
   if(dst == 0)
   {
      perror("Memory allocation error!");
      return(-1);
   }
   ExportImagePixels(nm_img, 0, 0, width, height, "RGBA", CharPixel, dst, exception);

   src = malloc(width * height * bpp);
   if(src == 0)
   {
      perror("Memory allocation error!");
      return(-1);
   }
   ExportImagePixels(img, 0, 0, width, height, "RGBA", CharPixel, src, exception);
   
   heights = calloc(width * height, sizeof(float));
   if(heights == 0)
   {
      perror("Memory allocation error!");
      return(-1);
   }

   if(!nmapvals.dudv && bpp == 4 && nmapvals.alpha == ALPHA_MAP && nmapvals.alphamap_id != 0)
   {
      //GimpDrawable *alphamap = gimp_drawable_get(nmapvals.alphamap_id);
      
      //amap_w = alphamap->width;
      //amap_h = alphamap->height;
      //amap = malloc(amap_w * amap_h);
      
      //gimp_pixel_rgn_init(&amap_rgn, alphamap, 0, 0, amap_w, amap_h, 0, 0);
      //gimp_pixel_rgn_get_rect(&amap_rgn, amap, 0, 0, amap_w, amap_h);
   }

   //gimp_pixel_rgn_init(&src_rgn, drawable, 0, 0, width, height, 0, 0);
   //gimp_pixel_rgn_get_rect(&src_rgn, src, 0, 0, width, height);
   
   switch(nmapvals.filter)
   {
      case FILTER_NONE:
         num_elements = 2;
         kernel_du = (kernel_element*)malloc(2 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(2 * sizeof(kernel_element));
         
         kernel_du[0].x = -1; kernel_du[0].y = 0; kernel_du[0].w = -0.5f;
         kernel_du[1].x =  1; kernel_du[1].y = 0; kernel_du[1].w =  0.5f;
         
         kernel_dv[0].x = 0; kernel_dv[0].y =  1; kernel_dv[0].w =  0.5f;
         kernel_dv[1].x = 0; kernel_dv[1].y = -1; kernel_dv[1].w = -0.5f;
         
         break;
      case FILTER_SOBEL_3x3:
         num_elements = 6;
         kernel_du = (kernel_element*)malloc(6 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(6 * sizeof(kernel_element));
      
         kernel_du[0].x = -1; kernel_du[0].y =  1; kernel_du[0].w = -1.0f;
         kernel_du[1].x = -1; kernel_du[1].y =  0; kernel_du[1].w = -2.0f;
         kernel_du[2].x = -1; kernel_du[2].y = -1; kernel_du[2].w = -1.0f;
         kernel_du[3].x =  1; kernel_du[3].y =  1; kernel_du[3].w =  1.0f;
         kernel_du[4].x =  1; kernel_du[4].y =  0; kernel_du[4].w =  2.0f;
         kernel_du[5].x =  1; kernel_du[5].y = -1; kernel_du[5].w =  1.0f;
         
         kernel_dv[0].x = -1; kernel_dv[0].y =  1; kernel_dv[0].w =  1.0f;
         kernel_dv[1].x =  0; kernel_dv[1].y =  1; kernel_dv[1].w =  2.0f;
         kernel_dv[2].x =  1; kernel_dv[2].y =  1; kernel_dv[2].w =  1.0f;
         kernel_dv[3].x = -1; kernel_dv[3].y = -1; kernel_dv[3].w = -1.0f;
         kernel_dv[4].x =  0; kernel_dv[4].y = -1; kernel_dv[4].w = -2.0f;
         kernel_dv[5].x =  1; kernel_dv[5].y = -1; kernel_dv[5].w = -1.0f;
         
         break;
      case FILTER_SOBEL_5x5:
         num_elements = 20;
         kernel_du = (kernel_element*)malloc(20 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(20 * sizeof(kernel_element));

         kernel_du[ 0].x = -2; kernel_du[ 0].y =  2; kernel_du[ 0].w =  -1.0f;
         kernel_du[ 1].x = -2; kernel_du[ 1].y =  1; kernel_du[ 1].w =  -4.0f;
         kernel_du[ 2].x = -2; kernel_du[ 2].y =  0; kernel_du[ 2].w =  -6.0f;
         kernel_du[ 3].x = -2; kernel_du[ 3].y = -1; kernel_du[ 3].w =  -4.0f;
         kernel_du[ 4].x = -2; kernel_du[ 4].y = -2; kernel_du[ 4].w =  -1.0f;
         kernel_du[ 5].x = -1; kernel_du[ 5].y =  2; kernel_du[ 5].w =  -2.0f;
         kernel_du[ 6].x = -1; kernel_du[ 6].y =  1; kernel_du[ 6].w =  -8.0f;
         kernel_du[ 7].x = -1; kernel_du[ 7].y =  0; kernel_du[ 7].w = -12.0f;
         kernel_du[ 8].x = -1; kernel_du[ 8].y = -1; kernel_du[ 8].w =  -8.0f;
         kernel_du[ 9].x = -1; kernel_du[ 9].y = -2; kernel_du[ 9].w =  -2.0f;
         kernel_du[10].x =  1; kernel_du[10].y =  2; kernel_du[10].w =   2.0f;
         kernel_du[11].x =  1; kernel_du[11].y =  1; kernel_du[11].w =   8.0f;
         kernel_du[12].x =  1; kernel_du[12].y =  0; kernel_du[12].w =  12.0f;
         kernel_du[13].x =  1; kernel_du[13].y = -1; kernel_du[13].w =   8.0f;
         kernel_du[14].x =  1; kernel_du[14].y = -2; kernel_du[14].w =   2.0f;
         kernel_du[15].x =  2; kernel_du[15].y =  2; kernel_du[15].w =   1.0f;
         kernel_du[16].x =  2; kernel_du[16].y =  1; kernel_du[16].w =   4.0f;
         kernel_du[17].x =  2; kernel_du[17].y =  0; kernel_du[17].w =   6.0f;
         kernel_du[18].x =  2; kernel_du[18].y = -1; kernel_du[18].w =   4.0f;
         kernel_du[19].x =  2; kernel_du[19].y = -2; kernel_du[19].w =   1.0f;
      
         kernel_dv[ 0].x = -2; kernel_dv[ 0].y =  2; kernel_dv[ 0].w =   1.0f;
         kernel_dv[ 1].x = -1; kernel_dv[ 1].y =  2; kernel_dv[ 1].w =   4.0f;
         kernel_dv[ 2].x =  0; kernel_dv[ 2].y =  2; kernel_dv[ 2].w =   6.0f;
         kernel_dv[ 3].x =  1; kernel_dv[ 3].y =  2; kernel_dv[ 3].w =   4.0f;
         kernel_dv[ 4].x =  2; kernel_dv[ 4].y =  2; kernel_dv[ 4].w =   1.0f;
         kernel_dv[ 5].x = -2; kernel_dv[ 5].y =  1; kernel_dv[ 5].w =   2.0f;
         kernel_dv[ 6].x = -1; kernel_dv[ 6].y =  1; kernel_dv[ 6].w =   8.0f;
         kernel_dv[ 7].x =  0; kernel_dv[ 7].y =  1; kernel_dv[ 7].w =  12.0f;
         kernel_dv[ 8].x =  1; kernel_dv[ 8].y =  1; kernel_dv[ 8].w =   8.0f;
         kernel_dv[ 9].x =  2; kernel_dv[ 9].y =  1; kernel_dv[ 9].w =   2.0f;
         kernel_dv[10].x = -2; kernel_dv[10].y = -1; kernel_dv[10].w =  -2.0f;
         kernel_dv[11].x = -1; kernel_dv[11].y = -1; kernel_dv[11].w =  -8.0f;
         kernel_dv[12].x =  0; kernel_dv[12].y = -1; kernel_dv[12].w = -12.0f;
         kernel_dv[13].x =  1; kernel_dv[13].y = -1; kernel_dv[13].w =  -8.0f;
         kernel_dv[14].x =  2; kernel_dv[14].y = -1; kernel_dv[14].w =  -2.0f;
         kernel_dv[15].x = -2; kernel_dv[15].y = -2; kernel_dv[15].w =  -1.0f;
         kernel_dv[16].x = -1; kernel_dv[16].y = -2; kernel_dv[16].w =  -4.0f;
         kernel_dv[17].x =  0; kernel_dv[17].y = -2; kernel_dv[17].w =  -6.0f;
         kernel_dv[18].x =  1; kernel_dv[18].y = -2; kernel_dv[18].w =  -4.0f;
         kernel_dv[19].x =  2; kernel_dv[19].y = -2; kernel_dv[19].w =  -1.0f;
      
         break;
      case FILTER_PREWITT_3x3:
         num_elements = 6;
         kernel_du = (kernel_element*)malloc(6 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(6 * sizeof(kernel_element));
      
         kernel_du[0].x = -1; kernel_du[0].y =  1; kernel_du[0].w = -1.0f;
         kernel_du[1].x = -1; kernel_du[1].y =  0; kernel_du[1].w = -1.0f;
         kernel_du[2].x = -1; kernel_du[2].y = -1; kernel_du[2].w = -1.0f;
         kernel_du[3].x =  1; kernel_du[3].y =  1; kernel_du[3].w =  1.0f;
         kernel_du[4].x =  1; kernel_du[4].y =  0; kernel_du[4].w =  1.0f;
         kernel_du[5].x =  1; kernel_du[5].y = -1; kernel_du[5].w =  1.0f;
         
         kernel_dv[0].x = -1; kernel_dv[0].y =  1; kernel_dv[0].w =  1.0f;
         kernel_dv[1].x =  0; kernel_dv[1].y =  1; kernel_dv[1].w =  1.0f;
         kernel_dv[2].x =  1; kernel_dv[2].y =  1; kernel_dv[2].w =  1.0f;
         kernel_dv[3].x = -1; kernel_dv[3].y = -1; kernel_dv[3].w = -1.0f;
         kernel_dv[4].x =  0; kernel_dv[4].y = -1; kernel_dv[4].w = -1.0f;
         kernel_dv[5].x =  1; kernel_dv[5].y = -1; kernel_dv[5].w = -1.0f;
         
         break;      
      case FILTER_PREWITT_5x5:
         num_elements = 20;
         kernel_du = (kernel_element*)malloc(20 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(20 * sizeof(kernel_element));

         kernel_du[ 0].x = -2; kernel_du[ 0].y =  2; kernel_du[ 0].w = -1.0f;
         kernel_du[ 1].x = -2; kernel_du[ 1].y =  1; kernel_du[ 1].w = -1.0f;
         kernel_du[ 2].x = -2; kernel_du[ 2].y =  0; kernel_du[ 2].w = -1.0f;
         kernel_du[ 3].x = -2; kernel_du[ 3].y = -1; kernel_du[ 3].w = -1.0f;
         kernel_du[ 4].x = -2; kernel_du[ 4].y = -2; kernel_du[ 4].w = -1.0f;
         kernel_du[ 5].x = -1; kernel_du[ 5].y =  2; kernel_du[ 5].w = -2.0f;
         kernel_du[ 6].x = -1; kernel_du[ 6].y =  1; kernel_du[ 6].w = -2.0f;
         kernel_du[ 7].x = -1; kernel_du[ 7].y =  0; kernel_du[ 7].w = -2.0f;
         kernel_du[ 8].x = -1; kernel_du[ 8].y = -1; kernel_du[ 8].w = -2.0f;
         kernel_du[ 9].x = -1; kernel_du[ 9].y = -2; kernel_du[ 9].w = -2.0f;
         kernel_du[10].x =  1; kernel_du[10].y =  2; kernel_du[10].w =  2.0f;
         kernel_du[11].x =  1; kernel_du[11].y =  1; kernel_du[11].w =  2.0f;
         kernel_du[12].x =  1; kernel_du[12].y =  0; kernel_du[12].w =  2.0f;
         kernel_du[13].x =  1; kernel_du[13].y = -1; kernel_du[13].w =  2.0f;
         kernel_du[14].x =  1; kernel_du[14].y = -2; kernel_du[14].w =  2.0f;
         kernel_du[15].x =  2; kernel_du[15].y =  2; kernel_du[15].w =  1.0f;
         kernel_du[16].x =  2; kernel_du[16].y =  1; kernel_du[16].w =  1.0f;
         kernel_du[17].x =  2; kernel_du[17].y =  0; kernel_du[17].w =  1.0f;
         kernel_du[18].x =  2; kernel_du[18].y = -1; kernel_du[18].w =  1.0f;
         kernel_du[19].x =  2; kernel_du[19].y = -2; kernel_du[19].w =  1.0f;
      
         kernel_dv[ 0].x = -2; kernel_dv[ 0].y =  2; kernel_dv[ 0].w =  1.0f;
         kernel_dv[ 1].x = -1; kernel_dv[ 1].y =  2; kernel_dv[ 1].w =  1.0f;
         kernel_dv[ 2].x =  0; kernel_dv[ 2].y =  2; kernel_dv[ 2].w =  1.0f;
         kernel_dv[ 3].x =  1; kernel_dv[ 3].y =  2; kernel_dv[ 3].w =  1.0f;
         kernel_dv[ 4].x =  2; kernel_dv[ 4].y =  2; kernel_dv[ 4].w =  1.0f;
         kernel_dv[ 5].x = -2; kernel_dv[ 5].y =  1; kernel_dv[ 5].w =  2.0f;
         kernel_dv[ 6].x = -1; kernel_dv[ 6].y =  1; kernel_dv[ 6].w =  2.0f;
         kernel_dv[ 7].x =  0; kernel_dv[ 7].y =  1; kernel_dv[ 7].w =  2.0f;
         kernel_dv[ 8].x =  1; kernel_dv[ 8].y =  1; kernel_dv[ 8].w =  2.0f;
         kernel_dv[ 9].x =  2; kernel_dv[ 9].y =  1; kernel_dv[ 9].w =  2.0f;
         kernel_dv[10].x = -2; kernel_dv[10].y = -1; kernel_dv[10].w = -2.0f;
         kernel_dv[11].x = -1; kernel_dv[11].y = -1; kernel_dv[11].w = -2.0f;
         kernel_dv[12].x =  0; kernel_dv[12].y = -1; kernel_dv[12].w = -2.0f;
         kernel_dv[13].x =  1; kernel_dv[13].y = -1; kernel_dv[13].w = -2.0f;
         kernel_dv[14].x =  2; kernel_dv[14].y = -1; kernel_dv[14].w = -2.0f;
         kernel_dv[15].x = -2; kernel_dv[15].y = -2; kernel_dv[15].w = -1.0f;
         kernel_dv[16].x = -1; kernel_dv[16].y = -2; kernel_dv[16].w = -1.0f;
         kernel_dv[17].x =  0; kernel_dv[17].y = -2; kernel_dv[17].w = -1.0f;
         kernel_dv[18].x =  1; kernel_dv[18].y = -2; kernel_dv[18].w = -1.0f;
         kernel_dv[19].x =  2; kernel_dv[19].y = -2; kernel_dv[19].w = -1.0f;
      
         break;
      case FILTER_3x3:
         num_elements = 6;
         kernel_du = (kernel_element*)malloc(6 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(6 * sizeof(kernel_element));
      
         weight = 1.0f / 6.0f;
      
         kernel_du[0].x = -1; kernel_du[0].y =  1; kernel_du[0].w = -weight;
         kernel_du[1].x = -1; kernel_du[1].y =  0; kernel_du[1].w = -weight;
         kernel_du[2].x = -1; kernel_du[2].y = -1; kernel_du[2].w = -weight;
         kernel_du[3].x =  1; kernel_du[3].y =  1; kernel_du[3].w =  weight;
         kernel_du[4].x =  1; kernel_du[4].y =  0; kernel_du[4].w =  weight;
         kernel_du[5].x =  1; kernel_du[5].y = -1; kernel_du[5].w =  weight;
         
         kernel_dv[0].x = -1; kernel_dv[0].y =  1; kernel_dv[0].w =  weight;
         kernel_dv[1].x =  0; kernel_dv[1].y =  1; kernel_dv[1].w =  weight;
         kernel_dv[2].x =  1; kernel_dv[2].y =  1; kernel_dv[2].w =  weight;
         kernel_dv[3].x = -1; kernel_dv[3].y = -1; kernel_dv[3].w = -weight;
         kernel_dv[4].x =  0; kernel_dv[4].y = -1; kernel_dv[4].w = -weight;
         kernel_dv[5].x =  1; kernel_dv[5].y = -1; kernel_dv[5].w = -weight;
         break;
      case FILTER_5x5:
      {
         int n;
         float usum = 0, vsum = 0;
         float wt22 = 1.0f / 16.0f;
         float wt12 = 1.0f / 10.0f;
         float wt02 = 1.0f / 8.0f;
         float wt11 = 1.0f / 2.8f;
         num_elements = 20;
         kernel_du = (kernel_element*)malloc(20 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(20 * sizeof(kernel_element));
         
         kernel_du[0 ].x = -2; kernel_du[0 ].y =  2; kernel_du[0 ].w = -wt22;
         kernel_du[1 ].x = -1; kernel_du[1 ].y =  2; kernel_du[1 ].w = -wt12;
         kernel_du[2 ].x =  1; kernel_du[2 ].y =  2; kernel_du[2 ].w =  wt12;
         kernel_du[3 ].x =  2; kernel_du[3 ].y =  2; kernel_du[3 ].w =  wt22;
         kernel_du[4 ].x = -2; kernel_du[4 ].y =  1; kernel_du[4 ].w = -wt12;
         kernel_du[5 ].x = -1; kernel_du[5 ].y =  1; kernel_du[5 ].w = -wt11;
         kernel_du[6 ].x =  1; kernel_du[6 ].y =  1; kernel_du[6 ].w =  wt11;
         kernel_du[7 ].x =  2; kernel_du[7 ].y =  1; kernel_du[7 ].w =  wt12;
         kernel_du[8 ].x = -2; kernel_du[8 ].y =  0; kernel_du[8 ].w = -wt02;
         kernel_du[9 ].x = -1; kernel_du[9 ].y =  0; kernel_du[9 ].w = -0.5f;
         kernel_du[10].x =  1; kernel_du[10].y =  0; kernel_du[10].w =  0.5f;
         kernel_du[11].x =  2; kernel_du[11].y =  0; kernel_du[11].w =  wt02;
         kernel_du[12].x = -2; kernel_du[12].y = -1; kernel_du[12].w = -wt12;
         kernel_du[13].x = -1; kernel_du[13].y = -1; kernel_du[13].w = -wt11;
         kernel_du[14].x =  1; kernel_du[14].y = -1; kernel_du[14].w =  wt11;
         kernel_du[15].x =  2; kernel_du[15].y = -1; kernel_du[15].w =  wt12;
         kernel_du[16].x = -2; kernel_du[16].y = -2; kernel_du[16].w = -wt22;
         kernel_du[17].x = -1; kernel_du[17].y = -2; kernel_du[17].w = -wt12;
         kernel_du[18].x =  1; kernel_du[18].y = -2; kernel_du[18].w =  wt12;
         kernel_du[19].x =  2; kernel_du[19].y = -2; kernel_du[19].w =  wt22;
         
         kernel_dv[0 ].x = -2; kernel_dv[0 ].y =  2; kernel_dv[0 ].w =  wt22;
         kernel_dv[1 ].x = -1; kernel_dv[1 ].y =  2; kernel_dv[1 ].w =  wt12;
         kernel_dv[2 ].x =  0; kernel_dv[2 ].y =  2; kernel_dv[2 ].w =  0.25f;
         kernel_dv[3 ].x =  1; kernel_dv[3 ].y =  2; kernel_dv[3 ].w =  wt12;
         kernel_dv[4 ].x =  2; kernel_dv[4 ].y =  2; kernel_dv[4 ].w =  wt22;
         kernel_dv[5 ].x = -2; kernel_dv[5 ].y =  1; kernel_dv[5 ].w =  wt12;
         kernel_dv[6 ].x = -1; kernel_dv[6 ].y =  1; kernel_dv[6 ].w =  wt11;
         kernel_dv[7 ].x =  0; kernel_dv[7 ].y =  1; kernel_dv[7 ].w =  0.5f;
         kernel_dv[8 ].x =  1; kernel_dv[8 ].y =  1; kernel_dv[8 ].w =  wt11;
         kernel_dv[9 ].x =  2; kernel_dv[9 ].y =  1; kernel_dv[9 ].w =  wt22;
         kernel_dv[10].x = -2; kernel_dv[10].y = -1; kernel_dv[10].w = -wt22;
         kernel_dv[11].x = -1; kernel_dv[11].y = -1; kernel_dv[11].w = -wt11;
         kernel_dv[12].x =  0; kernel_dv[12].y = -1; kernel_dv[12].w = -0.5f;
         kernel_dv[13].x =  1; kernel_dv[13].y = -1; kernel_dv[13].w = -wt11;
         kernel_dv[14].x =  2; kernel_dv[14].y = -1; kernel_dv[14].w = -wt12;
         kernel_dv[15].x = -2; kernel_dv[15].y = -2; kernel_dv[15].w = -wt22;
         kernel_dv[16].x = -1; kernel_dv[16].y = -2; kernel_dv[16].w = -wt12;
         kernel_dv[17].x =  0; kernel_dv[17].y = -2; kernel_dv[17].w = -0.25f;
         kernel_dv[18].x =  1; kernel_dv[18].y = -2; kernel_dv[18].w = -wt12;
         kernel_dv[19].x =  2; kernel_dv[19].y = -2; kernel_dv[19].w = -wt22;

         for(n = 0; n < 20; ++n)
         {
            usum += fabsf(kernel_du[n].w);
            vsum += fabsf(kernel_dv[n].w);
         }
         for(n = 0; n < 20; ++n)
         {
            kernel_du[n].w /= usum;
            kernel_dv[n].w /= vsum;
         }
         
         break;
      }
      case FILTER_7x7:
      {
         float du_weights[]=
         {
            -1, -2, -3, 0, 3, 2, 1,
            -2, -3, -4, 0, 4, 3, 2,
            -3, -4, -5, 0, 5, 4, 3,
            -4, -5, -6, 0, 6, 5, 4,
            -3, -4, -5, 0, 5, 4, 3,
            -2, -3, -4, 0, 4, 3, 2,
            -1, -2, -3, 0, 3, 2, 1   
         };
         float dv_weights[49];
         int n;
         float usum = 0, vsum = 0;
         
         num_elements = 49;
         kernel_du = (kernel_element*)malloc(49 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(49 * sizeof(kernel_element));
         
         make_kernel(kernel_du, du_weights, 7);
         rotate_array(dv_weights, du_weights, 7);
         make_kernel(kernel_dv, dv_weights, 7);
         
         for(n = 0; n < 49; ++n)
         {
            usum += fabsf(kernel_du[n].w);
            vsum += fabsf(kernel_dv[n].w);
         }
         for(n = 0; n < 49; ++n)
         {
            kernel_du[n].w /= usum;
            kernel_dv[n].w /= vsum;
         }
         
         break;
      }
      case FILTER_9x9:
      {
         float du_weights[]=
         {
            -1, -2, -3, -4, 0, 4, 3, 2, 1,
            -2, -3, -4, -5, 0, 5, 4, 3, 2,
            -3, -4, -5, -6, 0, 6, 5, 4, 3,
            -4, -5, -6, -7, 0, 7, 6, 5, 4,
            -5, -6, -7, -8, 0, 8, 7, 6, 5,
            -4, -5, -6, -7, 0, 7, 6, 5, 4,
            -3, -4, -5, -6, 0, 6, 5, 4, 3,
            -2, -3, -4, -5, 0, 5, 4, 3, 2,
            -1, -2, -3, -4, 0, 4, 3, 2, 1     
         };
         float dv_weights[81];
         int n;
         float usum = 0, vsum = 0;
         
         num_elements = 81;
         kernel_du = (kernel_element*)malloc(81 * sizeof(kernel_element));
         kernel_dv = (kernel_element*)malloc(81 * sizeof(kernel_element));
         
         make_kernel(kernel_du, du_weights, 9);
         rotate_array(dv_weights, du_weights, 9);
         make_kernel(kernel_dv, dv_weights, 9);
         
         for(n = 0; n < 81; ++n)
         {
            usum += fabsf(kernel_du[n].w);
            vsum += fabsf(kernel_dv[n].w);
         }
         for(n = 0; n < 81; ++n)
         {
            kernel_du[n].w /= usum;
            kernel_dv[n].w /= vsum;
         }
         
         break;
      }
   }

   if(nmapvals.conversion == CONVERT_BIASED_RGB)
   {
      /* approximated average color of the image
       * scale to 16x16, accumulate the pixels and average */
      unsigned int sum[3];

      tmp = malloc(16 * 16 * bpp);
      scale_pixels(tmp, 16, 16, src, width, height, bpp);

      sum[0] = sum[1] = sum[2] = 0;
      
      s = src;
      for(y = 0; y < 16; ++y)
      {
         for(x = 0; x < 16; ++x)
         {
            sum[0] += *s++;
            sum[1] += *s++;
            sum[2] += *s++;
            if(bpp == 4) s++;
         }
      }
         
      rgb_bias[0] = (float)sum[0] / 256.0f;
      rgb_bias[1] = (float)sum[1] / 256.0f;
      rgb_bias[2] = (float)sum[2] / 256.0f;
      
      free(tmp);
   }
   else
   {
      rgb_bias[0] = 0;
      rgb_bias[1] = 0;
      rgb_bias[2] = 0;
   }

   if(nmapvals.conversion != CONVERT_NORMALIZE_ONLY &&
      nmapvals.conversion != CONVERT_DUDV_TO_NORMAL &&
      nmapvals.conversion != CONVERT_HEIGHTMAP)
   {
      s = src;
      for(y = 0; y < height; ++y)
      {
         for(x = 0; x < width; ++x)
         {
            if(!nmapvals.height_source)
            {
               switch(nmapvals.conversion)
               {
                  case CONVERT_NONE:
                     val = (float)s[0] * 0.3f +
                           (float)s[1] * 0.59f +
                           (float)s[2] * 0.11f;
                     break;
                  case CONVERT_BIASED_RGB:
                     val = (((float)max(0, s[0] - rgb_bias[0])) * 0.3f ) +
                           (((float)max(0, s[1] - rgb_bias[1])) * 0.59f) +
                           (((float)max(0, s[2] - rgb_bias[2])) * 0.11f);
                     break;
                  case CONVERT_RED:
                     val = (float)s[0];
                     break;
                  case CONVERT_GREEN:
                     val = (float)s[1];
                     break;
                  case CONVERT_BLUE:
                     val = (float)s[2];
                     break;
                  case CONVERT_MAX_RGB:
                     val = (float)max(s[0], max(s[1], s[2]));
                     break;
                  case CONVERT_MIN_RGB:
                     val = (float)min(s[0], min(s[1], s[2]));
                     break;
                  case CONVERT_COLORSPACE:
                     val = (1.0f - ((1.0f - ((float)s[0] / 255.0f)) *
                                    (1.0f - ((float)s[1] / 255.0f)) *
                                    (1.0f - ((float)s[2] / 255.0f)))) * 255.0f;
                     break;
                  default:
                     val = 255.0f;
                     break;
               }
            }
            else
               val = (float)s[3];
         
            heights[x + y * width] = val * oneover255;
         
            s += bpp;
         }
      }
   }

#define HEIGHT(x,y) \
   (heights[(max(0, min(width - 1, (x)))) + (max(0, min(height - 1, (y)))) * width])
#define HEIGHT_WRAP(x,y) \
   (heights[((x) < 0 ? (width + (x)) : ((x) >= width ? ((x) - width) : (x)))+ \
            (((y) < 0 ? (height + (y)) : ((y) >= height ? ((y) - height) : (y))) * width)])

   for(y = 0; y < height; ++y)
   {
      for(x = 0; x < width; ++x)
      {
         d = dst + ((y * rowbytes) + (x * bpp));
         s = src + ((y * rowbytes) + (x * bpp));

         if(nmapvals.conversion == CONVERT_NORMALIZE_ONLY ||
            nmapvals.conversion == CONVERT_HEIGHTMAP)
         {
            n[0] = (((float)s[0] * oneover255) - 0.5f) * 2.0f;
            n[1] = (((float)s[1] * oneover255) - 0.5f) * 2.0f;
            n[2] = (((float)s[2] * oneover255) - 0.5f) * 2.0f;
            n[0] *= nmapvals.scale;
            n[1] *= nmapvals.scale;
         }
         else if(nmapvals.conversion == CONVERT_DUDV_TO_NORMAL)
         {
            n[0] = (((float)s[0] * oneover255) - 0.5f) * 2.0f;
            n[1] = (((float)s[1] * oneover255) - 0.5f) * 2.0f;
            n[2] = sqrtf(1.0f - (n[0] * n[0] - n[1] * n[1]));
            n[0] *= nmapvals.scale;
            n[1] *= nmapvals.scale;
         }
         else
         {
            du = 0; dv = 0;
            if(!nmapvals.wrap)
            {
               for(i = 0; i < num_elements; ++i)
                  du += HEIGHT(x + kernel_du[i].x,
                               y + kernel_du[i].y) * kernel_du[i].w;
               for(i = 0; i < num_elements; ++i)
                  dv += HEIGHT(x + kernel_dv[i].x,
                               y + kernel_dv[i].y) * kernel_dv[i].w;
            }
            else
            {
               for(i = 0; i < num_elements; ++i)
                  du += HEIGHT_WRAP(x + kernel_du[i].x,
                                    y + kernel_du[i].y) * kernel_du[i].w;
               for(i = 0; i < num_elements; ++i)
                  dv += HEIGHT_WRAP(x + kernel_dv[i].x,
                                    y + kernel_dv[i].y) * kernel_dv[i].w;
            }
            
            n[0] = -du * nmapvals.scale;
            n[1] = -dv * nmapvals.scale;
            n[2] = 1.0f;
         }
         
         NORMALIZE(n);
         
         if(n[2] < nmapvals.minz)
         {
            n[2] = nmapvals.minz;
            NORMALIZE(n);
         }
         
         if(nmapvals.xinvert) n[0] = -n[0];
         if(nmapvals.yinvert) n[1] = -n[1];
         if(nmapvals.swapRGB)
         {
            val = n[0];
            n[0] = n[2];
            n[2] = val;
         }
         
         if(!nmapvals.dudv)
         {
            *d++ = (unsigned char)((n[0] + 1.0f) * 127.5f);
            *d++ = (unsigned char)((n[1] + 1.0f) * 127.5f);
            *d++ = (unsigned char)((n[2] + 1.0f) * 127.5f);
         
            if(bpp == 4)
            {
               switch(nmapvals.alpha)
               {
                  case ALPHA_NONE:
                     *d++ = s[3]; break;
                  case ALPHA_HEIGHT:
                     *d++ = (unsigned char)(heights[x + y * width] * 255.0f); break;
                  case ALPHA_INVERSE_HEIGHT:
                     *d++ = 255 - (unsigned char)(heights[x + y * width] * 255.0f); break;
                  case ALPHA_ZERO:
                     *d++ = 0; break;
                  case ALPHA_ONE:
                     *d++ = 255; break;
                  case ALPHA_INVERT:
                     *d++ = 255 - s[3]; break;
                  case ALPHA_MAP:
                     *d++ = sample_alpha_map(amap, x, y, amap_w, amap_h,
                                             width, height); break;
                  default:
                     *d++ = s[3]; break;
               }
            }
         }
         else
         {
            if(nmapvals.dudv == DUDV_8BIT_SIGNED ||
               nmapvals.dudv == DUDV_8BIT_UNSIGNED)
            {
               if(nmapvals.dudv == DUDV_8BIT_UNSIGNED)
               {
                  n[0] += 1.0f;
                  n[1] += 1.0f;
               }
               *d++ = (unsigned char)(n[0] * 127.5f);
               *d++ = (unsigned char)(n[1] * 127.5f);
               *d++ = 0;
               //if(drawable->bpp == 4) *d++ = 255;
            }
            else if(nmapvals.dudv == DUDV_16BIT_SIGNED ||
                    nmapvals.dudv == DUDV_16BIT_UNSIGNED)
            {
               unsigned short *d16 = (unsigned short*)d;
               if(nmapvals.dudv == DUDV_16BIT_UNSIGNED)
               {
                  n[0] += 1.0f;
                  n[1] += 1.0f;
               }
               *d16++ = (unsigned short)(n[0] * 32767.5f);
               *d16++ = (unsigned short)(n[1] * 32767.5f);
            }
         }
      }

   }
   
   if(nmapvals.conversion == CONVERT_HEIGHTMAP)
      make_heightmap(dst, width, height, bpp, nmapvals);
   
#undef HEIGHT
#undef HEIGHT_WRAP

   ImportImagePixels(nm_img, 0, 0, width, height, "RGBA", CharPixel, dst);

   free(heights);
   free(src);
   free(dst);
   free(kernel_du);
   free(kernel_dv);
   if(amap) free(amap);
   
   return(0);
}
