#define __simplex_init
/*
    Initialize the system variables and set some defaults
    
    Copyright (C) 2015 Reuben Shea
    This program is NOT free software: However, you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
globalvar __simplex_g_seed,
          __simplex_g_hash
          __simplex_g_grad3,
          __simplex_g_grad4,
          __simplex_g_simplex4,
          __simplex_g_storedHashes;
__simplex_g_grad3 = __simplex_grad3();
__simplex_g_grad4 = __simplex_grad4();
__simplex_g_simplex4 = __simplex_simplex4();
__simplex_g_storedHashes = ds_map_create();
simplex_set_seed(random_get_seed());
#define __simplex_generate_hash
///simplex_generate_hash()
/*
    Generates a 512-size hash with 256 uniqe values,
    each repeated once (to avoid out-of-bounds.
    
    Returns -   Array of hash values
 */

if (ds_map_exists(__simplex_g_storedHashes, __simplex_g_seed))
    return __simplex_g_storedHashes[? __simplex_g_seed];
 
if (ds_map_size(__simplex_g_storedHashes) >= 16)
    ds_map_delete(__simplex_g_storedHashes, ds_map_find_first(__simplex_g_storedHashes)); 

var __final;
__final[512] = 0;

for (var i = 0; i < 256; ++i)
    __final[i] = i;

//Randomize hash by swapping:
//Use different seeds for different simplex results:
for (var i = 0; i < 256; ++i)
{
    var __j = irandom(255),
        __s = __final[i];
    __final[i] = __final[__j];
    __final[__j] = __s;
}

for (var i = 0; i < 256; ++i)
    __final[255 + i] = __final[i];

__simplex_g_storedHashes[? __simplex_g_seed] = __final;
    
return __final;

#define __simplex_raw2
/*
    Argument0   -   x position
    Argument1   -   y position
    Argument4   -   minimum range of final value
    Argument5   -   maximum range of final value

    Copyright (C) 2015 Reuben Shea
    This program is NOT free software: However, you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
var __result = 0,
    //Noise contributions from the three corners:
    __n0, __n1, __n2,
    //Skew input space to determine current simplex cell
    __f2 = 0.5 * (sqrt(3.0) - 1.0);
    //Hairy factor for 2D
var __s = (argument0 + argument1) * __f2,
    __i = floor(argument0 + __s), // Treat as int
    __j = floor(argument1 + __s), // Treat as int
    __g2 = (3.0 - sqrt(3.0)) / 6.0;

var __t = (__i + __j) * __g2,
//Unskew cell origin back to x / y
    __X0 = __i - __t,
    __Y0 = __j - __t;
    
//x / y distances from the cell origin
var __x0 = argument0 - __X0,
    __y0 = argument1 - __Y0;

var __i1, __j1;
if (__x0 > __y0)
{
    __i1 = 1;
    __j1 = 0;
}
else
{
    __i1 = 0;
    __j1 = 1;
}
    
var __x1 = __x0 - __i1 + __g2,
    __y1 = __y0 - __j1 + __g2,
    __x2 = __x0 - 1.0 + 2.0 * __g2,
    __y2 = __y0 - 1.0 + 2.0 * __g2;

//Work out hashed gradient indices of the three simplex corners:
var __ii = __i & 255,
    __jj = __j & 255;
var __gi0 = __simplex_g_hash[@ __ii + __simplex_g_hash[@ __jj]] % 12,
    __gi1 = __simplex_g_hash[@ __ii + __i1 + __simplex_g_hash[@ __jj + __j1]] % 12,
    __gi2 = __simplex_g_hash[@ __ii + 1 + __simplex_g_hash[@ __jj + 1]] % 12

//Calculate contribution of ea. corner:
var __t0 = 0.5 - sqr(__x0) - sqr(__y0);
if (__t0 < 0)
    __n0 = 0;
else
{
    __t0 *= __t0;
    __n0 = sqr(__t0) * __simplex_dot2(__simplex_g_grad3[@ __gi0], __x0, __y0);
}

var __t1 = 0.5 - sqr(__x1) - sqr(__y1);
if (__t1 < 0)
    __n1 = 0;
else
{
    __t1 *= __t1;
    __n1 = sqr(__t1) * __simplex_dot2(__simplex_g_grad3[@ __gi1], __x1, __y1);
}

var __t2 = 0.5 - sqr(__x2) - sqr(__y2);
if (__t2 < 0)
    __n2 = 0.0;
else
{
    __t2 *= __t2;
    __n2 = sqr(__t2) * __simplex_dot2(__simplex_g_grad3[@ __gi2], __x2, __y2);
}

//Scale result between [-1..1]
__result = 70 * (__n0 + __n1 + __n2);
 
//Scale between whatever we like:
return __result * (argument3 - argument2) / 2 + (argument3 + argument2) / 2;

#define __simplex_raw3
///simplex_raw3(x, y, z, min, max)

/*
    Calculates the simplex noise for a specified position.
    Assumes a size of 256!
    
    Argument0   -   x position
    Argument1   -   y position
    Argument2   -   z position
    Argument3   -   minimum range of final value
    Argument4   -   maximum range of final value
    
    Copyright (C) 2015 Reuben Shea
    This program is NOT free software: However, you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
var __n0, __n1, __n2, __n3; // Noise of the four corners

//Skew input space to determine which cell we are in:
var __F3 = 1.0 / 3.0;
var __s = (argument0 + argument1 + argument2) * __F3;
var __i = floor(argument0 + __s),
    __j = floor(argument1 + __s),
    __k = floor(argument2 + __s);

var __G3 = 1.0 / 6.0; // Unskew factor
var __t = (__i + __j + __k) * __G3;
var __X0 = __i - __t, //Unskey origin back to x, y, z
    __Y0 = __j - __t,
    __Z0 = __k - __t;
var __x0 = argument0 - __X0,
    __y0 = argument1 - __Y0,
    __z0 = argument2 - __Z0;

var __i1, __j1, __k1,
    __i2, __j2, __k2;
    
if (__x0 >= __y0)
{
    if (__y0 >= __z0)
    {__i1 = 1; __j1 = 0; __k1 = 0; __i2 = 1; __j2 = 1; __k2 = 0; }
    else if (__x0 >= __z0)
    {__i1 = 1; __j1 = 0; __k1 = 0; __i2 = 1; __j2 = 0; __k2 = 1; }
    else
    {__i1 = 0; __j1 = 0; __k1 = 1; __i2 = 1; __j2 = 0; __k2 = 1; }
}
else
{
    if (__y0 < __z0)
    {__i1 = 0; __j1 = 0; __k1 = 1; __i2 = 0; __j2 = 1; __k2 = 1; }
    else if (__x0 < __z0)
    {__i1 = 0; __j1 = 1; __k1 = 0; __i2 = 0; __j2 = 1; __k2 = 1; }
    else
    {__i1 = 0; __j1 = 1; __k1 = 0; __i2 = 1; __j2 = 1; __k2 = 0; }
}

var __x1 = __x0 - __i1 + __G3,
    __y1 = __y0 - __j1 + __G3,
    __z1 = __z0 - __k1 + __G3,
    __x2 = __x0 - __i2 + 2.0 * __G3,
    __y2 = __y0 - __j2 + 2.0 * __G3,
    __z2 = __z0 - __k2 + 2.0 * __G3,
    __x3 = __x0 - 1.0 + 3.0 * __G3,
    __y3 = __y0 - 1.0 + 3.0 * __G3,
    __z3 = __z0 - 1.0 + 3.0 * __G3;

var __ii = __i & 255,
    __jj = __j & 255, 
    __kk = __k & 255;
var __gi0 = __simplex_g_hash[@ __ii + __simplex_g_hash[@ __jj + __simplex_g_hash[@ __kk] ]] % 12,
    __gi1 = __simplex_g_hash[@ __ii + __i1 + __simplex_g_hash[@ __jj + __j1 + __simplex_g_hash[@ __kk + __k1] ]] % 12,
    __gi2 = __simplex_g_hash[@ __ii + __i2 + __simplex_g_hash[@ __jj + __j2 + __simplex_g_hash[@ __kk + __k2] ]] % 12,
    __gi3 = __simplex_g_hash[@ __ii + 1.0 + __simplex_g_hash[@ __jj + 1.0 + __simplex_g_hash[@ __kk + 1.0] ]] % 12;

var __t0 = 0.6 - sqr(__x0) - sqr(__y0) - sqr(__z0);
if (__t0 < 0)
    __n0= 0.0;
else
{
    __t0 *= __t0;
    __n0 = sqr(__t0) * __simplex_dot3(__simplex_g_grad3[@ __gi0], __x0, __y0, __z0);
}

var __t1 = 0.6 - sqr(__x1) - sqr(__y1) - sqr(__z1);
if (__t1 < 0)
    __n1= 0.0;
else
{
    __t1 *= __t1;
    __n1 = sqr(__t1) * __simplex_dot3(__simplex_g_grad3[@ __gi1], __x1, __y1, __z1);
}

var __t2 = 0.6 - sqr(__x2) - sqr(__y2) - sqr(__z2);
if (__t2 < 0)
    __n2= 0.0;
else
{
    __t2 *= __t2;
    __n2 = sqr(__t2) * __simplex_dot3(__simplex_g_grad3[@ __gi2], __x2, __y2, __z2);
}

var __t3 = 0.6 - sqr(__x3) - sqr(__y3) - sqr(__z3);
if (__t3 < 0)
    __n3= 0.0;
else
{
    __t3 *= __t3;
    __n3 = sqr(__t3) * __simplex_dot3(__simplex_g_grad3[@ __gi3], __x3, __y3, __z3);
}

var __result = 32.0 * (__n0 + __n1 + __n2 + __n3);
return  __result * (argument4 - argument3) / 2 + (argument4 + argument3) / 2;

#define __simplex_raw4
///simplex_raw4(x, y, z, w, [ARRAY 1D:INT] hash, [ARRAY 1D:1D:INT] gradient, [ARRAY 1D:1D:INT] simplex, min, max)
/*
    Calculates the simplex noise for a specified position.
    Assumes a size of 256!
    
    Argument0   -   x position
    Argument1   -   y position
    Argument2   -   z position
    Argument3   -   w position
    Argument4   -   minimum range of final value
    Argument5   -   maximum range of final value

    Copyright (C) 2015 Reuben Shea
    This program is NOT free software: However, you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


//Skewing / Unskewing factors:
var F4 = (sqrt(5.0) - 1.0) / 4.0,
    G4 = (5.0 - sqrt(5.0)) / 20.0;
var n0, n1, n2, n3, n4;

//Skew 4D space to determine which cell we are in:
var s = (argument0 + argument1 + argument2 + argument3) * F4;
var i = floor(argument0 + s),
    j = floor(argument1 + s),
    k = floor(argument2 + s),
    l = floor(argument3 + s);
var t = (i + j + k + l) * G4;
var X0 = i - t,
    Y0 = j - t,
    Z0 = k - t,
    W0 = l - t;
var x0 = argument0 - X0,
    y0 = argument1 - Y0,
    z0 = argument2 - Z0,
    w0 = argument3 - W0;
    

var c1 = (x0 > y0) * 32,
    c2 = (x0 > z0) * 16,
    c3 = (y0 > z0) * 8,
    c4 = (x0 > w0) * 4,
    c5 = (y0 > w0) * 2,
    c6 = (z0 > w0);
var c = c1 + c2 + c3 + c4 + c5 + c6;

var i1, j1, k1, l1, //Second simplex corner offset
    i2, j2, k2, l2, //Third     "       "
    i3, j3, k3, l3; //Fourth    "       "

var __a = __simplex_g_simplex4[@ c];
i1 = (__a[@ 0] >= 3);
j1 = (__a[@ 1] >= 3);
k1 = (__a[@ 2] >= 3);
l1 = (__a[@ 3] >= 3);

i2 = (__a[@ 0] >= 2);
j2 = (__a[@ 1] >= 2);
k2 = (__a[@ 2] >= 2);
l2 = (__a[@ 3] >= 2);

i3 = (__a[@ 0] >= 1);
j3 = (__a[@ 1] >= 1);
k3 = (__a[@ 2] >= 1);
l3 = (__a[@ 3] >= 1);

//Corner offsets in x,y,z,w coords:
var __G42 = G4 *2,
    __G43 = G4 *4;
var x1 = x0 - i1 + G4,
    y1 = y0 - j1 + G4,
    z1 = z0 - k1 + G4,
    w1 = w0 - l1 + G4,
    
    x2 = x0 - i2 + __G42,
    y2 = y0 - j2 + __G42,
    z2 = z0 - k2 + __G42,
    w2 = w0 - l2 + __G42,
    
    x3 = x0 - i3 + G4 * 3.0,
    y3 = y0 - j3 + G4 * 3.0,
    z3 = z0 - k3 + G4 * 3.0,
    w3 = w0 - l3 + G4 * 3.0,
    
    x4 = x0 - 1.0 + __G43,
    y4 = y0 - 1.0 + __G43,
    z4 = z0 - 1.0 + __G43,
    w4 = w0 - 1.0 + __G43;
//Hashed gradient indices:
var ii = i & 255,
    jj = j & 255,
    kk = k & 255,
    ll = l & 255;
    
var gi0 = __simplex_g_hash[@ ii + __simplex_g_hash[@ jj + __simplex_g_hash[@ kk + __simplex_g_hash[@ ll]]]] & 31,
    gi1 = __simplex_g_hash[@ ii + i1 + __simplex_g_hash[@ jj + j1 + __simplex_g_hash[@ kk + k1 + __simplex_g_hash[@ ll + l1]]]] & 31,
    gi2 = __simplex_g_hash[@ ii + i2 + __simplex_g_hash[@ jj + j2 + __simplex_g_hash[@ kk + k2 + __simplex_g_hash[@ ll + l2]]]] & 31,
    gi3 = __simplex_g_hash[@ ii + i3 + __simplex_g_hash[@ jj + j3 + __simplex_g_hash[@ kk + k3 + __simplex_g_hash[@ ll + l3]]]] & 31,
    gi4 = __simplex_g_hash[@ ii + 1.0 + __simplex_g_hash[@ jj + 1.0 + __simplex_g_hash[@ kk + 1.0 + __simplex_g_hash[@ ll + 1.0]]]] & 31;

//Calculate contribution for ea. corner:

t0 = 0.6 - sqr(x0) - sqr(y0) - sqr(z0) - sqr(w0);

if (t0 < 0)
    n0 = 0.0;
else
{
    t0 *= t0;
    n0 = sqr(t0) * __simplex_dot4(__simplex_g_grad4[@ gi0], x0, y0, z0, w0);
}

var t1 = 0.6 - sqr(x1) - sqr(y1) - sqr(z1) - sqr(w1);
if (t1 < 0)
    n1 = 0.0;
else
{
    t1*= t1;
    n1 = sqr(t1) * __simplex_dot4(__simplex_g_grad4[@ gi1], x1, y1, z1, w1);
}

var t2 = 0.6 - sqr(x2) - sqr(y2) - sqr(z2) - sqr(w2);
if (t2 < 0)
    n2 = 0.0;
else
{
    t2 *= t2;
    n2 = sqr(t2) * __simplex_dot4(__simplex_g_grad4[@ gi2], x2, y2, z2, w2);
}

var t3 = 0.6 - sqr(x3) - sqr(y3) - sqr(z3) - sqr(w3);
if (t3 < 0)
    n3 = 0.0;
else
{
    t3 *= t3;
    n3 = sqr(t3) * __simplex_dot4(__simplex_g_grad4[@ gi3], x3, y3, z3, w3);
}

var t4 = 0.6 - sqr(x4) - sqr(y4) - sqr(z4) - sqr(w4);
if (t4 < 0)
    n4 = 0.0;
else
{
    t4 *= t4;
    n4 = sqr(t4) * __simplex_dot4(__simplex_g_grad4[@ gi4], x4, y4, z4, w4);
}
    

return (27.0 * (n0 + n1 + n2 + n3 + n4)) * (argument5 - argument4) / 2 + (argument5 + argument4) / 2;
    

#define __simplex_grad3
/*
    Generates the 3D gradients
    
    Returns -   Array of arrays of points.
 */
var __a0, __a1, __a2, __a3, __a4, __a5,
    __a6, __a7, __a8, __a9, __a10, __a11;
    
__a0[0] = 1;    __a0[1] = 1;    __a0[2] = 0;
__a1[0] = -1;   __a1[1] = 1;    __a1[2] = 0;
__a2[0] = 1;    __a2[1] = -1;   __a2[2] = 0;
__a3[0] = -1;   __a3[1] = -1;   __a3[2] = 0;
__a4[0] = 1;    __a4[1] = 0;    __a4[2] = 1;
__a5[0] = -1;   __a5[1] = 0;    __a5[2] = 1;
__a6[0] = 1;    __a6[1] = 0;    __a6[2] = -1;
__a7[0] = -1;   __a7[1] = 0;    __a7[2] = -1;
__a8[0] = 0;    __a8[1] = 1;    __a8[2] = 1;
__a9[0] = 0;    __a9[1] = -1;   __a9[2] = 1;
__a10[0] = 0;   __a10[1] = 1;   __a10[2] = -1;
__a11[0] = 0;   __a11[1] = -1;  __a11[2] = -1;

var __return;
__return[11] = __a11;
__return[10] = __a10;
__return[9] = __a9;
__return[8] = __a8;
__return[7] = __a7;
__return[6] = __a6;
__return[5] = __a5;
__return[4] = __a4;
__return[3] = __a3;
__return[2] = __a2;
__return[1] = __a1;
__return[0] = __a0;

return __return;

#define __simplex_grad4
/*
    Generates the 4D gradients
    
    Returns -   Array of arrays of points.
 */
var __a0, __a1, __a2, __a3, __a4, __a5,
    __a6, __a7, __a8, __a9, __a10, __a11,
    __a12, __a13, __a14, __a15, __a16, __a17,
    __a18, __a19, __a20, __a21, __a22, __a23,
    __a24, __a25, __a26, __a27, __a28, __a29,
    __a30, __a31;
    
__a0[0] = 0;     __a0[1] = 1;     __a0[2] = 1;     __a0[3] = 1;
__a1[0] = 0;     __a1[1] = 1;     __a1[2] = 1;     __a1[3] = -1;
__a2[0] = 0;     __a2[1] = 1;     __a2[2] = -1;    __a2[3] = 1;
__a3[0] = 0;     __a3[1] = 1;     __a3[2] = -1;    __a3[3] = -1;
__a4[0] = 0;     __a4[1] = -1;    __a4[2] = 1;     __a4[3] = 1;
__a5[0] = 0;     __a5[1] = -1;    __a5[2] = 1;     __a5[3] = -1;
__a6[0] = 0;     __a6[1] = -1;    __a6[2] = -1;    __a6[3] = 1;
__a7[0] = 0;     __a7[1] = -1;    __a7[2] = -1;    __a7[3] = -1;
__a8[0] = 1;     __a8[1] = 0;     __a8[2] = 1;     __a8[3] = 1;
__a9[0] = 1;     __a9[1] = 0;     __a9[2] = 1;     __a9[3] = -1;
__a10[0] = 1;    __a10[1] = 0;    __a10[2] = -1;   __a10[3] = 1;
__a11[0] = 1;    __a11[1] = 0;    __a11[2] = -1;   __a11[3] = -1;
__a12[0] = -1;   __a12[1] = 0;    __a12[2] = 1;    __a12[3] = 1;
__a13[0] = -1;   __a13[1] = 0;    __a13[2] = 1;    __a13[3] = -1;
__a14[0] = -1;   __a14[1] = 0;    __a14[2] = -1;   __a14[3] = 1;
__a15[0] = -1;   __a15[1] = 0;    __a15[2] = -1;   __a15[3] = -1;
__a16[0] = 1;    __a16[1] = 1;    __a16[2] = 0;    __a16[3] = 1;
__a17[0] = 1;    __a17[1] = 1;    __a17[2] = 0;    __a17[3] = -1;
__a18[0] = 1;    __a18[1] = -1;   __a18[2] = 0;    __a18[3] = 1;
__a19[0] = 1;    __a19[1] = -1;   __a19[2] = 0;    __a19[3] = -1;
__a20[0] = -1;   __a20[1] = 1;    __a20[2] = 0;    __a20[3] = 1;
__a21[0] = -1;   __a21[1] = 1;    __a21[2] = 0;    __a21[3] = -1;
__a22[0] = -1;   __a22[1] = -1;   __a22[2] = 0;    __a22[3] = 1;
__a23[0] = -1;   __a23[1] = -1;   __a23[2] = 0;    __a23[3] = -1;
__a24[0] = 1;    __a24[1] = 1;    __a24[2] = 1;    __a24[3] = 0;
__a25[0] = 1;    __a25[1] = 1;    __a25[2] = -1;   __a25[3] = 0;
__a26[0] = 1;    __a26[1] = -1;   __a26[2] = 1;    __a26[3] = 0;
__a27[0] = 1;    __a27[1] = -1;   __a27[2] = -1;   __a27[3] = 0;
__a28[0] = -1;   __a28[1] = 1;    __a28[2] = 1;    __a28[3] = 0;
__a29[0] = -1;   __a29[1] = 1;    __a29[2] = -1;   __a29[3] = 0;
__a30[0] = -1;   __a30[1] = -1;   __a30[2] = 1;    __a30[3] = 0;
__a31[0] = -1;   __a31[1] = -1;   __a31[2] = -1;   __a31[3] = 0;


var __return;
__return[31] = __a31;
__return[30] = __a30;
__return[29] = __a29;
__return[28] = __a28;
__return[27] = __a27;
__return[26] = __a26;
__return[25] = __a25;
__return[24] = __a24;
__return[23] = __a23;
__return[22] = __a22;
__return[21] = __a21;
__return[20] = __a20;
__return[19] = __a19;
__return[18] = __a18;
__return[17] = __a17;
__return[16] = __a16;
__return[15] = __a15;
__return[14] = __a14;
__return[13] = __a13;
__return[12] = __a12;
__return[11] = __a11;
__return[10] = __a10;
__return[9] = __a9;
__return[8] = __a8;
__return[7] = __a7;
__return[6] = __a6;
__return[5] = __a5;
__return[4] = __a4;
__return[3] = __a3;
__return[2] = __a2;
__return[1] = __a1;
__return[0] = __a0;

return __return;

#define __simplex_simplex4
/*
    Returns -   Array of arrays of points.
 */
var __a0, __a1, __a2, __a3, __a4, __a5,
    __a6, __a7, __a8, __a9, __a10, __a11,
    __a12, __a13, __a14, __a15, __a16, __a17,
    __a18, __a19, __a20, __a21, __a22, __a23,
    __a24, __a25, __a26, __a27, __a28, __a29,
    __a30, __a31,
    __a32, __a33, __a34, __a35, __a36, __a37,
    __a38, __a39, __a40, __a41, __a42, __a43,
    __a44, __a45, __a46, __a47, __a48, __a49,
    __a50, __a51, __a52, __a53, __a54, __a55,
    __a56, __a57, __a58, __a59, __a60, __a61,
    __a62, __a63;
    

__a0 [0] = 0;   __a0 [1] = 1;   __a0 [2] = 2;   __a0[3] = 3;
__a1 [0] = 0;   __a1 [1] = 1;   __a1 [2] = 3;   __a1[3] = 2;
__a2 [0] = 0;   __a2 [1] = 0;   __a2 [2] = 0;   __a2[3] = 0;
__a3 [0] = 0;   __a3 [1] = 2;   __a3 [2] = 3;   __a3[3] = 1;
__a4 [0] = 0;   __a4 [1] = 0;   __a4 [2] = 0;   __a4[3] = 0;
__a5 [0] = 0;   __a5 [1] = 0;   __a5 [2] = 0;   __a5[3] = 0;
__a6 [0] = 0;   __a6 [1] = 0;   __a6 [2] = 0;   __a6[3] = 0;
__a7 [0] = 1;   __a7 [1] = 2;   __a7 [2] = 3;   __a7[3] = 0;

__a8 [0] = 0;   __a8 [1] = 2;   __a8 [2] = 1;   __a8[3] = 3;
__a9 [0] = 0;   __a9 [1] = 0;   __a9 [2] = 0;   __a9[3] = 0;
__a10 [0] = 0;   __a10 [1] = 3;   __a10 [2] = 1;   __a10[3] = 2;
__a11 [0] = 0;   __a11 [1] = 3;   __a11 [2] = 2;   __a11[3] = 1;
__a12 [0] = 0;   __a12 [1] = 0;   __a12 [2] = 0;   __a12[3] = 0;
__a13 [0] = 0;   __a13 [1] = 0;   __a13 [2] = 0;   __a13[3] = 0;
__a14 [0] = 0;   __a14 [1] = 0;   __a14 [2] = 0;   __a14[3] = 0;
__a15 [0] = 1;   __a15 [1] = 3;   __a15 [2] = 2;   __a15[3] = 0;

__a16 [0] = 0;   __a16 [1] = 0;   __a16 [2] = 0;   __a16[3] = 0;
__a17 [0] = 0;   __a17 [1] = 0;   __a17 [2] = 0;   __a17[3] = 0;
__a18 [0] = 0;   __a18 [1] = 0;   __a18 [2] = 0;   __a18[3] = 0;
__a19 [0] = 0;   __a19 [1] = 0;   __a19 [2] = 0;   __a19[3] = 0;
__a20 [0] = 0;   __a20 [1] = 0;   __a20 [2] = 0;   __a20[3] = 0;
__a21 [0] = 0;   __a21 [1] = 0;   __a21 [2] = 0;   __a21[3] = 0;
__a22 [0] = 0;   __a22 [1] = 0;   __a22 [2] = 0;   __a22[3] = 0;
__a23 [0] = 0;   __a23 [1] = 0;   __a23 [2] = 0;   __a23[3] = 0;

__a24 [0] = 1;   __a24 [1] = 2;   __a24 [2] = 0;   __a24[3] = 3;
__a25 [0] = 0;   __a25 [1] = 0;   __a25 [2] = 0;   __a25[3] = 0;
__a26 [0] = 1;   __a26 [1] = 3;   __a26 [2] = 0;   __a26[3] = 2;
__a27 [0] = 0;   __a27 [1] = 0;   __a27 [2] = 0;   __a27[3] = 0;
__a28 [0] = 0;   __a28 [1] = 0;   __a28 [2] = 0;   __a28[3] = 0;
__a29 [0] = 0;   __a29 [1] = 0;   __a29 [2] = 0;   __a29[3] = 0;
__a30 [0] = 2;   __a30 [1] = 3;   __a30 [2] = 0;   __a30[3] = 1;
__a31 [0] = 2;   __a31 [1] = 3;   __a31 [2] = 1;   __a31[3] = 0;

__a32 [0] = 1;   __a32 [1] = 0;   __a32 [2] = 2;   __a32[3] = 3;
__a33 [0] = 1;   __a33 [1] = 0;   __a33 [2] = 3;   __a33[3] = 2;
__a34 [0] = 0;   __a34 [1] = 0;   __a34 [2] = 0;   __a34[3] = 0;
__a35 [0] = 0;   __a35 [1] = 0;   __a35 [2] = 0;   __a35[3] = 0;
__a36 [0] = 0;   __a36 [1] = 0;   __a36 [2] = 0;   __a36[3] = 0;
__a37 [0] = 2;   __a37 [1] = 0;   __a37 [2] = 3;   __a37[3] = 1;
__a38 [0] = 0;   __a38 [1] = 0;   __a38 [2] = 0;   __a38[3] = 0;
__a39 [0] = 2;   __a39 [1] = 1;   __a39 [2] = 3;   __a39[3] = 0;

__a40 [0] = 0;   __a40 [1] = 0;   __a40 [2] = 0;   __a40[3] = 0;
__a41 [0] = 0;   __a41 [1] = 0;   __a41 [2] = 0;   __a41[3] = 0;
__a42 [0] = 0;   __a42 [1] = 0;   __a42 [2] = 0;   __a42[3] = 0;
__a43 [0] = 0;   __a43 [1] = 0;   __a43 [2] = 0;   __a43[3] = 0;
__a44 [0] = 0;   __a44 [1] = 0;   __a44 [2] = 0;   __a44[3] = 0;
__a45 [0] = 0;   __a45 [1] = 0;   __a45 [2] = 0;   __a45[3] = 0;
__a46 [0] = 0;   __a46 [1] = 0;   __a46 [2] = 0;   __a46[3] = 0;
__a47 [0] = 0;   __a47 [1] = 0;   __a47 [2] = 0;   __a47[3] = 0;

__a48 [0] = 2;   __a48 [1] = 0;   __a48 [2] = 1;   __a48[3] = 3;
__a49 [0] = 0;   __a49 [1] = 0;   __a49 [2] = 0;   __a49[3] = 0;
__a50 [0] = 0;   __a50 [1] = 0;   __a50 [2] = 0;   __a50[3] = 0;
__a51 [0] = 0;   __a51 [1] = 0;   __a51 [2] = 0;   __a51[3] = 0;
__a52 [0] = 3;   __a52 [1] = 0;   __a52 [2] = 1;   __a52[3] = 2;
__a53 [0] = 3;   __a53 [1] = 0;   __a53 [2] = 2;   __a53[3] = 1;
__a54 [0] = 0;   __a54 [1] = 0;   __a54 [2] = 0;   __a54[3] = 0;
__a55 [0] = 3;   __a55 [1] = 1;   __a55 [2] = 2;   __a55[3] = 0;

__a56 [0] = 2;   __a56 [1] = 1;   __a56 [2] = 0;   __a56[3] = 3;
__a57 [0] = 0;   __a57 [1] = 0;   __a57 [2] = 0;   __a57[3] = 0;
__a58 [0] = 0;   __a58 [1] = 0;   __a58 [2] = 0;   __a58[3] = 0;
__a59 [0] = 0;   __a59 [1] = 0;   __a59 [2] = 0;   __a59[3] = 0;
__a60 [0] = 3;   __a60 [1] = 1;   __a60 [2] = 0;   __a60[3] = 2;
__a61 [0] = 0;   __a61 [1] = 0;   __a61 [2] = 0;   __a61[3] = 0;
__a62 [0] = 3;   __a62 [1] = 2;   __a62 [2] = 0;   __a62[3] = 1;
__a63 [0] = 3;   __a63 [1] = 2;   __a63 [2] = 1;   __a63[3] = 0;

var __return;
__return[63] = __a63;
__return[62] = __a62;
__return[61] = __a61;
__return[60] = __a60;
__return[59] = __a59;
__return[58] = __a58;
__return[57] = __a57;
__return[56] = __a56;
__return[55] = __a55;
__return[54] = __a54;
__return[53] = __a53;
__return[52] = __a52;
__return[51] = __a51;
__return[50] = __a50;
__return[49] = __a49;
__return[48] = __a48;
__return[47] = __a47;
__return[46] = __a46;
__return[45] = __a45;
__return[44] = __a44;
__return[43] = __a43;
__return[42] = __a42;
__return[41] = __a41;
__return[40] = __a40;
__return[39] = __a39;
__return[38] = __a38;
__return[37] = __a37;
__return[36] = __a36;
__return[35] = __a35;
__return[34] = __a34;
__return[33] = __a33;
__return[32] = __a32;
__return[31] = __a31;
__return[30] = __a30;
__return[29] = __a29;
__return[28] = __a28;
__return[27] = __a27;
__return[26] = __a26;
__return[25] = __a25;
__return[24] = __a24;
__return[23] = __a23;
__return[22] = __a22;
__return[21] = __a21;
__return[20] = __a20;
__return[19] = __a19;
__return[18] = __a18;
__return[17] = __a17;
__return[16] = __a16;
__return[15] = __a15;
__return[14] = __a14;
__return[13] = __a13;
__return[12] = __a12;
__return[11] = __a11;
__return[10] = __a10;
__return[9] = __a9;
__return[8] = __a8;
__return[7] = __a7;
__return[6] = __a6;
__return[5] = __a5;
__return[4] = __a4;
__return[3] = __a3;
__return[2] = __a2;
__return[1] = __a1;
__return[0] = __a0;

return __return;

#define __simplex_dot2
/*
    2D dot product.
    
    Argument0   -   1D array of 2 vector values
    Argument1   -   x val
    Argument2   -   y val
    
    Returns     -   Result
 */

return (argument0[@ 0] * argument1) + (argument0[@ 1] * argument2);

#define __simplex_dot3
/*
    3D dot product.
    
    Argument0   -   1D array of 3 vector values
    Argument1   -   x val
    Argument2   -   y val
    Argument3   -   z val
    
    Returns     -   Result
 */

return (argument0[@ 0] * argument1) + (argument0[@ 1] * argument2) + (argument0[@ 2] * argument3);

#define __simplex_dot4
/*
    4D dot product.
    
    Argument0   -   1D array of 4 vector values
    Argument1   -   x val
    Argument2   -   y val
    Argument3   -   z val
    Argument4   -   w val
    
    Returns     -   Result
 */

return (argument0[@ 0] * argument1) + (argument0[@ 1] * argument2) + (argument0[@ 2] * argument3) + (argument0[@ 3] * argument4);

#define simplex_get_seed
///simplex_get_seed()
/*
    Returns the current seed being used by the simplex noise
    algorithm.
 */

return __simplex_g_seed;

#define simplex_set_seed
///simplex_set_seed(seed)
/*
    sets the current seed for the simplex algorithm.
    
    Argument0   -   (int) the seed to use
    Returns     -   Argument0
 */
var __seed = random_get_seed();
random_set_seed(argument0)
__simplex_g_seed = argument0;
__simplex_g_hash = __simplex_generate_hash();
random_set_seed(__seed);

#define simplex_calculate_2d
///simplex_calculate_2d(x, y, min, max, octaves, persistence, scale)
/*
    Generates fractal simplex noise at the specified position.
    
    Argument0   -   x position
    Argument1   -   y position
    Argument2   -   minimum range of final value
    Argument3   -   maximum range of final value
    Argument4   -   number of samples
    Argument5   -   delta octave intensity [0..1]
    Argument6   -   scale of deltas
    
    Returns     -   Calculated result 
 
    Copyright (C) 2015 Reuben Shea
    This program is NOT free software: However, you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

var __total = 0,
    __freq = argument6,
    __amp = 1,
    __maxAmp = 0; // Will keep things between [-1..1]

for (var i = 0; i < argument4; ++i)
{
    __total += __simplex_raw2(argument0 * __freq, argument1 * __freq, -1, 1) * __amp;
    
    __freq *= 2;
    __maxAmp += __amp;
    __amp *= argument5;
}

return (__total / __maxAmp) * (argument3 - argument2) / 2 + (argument3 + argument2) / 2 ;

#define simplex_calculate_2dr
///simplex_calculate_2dr(x, y, min, max, octaves, persistence, scale)
/*
    Generates ridged fractal simplex noise at the specified position.
    
    Argument0   -   x position
    Argument1   -   y position
    Argument2   -   minimum range of final value
    Argument3   -   maximum range of final value
    Argument4   -   number of samples
    Argument5   -   delta octave intensity [0..1]
    Argument6   -   scale of deltas
    
    Returns     -   Calculated result 
 
    Copyright (C) 2015 Reuben Shea
 */

var __total = 0,
    __freq = argument6,
    __amp = 1,
    __maxAmp = 0; // Will keep things between [-1..1]

for (var i = 0; i < argument4; ++i)
{
    __total += (1 - abs(__simplex_raw2(argument0 * __freq, argument1 * __freq, -1, 1)) * __amp);
    
    __freq *= 2;
    __maxAmp += 1;
    __amp *= argument5;
}

return (__total / __maxAmp) * (argument3 - argument2) / 2 + (argument3 + argument2) / 2 ;

#define simplex_calculate_3d
///simplex_calculate_3d(x, y, z, min, max, octaves, persistence, scale)
/*
    Generates fractal simplex noise at the specified position.
    
    Argument0   -   x position
    Argument1   -   y position
    Argument2   -   z position
    Argument3   -   minimum range of final value
    Argument4   -   maximum range of final value
    Argument5   -   number of samples
    Argument6   -   delta octave intensity [0..1]
    Argument7   -   scale of deltas
    
    Returns     -   Calculated result 
   
    Copyright (C) 2015 Reuben Shea
    This program is NOT free software: However, you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

var __total = 0,
    __freq = argument7,
    __amp = 1,
    __maxAmp = 0; // Will keep things between [-1..1]

for (var i = 0; i < argument5; ++i)
{
    __total += __simplex_raw3(argument0 * __freq, argument1 * __freq, argument2 * __freq, -1, 1) * __amp;
    
    __freq *= 2;
    __maxAmp += __amp;
    __amp *= argument6;
}

return (__total / __maxAmp) * (argument4 - argument3) / 2 + (argument4 + argument3) / 2 ;

#define simplex_calculate_3dr
///simplex_calculate_3dr(x, y, z, min, max, octaves, persistence, scale)
/*
    Generates ridged fractal simplex noise at the specified position.
    
    Argument0   -   x position
    Argument1   -   y position
    Argument2   -   z position
    Argument3   -   minimum range of final value
    Argument4   -   maximum range of final value
    Argument5   -   number of samples
    Argument6   -   delta octave intensity [0..1]
    Argument7   -   scale of deltas
    
    Returns     -   Calculated result 
   
    Copyright (C) 2015 Reuben Shea
 */

var __total = 0,
    __freq = argument7,
    __amp = 1,
    __maxAmp = 0; // Will keep things between [-1..1]

for (var i = 0; i < argument5; ++i)
{
    __total += (1 - abs(__simplex_raw3(argument0 * __freq, argument1 * __freq, argument2 * __freq, -1, 1)) * __amp);
    
    __freq *= 2;
    __maxAmp += 1;
    __amp *= argument6;
}

return (__total / __maxAmp) * (argument4 - argument3) / 2 + (argument4 + argument3) / 2 ;

#define simplex_calculate_4d
///simplex_calculate_4d(x, y, z, w, min, max, octaves, persistence, scale)
/*
    Generates fractal simplex noise at the specified position.
    
    Argument0   -   x position
    Argument1   -   y position
    Argument2   -   z position
    Argument3   -   w position
    Argument4   -   minimum range of final value
    Argument5   -   maximum range of final value
    Argument6   -   number of samples
    Argument7   -   delta octave intensity [0..1]
    Argument8   -   scale of deltas
    
    Returns     -   Calculated result  


    Copyright (C) 2015 Reuben Shea
    This program is NOT free software: However, you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

var __total = 0,
    __freq = argument8,
    __amp = 1,
    __maxAmp = 0; // Will keep things between [-1..1]

for (var i = 0; i < argument6; ++i)
{
    __total += __simplex_raw4(argument0 * __freq, argument1 * __freq, 
                              argument2 * __freq, argument3 * __freq, 
                              -1, 1) * __amp;
    
    __freq *= 2;
    __maxAmp += __amp;
    __amp *= argument7;
}

return (__total / __maxAmp) * (argument5 - argument4) / 2 + (argument5 + argument4) / 2 ;

#define simplex_calculate_4dr
///simplex_calculate_4dr(x, y, z, w, min, max, octaves, persistence, scale)
/*
    Generates ridged fractal simplex noise at the specified position.
    
    Argument0   -   x position
    Argument1   -   y position
    Argument2   -   z position
    Argument3   -   w position
    Argument4   -   minimum range of final value
    Argument5   -   maximum range of final value
    Argument6   -   number of samples
    Argument7   -   delta octave intensity [0..1]
    Argument8   -   scale of deltas
    
    Returns     -   Calculated result  


    Copyright (C) 2015 Reuben Shea
 */

var __total = 0,
    __freq = argument8,
    __amp = 1,
    __maxAmp = 0; // Will keep things between [-1..1]

for (var i = 0; i < argument6; ++i)
{
    __total += (1 - abs(__simplex_raw4(argument0 * __freq, argument1 * __freq, 
                              argument2 * __freq, argument3 * __freq, 
                              -1, 1)) * __amp);
    
    __freq *= 2;
    __maxAmp += 1;
    __amp *= argument7;
}

return (__total / __maxAmp) * (argument5 - argument4) / 2 + (argument5 + argument4) / 2 ;

