//Generates the perlin noise
if (surface_exists(global.surPN)) surface_free(global.surPN); // memleak fix
if (surface_exists(global.surSN)) surface_free(global.surSN); // memleak fix
if (surface_exists(global.surIN)) surface_free(global.surIN); // memleak fix
if (surface_exists(global.surNN)) surface_free(global.surNN); // memleak fix


// perlin noise graph
global.surPN = surface_create(global.graphWidth, 100);
surface_set_target(global.surPN);

draw_set_colour(c_white);
draw_set_alpha(1);
global.arrayPN[0] = 0;
for (var i = 0; i< global.graphWidth; i++) {
    // calc starts
    draw_set_colour(c_white);
    var yy = 50+50*PN_1D_perlinNoise(i+x,
        global.seed,
        global.persistence,
        global.octaves, 
        global.wavelength,
        global.scale);
    global.arrayPN[@ i] = yy;
    draw_point(i, yy);
    // end of the calculation
}
surface_reset_target();

// smoothed noise graph
global.surSN = surface_create(global.graphWidth, 100);
surface_set_target(global.surSN);

draw_set_colour(c_white);
draw_set_alpha(.5);

for (var i = 0; i< global.graphWidth; i++) {
    // calc starts
    draw_set_colour(c_white);
    var yy = 50+50*PN_1D_smoothedNoise(i+x,
        global.seed,
        global.wavelength,
        global.scale);
    draw_point(i, yy);
    // end of the calculation
}
surface_reset_target();

// interpolated noise graph
global.surIN = surface_create(global.graphWidth, 100);
surface_set_target(global.surIN);

draw_set_colour(c_white);
draw_set_alpha(.5);

for (var i = 0; i< global.graphWidth; i++) {
    // calc starts
    draw_set_colour(c_white);
    var yy = 50+50*PN_1D_interpolatedNoise(i+x,
        global.seed,
        global.wavelength,
        global.scale);
    draw_point(i, yy);
    // end of the calculation
}
surface_reset_target();

// noise graph
global.surNN = surface_create(global.graphWidth, 100);
surface_set_target(global.surNN);

draw_set_colour(c_white);
draw_set_alpha(1);
for (var i = 0; i< global.graphWidth; i++) {
    // calc starts
    draw_set_colour(c_white);
    var yy = 50+50*PN_1D_noise(i+x, global.seed);
    draw_point(i, yy);
    // end of the calculation
}

surface_reset_target();

