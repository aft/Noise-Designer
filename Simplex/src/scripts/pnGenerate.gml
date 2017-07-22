//Generates the perlin noise
if (surface_exists(global.surSN1d)) surface_free(global.surSN1d); // memleak fix
if (surface_exists(global.surSN2d)) surface_free(global.surSN2d); // memleak fix



// perlin noise graph
global.surSN1d = surface_create(global.graphWidth, 100);
surface_set_target(global.surSN1d);

draw_set_colour(c_white);
draw_set_alpha(1);
global.arrayPN[0] = 0;
for (var i = 0; i< global.graphWidth; i++) {
    // calc starts
    draw_set_colour(c_white);
    var yy = 100-simplex_calculate_3d(i+controller.x, 0, controller.z,
        global.minVal,
        global.maxVal,
        global.octaves, 
        global.persistence,
        global.scale/100);
    global.arrayPN[@ i] = 100-yy;
    draw_point(i, yy);
    // end of the calculation
}
surface_reset_target();

// smoothed noise graph
global.surSN2d = surface_create(global.graphWidth, global.graphHeight);
surface_set_target(global.surSN2d);

draw_set_colour(c_white);
draw_set_alpha(1);

for (var i = 0; i< global.graphWidth; i+=global.resolution) {
    // calc starts
    for (var j=0; j < global.graphHeight; j+=global.resolution) {
        var c = simplex_calculate_3d(i+controller.x, j, controller.z, 
                            global.minVal,
                            global.maxVal,
                            global.octaves, 
                            global.persistence,
                            global.scale/100);
        c = inverseLerp(global.minVal, global.maxVal, c);
        c = (c*255);
        draw_set_colour(make_colour_rgb(c, c, c));
        draw_rectangle(i, j, i+global.resolution, j+global.resolution, false);
    }

    // end of the calculation
}
surface_reset_target();
