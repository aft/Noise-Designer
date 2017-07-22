var c = instance_find(controller, 0);

if (instance_exists(c)) {
    c.objSeed.val = real(c.seedInput);
    global.seed = c.objSeed.val;
    simplex_set_seed(global.seed);
}
