//resets the global seed
randomize();
global.seed = irandom_range(1000000, 10000000);
random_set_seed(global.seed);
show_debug_message(string(global.seed));
