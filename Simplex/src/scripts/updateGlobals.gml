var c = instance_find(controller, 0);

if (instance_exists(c)) {
    global.seed = c.objSeed.val;
    global.persistence = c.objPersistence.val;
    global.octaves = c.objOctaves.val;
    global.minVal = c.objMin.val;
    global.maxVal = c.objMax.val;
    global.scale = c.objScale.val;
}
