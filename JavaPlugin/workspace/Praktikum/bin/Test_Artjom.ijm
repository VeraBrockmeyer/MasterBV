open("D:\\Documents\\GitHub\\MasterBV\\JavaPlugin\\workspace\\18mm_11-1.tif");
open("D:\\Documents\\GitHub\\MasterBV\\JavaPlugin\\workspace\\Praktikum\\Landmarks_complete.txt");
run("point grid radial affin distor ");
run("Images to Stack", "name=Stack title=[] use");