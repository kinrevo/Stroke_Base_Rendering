gcc -owater -lm ImageIO/png.c ImageIO/image.c ImageIO/jpeg.c other_func.c sbr.c -lpng -ljpeg -O2 water.c water_main.c