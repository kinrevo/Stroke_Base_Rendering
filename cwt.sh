gcc -owater ImageIO/png.c ImageIO/image.c ImageIO/jpeg.c other_func.c sbr.c water.c -lm -lpng -ljpeg -O2 -Wall -Wno-unknown-pragmass -Wno-unused-result water_test.c