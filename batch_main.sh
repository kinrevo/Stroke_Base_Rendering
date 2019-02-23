gcc -Wall -omake_batch make_batch.c
gcc -O2 -omulti_file_CIB -lm ImageIO/png.c ImageIO/image.c ImageIO/jpeg.c other_func.c sbr.c sbr_main.c -lpng -ljpeg
./make_batch $1
bash make_dataset_batch.sh