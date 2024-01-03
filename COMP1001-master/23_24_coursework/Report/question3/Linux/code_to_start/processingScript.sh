echo Script running
gcc image_processing.c   -o p -O3  -lm
./image_processingCOMP ./input_images/a15.pgm ./output_images/blurred.pgm ./output_images/edge_detection.pgm

