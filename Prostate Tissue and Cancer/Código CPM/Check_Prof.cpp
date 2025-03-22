int CPM::tumor_area () {    // calculate tumor surface area
    int pixel_index;
    int tumor_area = 0;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                pixel_index = i*Ny*Nz + j*Nz + k;
                if (PIXELS[pixel_index].tag->tau->index == 5) {
                    int neighbours_non_tumoral = 0;
                    if (i > 0 && PIXELS[pixel_index - Ny*Nz].tag->tau->index != 5) {
						neighbours_non_tumoral++;
                    }
                    if (i < Nx-1 && PIXELS[pixel_index + Ny*Nz].tag->tau->index != 5) {
                        neighbours_non_tumoral++;
                    }
                    if (j > 0 && PIXELS[pixel_index - Nz].tag->tau->index != 5) {
                        neighbours_non_tumoral++;
                    }
                    if (j < Ny-1 && PIXELS[pixel_index + Nz].tag->tau->index != 5) {
                        neighbours_non_tumoral++;
                    }
                    if (k > 0 && PIXELS[pixel_index - 1].tag->tau->index != 5) {
                        neighbours_non_tumoral++;
                    }
                    if (k < Nz-1 && PIXELS[pixel_index + 1].tag->tau->index != 5) {
                        neighbours_non_tumoral++;
                    }
					if (neighbours_non_tumoral != 0) {
						tumor_area++;
					}
                }
            }
        }
    }
    return tumor_area;
}