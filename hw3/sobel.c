#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <math.h>

#include <jpeglib.h>

#include <mpi.h>

/* Filters */
static int sobel_x_filter[9]={ -1, 0,+1  ,-2, 0,+2  ,-1, 0,+1};
static int sobel_y_filter[9]={ -1,-2,-1  , 0, 0, 0  , 1, 2,+1};

/* Structure describing the image */
struct image_t {
	int xsize;
	int ysize;
	int depth;	/* bytes */
	unsigned char *pixels;	/*flattened*/
};

struct convolve_data_t {
	struct image_t *old;
	struct image_t *new;
	int *filter;
	int ystart;
	int yend;
};


// /* very inefficient convolve code */
// static void *generic_convolve(void *argument) {

// 	int x,y,k,l,d;
// 	uint32_t color;
// 	int sum,depth,width;

// 	struct image_t *old;
// 	struct image_t *new;
// 	int *filter;
// 	struct convolve_data_t *data;
// 	int ystart, yend;

// 	/* Convert from void pointer to the actual data type */
// 	data=(struct convolve_data_t *)argument;
// 	old=data->old;
// 	new=data->new;
// 	filter=data->filter;

// 	ystart=data->ystart;
// 	yend=data->yend;

// 	depth=old->depth;
// 	width=old->xsize*old->depth;

// 	/* handle border*/
// 	//3*3 kernel, shrinking size by 1 at each direction.
// 	if (ystart==0) ystart=1;
// 	if (yend==old->ysize) yend=old->ysize-1;

// 	for(d=0;d<3;d++) {
// 	   for(x=1;x<old->xsize-1;x++) {
// 	     for(y=ystart;y<yend;y++) {
// 		sum=0;
// 		for(k=-1;k<2;k++) {
// 		   for(l=-1;l<2;l++) {
// 			// convolve on elements of each color
// 			color=old->pixels[((y+l)*width)+(x*depth+d+k*depth)];
// 			sum+=color * filter[(l+1)*3+(k+1)];
// 		   }
// 		}

// 		if (sum<0) sum=0;
// 		if (sum>255) sum=255;

// 		new->pixels[(y*width)+x*depth+d]=sum;
// 	     }
// 	   }
// 	}

// 	return NULL;
// }

static void *generic_convolve(void *argument) {

	int x,y,k,l,d;
	uint32_t color;
	int sum,depth,width;

	struct image_t *old;
	struct image_t *new;
	int *filter;
	struct convolve_data_t *data;
	int ystart, yend;

	/* Convert from void pointer to the actual data type */
	data=(struct convolve_data_t *)argument;
	old=data->old;
	new=data->new;
	filter=data->filter;

	ystart=data->ystart;
	yend=data->yend;

	depth=old->depth;
	width=old->xsize*old->depth;

	/* handle border*/
	//3*3 kernel, shrinking size by 1 at each direction.
	if (ystart==0) ystart=1;
	if (yend==old->ysize) yend=old->ysize-1;

	int current_rank, num_ranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	
	//get step for current process
	int ystep = (yend-ystart) / num_ranks;
	if((yend-ystart) % num_ranks != 0) ystep++;

	//get upper bound & lower bound for current process
	int upper_bound;
	if(yend > ystart + (current_rank+1)*ystep) upper_bound = ystart + (current_rank+1)*ystep;
	else upper_bound = yend;
	int lower_bound = ystart + current_rank * ystep;
	
	// for each color
	for(d=0;d<3;d++) {
	   for(x=1;x<old->xsize-1;x++) {
	     for(y=lower_bound;y<upper_bound;y++) {
		
		sum=0;
		for(k=-1;k<2;k++) {
		   for(l=-1;l<2;l++) {
			// convolve on elements of each color
			color=old->pixels[((y+l)*width)+(x*depth+d+k*depth)];
			sum+=color * filter[(l+1)*3+(k+1)];
			
		   }
		}
		
		if (sum<0) sum=0;
		if (sum>255) sum=255;
		
		new->pixels[(y*width)+x*depth+d]=sum;
		
	     }
	   }
	}
	
	unsigned char* buf = (unsigned char*)malloc(width*ystep*num_ranks);

	
	MPI_Allgather(new->pixels + lower_bound*width, ystep*width, MPI_UNSIGNED_CHAR,
					buf, ystep*width, MPI_UNSIGNED_CHAR,
					MPI_COMM_WORLD);
	
	memcpy(new->pixels + width, buf, width * (old->ysize-2));
	return NULL;
}

//original ver.
// static int combine(struct image_t *s_x,
// 			struct image_t *s_y,
// 			struct image_t *new) {
// 	int i;
// 	int out;

// 	for(i=0;i<( s_x->depth * s_x->xsize * s_x->ysize );i++) {

// 		out=sqrt(
// 			(s_x->pixels[i]*s_x->pixels[i])+
// 			(s_y->pixels[i]*s_y->pixels[i])
// 			);
// 		if (out>255) out=255;
// 		if (out<0) out=0;
// 		new->pixels[i]=out;
// 	}

// 	return 0;
// }

static int combine(struct image_t *s_x,
			struct image_t *s_y,
			struct image_t *new) {
	int i;
	int out;

	int current_rank, num_ranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	
	int total = s_x->depth * s_x->xsize * s_x->ysize;
	int width = s_x->depth * s_x->xsize;
	int step = s_x->ysize / num_ranks;
	if(total % num_ranks != 0)step ++;
	int lower_bound = current_rank * step * width;
	int upper_bound;
	
	if(total > lower_bound + step * width)upper_bound = lower_bound + step * width;
	else upper_bound = total;

	
	for(i = lower_bound; i < upper_bound;i++){
		out=sqrt(
			(s_x->pixels[i]*s_x->pixels[i])+
			(s_y->pixels[i]*s_y->pixels[i])
			);
		if (out>255) out=255;
		if (out<0) out=0;
		new->pixels[i]=out;
	}
	unsigned char* buf = (unsigned char*)malloc(width * step * num_ranks);
	
	MPI_Allgather(new->pixels+lower_bound, step*width, MPI_UNSIGNED_CHAR,
				buf, step * width, MPI_UNSIGNED_CHAR,
				MPI_COMM_WORLD);

	memcpy(new->pixels, buf, total);
	

	return 0;
}


//original ver.
static int load_jpeg(char *filename, struct image_t *image) {

	FILE *fff;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW output_data;
	unsigned int scanline_len;
	int scanline_count=0;

	fff=fopen(filename,"rb");
	if (fff==NULL) {
		fprintf(stderr, "Could not load %s: %s\n",
			filename, strerror(errno));
		return -1;
	}

	/* set up jpeg error routines */
	cinfo.err = jpeg_std_error(&jerr);

	/* Initialize cinfo */
	jpeg_create_decompress(&cinfo);

	/* Set input file */
	jpeg_stdio_src(&cinfo, fff);

	/* read header */
	jpeg_read_header(&cinfo, TRUE);

	/* Start decompressor */
	jpeg_start_decompress(&cinfo);

	printf("output_width=%d, output_height=%d, output_components=%d\n",
		cinfo.output_width,
		cinfo.output_height,
		cinfo.output_components);

	image->xsize=cinfo.output_width;
	image->ysize=cinfo.output_height;
	image->depth=cinfo.output_components;

	scanline_len = cinfo.output_width * cinfo.output_components;
	image->pixels=malloc(cinfo.output_width * cinfo.output_height * cinfo.output_components);

	while (scanline_count < cinfo.output_height) {
		output_data = (image->pixels + (scanline_count * scanline_len));
		jpeg_read_scanlines(&cinfo, &output_data, 1);
		scanline_count++;
	}

	/* Finish decompressing */
	jpeg_finish_decompress(&cinfo);

	jpeg_destroy_decompress(&cinfo);

	fclose(fff);

	return 0;
}

// 1 loading + broadcasting
// static int load_jpeg(char *filename, struct image_t *image) {

// 	FILE *fff;
// 	struct jpeg_decompress_struct cinfo;
// 	struct jpeg_error_mgr jerr;
// 	JSAMPROW output_data;
// 	unsigned int scanline_len;
// 	int scanline_count=0;

// 	fff=fopen(filename,"rb");
// 	if (fff==NULL) {
// 		fprintf(stderr, "Could not load %s: %s\n",
// 			filename, strerror(errno));
// 		return -1;
// 	}

// 	int current_rank;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
// 	int metadata[3];

// 	if(current_rank == 0){
// 		/* set up jpeg error routines */
// 		cinfo.err = jpeg_std_error(&jerr);

// 		/* Initialize cinfo */
// 		jpeg_create_decompress(&cinfo);

// 		/* Set input file */
// 		jpeg_stdio_src(&cinfo, fff);

// 		/* read header */
// 		jpeg_read_header(&cinfo, TRUE);

// 		/* Start decompressor */
// 		jpeg_start_decompress(&cinfo);

// 		printf("output_width=%d, output_height=%d, output_components=%d\n",
// 			cinfo.output_width,
// 			cinfo.output_height,
// 			cinfo.output_components);

		
// 		metadata[0] = cinfo.output_width; 
// 		metadata[1] = cinfo.output_height;
// 		metadata[2] = cinfo.output_components;

		

// 		image->xsize=cinfo.output_width;
// 		image->ysize=cinfo.output_height;
// 		image->depth=cinfo.output_components;

// 		scanline_len = cinfo.output_width * cinfo.output_components;
// 		image->pixels=malloc(cinfo.output_width * cinfo.output_height * cinfo.output_components);
		

// 		while (scanline_count < cinfo.output_height) {
// 			output_data = (image->pixels + (scanline_count * scanline_len));
// 			jpeg_read_scanlines(&cinfo, &output_data, 1);
// 			scanline_count++;
// 		}

// 		/* Finish decompressing */
// 		jpeg_finish_decompress(&cinfo);

// 		jpeg_destroy_decompress(&cinfo);

// 		fclose(fff);
// 	}
// 	// broadcasting
// 	MPI_Bcast(metadata, 3, MPI_INT, 0, MPI_COMM_WORLD);

// 	if(current_rank){
// 		image->xsize=metadata[0];
// 		image->ysize=metadata[1];
// 		image->depth=metadata[2];
// 		image->pixels = malloc(metadata[0] * metadata[1] * metadata[2]);
// 	}
	

// 	MPI_Bcast(image->pixels, metadata[0] * metadata[1] * metadata[2], MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
	
// 	return 0;
	
// }

static int store_jpeg(char *filename, struct image_t *image) {

	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	int quality=90; /* % */
	int i;

	FILE *fff;

	JSAMPROW row_pointer[1];
	int row_stride;

	/* setup error handler */
	cinfo.err = jpeg_std_error(&jerr);

	/* initialize jpeg compression object */
	jpeg_create_compress(&cinfo);

	/* Open file */
	fff = fopen(filename, "wb");
	if (fff==NULL) {
		fprintf(stderr, "can't open %s: %s\n",
			filename,strerror(errno));
		return -1;
	}

	jpeg_stdio_dest(&cinfo, fff);

	/* Set compression parameters */
	cinfo.image_width = image->xsize;
	cinfo.image_height = image->ysize;
	cinfo.input_components = image->depth;
	cinfo.in_color_space = JCS_RGB;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	/* start compressing */
	jpeg_start_compress(&cinfo, TRUE);

	row_stride=image->xsize*image->depth;

	for(i=0;i<image->ysize;i++) {
		row_pointer[0] = & image->pixels[i * row_stride];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	/* finish compressing */
	jpeg_finish_compress(&cinfo);

	/* close file */
	fclose(fff);

	/* clean up */
	jpeg_destroy_compress(&cinfo);

	return 0;
}

int main(int argc, char **argv) {

	struct image_t image,sobel_x,sobel_y,new_image;
	struct convolve_data_t sobel_data[2];
	double start_time,load_time,store_time,convolve_time,combine_time;
	int result;
	int num_ranks,rank;

	/* Check command line usage */
	if (argc<2) {
		fprintf(stderr,"Usage: %s image_file\n",argv[0]);
		return -1;
	}

	/* Initialize MPI */
	result = MPI_Init(&argc,&argv);
	if (result != MPI_SUCCESS) {
		fprintf (stderr, "Error starting MPI program!.\n");
		MPI_Abort(MPI_COMM_WORLD, result);
	}

	/* Initial time */
	start_time=MPI_Wtime();

	/* Get number of tasks and our process number (rank) */
	MPI_Comm_size(MPI_COMM_WORLD,&num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (rank==0) {
		printf("Number of ranks= %d My rank= %d\n",num_ranks,rank);
	}

	/* Load an image */
	load_jpeg(argv[1],&image);

	load_time=MPI_Wtime();

	/* Allocate space for sobel_x output */
	sobel_x.xsize=image.xsize;
	sobel_x.ysize=image.ysize;
	sobel_x.depth=image.depth;
	sobel_x.pixels=calloc(image.xsize*image.ysize*image.depth,sizeof(char));

	/* Allocate space for sobel_y output */
	sobel_y.xsize=image.xsize;
	sobel_y.ysize=image.ysize;
	sobel_y.depth=image.depth;
	sobel_y.pixels=calloc(image.xsize*image.ysize*image.depth,sizeof(char));

	/* Allocate space for output image */
	new_image.xsize=image.xsize;
	new_image.ysize=image.ysize;
	new_image.depth=image.depth;
	new_image.pixels=calloc(image.xsize*image.ysize*image.depth,sizeof(char));

	/* convolution */
	sobel_data[0].old=&image;
	sobel_data[0].new=&sobel_x;
	sobel_data[0].filter=sobel_x_filter;
	sobel_data[0].ystart=0;
	sobel_data[0].yend=image.ysize;
	generic_convolve((void *)&sobel_data[0]);

	/* Write out sobel_x output (for debugging) */
	// store_jpeg("outx.jpg",&sobel_x);

	sobel_data[1].old=&image;
	sobel_data[1].new=&sobel_y;
	sobel_data[1].filter=sobel_y_filter;
	sobel_data[1].ystart=0;
	sobel_data[1].yend=image.ysize;
	generic_convolve((void *)&sobel_data[1]);

	/* Write out sobel_y output (for debugging) */
	// store_jpeg("outy.jpg",&sobel_y);

	convolve_time=MPI_Wtime();

	/* Combine to form output */
	combine(&sobel_x,&sobel_y,&new_image);

	combine_time=MPI_Wtime();

	/* Write data back out to disk */
	if(rank == 0)store_jpeg("out.jpg",&new_image);

	store_time=MPI_Wtime();

	if (rank==0) {
		printf("Load time: %lf\n",load_time-start_time);
		printf("Convolve time: %lf\n",convolve_time-load_time);
		printf("Combine time: %lf\n",combine_time-convolve_time);
		printf("Store time: %lf\n",store_time-combine_time);
		printf("Total time = %lf\n",store_time-start_time);
	}

	MPI_Finalize();

	return 0;
}
