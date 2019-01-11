#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

/* Opens and reads a binary tile file and returns its contents into a 2x2 
   array.
 
    interface
        function c_get_tile(file, dx, dy, x_offset, y_offset, word_size, tile) bind(C)
            use iso_c_binding, only : c_int, c_char, c_float
            character (c_char), intent(in) :: file
            integer (c_int), intent(in), value :: dx
            integer (c_int), intent(in), value :: dy
            integer (c_int), intent(in), value :: x_offset
            integer (c_int), intent(in), value :: y_offset
            integer (c_int), intent(in), value :: word_size
            real (c_float) :: tile(dx, dy)
        end function c_get_tile
    end interface
 
    integer (c_int) :: dx, dy, x_offset, y_offset, word_size
    integer (c_int) :: tile(dx, dy)

*/

#define GEOG_BIG_ENDIAN 0
#define GEOG_LITTLE_ENDIAN 1

/* int c_get_tile(char *file, int dx, int dy, int *tile[dx][dy])
 *
 * char *file        - The path, realtive or absolute to the file
 * int dx            - The size of the x direction of the tile
 * int dy            - The size of the y direction of the tile
 * int *tile[dx][dy] - The array to hold the tile values on return
 *
 * returns - 1 if success and -1 if file does not exist
 */

int c_get_tile(char *file, 
               int dx, 
               int dy, 
               int x_offset, 
               int y_offset, 
               int word_size, 
               float* tile)
{
    int fd;                /* File Descriptor */
    int i, j;
    unsigned char *buf;
	 size_t numBytes = 0;
    int value = 0;

    int narray = (dx + x_offset) * (dy + y_offset); /* Extent of bytes we have to read */
    float tile_1d[narray];             

    /* Allocate enough space for the entire file */
    buf = (unsigned char *) malloc(sizeof(unsigned char)*(word_size) * narray);

    fd = open(file, O_RDONLY);
    if (fd == -1){
		perror("OPEN ERROR: ");
        printf("C filename was: %s\n", file);
        return -1;
    }

    numBytes = read(fd, buf, word_size * narray); // Read in all the values at once
    if (numBytes == -1){
        perror("READ ERROR: ");
        return -1;
    }
    if(close(fd) == -1){
        perror("CLOSE ERROR: ");
        return -1;
    }


    for(i=0; i < narray ; i++){
        switch(word_size) {
            case 2: 
                /* Shift the first byte read by 8 bytes to make room for the 
                 * 2nd byte we have read.
                 */

                value = (int) (buf[ word_size * i ] << 8 | buf[ word_size * i + 1 ]);

                /* Special case for a negative value. Our sign bit is currently
                 * at the most significant bit (MSB) for an 8-bit value, but it 
                 * needs to bet at the MSB for a 16-bit value.
                 */
                if(buf[word_size * i] >> 7 == 1)
                    value -= 1 << ( 8 * word_size);
                tile_1d[i] = (float) value;
                break;
        }
    }

    for(j=0; j < dx + x_offset; j++){
        for(i=0; i < dy + y_offset; i++){
            /* Place the values into the fortran interoperable array and return */
            tile[i * dx + j] = tile_1d[ (dx + x_offset) * j + i ];
        }
    }

    free(buf);

    return 1;
}
