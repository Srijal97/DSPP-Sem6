
/*
 * Program to perform Fast Fourier Transform (DFT) using Radix 2 DITFFT Algorithm
 *
 * Author: Srijal Poojari
 * Contact: srijal97@gmail.com
 * Description: Provides functions to calculate FFT of any signal of length N (radix 2)
 *
 */

#include<stdio.h>
#include<math.h>


//----------------------------------------------------------------------------------//

struct complex_number {   // structure to store and represent complex numbers
    double real;
    double imag;
};

struct complex_number twiddle_factor(int signal_length, int degree){  // calculates and returns the twiddle factor
    double theta = -2 * M_PI * degree / signal_length;

    // We know, twiddle_factor = exp(-j * 2*pi * n*k / N)
    // exp(theta) = cos(theta) + i*sin(theta)

    struct complex_number Wn = {cos(theta), sin(theta)};

    return Wn;
};

struct complex_number multiply_complex(struct complex_number a, struct complex_number b){  // multiplies 2 complex numbers
    struct complex_number result = {(a.real * b.real - a.imag * b.imag),
                                    (a.real * b.imag + a.imag * b.real)};

    return result;
};

struct complex_number add_complex(struct complex_number a, struct complex_number b){  // adds 2 complex numbers
    struct complex_number result = {(a.real + b.real),
                                    (a.imag + b.imag)};

    return result;
};

struct complex_number sub_complex(struct complex_number a, struct complex_number b){  // subtracts 2 complex numbers
    struct complex_number result = {(a.real - b.real),
                                    (a.imag - b.imag)};

    return result;
};

//----------------------------------------------------------------------------------//

int radix2_greater_than_equal_to(int n){
    int result = 1;

    while( result < n) {
        result = result << 1;
    }

    return result;
}

void fft(int N, struct complex_number a[], struct complex_number A[]) {

    int k, n;

    if (N == 2){
        A[0] = add_complex(a[0], a[1]);
        A[1] = sub_complex(a[0], a[1]);
        return;
    }
    else {
        struct complex_number h[N/2];
        struct complex_number H[N/2];

        struct complex_number g[N/2];
        struct complex_number G[N/2];

        for (n = 0; n < N/2; n++){
            g[n] = a[2*n];
            h[n] = a[2*n + 1];
        }

        fft(N/2, g, G);
        fft(N/2, h, H);

        for (k = 0; k < N; k++){
            A[k] = add_complex(G[k % (N/2)], multiply_complex(twiddle_factor(N, k), H[k % (N/2)]));
        }

        return;
    }
}

void ifft(int N, struct complex_number A[], struct complex_number a[]) {

    // x[n] = (1/N) * FFT(X*[k])*

    int i;

	struct complex_number A_conj[N];  // Conjugate of A[k]
	struct complex_number a_conj[N];  // Conjugate of a[n]

	for(i = 0; i < N; i++) {
        A_conj[i].real = A[i].real;
        A_conj[i].imag = -1 * A[i].imag;
	}

	fft(N, A_conj, a_conj);

	for(i = 0; i < N; i++) {
        a[i].real = a_conj[i].real / N;
        a[i].imag = -1 * a_conj[i].imag / N;
	}

	return;
}

void fast_circular_convolve(int N, struct complex_number x[], struct complex_number h[], struct complex_number y[]) {
    // Find circular convolution of 2 signals of length N. Result will be of length N.

    int N_rad2 = radix2_greater_than_equal_to(N);  // N for radix 2 fft algorithm

    struct complex_number x_new[N_rad2];
    struct complex_number h_new[N_rad2];
    struct complex_number y_new[N_rad2];

    int i;

    for(i = 0; i < N_rad2; i++) {  // Zero padding
        if(i < N) {
            x_new[i].real = x[i].real;
            x_new[i].imag = x[i].imag;

            h_new[i].real = h[i].real;
            h_new[i].imag = h[i].imag;
        }
        else {
            x_new[i].real = 0;
            x_new[i].imag = 0;

            h_new[i].real = 0;
            h_new[i].imag = 0;
        }
    }

    struct complex_number X[N_rad2];
    struct complex_number H[N_rad2];
    struct complex_number Y[N_rad2];

    fft(N_rad2, x_new, X);
    fft(N_rad2, h_new, H);

    for(i = 0; i < N_rad2; i++) {
        Y[i] = multiply_complex(X[i], H[i]);
    }

    ifft(N_rad2, Y, y_new);

    for(i = 0; i < N; i++) {  // copy result to y[n]
        y[i].real = y_new[i].real;
        y[i].imag = y_new[i].imag;
    }

    return;
}

void fast_linear_convolve(int len_x, int len_h, struct complex_number x[], struct complex_number h[], struct complex_number y[]) {
    // Find linear convolution of 2 signals of lengths len_x and len_h. Result will be of length (len_x + len_h - 1).

    int N = len_x + len_h - 1;  // length of result

    struct complex_number x_new[N];
    struct complex_number h_new[N];
    struct complex_number y_new[N];

    int i;

    for(i = 0; i < N; i++) {  // Zero padding
        if(i < len_x) {
            x_new[i].real = x[i].real;
            x_new[i].imag = x[i].imag;
        }
        else {
            x_new[i].real = 0;
            x_new[i].imag = 0;
        }

        if(i < len_h) {
            h_new[i].real = h[i].real;
            h_new[i].imag = h[i].imag;
        }
        else {
            h_new[i].real = 0;
            h_new[i].imag = 0;
        }
    }

    fast_circular_convolve(N, x_new, h_new, y_new);

    for(i = 0; i < N; i++) {  // copy result to y[n]
        y[i].real = y_new[i].real;
        y[i].imag = y_new[i].imag;
    }

    return;
}

//----------------------------------------------------------------------------------//

int main()
{

	int i, j;
	int signal_length, M, N, L;

	printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
	scanf("%d", &i);

	printf( "Length of x[n] = ");
	scanf("%d", &signal_length);

	struct complex_number x[signal_length];

	printf( "Enter the values of x[n] : \n");

	if (i == 1){
        for(i = 0; i < signal_length; i++) {
		    scanf("%lf", &x[i].real);
		    x[i].imag = 0;
    	}
	}
	else{
        for(i = 0; i < signal_length; i++) {
            printf("Real: ");
		    scanf("%lf", &x[i].real);

		    printf("Imaginary: ");
		    scanf("%lf", &x[i].imag);
    	}
	}

	printf( "Is h[n] a real valued signal? (1: Yes, 0: No): ");
	scanf("%d", &i);

	printf( "Length of h[n] = ");
	scanf("%d", &M);

	struct complex_number h[M];

	printf( "Enter the values of h[n] : \n");

	if (i == 1){
        for(i = 0; i < M; i++) {
		    scanf("%lf", &h[i].real);
		    h[i].imag = 0;
    	}
	}
	else{
        for(i = 0; i < M; i++) {
            printf("Real: ");
		    scanf("%lf", &h[i].real);

		    printf("Imaginary: ");
		    scanf("%lf", &h[i].imag);
    	}
	}

    N = 1; // Initial Assumption
    do {
        N = N << 1; // get radix-2 value of N
        L = N - M + 1; // Length of decomposed x[n]
    } while( L < 2);

	int num_of_decompositions = (signal_length / L) + 2;

	printf("Number of decompositions: %d \n", num_of_decompositions);

	int result_length = (num_of_decompositions)*L + M - 1;
	int overlap_length = N-L;

	struct complex_number y[result_length];

	for(i = 0; i < result_length; i++){
        y[i].real = 0;
        y[i].imag = 0;
	}

	for(i = 0; i < num_of_decompositions; i++){
	    struct complex_number x_decom[N];  // decomposed x[n]
        struct complex_number y_decom[N];  // decomposed y[n]

        for(j = 0; j < L; j++) {
            if ((L*i + j) < signal_length) {
                x_decom[overlap_length + j] = x[L*i + j];
            }
            else {
                x_decom[overlap_length + j].real = 0;
                x_decom[overlap_length + j].imag = 0;
            }
        }

        if (i > 0) {
            for(j = 0; j < overlap_length; j++) {
                if ((L*i - overlap_length + j) < signal_length) {
                    x_decom[j] = x[L*i - overlap_length + j];
                }
                else {
                    x_decom[j].real = 0;
                    x_decom[j].imag = 0;
                }
            }
        }

        struct complex_number h_padded[N];

        for(j = 0; j < N; j++) {  // Zero padding
            if(j < M) {
                h_padded[j].real = h[j].real;
                h_padded[j].imag = h[j].imag;
            }
            else {
                h_padded[j].real = 0;
                h_padded[j].imag = 0;
            }
        }

        fast_circular_convolve(N, x_decom, h_padded, y_decom);

        for(j = 0; j < L; j++) {
            y[i*L + j] = y_decom[overlap_length + j];
        }
	}

	printf("OSM result, y[n]: \n");
	for(i = 0; i < result_length; i++) {
        printf("y[%d] = %lf \n", i, y[i].real);

	}

	return 0;

}
/*------------------------------------------------------------------------*/





