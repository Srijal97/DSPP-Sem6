
/*
 * Program to perform Discrete Fourier Transform (DFT) and Inverse DFT (IDFT)
 *
 * Author: Srijal Poojari
 * Contact: srijal97@gmail.com
 * Description: Provides functions for the calculation of DFT and IDFT of a given
 *              signal of any length.
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

//----------------------------------------------------------------------------------//

void dft(int N, struct complex_number x[], struct complex_number X[]) {

    int k, n;

    for(k = 0; k < N; k++) {
        struct complex_number summation;
        summation.real = 0;
        summation.imag = 0;

        for(n = 0; n < N; n++) {
             struct complex_number term = multiply_complex(x[n], twiddle_factor(N, n*k));
             summation = add_complex(summation, term);
        }

        X[k] = summation;

    }
}

void idft(int N, struct complex_number X[], struct complex_number x[]) {

    int n, k;

    for(n = 0; n < N; n++) {
        struct complex_number summation;
        summation.real = 0;
        summation.imag = 0;

        for(k = 0; k < N; k++) {
             struct complex_number term = multiply_complex(X[k], twiddle_factor(N, -n*k));
             summation = add_complex(summation, term);
        }

        summation.real = summation.real / N;
        summation.imag = summation.imag / N;

        x[n] = summation;

    }
}

//----------------------------------------------------------------------------------//

int main()
{

	int i, N;
	
	printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
	scanf("%d", &i);

	printf( "Enter the length of x[n] i.e. N = ");
	scanf("%d", &N);

	struct complex_number x[N];
	struct complex_number X[N];

	printf( "Enter the values of x[n] : \n");

	if (i == 1){
        for(i = 0; i < N; i++) {
		    scanf("%lf", &x[i].real);
		    x[i].imag = 0;
    	}
	}
	else{
        for(i = 0; i < N; i++) {
            printf("Real: ");
		    scanf("%lf", &x[i].real);

		    printf("Imaginary: ");
		    scanf("%lf", &x[i].imag);
    	}
	}

	dft(N, x, X);

	printf("DFT result, X[k]: \n");
	for(i = 0; i < N; i++) {
        printf("X[%d] = %lf + j%lf \n", i, X[i].real, X[i].imag);

	}

	idft(N, X, x);

	printf("Inverse DFT result, x[n]: \n");
	for(i = 0; i < N; i++) {
        printf("x[%d] = %lf + j%lf \n", i, x[i].real, x[i].imag);

	}

	return 0;

}
/*------------------------------------------------------------------------*/





