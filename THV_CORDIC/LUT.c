
#include <stdio.h>
#include <math.h>

#define N_ITERATIONS 16

int main() {
    int j;
    double val_atanh;
    double val_atanh_div_ln2;

//    printf(" j |  2^-j           | atanh(2^-j)           | atanh(2^-j) / ln(2)\n");
    printf(" j | atanh(2^-j) / ln(2)\n");
//    printf("--------------------------------------------------------------------------\n");
    printf("--------------------------\n");

    for (j = 1; j <= N_ITERATIONS; j++) {
        double pow_val = pow(2.0, -j);        // 2^-j
        val_atanh = atanh(pow_val);           // atanh(2^-j)
        val_atanh_div_ln2 = val_atanh / log(2.0);  // normalization
//        printf("%2d | %.12f | %.12f | %.12f\n",
//               j, pow_val, val_atanh, val_atanh_div_ln2);

        if(j<11)      printf("%.15lf,     // j=%d\n", val_atanh_div_ln2, j);
        else if(j<16) printf("%.15e, // j=%d\n", val_atanh_div_ln2, j);
        else          printf("%.15e  // j=%d\n", val_atanh_div_ln2, j);
//         printf("%2d | %.19lf\n", j, val_atanh_div_ln2);
    }

    return 0;
}
