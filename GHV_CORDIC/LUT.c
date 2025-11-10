
#include <stdio.h>
#include <math.h>

#define N 24

int main() {
    printf("Pre-calculated Look-Up Table (LUT) for GHV CORDIC (base 2)\n");
    printf("GHV_LUT[i-1] = atanh(pow(2, -i)) / log(2)\n\n");
    printf("const double GHV_LUT[%d] = {\n", N);

    for (int i = 1; i <= N; i++) {
        double value = atanh(pow(2.0, -i)) / log(2.0);

        printf("    %.16f, // i=%d", value, i);

        // Mark repeated rounds (hyperbolic CORDIC)
//        if (i == 4 || i == 13) {
//            printf(" (this round is repeated)");
//        } 

        printf("\n");

    }

    printf("};\n");
    return 0;
}
