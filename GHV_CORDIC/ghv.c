#include <stdio.h>
#include <math.h> // Required for pow(), log2(), log(), and atanh()

// Define the number of iterations (Rounds)
// As per the paper, N=24 is needed for 23-bit precision
#define N 24

// Pre-calculated Look-Up Table (LUT) for GHV CORDIC (base 2)
// The value in the table is GHV_LUT[i-1] = atanh(pow(2, -i)) / log(2)
const double GHV_LUT[N] = {

    0.7924812503605780, // i=1
    0.3684827970831031, // i=2
    0.1812850396923541, // i=3
    0.0902861228209104, // i=4 (this round is repeated)
    0.0450989044857891, // i=5
    0.0225439447642690, // i=6
    0.0112712843255441, // i=7
    0.0056355561675101, // i=8
    0.0028177673347164, // i=9
    0.0014088823237399, // i=10
    0.0007044409939180, // i=11
    0.0003522204759650, // i=12
    0.0001761102353583, // i=13 (this round is repeated)
    0.0000880551173511, // i=14
    0.0000440275586345, // i=15
    0.0000220137793121, // i=16
    0.0000110068896554, // i=17
    0.0000055034448276, // i=18
    0.0000027517224138, // i=19
    0.0000013758612069, // i=20
    0.0000006879306035, // i=21
    0.0000003439653017, // i=22
    0.0000001719826509, // i=23
    0.0000000859913254, // i=24
};

/*
 * @brief Performs a single iteration of the GHV CORDIC algorithm.
 * @param i The current iteration index (starts from 1).
 * @param x Pointer to the x-register.
 * @param y Pointer to the y-register.
 * @param z Pointer to the z-register (the accumulator).
 */
void ghv_iteration(int i, double* x, double* y, double* z) {
    // 1. Determine the direction 'd' from the sign of y
    // d = sign(y^i)
    double d = (*y >= 0) ? 1.0 : -1.0;
    
    // 2. Calculate the shift value (2^-i)
    // In hardware, this is a bit-shift.
    // In software (using double), we use pow().
    double shift = pow(2.0, -(double)i);
    
    // 3. Get the pre-calculated constant from the LUT
    // (C arrays are 0-indexed, but our iteration 'i' starts at 1)
    double lut_val = GHV_LUT[i - 1];
    
    // 4. Store the current x-value (x^i)
    // This is needed to calculate y^(i+1)
    double x_old = *x;
    
    // 5. Calculate the CORDIC equations
    // (Eq. 18a): x_i+1 = x_i - d * (2^-i * y_i)
    *x = *x - d * shift * (*y);
    
    // (Eq. 18b): y_i+1 = y_i - d * (2^-i * x_i)
    *y = *y - d * shift * x_old;
    
    // (Eq. 18c): z_i+1 = z_i + d * (atanh(2^-i) / ln(2))
    *z = *z + d * lut_val;
}

int main() {
    // --- 1. Set up initial values ---
    // We will test by calculating log2(M)
    // Let's use M = 1.5 (we expect log2(1.5) ≈ 0.58496)
    double M = 1.5;

    printf("--- Convergence Range of GHV ---\n");
    printf("|tanh-1(y0/x0)| <= 1.1178\n");
    printf("rage of M = [1/9.36, 9.36] = [0.1068, 9.36]\n");
    printf("input M :");
    scanf(" %lf", &M);

    // Initial inputs for GHV CORDIC (as per the paper)
    // x0 = M + 1, y0 = M - 1, z0 = 0
    double x = M + 1.0;
    double y = M - 1.0;
    double z = 0.0;

    printf("--- Starting GHV CORDIC Calculation ---\n");
    printf("Calculating log2(%.4f)\n", M);
    printf("Initial values: x0=%.4f, y0=%.4f, z0=%.4f\n\n", x, y, z);
    double x0 = x;
    double y0 = y; 
    // --- 2. Run the CORDIC iterations ---
    for (int i = 1; i <= N; i++) {
        
        ghv_iteration(i, &x, &y, &z);
        printf("i=%2d: x=%.8f, y=%.8f, z=%.8f\n", i, x, y, z);
        
        // ** CRITICAL: Repeat iteration at i = 4 and i = 13 **
        // This is required for hyperbolic CORDIC to converge.
        if (i == 4 || i == 13) {
            // Repeat the iteration using the *same* index 'i'
            ghv_iteration(i, &x, &y, &z);
            printf("i=%2d (repeat): x=%.8f, y=%.8f, z=%.8f\n", i, x, y, z);
        }
    }
    
    // --- 3. Finalize and display results ---
    
    // NOTE: Based on the mathematical derivation (using atanh),
    // the raw output z_n from this specific input setup (M+1, M-1)
    // converges to (1/2) * log2(M).
    // Therefore, we must multiply the final z_n by 2.0.
    
    double final_zn = z;
    double log2_M_result = 2.0 * final_zn;
    
    printf("\n--- Final Results ---\n");
    printf("Calculated z_n (raw output)   : %.10f\n", final_zn);
    printf("Calculated log2(M) (z_n * 2)  : %.10f\n", log2_M_result);
    
    // Compare against the standard math library function
    double expected_log2_M = log2(M);
    printf("\nExpected log2(M) (from math.h): %.10f\n", expected_log2_M);
   
    printf("\n--- Error Caluculation ---\n");
    printf("Calculation diff Error: %.10e\n", log2_M_result - expected_log2_M);
    printf("Calculation     %%Error: %.10f\n", (log2_M_result - expected_log2_M)/expected_log2_M*100);

    return 0;
}
