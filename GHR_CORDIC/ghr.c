#include <stdio.h>
#include <math.h> // Required for pow(), log2(), log(), atanh(), and sqrt()

// Define the number of iterations (Rounds)
#define N 24

// Pre-calculated Look-Up Table (LUT) for GH CORDIC (base 2)
// This is the *same* LUT as GHV
// value = atanh(pow(2, -i)) / log(2)
const double GH_LUT[N] = {
    0.7924812503605780, // i=1
    0.3684827970831031, // i=2
    0.1812850396923541, // i=3
    0.0902861228209104, // i=4 (repeated)
    0.0450989044857891, // i=5
    0.0225439447642690, // i=6
    0.0112712843255441, // i=7
    0.0056355561675101, // i=8
    0.0028177673347164, // i=9
    0.0014088823237399, // i=10
    0.0007044409939180, // i=11
    0.0003522204759650, // i=12
    0.0001761102353583, // i=13 (repeated)
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
 * @brief Calculates the hyperbolic scale factor (Kh)
 * @param iterations The total number of iterations (N).
 * @return The value of Kh.
 */
double calculate_Kh(int iterations) {
    double Kh = 1.0;
    for (int i = 1; i <= iterations; i++) {
        // Kh = product( sqrt(1 - 2^(-2i)) ) 
        Kh *= sqrt(1.0 - pow(2.0, -2.0 * i));
        
        // Account for repeated iterations
        if (i == 4 || i == 13) {
            Kh *= sqrt(1.0 - pow(2.0, -2.0 * i));
        }
    }
    return Kh;
}

/*
 * @brief Performs a single iteration of the GHR CORDIC algorithm.
 * @param i The current iteration index (starts from 1).
 * @param x Pointer to the x-register.
 * @param y Pointer to the y-register.
 * @param z Pointer to the z-register (the angle).
 */
void ghr_iteration(int i, double* x, double* y, double* z) {
    // 1. Determine the direction 'd' from the sign of z
    // d = sign(z^i)
    double d = (*z >= 0) ? 1.0 : -1.0;
    
    // 2. Calculate the shift value (2^-i)
    double shift = pow(2.0, -(double)i);
    
    // 3. Get the pre-calculated constant from the LUT
    double lut_val = GH_LUT[i - 1];
    
    // 4. Store the current x and y values (x^i, y^i)
    double x_old = *x;
    double y_old = *y;
    
    // 5. Calculate the CORDIC equations
    // (Eq. 19a): x_i+1 = x_i + d * (2^-i * y_i)
    *x = x_old + d * shift * y_old;
    
    // (Eq. 19b): y_i+1 = y_i + d * (2^-i * x_i)
    *y = y_old + d * shift * x_old;
    
    // (Eq. 19c): z_i+1 = z_i - d * (atanh(2^-i) / ln(2))
    *z = *z - d * lut_val;
}

int main() {
    // --- 1. Set up initial values ---
    // Let's test by calculating 2^z0
    // We'll use the z0 we got from the last example: log2(1.5) ≈ 0.58496
    double z0 = 0.5849625007211562;

    printf("--- Convergence Range of GHV ---\n");
    printf("|z0| <= 1.6126\n");
    printf("input z0 : ");
    scanf("%lf",&z0);
    // Calculate the inverse scale factor, 1/Kh [cite: 321]
    double Kh = calculate_Kh(N);
    double K_inv = 1.0 / Kh;

    // Initial inputs for GHR CORDIC [cite: 321]
    double x = K_inv; // x0 = 1/Kh
    double y = 0.0;   // y0 = 0
    double z = z0;    // z0 = target exponent

    printf("--- Starting GHR CORDIC Calculation ---\n");
    printf("Calculating 2^(%.10f)\n", z0);
    printf("Initial values: x0=1/Kh=%.10f, y0=%.4f, z0=%.10f\n\n", x, y, z);
    
    // --- 2. Run the CORDIC iterations ---
    for (int i = 1; i <= N; i++) {
        
        ghr_iteration(i, &x, &y, &z);
        printf("i=%2d: x=%.8f, y=%.8f, z=%.8f\n", i, x, y, z);
        
        // ** CRITICAL: Repeat iteration at i = 4 and i = 13 **
        if (i == 4 || i == 13) {
            ghr_iteration(i, &x, &y, &z);
            printf("i=%2d (repeat): x=%.8f, y=%.8f, z=%.8f\n", i, x, y, z);
        }
    }
    
    // --- 3. Finalize and display results ---
    
    // The final M_root is the sum of x_n and y_n 
    double M_root = x + y; 
    
    printf("\n--- Final Results ---\n");
    printf("Calculated x_n          : %.10f\n", x);
    printf("Calculated y_n          : %.10f\n", y);
    printf("Final M_root (x_n + y_n): %.10f\n", M_root);
    
    // Compare against the standard math library function
    double expected_M_root = pow(2.0, z0);
    printf("\nExpected 2^z0 (from math.h): %.10f\n", expected_M_root);

    printf("\n--- Calculation Error ---\n");
    printf("Calculation Diff Error: %.10e\n", M_root - expected_M_root);
    printf("Calculation     %%Error: %.10f\n", (M_root - expected_M_root)/expected_M_root*100);
    return 0;
}
