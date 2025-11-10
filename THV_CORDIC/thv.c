/**
 * C implementation of the THV-CORDIC log2(R) algorithm described in
 * "Hyperbolic CORDIC-Based Architecture for Computing Logarithm"
 * (IEEE TRANSACTIONS ON CIRCUITS AND SYSTEMS-II: EXPRESS BRIEFS, VOL. 67, NO. 11, NOVEMBER 2020)
 */

#include <stdio.h>
#include <math.h> // For log2 (true value), frexp, pow, fabs

// Define the number of iterations, e.g., 16 as in the paper's specific case
#define N_ITERATIONS 16

/**
 * Constants array for theta_new_j = tanh_dagger_inv(2^-j)
 * These values are from Table I [cite: 119]
 * We use a 0-based index, so THETA_NEW_CONSTANTS[0] corresponds to j=1
 */
const double THETA_NEW_CONSTANTS[N_ITERATIONS] = {
    0.792481250360578,     // j=1
    0.368482797083103,     // j=2
    0.181285039692354,     // j=3
    0.090286122820910,     // j=4
    0.045098904485789,     // j=5
    0.022543944764269,     // j=6
    0.011271284325544,     // j=7
    0.005635556167510,     // j=8
    0.002817767334716,     // j=9
    0.001408882323740,     // j=10
    7.044409939180080e-04, // j=11
    3.522204759650244e-04, // j=12
    1.761102353582651e-04, // j=13
    8.805511735110164e-05, // j=14
    4.402755863454696e-05, // j=15
    2.201377931214800e-05  // j=16
};

/**
 * Checks if a specific iteration j needs to be repeated.
 * According to the paper, iterations j = 4, 13, 40, ... must be repeated[cite: 108, 122].
 * @param j The current iteration number (1-based).
 * @return 1 if iteration should be repeated, 0 otherwise.
 */
int is_repeated_iteration(int j) {
    if (j == 4 || j == 13 || j == 40) {
        return 1;
    }
    return 0;
}

/**
 * Computes log2(V) using the THV-CORDIC algorithm (Equation 7) .
 * This function is designed to compute log2 of the mantissa V,
 * where V is in the range [1.0, 2.0).
 *
 * @param V The input value (1.M_in), expected to be in [1.0, 2.0).
 * @param n_iterations The total number of iterations to perform.
 * @return The CORDIC approximation of log2(V).
 */
double thv_cordic_log2(double V, int n_iterations) {
    // 1. Initialization
    // As per Section III-A, for log calculation[cite: 76, 110]:
    double x = V + 1.0;
    double y = V - 1.0;
    double z = 0.0;

    // 2. Iteration Loop
    for (int j = 1; j <= n_iterations; j++) {
        
        // Get pre-computed constant for this iteration
        // (Note: j is 1-based, array is 0-based)
        double theta_j = THETA_NEW_CONSTANTS[j - 1];
        
        // d is sign(y_j)
        double d = (y >= 0.0) ? 1.0 : -1.0;
        
        // 2^(-j) term, implemented with pow() for floating-point simulation
        // In hardware, this would be a simple bit-shift.
        double shift_val = pow(2.0, -j);

        // 3. Apply THV-CORDIC iterative formulas (Equation 7) [cite: 96, 97]
        // We use temporary variables to ensure simultaneous updates
        double x_new = x - d * shift_val * y;
        double y_new = y - d * shift_val * x;
        double z_new = z + d * theta_j;

        // Update the values for the next iteration
        x = x_new;
        y = y_new;
        z = z_new;

        // 4. Handle Repeated Iterations [cite: 108, 122]
        if (is_repeated_iteration(j)) {
            // Perform the iteration again with the *same* j values
            d = (y >= 0.0) ? 1.0 : -1.0;
            // shift_val and theta_j are the same as before

            x_new = x - d * shift_val * y;
            y_new = y - d * shift_val * x;
            z_new = z + d * theta_j;

            x = x_new;
            y = y_new;
            z = z_new;
            printf("j = %3d : x = %13.10lf, y = %13.10lf, z = %13.10lf (repeat)\n", j, x, y, z);
        }
        printf("j = %3d : x = %13.10lf, y = %13.10lf, z = %13.10lf\n", j, x, y, z);
    }

    // 5. Final Result
    // The paper states z_n approximates (1/2)log2(V)[cite: 88, 94].
    // We must multiply by 2 ("Shift left one bit")[cite: 84, 193].
    return z * 2.0;
}

/**
 * Computes log2(R) for an arbitrary single-precision floating-point number R.
 * This implements the full method from Section III-B.
 * log2(R) = log2(1.M_in * 2^E_in) = E_in + log2(1.M_in) [cite: 131]
 *
 * @param R The input floating-point number.
 * @param n_iterations The number of iterations to use for the CORDIC module.
 * @return The CORDIC approximation of log2(R).
 */
double compute_log2_R(double R, int n_iterations) {
    // Handle edge case: log of non-positive number is undefined
    if (R <= 0.0) {
        return -INFINITY; // Or NAN
    }

    // 1. Decompose R into Mantissa (M) and Exponent (E)
    // R = M * 2^E
    // C's frexp(R, &E) returns M in [0.5, 1.0) and integer E
    int E_raw;
    double M_raw = frexp(R, &E_raw);

    // 2. Adjust to match the paper's format
    // The paper requires V = (1.M_in) in the range [1.0, 2.0)[cite: 130].
    // We convert R = (M_raw * 2) * 2^(E_raw - 1)
    double V_mantissa = M_raw * 2.0; // This is (1.M_in) in [1.0, 2.0)
    int E_in = E_raw - 1;            // This is the corresponding exponent E_in

    // 3. Compute log2(1.M_in) using the THV-CORDIC module
    double log2_of_mantissa = thv_cordic_log2(V_mantissa, n_iterations);

    // 4. Combine results: log2(R) = E_in + log2(1.M_in) [cite: 131]
    return (double)E_in + log2_of_mantissa;
}

/**
 * Main function to test the THV-CORDIC implementation.
 */
int main() {
    // Array of test values for R
    double R_values;
    int n_iter = N_ITERATIONS;

    printf("--- Testing log2(R) with %d iterations ---\n", n_iter);
    printf("--- (Comparing CORDIC vs. math.h log2) ---\n");

    printf("input R_values : ");
    scanf("%lf", &R_values); 
    double R = R_values;
    // Calculate with our CORDIC function
    double cordic_result = compute_log2_R(R, n_iter);
    
    // Calculate the "true" value using the standard math library
    double true_result = log2(R);
    
    // Calculate the relative error
    double rel_error = fabs(cordic_result - true_result) / fabs(true_result);
    
    printf("\nInput R = %.4f\n", R);
    printf("  True log2(R):   %12.10f\n", true_result);
    printf("  CORDIC log2(R): %12.10f\n", cordic_result);
    printf("  Relative  Error: %9.2e\n", rel_error);
    printf("  Relative %%Error: %.2lf %%\n", rel_error*100);
    
    printf("  (Paper's ARE for n=16 is 2.09e-6) \n");

    return 0;
}
