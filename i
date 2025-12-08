// newton_vdw_si.c
#include <stdio.h>
#include <math.h>

#define R 8.314      // J/(mol·K)  (SI)
#define PA_PER_ATM 101325.0

/* Compute van der Waals 'a' and 'b' using critical properties (SI adapted) */
double compute_a(double Tc, double Pc_atm) {
    double Pc = Pc_atm * PA_PER_ATM; // convert atm -> Pa
    // a = 27 * R^2 * Tc^2 / (64 * Pc)
    return (27.0 * R * R * Tc * Tc) / (64.0 * Pc);
}
double compute_b(double Tc, double Pc_atm) {
    double Pc = Pc_atm * PA_PER_ATM;
    // b = R * Tc / (8 * Pc)
    return (R * Tc) / (8.0 * Pc);
}

/* f(V) = (P + a/V^2)*(V - b) - R*T ; P in Pa, V in m^3/mol, a in (J·m^3)/mol^2? units consistent if Pc in Pa */
double f(double P_pa, double T, double V, double a, double b) {
    return (P_pa + a/(V*V)) * (V - b) - R * T;
}

/* derivative f'(V) */
double fprime(double P_pa, double T, double V, double a, double b) {
    // derivative = (P + a/V^2) - (2*a*(V-b))/V^3
    return (P_pa + a/(V*V)) - (2.0 * a * (V - b)) / (V*V*V);
}

int main() {
    char gas_name[100];
    double P_atm, T, Tc, Pc_atm, accuracy;
    int max_iterations;

    printf("Enter gas name: ");
    scanf("%99s", gas_name);
    printf("Enter pressure P (in atm): ");
    scanf("%lf", &P_atm);
    printf("Enter temperature T (in K): ");
    scanf("%lf", &T);
    printf("Enter critical temperature Tc (in K): ");
    scanf("%lf", &Tc);
    printf("Enter critical pressure Pc (in atm): ");
    scanf("%lf", &Pc_atm);
    printf("Enter accuracy (e.g. 1e-6): ");
    scanf("%lf", &accuracy);
    printf("Enter maximum iterations: ");
    scanf("%d", &max_iterations);

    double P_pa = P_atm * PA_PER_ATM;   // convert atm -> Pa
    double a = compute_a(Tc, Pc_atm);
    double b = compute_b(Tc, Pc_atm);

    // Initial guess: ideal gas law in SI (V in m^3/mol)
    double V0 = (R * T) / P_pa;
    double V1 = V0;
    double eps = 1.0;
    int n = 0;

    // Newton-Raphson loop (do-while to match handwritten style)
    do {
        double fv = f(P_pa, T, V0, a, b);
        double fvp = fprime(P_pa, T, V0, a, b);

        if (fabs(fvp) < 1e-16) {
            printf("Derivative too small. Stopping to avoid division by zero.\n");
            return 1;
        }

        V1 = V0 - fv / fvp;               // Newton update
        eps = fabs((V1 - V0) / V0);       // relative change
        V0 = V1;
        n++;
    } while (eps > accuracy && n < max_iterations);

    if (n >= max_iterations && eps > accuracy) {
        printf("Warning: did not converge within %d iterations. eps = %.6e\n", max_iterations, eps);
    } else {
        printf("Converged in %d iterations. Relative error eps = %.6e\n", n, eps);
    }

    printf("Gas: %s\n", gas_name);
    printf("Pressure = %.6f atm (%.3f Pa)\n", P_atm, P_pa);
    printf("Temperature = %.6f K\n", T);
    printf("Molar volume V = %.8e m^3/mol (%.6f L/mol)\n", V1, V1 * 1000.0); // 1 m^3 = 1000 L
    return 0;
}
