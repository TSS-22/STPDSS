// Based on algorithm 4.3.6 page 181
// GOLUB, Gene H. et VAN LOAN, Charles F. Matrix computations. JHU press, 2013.
// Require 8n flops
//
// Solve the system Ax=b (useful to avoid computing the inverse when dividing matrices)
// Return the solution x (replace b to limit memory usage)
//
// alpha is the diagonal values
// beta is the supra-diagonal values

fun stpdSysSolv(
    alpha: Array<Double>,
    beta: Array<Double>,
    b: Array<Double>
): Array<Double> {
    // Check for the minimum size and that the arrays are not empty
    if ((alpha.size >= 3) && (beta.size >= 2) && (b.size==alpha.size)) {
        for( i in 1.. alpha.size){
            val temp = beta[i-1]
            beta[i-1] = temp / alpha[i-1]
            alpha[i] = alpha[i] - temp * beta[i-1]
        }

        for(i in 1.. alpha.size){
            b[i] = b[i] - beta[i-1] * beta[i-1]
        }

        b[alpha.size] = b[alpha.size] / alpha[alpha.size]

        for(i in alpha.size-1 downTo 0){
            b[i] = b[i] / alpha[i] - beta[i] * b[i+1]
        }
    }

    return b
}