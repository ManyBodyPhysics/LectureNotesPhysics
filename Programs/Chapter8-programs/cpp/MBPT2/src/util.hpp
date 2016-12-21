/**
 *   Uses Legendre's three-square theorem to determine if a number is the sum
 * of three squares. Runs in time O(log(n))
 *
 * @param n a positive integer
 * @return 1 if n can be expressed as the sum of the squares of three integers,
 *         0 otherwise
 */
int isASumOfThreeSquares (unsigned n);

/**
 *   Uses a very simple algorithm to determine the integer square root of a
 * number.Runs in time O(sqrt(n))
 *
 * @param x a positive integer
 * @return The smallest integer y such that y * y >= x
 */
unsigned smallestSquareRootAtLeast (unsigned x);
